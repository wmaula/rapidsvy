// Single-pass Taylor linearisation variance for stratified one-stage
// (with-replacement) cluster designs, for many domains at once.
//
// Rows are pre-sorted by stratum and PSU (see R/design.R), so one sweep over
// the data gives the PSU totals of every domain. Strata are split across
// threads; each thread keeps its own accumulators and they are summed at the
// end. Formulas follow survey:::onestrat()/onestage() exactly, including the
// lonely PSU rules and the domain conventions of subset.survey.design.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

using namespace Rcpp;

namespace {

enum Lonely { FAIL = 0, REMOVE = 1, CERTAINTY = 2, ADJUST = 3, AVERAGE = 4 };

struct DesignPtr {
  int H;
  const int* strat_ptr;  // H + 1 offsets into the PSU list
  const int* psu_ptr;    // C + 1 offsets into rows
  const int* rows;       // 0-based row ids, grouped by PSU
  const double* fpc_f;   // (N_h - n_h) / N_h, 1 without fpc
};

// Score of row i added into the K-vector acc. Z is n x K column-major and,
// when mult is given, row i of Z is multiplied by mult[i] (used for GLM
// estimating functions so the n x p product never needs to be stored).
// zrow maps a data row to its row in Z when Z holds only domain rows.
struct NumScorer {
  const double* Z;
  const double* mult;
  const int* zrow;
  R_xlen_t n;
  int K;
  inline void add(int i, int, double* acc) const {
    const double* zi = Z + (zrow ? zrow[i] : i);
    if (mult) {
      const double m = mult[i];
      for (int k = 0; k < K; ++k) acc[k] += zi[(R_xlen_t)k * n] * m;
    } else {
      for (int k = 0; k < K; ++k) acc[k] += zi[(R_xlen_t)k * n];
    }
  }
};

// Proportions of a categorical y within domain g:
// score_ik = w_i / W_g * (1{y_i = k} - p_gk); p is G x K row-major.
struct CatScorer {
  const int* y;
  const double* w;
  const double* p;
  const double* invW;
  int K;
  inline void add(int i, int g, double* acc) const {
    const double a = w[i] * invW[g];
    const double* pg = p + (size_t)g * K;
    for (int k = 0; k < K; ++k) acc[k] -= a * pg[k];
    acc[y[i]] += a;
  }
};

// Stratum ranges with roughly equal numbers of rows per thread
std::vector<int> partition(const DesignPtr& d, int T) {
  std::vector<int> b(T + 1, d.H);
  b[0] = 0;
  const double total = d.psu_ptr[d.strat_ptr[d.H]];
  int h = 0;
  for (int t = 1; t < T; ++t) {
    const double target = total * t / T;
    while (h < d.H && d.psu_ptr[d.strat_ptr[h]] < target) ++h;
    b[t] = h;
  }
  return b;
}

template <class Scorer>
struct RecenterWorker {
  const DesignPtr* d;
  const Scorer* sc;
  const int* g;
  int G, K;
  std::vector<double> tot;   // G x K
  std::vector<double> nhsum; // sum of n_h over strata present in domain
  void run(int h0, int h1) {
    tot.assign((size_t)G * K, 0.0);
    nhsum.assign(G, 0.0);
    std::vector<int> last(G, -1);
    for (int h = h0; h < h1; ++h) {
      const int c0 = d->strat_ptr[h], c1 = d->strat_ptr[h + 1];
      for (int r = d->psu_ptr[c0]; r < d->psu_ptr[c1]; ++r) {
        const int i = d->rows[r];
        const int gi = g[i];
        if (gi < 0) continue;
        sc->add(i, gi, &tot[(size_t)gi * K]);
        if (last[gi] != h) {
          last[gi] = h;
          nhsum[gi] += c1 - c0;
        }
      }
    }
  }
};

template <class Scorer>
struct VarWorker {
  const DesignPtr* d;
  const Scorer* sc;
  const int* g;
  int G, K;
  bool full;
  int lonely;
  const double* recenter;  // G x K row-major, only for ADJUST

  std::vector<double> V;   // G x KK, KK = K*K (full) or K (diagonal)
  std::vector<int> npsu, nstr, nok;
  int fail_stratum = -1;

  void run(int h0, int h1) {
    const int KK = full ? K * K : K;
    V.assign((size_t)G * KK, 0.0);
    npsu.assign(G, 0);
    nstr.assign(G, 0);
    nok.assign(G, 0);

    std::vector<int> slot(G, -1), sslot(G, -1);
    std::vector<int> ptg, stg, snsub, cg;
    std::vector<double> pcell, ssum, cS, center, scale, e(K);
    std::vector<char> skip;

    for (int h = h0; h < h1; ++h) {
      const int c0 = d->strat_ptr[h], c1 = d->strat_ptr[h + 1];
      const int nh = c1 - c0;
      stg.clear(); ssum.clear(); snsub.clear(); cg.clear(); cS.clear();

      // PSU totals of each domain present in this stratum
      for (int c = c0; c < c1; ++c) {
        ptg.clear(); pcell.clear();
        for (int r = d->psu_ptr[c]; r < d->psu_ptr[c + 1]; ++r) {
          const int i = d->rows[r];
          const int gi = g[i];
          if (gi < 0) continue;
          int s = slot[gi];
          if (s < 0) {
            s = slot[gi] = (int)ptg.size();
            ptg.push_back(gi);
            pcell.resize(pcell.size() + K, 0.0);
          }
          sc->add(i, gi, &pcell[(size_t)s * K]);
        }
        for (size_t s = 0; s < ptg.size(); ++s) {
          const int gi = ptg[s];
          slot[gi] = -1;
          int ss = sslot[gi];
          if (ss < 0) {
            ss = sslot[gi] = (int)stg.size();
            stg.push_back(gi);
            ssum.resize(ssum.size() + K, 0.0);
            snsub.push_back(0);
          }
          const double* src = &pcell[s * K];
          for (int k = 0; k < K; ++k) ssum[(size_t)ss * K + k] += src[k];
          snsub[ss]++;
          cg.push_back(ss);
          cS.insert(cS.end(), src, src + K);
        }
      }

      const size_t ns = stg.size();
      if (ns == 0) continue;
      center.assign(ns * K, 0.0);
      scale.assign(ns, 0.0);
      skip.assign(ns, 0);
      const double f = d->fpc_f[h];

      for (size_t ss = 0; ss < ns; ++ss) {
        const int gi = stg[ss];
        nstr[gi]++;
        npsu[gi] += snsub[ss];
        if (nh > 1) {
          // centre at the mean over all n_h PSUs (zeros for PSUs outside domain)
          for (int k = 0; k < K; ++k) center[ss * K + k] = ssum[ss * K + k] / nh;
          scale[ss] = (f < 1e-7) ? 0.0 : f * nh / (nh - 1.0);
          nok[gi]++;
        } else if (f < 1e-7) {
          scale[ss] = 0.0;      // census stratum: survey returns zero first
          nok[gi]++;
        } else {
          switch (lonely) {
            case FAIL:
              fail_stratum = h;
              skip[ss] = 1;
              break;
            case REMOVE:
            case CERTAINTY:
              scale[ss] = 0.0;  // centred on itself: contributes zero
              nok[gi]++;
              break;
            case ADJUST:
              for (int k = 0; k < K; ++k) center[ss * K + k] = recenter[(size_t)gi * K + k];
              scale[ss] = (f < 1e-7) ? 0.0 : f;
              nok[gi]++;
              break;
            case AVERAGE:
              skip[ss] = 1;     // NA in survey; rescaled afterwards
              break;
          }
        }
      }

      for (size_t q = 0; q < cg.size(); ++q) {
        const int ss = cg[q];
        if (skip[ss] || scale[ss] == 0.0) continue;
        const double* S = &cS[q * K];
        const double* m = &center[(size_t)ss * K];
        const double sc_ = scale[ss];
        double* Vg = &V[(size_t)stg[ss] * KK];
        if (!full) {
          for (int k = 0; k < K; ++k) {
            const double x = S[k] - m[k];
            Vg[k] += sc_ * x * x;
          }
        } else {
          // upper triangle only; mirrored once at the end
          for (int k = 0; k < K; ++k) e[k] = S[k] - m[k];
          for (int a = 0; a < K; ++a) {
            const double ea = sc_ * e[a];
            double* Va = Vg + a * K;
            for (int b = a; b < K; ++b) Va[b] += ea * e[b];
          }
        }
      }

      // PSUs of the stratum with no member of the domain have total zero
      for (size_t ss = 0; ss < ns; ++ss) {
        sslot[stg[ss]] = -1;
        const int nzero = nh - snsub[ss];
        if (skip[ss] || scale[ss] == 0.0 || nzero <= 0) continue;
        const double* m = &center[ss * K];
        const double sc_ = scale[ss] * nzero;
        double* Vg = &V[(size_t)stg[ss] * KK];
        if (!full) {
          for (int k = 0; k < K; ++k) Vg[k] += sc_ * m[k] * m[k];
        } else {
          for (int a = 0; a < K; ++a)
            for (int b = a; b < K; ++b) Vg[a * K + b] += sc_ * m[a] * m[b];
        }
      }
    }
  }
};

template <class W>
void run_threads(std::vector<W>& workers, const std::vector<int>& b) {
  const int T = (int)workers.size();
  if (T == 1) {
    workers[0].run(b[0], b[1]);
    return;
  }
  std::vector<std::thread> th;
  th.reserve(T);
  for (int t = 0; t < T; ++t)
    th.emplace_back([&workers, &b, t]() { workers[t].run(b[t], b[t + 1]); });
  for (auto& x : th) x.join();
}

template <class Scorer>
List variance_engine(const DesignPtr& d, const Scorer& sc, const int* g, int G,
                     int K, bool full, int lonely, bool any_lonely, int nthreads) {
  int T = std::max(1, std::min(nthreads, d.H));
  std::vector<int> b = partition(d, T);

  std::vector<double> recenter;
  if (lonely == ADJUST && any_lonely) {
    std::vector<RecenterWorker<Scorer>> rw(T);
    for (auto& w : rw) { w.d = &d; w.sc = &sc; w.g = g; w.G = G; w.K = K; }
    run_threads(rw, b);
    std::vector<double> tot((size_t)G * K, 0.0), nhs(G, 0.0);
    for (auto& w : rw) {
      for (size_t j = 0; j < tot.size(); ++j) tot[j] += w.tot[j];
      for (int j = 0; j < G; ++j) nhs[j] += w.nhsum[j];
    }
    recenter.assign((size_t)G * K, 0.0);
    for (int gi = 0; gi < G; ++gi)
      if (nhs[gi] > 0)
        for (int k = 0; k < K; ++k) recenter[(size_t)gi * K + k] = tot[(size_t)gi * K + k] / nhs[gi];
  }

  std::vector<VarWorker<Scorer>> vw(T);
  for (auto& w : vw) {
    w.d = &d; w.sc = &sc; w.g = g; w.G = G; w.K = K; w.full = full;
    w.lonely = lonely; w.recenter = recenter.empty() ? nullptr : recenter.data();
  }
  run_threads(vw, b);

  const int KK = full ? K * K : K;
  NumericMatrix V(G, KK);
  IntegerVector npsu(G), nstr(G), nok(G);
  int fail = -1;
  for (auto& w : vw) {
    if (w.fail_stratum >= 0 && fail < 0) fail = w.fail_stratum;
    for (int gi = 0; gi < G; ++gi) {
      for (int k = 0; k < KK; ++k) V(gi, k) += w.V[(size_t)gi * KK + k];
      npsu[gi] += w.npsu[gi];
      nstr[gi] += w.nstr[gi];
      nok[gi] += w.nok[gi];
    }
  }
  if (full) {
    for (int gi = 0; gi < G; ++gi)
      for (int a = 0; a < K; ++a)
        for (int b = 0; b < a; ++b) V(gi, a * K + b) = V(gi, b * K + a);
  }
  if (lonely == AVERAGE) {
    for (int gi = 0; gi < G; ++gi) {
      const double r = nok[gi] > 0 ? (double)nstr[gi] / nok[gi] : NA_REAL;
      for (int k = 0; k < KK; ++k) V(gi, k) *= r;
    }
  }
  return List::create(_["V"] = V, _["npsu"] = npsu, _["nstrata"] = nstr,
                      _["fail_stratum"] = fail + 1);
}

DesignPtr make_design(const IntegerVector& strat_ptr, const IntegerVector& psu_ptr,
                      const IntegerVector& rows, const NumericVector& fpc_f) {
  DesignPtr d;
  d.H = strat_ptr.size() - 1;
  d.strat_ptr = strat_ptr.begin();
  d.psu_ptr = psu_ptr.begin();
  d.rows = rows.begin();
  d.fpc_f = fpc_f.begin();
  return d;
}

bool has_lonely(const DesignPtr& d) {
  for (int h = 0; h < d.H; ++h)
    if (d.strat_ptr[h + 1] - d.strat_ptr[h] == 1) return true;
  return false;
}

}  // namespace

// [[Rcpp::export]]
List fs_var_num(NumericMatrix Z, Nullable<NumericVector> mult, Nullable<IntegerVector> zrow,
                IntegerVector g, int G,
                IntegerVector strat_ptr, IntegerVector psu_ptr, IntegerVector rows,
                NumericVector fpc_f, int lonely, bool full, int nthreads) {
  DesignPtr d = make_design(strat_ptr, psu_ptr, rows, fpc_f);
  NumScorer sc;
  sc.Z = Z.begin();
  sc.n = Z.nrow();
  sc.K = Z.ncol();
  NumericVector m;
  IntegerVector zr;
  sc.mult = nullptr;
  sc.zrow = nullptr;
  if (mult.isNotNull()) {
    m = NumericVector(mult);
    sc.mult = m.begin();
  }
  if (zrow.isNotNull()) {
    zr = IntegerVector(zrow);
    sc.zrow = zr.begin();
  }
  return variance_engine(d, sc, g.begin(), G, sc.K, full, lonely, has_lonely(d), nthreads);
}

// [[Rcpp::export]]
List fs_var_cat(IntegerVector y, NumericVector w, IntegerVector g, int G, int K,
                NumericMatrix p, NumericVector Wg, IntegerVector strat_ptr,
                IntegerVector psu_ptr, IntegerVector rows, NumericVector fpc_f,
                int lonely, bool full, int nthreads) {
  DesignPtr d = make_design(strat_ptr, psu_ptr, rows, fpc_f);
  // row-major copy of p and inverse domain weight totals
  std::vector<double> prm((size_t)G * K), invW(G);
  for (int gi = 0; gi < G; ++gi) {
    invW[gi] = Wg[gi] > 0 ? 1.0 / Wg[gi] : 0.0;
    for (int k = 0; k < K; ++k) prm[(size_t)gi * K + k] = p(gi, k);
  }
  CatScorer sc;
  sc.y = y.begin();
  sc.w = w.begin();
  sc.p = prm.data();
  sc.invW = invW.data();
  sc.K = K;
  return variance_engine(d, sc, g.begin(), G, K, full, lonely, has_lonely(d), nthreads);
}

// Grouped column sums (long double accumulation, as in colSums)
// [[Rcpp::export]]
NumericMatrix fs_grp_sum(NumericMatrix X, IntegerVector g, int G) {
  const R_xlen_t n = X.nrow();
  const int K = X.ncol();
  std::vector<long double> acc((size_t)G * K, 0.0L);
  for (int k = 0; k < K; ++k) {
    const double* xk = X.begin() + (R_xlen_t)k * n;
    long double* ak = acc.data() + (size_t)k * G;
    for (R_xlen_t i = 0; i < n; ++i) {
      const int gi = g[i];
      if (gi >= 0) ak[gi] += xk[i];
    }
  }
  NumericMatrix out(G, K);
  for (int k = 0; k < K; ++k)
    for (int gi = 0; gi < G; ++gi) out(gi, k) = (double)acc[(size_t)k * G + gi];
  return out;
}

// Weighted quantiles with survey's qrule_math: x and w sorted by x within
// groups (ptr gives G + 1 offsets), P is G x Q matrix of probabilities.
// [[Rcpp::export]]
NumericMatrix fs_wquantile(NumericVector x, NumericVector w, IntegerVector ptr,
                           NumericMatrix P) {
  const int G = ptr.size() - 1, Q = P.ncol();
  NumericMatrix out(G, Q);
  std::vector<double> cumw;
  for (int gi = 0; gi < G; ++gi) {
    const int a = ptr[gi], b = ptr[gi + 1], n = b - a;
    if (n == 0) {
      for (int q = 0; q < Q; ++q) out(gi, q) = NA_REAL;
      continue;
    }
    cumw.resize(n);
    long double s = 0.0L;
    for (int j = 0; j < n; ++j) { s += w[a + j]; cumw[j] = (double)s; }
    const double sw = (double)s;
    for (int q = 0; q < Q; ++q) {
      const double p = P(gi, q);
      if (ISNAN(p)) { out(gi, q) = NA_REAL; continue; }
      const double thr = p * sw;
      // last position with cumw <= p * sum(w); 0 if none (survey:::last)
      int pos = (int)(std::upper_bound(cumw.begin(), cumw.end(), thr) - cumw.begin()) - 1;
      if (pos < 0) pos = 0;
      const int nxt = (pos == n - 1) ? pos : pos + 1;
      const double wlow = p - cumw[pos] / sw;
      out(gi, q) = (wlow <= 0) ? x[a + pos] : x[a + nxt];
    }
  }
  return out;
}

// X' W X and X' W z for IRLS, rows split across threads
// [[Rcpp::export]]
List fs_xtwx(NumericMatrix X, NumericVector w, NumericVector z, int nthreads) {
  const R_xlen_t n = X.nrow();
  const int p = X.ncol();
  const double* Xp = X.begin();
  const double* wp = w.begin();
  const double* zp = z.begin();
  int T = (int)std::max<R_xlen_t>(1, std::min<R_xlen_t>(nthreads, n / 20000 + 1));
  std::vector<std::vector<double>> A(T, std::vector<double>((size_t)p * p, 0.0));
  std::vector<std::vector<double>> B(T, std::vector<double>(p, 0.0));

  auto work = [&](int t) {
    const R_xlen_t i0 = n * t / T, i1 = n * (t + 1) / T;
    const R_xlen_t BS = 1024;
    std::vector<double> blk((size_t)BS * p), zb(BS);
    std::vector<double>& At = A[t];
    std::vector<double>& Bt = B[t];
    for (R_xlen_t s = i0; s < i1; s += BS) {
      const R_xlen_t m = std::min(BS, i1 - s);
      for (R_xlen_t r = 0; r < m; ++r) zb[r] = std::sqrt(wp[s + r]);
      for (int a = 0; a < p; ++a) {
        const double* xa = Xp + (R_xlen_t)a * n + s;
        double* ba = &blk[(size_t)a * BS];
        for (R_xlen_t r = 0; r < m; ++r) ba[r] = xa[r] * zb[r];
      }
      for (R_xlen_t r = 0; r < m; ++r) zb[r] *= zp[s + r];
      for (int a = 0; a < p; ++a) {
        const double* ba = &blk[(size_t)a * BS];
        for (int c = a; c < p; ++c) {
          const double* bc = &blk[(size_t)c * BS];
          double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
          R_xlen_t r = 0;
          for (; r + 3 < m; r += 4) {
            s0 += ba[r] * bc[r]; s1 += ba[r + 1] * bc[r + 1];
            s2 += ba[r + 2] * bc[r + 2]; s3 += ba[r + 3] * bc[r + 3];
          }
          for (; r < m; ++r) s0 += ba[r] * bc[r];
          At[(size_t)a * p + c] += (s0 + s1) + (s2 + s3);
        }
        double t0 = 0;
        for (R_xlen_t r = 0; r < m; ++r) t0 += ba[r] * zb[r];
        Bt[a] += t0;
      }
    }
  };
  if (T == 1) work(0);
  else {
    std::vector<std::thread> th;
    for (int t = 0; t < T; ++t) th.emplace_back(work, t);
    for (auto& x : th) x.join();
  }
  NumericMatrix XtWX(p, p);
  NumericVector XtWz(p);
  for (int t = 0; t < T; ++t) {
    for (int a = 0; a < p; ++a) {
      XtWz[a] += B[t][a];
      for (int c = a; c < p; ++c) XtWX(a, c) += A[t][(size_t)a * p + c];
    }
  }
  for (int a = 0; a < p; ++a)
    for (int c = 0; c < a; ++c) XtWX(a, c) = XtWX(c, a);
  return List::create(_["XtWX"] = XtWX, _["XtWz"] = XtWz);
}
