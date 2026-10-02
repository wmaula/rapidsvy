# rapidsvy

<!-- badges: start -->
[![R-CMD-check](https://github.com/wmaula/rapidsvy/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/wmaula/rapidsvy/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

*English:* fast Taylor linearisation for stratified one-stage cluster
survey designs, with results identical to the `survey` package and a
`srvyr`/dplyr interface. All `group_by()` domains are estimated in one
multithreaded C++ pass.

Estimasi berbasis desain survei kompleks yang cepat untuk R, dengan antarmuka
gaya `srvyr`/dplyr. Hasilnya identik dengan package `survey` (diuji sampai
toleransi 1e-8) untuk desain stratifikasi klaster satu tahap (PSU dengan
pengembalian, fpc opsional), yaitu desain yang biasa dipakai untuk data
sampel BPJS Kesehatan.

## Mengapa `survey`/`srvyr` lambat

Masalah utamanya bukan karena R hanya memakai satu core, melainkan algoritma:

* `svyby()` dan `group_by() |> summarise()` di `srvyr` memotong desain per
  domain lalu menghitung ulang variansi **seluruh desain** untuk tiap domain.
  Biaya total sekitar (jumlah domain) x (ukuran data). Tabel 34 provinsi x 5
  segmen berarti ratusan kali pemrosesan 2 juta baris.
* `survey_prop()` di `srvyr` memanggil `svyciprop()`, yang menjalankan
  `svyglm()` untuk **setiap sel** tabel.
* `svydesign()` sendiri menyalin dan memeriksa data berulang kali.

`rapidsvy` menghitung total per PSU untuk **semua domain sekaligus** dalam satu
kali lintasan data di C++, lalu membagi strata ke beberapa thread. Biaya tidak
lagi bertambah dengan jumlah domain. Paralelisme menambah percepatan, tetapi
percepatan terbesar datang dari algoritma.

## Benchmark pada data sampel BPJS 2020

Data `01_kepesertaan` (1.971.744 peserta, 747.221 PSU keluarga, 2.726 strata
kab/kota faskes x segmen, `nest = TRUE`, `lonely_psu = "adjust"`), Apple M1 Pro
8 core, 16 GB RAM, 7 thread. Waktu dalam detik.

| Tugas | survey / srvyr | rapidsvy | Selisih hasil |
|---|---|---|---|
| Membuat desain | 11,5 | 0,25 | |
| Proporsi segmen per provinsi (34 x 6) | > 900 (dihentikan); 88 per provinsi | 0,25 | SE identik (beda < 1e-16) |
| Sama, dengan `srvyr` | > 900 (dihentikan) | 0,25 | |
| Median umur per provinsi | > 900 (dihentikan); 88 per provinsi | 0,24 | identik |
| Rerata umur per kab/kota (515 domain) | > 900 (dihentikan) | 0,16 | |
| `svychisq` segmen x jenis kelamin | 434 | 0,26 | F dan df identik |
| `svyglm` logistik, 36 parameter | 460 | 2,9 | SE beda relatif 2e-12 |

Angka "per provinsi" adalah waktu satu domain di `survey`; dikalikan 34
provinsi berarti sekitar 50 menit per tabel (perkiraan, bukan pengukuran
langsung). Mesin juga sedang memakai swap, yang memperlambat kedua package.

## Instalasi

```r
# install.packages("pak")
pak::pak("wmaula/rapidsvy")

# atau
# install.packages("remotes")
remotes::install_github("wmaula/rapidsvy")
```

Package dikompilasi saat instalasi, jadi butuh compiler C++17: Xcode Command
Line Tools di macOS (`xcode-select --install`), Rtools di Windows. Tidak butuh
OpenMP atau gfortran.

## Contoh dengan data sampel BPJS

```r
library(dplyr)
library(haven)
library(rapidsvy)

peserta <- read_dta("01_kepesertaan.dta") |>
  mutate(prov = as_factor(PSTV09), seg = as_factor(PSTV08),
         sex = as_factor(PSTV05), kelas3 = PSTV07 == 3,
         umur = as.numeric(as.Date("2016-12-31") - PSTV03) / 365.25)

# keluarga sebagai klaster, strata kab/kota faskes x segmen
des <- as_fsvy(peserta, ids = PSTV02, strata = c(PSTV14, seg),
               weights = PSTV15, nest = TRUE, lonely_psu = "adjust")

# proporsi segmen kepesertaan per provinsi (persentase, SE, CI logit)
fs_tab(des, seg, by = prov, percent = TRUE)

# hal yang sama dengan gaya srvyr
des |>
  group_by(prov, seg) |>
  summarise(n = fs_n(), p = fs_prop(vartype = c("se", "ci")))

# rerata, total, rasio, median per domain; filter() = analisis subpopulasi
des |>
  filter(sex == "PEREMPUAN") |>
  group_by(prov) |>
  summarise(kelas3 = fs_mean(kelas3, vartype = "ci", deff = TRUE),
            N = fs_total(),
            umur_median = fs_median(umur))

# uji Rao Scott (svychisq) dan regresi (svyglm)
fs_chisq(des, seg, sex)
fit <- fs_glm(des, kelas3 ~ umur + sex + prov, family = quasibinomial())
tidy(fit, conf.int = TRUE, exponentiate = TRUE)
```

**Perhatian soal strata.** Di `rapidsvy`, `strata = c(PSTV14, seg)` berarti
strata silang (kab/kota x segmen). Di `survey`, `strata = ~PSTV14 + seg`
**tidak** berarti silang: kolom kedua dibaca sebagai strata tahap kedua,
sehingga variansi hanya memakai PSTV14. Padanan yang benar di `survey` adalah
`strata = ~interaction(PSTV14, seg)`.

Desain dari `survey` atau `srvyr` yang sudah ada bisa dikonversi langsung:
`as_fsvy(svydesign_object)` atau `as_fsvy(tbl_svy_object)`.

## Fungsi

| rapidsvy | Padanan survey / srvyr |
|---|---|
| `as_fsvy()` | `svydesign()` / `as_survey_design()` |
| `fs_mean()`, `fs_total()`, `fs_ratio()` | `svymean`, `svytotal`, `svyratio` / `survey_mean`, `survey_total`, `survey_ratio` |
| `fs_prop()` | `svymean` pada faktor, `svyciprop` / `survey_prop` |
| `fs_quantile()`, `fs_median()` | `svyquantile` (qrule "math", interval Woodruff) |
| `fs_n()` | `unweighted(n())` |
| `fs_tab()` | `svytable` + `svyby(svymean)` dalam satu tabel |
| `fs_chisq()` | `svychisq` (statistik "F" dan "Chisq") |
| `fs_glm()` + `tidy()` | `svyglm` |
| `fs_degf()` | `degf` |

Argumen `vartype` ("se", "ci", "var", "cv"), `level`, `deff`, `df`, dan
`na.rm` bekerja seperti di `srvyr`. Kolom keluaran juga dinamai seperti `srvyr`
(`p`, `p_se`, `p_low`, `p_upp`, `p_deff`).

## Aturan statistik yang perlu diketahui

* **Lonely PSU**: bawaan `"fail"`, sama dengan `survey`. Pilihan `"adjust"`,
  `"average"`, `"remove"`, `"certainty"` memberi hasil identik dengan
  `options(survey.lonely.psu = ...)`.
* **Domain**: `filter()` tidak membuang baris dari desain; baris di luar
  domain tetap dihitung sebagai PSU bernilai nol, seperti `subset()` di
  `survey`. Jangan memfilter data sebelum `as_fsvy()` bila tujuannya analisis
  subpopulasi.
* **Derajat bebas CI**: bawaan `fs_degf()` desain (jumlah PSU dikurangi jumlah
  strata), seperti `srvyr`. Pakai `df = Inf` untuk interval normal.
* **Proporsi 0 atau 1**: CI logit tidak terdefinisi; `rapidsvy` mengembalikan
  interval selebar nol pada nilai titik.
* `svyglm` di `survey` mengembalikan `NA` untuk koefisien yang alias;
  `fs_glm` mengikuti aturan yang sama.

## Belum didukung

Desain multitahap dengan variansi per tahap (hanya tahap pertama yang dipakai,
yaitu pendekatan ultimate cluster), bobot replikasi (bootstrap, jackknife, BRR),
poststratifikasi dan kalibrasi (`postStratify`, `calibrate`, `rake`), serta
`svyby` untuk model regresi per kelompok.

## Jumlah thread

`options(rapidsvy.threads = 4)`. Bawaan: jumlah core fisik dikurangi satu.
Hasil tidak bergantung pada jumlah thread.

## Validasi

`tests/testthat/` membandingkan setiap fungsi dengan `survey` pada data
simulasi yang memuat strata dengan satu PSU, nilai hilang, domain, fpc, dan
`nest = TRUE`. Benchmark data BPJS ada di `inst/bench/bench_bpjs.R`.
