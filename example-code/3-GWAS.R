library(bigsnpr)
obj.bigsnp <- snp_attach("dosage_ipsych2015.rds")
G <- obj.bigsnp$genotypes

table(y <- obj.bigsnp$fam$affection - 1, exclude = NULL)
table(sex <- obj.bigsnp$fam$sex - 1, exclude = NULL)
sex[sex == -1] <- NA
table(is_2012 <- obj.bigsnp$fam$is_2012, exclude = NULL)

PC <- readRDS("PC.rds")
hist(log_dist <- log(bigutilsr::dist_ogk(PC)))
is_homogeneous <- log_dist < 4.5

str(COVAR <- cbind(PC, sex, is_2012))
mean(has_no_na <- complete.cases(COVAR)) # 98.8%

rel <- readRDS("rel.rds")
is_rel <- obj.bigsnp$fam$sample.ID %in% subset(rel, KINSHIP > 2^-3.5)$IID2
mean(is_rel)

mean(KEEP <- !is_rel & has_no_na & is_homogeneous) # 80.2%
ind_keep <- which(KEEP)

# intervals <- bigparallelr::split_len(ncol(G), nb_split = 20)
intervals <- bigparallelr::split_len(ncol(G), block_len = 500e3)

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("gwas-results")

res_files <- paste0("gwas-results/part", rows_along(intervals), ".rds")
mean(already_done <- file.exists(res_files))

future.apply::future_lapply(which(!already_done), function(ic) {

  ind <- seq(intervals[ic, "lower"], intervals[ic, "upper"])

  # Lin is faster but you might want to use Log instead
  gwas <- bigstatsr::big_univLinReg(
    G, y[ind_keep], ind.train = ind_keep, ind.col = ind,
    covar.train = COVAR[ind_keep, ], ncores = NCORES)

  saveRDS(gwas, res_files[ic])
}) # ~20 min for each part

gwas <- do.call("rbind", lapply(res_files, readRDS))

plot(gwas)
snp_qq(gwas) + ggplot2::xlim(3, NA)

CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS
snp_manhattan(gwas, CHR, POS, npoints = 50e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)

lpval <- -predict(gwas)
ind_top <- snp_clumping(G, infos.chr = CHR, infos.pos = POS,
                        thr.r2 = 0.01, size = 10e3,
                        S = lpval, ind.row = ind_keep,
                        exclude = which(lpval < -log10(5e-8)))
obj.bigsnp$map[ind_top, c(1:3, 7:8, 5, 10, 6, 11)]
#         CHR        SNP       POS a1 a2 info_2012 info_2015 freq_2012 freq_2015
# 1791856   3 rs55687416 117881468  G  A    0.9820    0.9953    0.8173    0.8167
# 2840893   5  rs6451675  43110855  C  G    0.9818    1.0000    0.3185    0.3184
# 2962460   5  rs4916723  87854395  A  C    0.9280    0.9607    0.5693    0.5709
# 3139582   5  rs4256374 144512126  C  T    0.9867    0.9764    0.4831    0.4829
# 5569246  10  rs1021363 106610839  A  G    0.9995    0.9961    0.3646    0.3648
# 8421756  20   rs809220  21414037  T  C    1.0050    0.9766    0.7578    0.7599
