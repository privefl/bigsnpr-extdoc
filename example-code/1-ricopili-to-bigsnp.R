library(bigreadr)
library(tidyverse)
library(bigsnpr)

# file.symlink("~/iPSYCH2015/HRC_Imputed/Ricopili/Autosomes/iPSYCH2012_HRC_2020-05-14/dasuqc1_iPSYCH2012_HRC_2020-05-14.hg19.ch.fl", "ipsych_2012")
# file.symlink("~/iPSYCH2015/HRC_Imputed/Ricopili/Autosomes/iPSYCH2015i_HRC_2020-02-07/dasuqc1_iPSYCH2015i_HRC_2020-02-07.hg19.ch.fl", "ipsych_2015")


#### GET SIZES ####

fam_2012 <- fread2("ipsych_2012/qc1/dos_iPSYCH2012_HRC_2020-05-14.hg19.ch.fl.chr1_000_023.out.dosage.fam")
fam_2015 <- fread2("ipsych_2015/qc1/dos_iPSYCH2015i_HRC_2020-02-07.hg19.ch.fl.chr1_000_023.out.dosage.fam")

ind_info <- tibble(
  info_file_2012 = list.files("ipsych_2012/info", full.names = TRUE),
  info_file_2015 = list.files("ipsych_2015/info", full.names = TRUE)
) %>%
  mutate(
    ind_2012 = map(info_file_2012, ~ {
      cat(".")
      info <- bigreadr::fread2(.)
      which(info$pass == 1)
    }),
    ind_2015 = map(info_file_2015, ~ {
      cat(".")
      info <- bigreadr::fread2(.)
      which(info$pass == 1)
    })
  )
ind_info <- ind_info[gtools::mixedorder(ind_info$info_file_2012), ]

len <- pmap_int(ind_info[3:4], function(ind_2012, ind_2015) {
  length(intersect(ind_2012, ind_2015))
})
sum(len)  # 8,785,478
(offset <- cumsum(c(0, len)))


#### CREATE ACTUAL DATA BASED ON SIZES ####

fam <- rbind(fam_2012, fam_2015) %>% setNames(bigsnpr:::NAMES.FAM)
str(fam)
(N <- nrow(fam))  # 134,677
sum(len) * 1 * N / 1024^3  # 1102 GB

# Prepare Filebacked Big Matrix (FBM)
G <- FBM.code256(
  nrow = N,
  ncol = sum(len),
  code = CODE_DOSAGE,
  backingfile = "dosage_ipsych2015",
  create_bk = !file.exists("dosage_ipsych2015.bk")
)

# Write dosages to FBM
plink2 <- download_plink2("tmp-data", AVX2 = FALSE)
bigassertr::assert_dir("log")
bigassertr::assert_dir("map")

get_dosage <- function(info_file, ind, NCORES) {

  tmp <- paste(tempfile(tmpdir = "tmp-data"), Sys.getpid(), sep = "_")
  on.exit(file.remove(paste0(tmp, c(".log", ".raw"))))

  gz_file <- sub("^(.*/)info(.*)info$", "\\1qc1\\2gz", info_file)
  system(glue::glue(
    "{plink2} --import-dosage {gz_file} format=2",
    " --fam {sub('\\\\.gz$', '.fam', gz_file)}",
    " --recode A --out {tmp}",
    " --threads {NCORES} --memory {round(128 / 16 * NCORES * 0.9 * 1000)}"
  ), ignore.stdout = TRUE)

  dosage <- bigreadr::fread2(paste0(tmp, ".raw"),
                             select = ind + 6L, colClasses = "double",
                             nThread = NCORES, showProgress = FALSE)

  for (j in seq_along(dosage)) {
    dosage[[j]] <- as.raw(207 - round(dosage[[j]] * 100))  # needed to reverse
  }

  dosage
}

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))

map_files <- paste0("map/map", rows_along(ind_info), ".rds")
mean(already_done <- file.exists(map_files))

# future.apply::future_lapply(order(len)[1:4], function(ic) {
future.apply::future_lapply(sample(which(!already_done)), function(ic) {

  plink2 <- plink2

  # ic <- which.min(len)
  ind_2012 <- ind_info$ind_2012[[ic]]
  ind_2015 <- ind_info$ind_2015[[ic]]
  info_file_2012 <- ind_info$info_file_2012[ic]
  info_file_2015 <- ind_info$info_file_2015[ic]

  dosage_2012 <- get_dosage(info_file_2012, which(ind_2012 %in% ind_2015), NCORES)
  stopifnot(nrow(dosage_2012) == nrow(fam_2012))
  stopifnot(ncol(dosage_2012) == len[ic])
  G[rows_along(fam_2012), offset[ic] + cols_along(dosage_2012)] <- dosage_2012
  rm(dosage_2012); gc()

  dosage_2015 <- get_dosage(info_file_2015, which(ind_2015 %in% ind_2012), NCORES)
  stopifnot(nrow(dosage_2015) == nrow(fam_2015))
  stopifnot(ncol(dosage_2015) == len[ic])
  G[nrow(fam_2012) + rows_along(fam_2015),
    offset[ic] + cols_along(dosage_2015)] <- dosage_2015
  rm(dosage_2015); gc()

  info_2012 <- dplyr::filter(bigreadr::fread2(info_file_2012), pass == 1)
  info_2015 <- dplyr::filter(bigreadr::fread2(info_file_2015), pass == 1)
  map <- dplyr::full_join(info_2012[which(ind_2012 %in% ind_2015), ],
                          info_2015[which(ind_2015 %in% ind_2012), ],
                          by = c("CHR", "SNP", "POS", "a1", "a2"),
                          suffix = c("_2012", "_2015"))
  saveRDS(map, map_files[ic])
})

map <- do.call("rbind", lapply(map_files, readRDS))
all(map$pass_2012 == 1)
all(map$pass_2015 == 1)
map$pass_2012 <- map$pass_2015 <- NULL
fam$is_2012 <- rep(1:0, c(nrow(fam_2012), nrow(fam_2015)))

# Verif freq (whether reversed mostly)
ind <- sample(ncol(G), 100)
freq <- colMeans(G[, ind]) / 2
plot(freq, map$freq_2012[ind], pch = 20)
points(freq, map$freq_2015[ind], col = "red")
abline(0, 1, col = "chartreuse")

# Create the bigSNP object and save it
G$is_read_only <- TRUE
obj.bigsnp <- structure(
  list(genotypes = G, map = map, fam = fam),
  class = "bigSNP")
str(obj.bigsnp)
saveRDS(obj.bigsnp, file = sub_bk(G$bk, ".rds"))

# Make the files read-only
Sys.chmod("dosage_ipsych2015.bk",  "0444")
Sys.chmod("dosage_ipsych2015.rds", "0444")
