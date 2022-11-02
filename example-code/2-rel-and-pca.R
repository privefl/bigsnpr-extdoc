library(bigsnpr)

obj.bigsnp <- snp_attach("dosage_ipsych2015.rds")
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS
(NCORES <- nb_cores())

# Restrict to variants on both chips
ind_chip <- with(obj.bigsnp$map,
                 which(genotyped_2012 == 1 & genotyped_2015 == 1))

# First quick round of PCA
system.time(
  obj.svd <- snp_autoSVD(G, CHR, POS, ind.col = ind_chip, k = 5,
                         max.iter = 1, ncores = NCORES)
)
plot(obj.svd)
plot(obj.svd, type = "scores")
plot(obj.svd, type = "scores", scores = 4:5)

# Get variants associated with pop struct, discard them and write bed file
obj.pcadapt <- snp_pcadapt(G, obj.svd$u, ind.col = ind_chip, ncores = NCORES)
plot(obj.pcadapt, type = "Manhattan")

length(ind_keep <- ind_chip[predict(obj.pcadapt, log10 = FALSE) > 0.05])
obj.bigsnp2 <- obj.bigsnp
obj.bigsnp2$map <- dplyr::transmute(
  obj.bigsnp$map, chromosome = CHR, marker.ID = SNP, genetic.dist = 0,
  physical.pos = POS, allele1 = a1, allele2 = a2)
snp_writeBed(obj.bigsnp2, bedfile = "tmp-data/ipsych2015.bed", ind.col = ind_keep)

# Compute the relatedness
plink2 <- download_plink2("tmp-data")
rel <- runonce::save_run(
  snp_plinkKINGQC(plink2, bedfile.in = "tmp-data/ipsych2015.bed",
                  thr.king = 2^-4.5, make.bed = FALSE, ncores = NCORES),
  file = "rel.rds"
)
dim(rel)  # 189,742 pairs
hist(rel$KINSHIP, breaks = 100); abline(v = 2^-(1.5:4.5), col = "red")

library(dplyr)
rel2 <- rbind(data.frame(IID = rel$IID1, K = rel$KINSHIP),
              data.frame(IID = rel$IID2, K = rel$KINSHIP)) %>%
  group_by(IID) %>%
  summarise(sum_K = sum(K))

rel3 <- subset(rel2, sum_K > 2^-3.5)
hist(rel2$sum_K, breaks = 100)
is_rel <- obj.bigsnp$fam$sample.ID %in% rel3$IID
sum(!is_rel)  # 97,013

# Info on individuals (country of origin of parents)
info <- bigreadr::fread2("~/Register/FromCrome/stam2012f.csv")
dim(info)
str(info)

region_father <- info$region_birth_f[match(obj.bigsnp$fam$family.ID, info$pid)]
region_mother <- info$region_birth_m[match(obj.bigsnp$fam$family.ID, info$pid)]
table(region_father, exclude = NULL)
parents_from_DK <- region_father == "Denmark" & region_mother == "Denmark"

# Sub-sampling trick, see https://doi.org/10.1093/bioinformatics/btaa520
ind_row <- c(sample(intersect(which(!is_rel), which(parents_from_DK)), 10e3),
             setdiff(which(!is_rel), which(parents_from_DK)))
length(ind_row)  # 47,550

system.time(
  obj.svd2 <- snp_autoSVD(G, CHR, POS,
                          ind.row = ind_row,
                          ind.col = ind_chip,
                          k = 30, ncores = NCORES)
) # 38 min

plot(obj.svd2)

library(ggplot2)
source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(plotlist = lapply(1:12, function(k) {
  plot(obj.svd2, type = "scores", scores = 2 * k - 1:0, coeff = 0.5) +
    aes(color = region_father[ind_row],
        alpha = I(ifelse(is.na(region_father[ind_row]) |
                           region_father[ind_row] == "Unknown", 0, 1))) +
    labs(color = "Region of birth of father")
}))
plot(obj.svd2, type = "loadings", loadings = 11:20, coeff = 0.5)


# Project remaining individuals
proj <- big_apply(G, function(X, ind, obj.svd) {
  X_scaled <- scale(X[ind, attr(obj.svd, "subset")],
                    center = obj.svd$center, obj.svd$scale)
  bigutilsr::pca_OADP_proj(X = X_scaled, loadings = obj.svd$v[, 1:20],
                           sval = obj.svd$d)$OADP_proj
}, ind = rows_along(G)[-ind_row], obj.svd = obj.svd2, a.combine = "rbind",
ncores = NCORES, block.size = 1000)

PC <- matrix(0, nrow(G), 20)
PC[ind_row, ] <- predict(obj.svd2)[, 1:20]
PC[-ind_row, ] <- proj

# Verif projection is okay
plot_grid2(plotlist = lapply(1:10, function(k) {
  qplot(PC[, 2 * k - 1], PC[, 2 * k]) +
    theme_bigstatsr() +
    aes(color = !rows_along(G) %in% ind_row,
        alpha = I(ifelse(rows_along(G) %in% ind_row, 0.3, 1))) +
    labs(color = "Is projected?") +
    scale_colour_viridis_d(direction = -1)
}))

saveRDS(PC, "PC.rds")
