library("ballgown")
cas9 <- list.files("~/m4_cas9_bg/annotated",pattern = "cas9.*", full.names = T)
m4 <- list.files("~/m4_cas9_bg/annotated",pattern = "m4.*", full.names = T)
BG=ballgown(samples = c(cas9, m4), meas = "all")
pData(BG) = data.frame(id=sampleNames(BG), group=rep(c("cas9","m4"), each=3))

eTOt <- data.frame(indexes(BG)$e2t)
tTOg <- data.frame(indexes(BG)$t2g)
eTOg <- merge(eTOt, tTOg, by = "t_id")
eTOg <- eTOg[,-c(1)]
eTOg <- unique(eTOg)
write.table(tTOg, "~/m4_cas9_bg/data/t2g_annotated.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(eTOg, "~/m4_cas9_bg/data/e2g_annotated.txt", quote = F, sep = "\t", row.names = F, col.names = T)

df_gene <- stattest(BG, feature = "gene", meas='FPKM', covariate = 'group', getFC = T)
df_gene[is.na(df_gene)] <- 1
df_gene <- df_gene[,-c(1)]

df_trans = stattest(BG, feature = "transcript", meas = "FPKM", covariate = "group", getFC = T)
df_trans[is.na(df_trans)] <- 1
df_trans <- df_trans[,-c(1)]

df_exon = stattest(BG, feature = "exon", covariate = "group", getFC = T)
df_exon[is.na(df_exon)] <- 1
df_exon <- df_exon[,-c(1)]

df_gene_exp = data.frame(gexpr(BG))
df_gene_exp$cas9_FPKM <- rowMeans(df_gene_exp[, c(1:3)])
df_gene_exp$m4_FPKM <- rowMeans(df_gene_exp[, c(4:6)])
df_gene_exp$id <- rownames(df_gene_exp)
rownames(df_gene_exp) <- NULL
df_gene_exp <- df_gene_exp[,c(ncol(df_gene_exp), 1:(ncol(df_gene_exp) - 1))]
df_gene <- merge(df_gene_exp, df_gene, by = "id")

df_trans_exp = data.frame(texpr(BG, meas = "all"))
colnames(df_trans_exp)[1] <- "id"
df_trans_exp$cas9_FPKM <- rowMeans(df_trans_exp[, c(12,14,16)])
df_trans_exp$m4_FPKM <- rowMeans(df_trans_exp[, c(18,20,22)])
df_trans <- merge(df_trans_exp, df_trans, by = "id")

df_exon_exp = data.frame(eexpr(BG, meas = "all"))
colnames(df_exon_exp)[1] <- "id"
df_exon_exp$cas9_FPKM <- rowMeans(df_exon_exp[, c(8,15,22)])
df_exon_exp$m4_FPKM <- rowMeans(df_exon_exp[, c(29,36,43)])
df_exon <- merge(df_exon_exp, df_exon, by = "id")

gene_table <- read.delim("~/annotations/Dmel_gene_names.txt", header = T)
gene_table <- gene_table[, c(1,2)]

df_gene <- merge(gene_table, df_gene, by = "id")

colnames(df_trans)[9] <- "id"
colnames(df_trans)[1] <- "t_id"
df_trans <- merge(gene_table, df_trans, by = "id")

colnames(df_exon)[1] <- "e_id"
df_exon <- merge(eTOg, df_exon, by = "e_id")
colnames(df_exon)[2] <- "id"
df_exon <- merge(gene_table, df_exon, by = "id")

df_gene_sig <- df_gene[df_gene$qval < .05,]
df_trans_sig <- df_trans[df_trans$qval < .05,]
df_exon_sig <- df_exon[df_exon$qval < .05,]


write.table(df_gene, "~/m4_cas9_bg/data/gene_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df_trans, "~/m4_cas9_bg/data/transcript_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df_exon, "~/m4_cas9_bg/data/exon_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df_gene_sig, "~/m4_cas9_bg/data/significnat_gene_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df_trans_sig, "~/m4_cas9_bg/data/significnat_transcript_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df_exon_sig, "~/m4_cas9_bg/data/significnat_exon_level_differential_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)


