# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ballgown")
# BiocManager::install("genefilter")
# BiocManager::install("dplyr")
# BiocManager::install("devtools")
# devtools::install_github('alyssafrazee/RSkittleBrewer')

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(RSkittleBrewer)

setwd("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/")

# read in sample treatment key
pheno = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/pheno_data.csv")
# arrange pheno data sample order to  bg object sample order
pheno_sort = pheno[c(5:8,1:4,13:16,9:12),]
pheno_data = data.frame(id = paste(pheno_sort$Sample,"_abund",sep=""), trt = pheno_sort$Trt,day = pheno_sort$Day,rep = pheno_sort$Rep,sample.id =pheno_sort$Sample)


# abund_dir = "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/"
# 
# bg = ballgown(dataDir=abund_dir, samplePattern='S*_abund', meas='all', pData=pheno_data)
# # examine ballgown object:
# bg

# add phenodata
# pData(bg) <- data.frame(treatment=pheno$Trt[c(5:8,1:4,13:16,9:12)], day=pheno$Day[c(5:8,1:4,13:16,9:12)])
# pData(bg) <- data.frame(pData(bg), rep=pheno$Rep[c(5:8,1:4,13:16,9:12)], id=pheno$Sample[c(5:8,1:4,13:16,9:12)])

# check which rows correspond to the 5d treatment
pheno_data[1:8,]

five = c("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S14_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S15_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S16_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S17_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S18_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S19_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S20_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S21_abund/")
bg_5d = ballgown(samples=five, meas='all', pData=pheno_data[1:8,])


# check which rows correspond to the 7d treatment
pheno_data[9:16,]

seven = c("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S22_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S23_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S24_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S25_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S26_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S27_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S28_abund/",
         "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/S29_abund/")

bg_7d = ballgown(samples=seven, meas='all', pData=pheno_data[9:16,])


# save object to run interactively
save(bg_5d, file='bg_5d.with.pheno.rda')
save(bg_7d, file='bg_5d.with.pheno.rda')

# A ballgown object has six slots: structure, expr, indexes, dirs, mergedDate, and meas.
# structure(bg)$exon

# All of the following are valid ways to extract expression data from the bg ballgown object:
# transcript_fpkm = texpr(bg, 'FPKM')
# transcript_tpm = texpr(bg, 'TPM')
# transcript_cov = texpr(bg, 'cov')
transcript_fpkm_5d = texpr(bg_5d, 'FPKM')
transcript_fpkm_7d = texpr(bg_7d, 'FPKM')

# whole_tx_table = texpr(bg, 'all')
# exon_mcov = eexpr(bg, 'mcov')
# exon_rcount = eexpr(bg, 'rcount')
# junction_rcount = iexpr(bg, 'rcount')
# whole_intron_table = iexpr(bg, 'all')
# gene_expression = gexpr(bg)

# save raw results
write.csv(transcript_fpkm_5d,"Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/transcript_fpkm_5d.csv")
write.csv(transcript_fpkm_7d,"Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/transcript_fpkm_7d.csv")



texpr_5d_filt = subset(bg_5d,"rowVars(texpr(bg_5d)) >1",genomesubset=TRUE)
texpr_7d_filt = subset(bg_7d,"rowVars(texpr(bg_7d)) >1",genomesubset=TRUE)

results5d_transcripts_trt_adjrep_FC_FPKM = stattest(bg_5d, feature="transcript",covariate="trt",adjustvars = c("rep"), getFC=TRUE, meas="FPKM")
results7d_transcripts_trt_adjrep_FC_FPKM = stattest(bg_7d, feature="transcript",covariate="trt",adjustvars = c("rep"), getFC=TRUE, meas="FPKM")

results5d_transcripts_trt_adjrep_FC_FPKM_sorted = arrange(results5d_transcripts_trt_adjrep_FC_FPKM,qval)
results7d_transcripts_trt_adjrep_FC_FPKM_sorted = arrange(results7d_transcripts_trt_adjrep_FC_FPKM,qval)

results5d_transcripts_trt_adjrep_FC_FPKM_sorted_df = data.frame(geneIDs=ballgown::geneIDs(bg_5d), results5d_transcripts_trt_adjrep_FC_FPKM_sorted)
results7d_transcripts_trt_adjrep_FC_FPKM_sorted_df = data.frame(geneIDs=ballgown::geneIDs(bg_7d), results7d_transcripts_trt_adjrep_FC_FPKM_sorted)


# results5d_filt_transcripts_trt_adjrep_FC_FPKM = stattest(subset(bg_5d,"rowVars(texpr(bg_5d)) >1",genomesubset=TRUE), feature="transcript",covariate="trt",adjustvars = c("rep"), getFC=TRUE, meas="FPKM")
# results7d_filt_transcripts_trt_adjrep_FC_FPKM = stattest(subset(bg_7d,"rowVars(texpr(bg_7d)) >1",genomesubset=TRUE), feature="transcript",covariate="trt",adjustvars = c("rep"), getFC=TRUE, meas="FPKM")
# 
# results5d_filt_transcripts_trt_adjrep_FC_FPKM_sorted = arrange(results5d_filt_transcripts_trt_adjrep_FC_FPKM,qval)
# results7d_filt_transcripts_trt_adjrep_FC_FPKM_sorted = arrange(results7d_filt_transcripts_trt_adjrep_FC_FPKM,qval)
# 
# results5d_transcripts_trt_adjrep_FC_FPKM_sort = data.frame(geneIDs=ballgown::geneIDs(bg_5d), results5d_filt_transcripts_trt_adjrep_FC_FPKM)



DE.res_5d = results5d_transcripts_trt_adjrep_FC_FPKM_sorted_df
DE.res_7d = results7d_transcripts_trt_adjrep_FC_FPKM_sorted_df


# read in latest summary file and add in DE results for 5- and 7-day analysis
# texpr_table = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/DE_texpression_results_adj.day.rep_texpr_matrix_annotated_with_tclass_FC.csv")

texpr_table$DE.5d.geneIDs = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$geneIDs
texpr_table$DE.5d.feature = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$feature
texpr_table$DE.5d.id = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$id
texpr_table$DE.5d.fc = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$fc
texpr_table$DE.5d.pval = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$pval
texpr_table$DE.5d.qval = DE.res_5d[match(texpr_table$gene_id, DE.res_5d$geneIDs),]$qval

texpr_table$DE.7d.geneIDs = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$geneIDs
texpr_table$DE.7d.feature = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$feature
texpr_table$DE.7d.id = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$id
texpr_table$DE.7d.fc = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$fc
texpr_table$DE.7d.pval = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$pval
texpr_table$DE.7d.qval = DE.res_7d[match(texpr_table$gene_id, DE.res_7d$geneIDs),]$qval
# texpr_table$Cs.id = Cs_anno[match(texpr_table$gene_id, Cs_anno$Cs_gene),]$Cs_id
# texpr_table$At.id = Cs_anno[match(texpr_table$gene_id, Cs_anno$Cs_gene),]$At_id
# texpr_table$Cs.description = Cs_anno[match(texpr_table$gene_id, Cs_anno$Cs_gene),]$Description
# texpr_table$Cs_t_name = Cs_anno[match(texpr_table$t_name, Cs_anno$Cs_id),]$Cs_id
# texpr_table$tmap_to_Cs = tmap[match(texpr_table$t_name, tmap$qry_id),]$ref_id
# texpr_table$Cs_id_merge = ifelse(is.na(texpr_table$Cs_t_name),texpr_table$tmap_to_Cs,texpr_table$Cs_t_name)
# texpr_table$Cs.description = Cs_anno[match(texpr_table$Cs_id_merge, Cs_anno$Cs_id),]$Description
# texpr_table$t_class = tmap[match(texpr_table$t_name, tmap$qry_id),]$class_code
# texpr_table$Cs.id.new = Cs_anno[match(texpr_table$Cs_id_merge, Cs_anno$Cs_id),]$Cs_id
# texpr_table$At.id.new = Cs_anno[match(texpr_table$Cs_id_merge, Cs_anno$Cs_id),]$At_id
# texpr_table$Description.new = Cs_anno[match(texpr_table$Cs_id_merge, Cs_anno$Cs_id),]$Description

write.csv(texpr_table, "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/DE_texpression_results_adj.day.rep_texpr_matrix_annotated_with_tclass_FC_with_5d.7d.DE.results.csv")

# create a column to mark genes with a circRNA
texpr_table$circRNA_identified = NA
texpr_table[texpr_table$Cs_id_merge %in% clear_res$isoformName,]$circRNA_identified <- "X"

# get rid of columns I don't want to keep in the final output
texpr_table$gene_name = NULL
texpr_table$cov.S14_abund = NULL
texpr_table$cov.S15_abund = NULL
texpr_table$cov.S16_abund = NULL
texpr_table$cov.S17_abund = NULL
texpr_table$cov.S18_abund = NULL
texpr_table$cov.S19_abund = NULL
texpr_table$cov.S20_abund = NULL
texpr_table$cov.S21_abund = NULL
texpr_table$cov.S22_abund = NULL
texpr_table$cov.S23_abund = NULL
texpr_table$cov.S24_abund = NULL
texpr_table$cov.S25_abund = NULL
texpr_table$cov.S26_abund = NULL
texpr_table$cov.S27_abund = NULL
texpr_table$cov.S28_abund = NULL
texpr_table$cov.S29_abund = NULL

texpr_table$Cs.id = NULL
texpr_table$At.id = NULL
texpr_table$Cs.description = NULL
texpr_table$Cs_t_name = NULL

write.csv(texpr_table, "Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/DE_texpression_results_adj.day.rep_texpr_matrix_annotated_with_tclass_FC_with_5d.7d.DE.results_final.csv")



# common_alltrts = c("MSTRG.13110",	"MSTRG.13110",	"MSTRG.13110",	"MSTRG.13110",	"MSTRG.13110",	"MSTRG.13110",	"Csa13g018960",	"Csa13g018960",	"MSTRG.11957",	"MSTRG.11957",	"MSTRG.11957",	"MSTRG.11957",	"Csa12g084210",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"Csa18g019070",	"Csa18g019070",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.26960",	"MSTRG.35161",	"MSTRG.35161",	"MSTRG.35161",	"MSTRG.35161",	"MSTRG.35161",	"Csa20g052830",	"MSTRG.35161",	"MSTRG.35162",	"MSTRG.35163",	"MSTRG.35162",	"MSTRG.42678",	"MSTRG.42678",	"Csa05g018230",	"Csa05g018230",	"MSTRG.42678",	"MSTRG.42678",	"MSTRG.42678",	"MSTRG.42678",	"MSTRG.42678",	"MSTRG.15414",	"Csa14g015000",	"Csa14g015000",	"MSTRG.6861",	"Csa11g034060",	"Csa11g034060",	"MSTRG.6861",	"Csa03g011570",	"MSTRG.36359",	"MSTRG.36359",	"Csa03g011570",	"MSTRG.2101",	"MSTRG.2101",	"MSTRG.2101",	"Csa01g038070",	"MSTRG.2101",	"MSTRG.2101",	"MSTRG.2101",	"Csa15g064470",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"MSTRG.19720",	"Csa07g062760",	"Csa10g008100",	"Csa10g008100",	"MSTRG.2820",	"MSTRG.2820",	"MSTRG.2820",	"MSTRG.2820",	"MSTRG.2820",	"MSTRG.2820",	"Csa10g008100",	"MSTRG.2820",	"MSTRG.39353",	"MSTRG.39353",	"MSTRG.39353",	"Csa04g016640",	"MSTRG.39353")
genes_common_alltrts = c("MSTRG.13110",	"MSTRG.11957",	"MSTRG.26960",	
                         "MSTRG.35161",	"MSTRG.35162",	"MSTRG.35163",	
                         "MSTRG.42678",	"MSTRG.15414",	"MSTRG.6861",	
                         "MSTRG.36359",	"MSTRG.2101",	"MSTRG.19720",	
                         "MSTRG.50028",	"MSTRG.2820",	"MSTRG.39353")

genes_common_alltrts.reduced = genes_common_alltrts[-c(5,6,13)]

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/LipaseClass3_assembled_transcripts_bySample.pdf", width=25, height=15)
plotTranscripts("MSTRG.36174", bg, 
                samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                          'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                          'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                          'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                meas='FPKM', colorby='transcript')
dev.off()


pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/DNALigaseI_assembled_transcripts_bySample.pdf", width=25, height=15)
plotTranscripts("MSTRG.36359", bg, 
                samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                          'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                          'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                          'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                meas='FPKM', colorby='transcript')
dev.off()

for(i in 1:length(genes_common_alltrts)){
  plotTranscripts(genes_common_alltrts[i], bg, 
                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                  meas='FPKM', colorby='transcript')
}

dev.off()


pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/DE.circ.genes_common_in_all_samples_width25_latent.clusters.pdf", width=25, height=15)
for(i in 1:length(genes_common_alltrts.reduced)){
plotLatentTranscripts(genes_common_alltrts.reduced[i], bg, 
                      k=NULL, method='kmeans', choosek = 'var90',
                      returncluster=FALSE)
}
dev.off()




## now overlay the circRNA start & end coordinates to the plots
# circRNA_coords = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/common_to_all.csv")

gene_key = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/transcript_gene_key.csv")
clear_res = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/CLEAR/data_summaries/CLEAR_output_allsamples.merged.csv")
clear_res$circID = paste(clear_res$chrom,":",clear_res$start, "-",clear_res$end, sep='')
clear_res$Sample = paste(clear_res$Sample,"_abund",sep='')


hk_genes = c("Csa13g018960",	"Csa12g084210",	"Csa11g025660",	"Csa18g019070",	"Csa20g052830",	"Csa05g018230",	"Csa14g015000",	"Csa11g034060",	"Csa03g011570",	"Csa01g038070",	"Csa15g064470",	"Csa07g062760",	"Csa10g008100",	"Csa04g016640")
hk_df = clear_res[clear_res$geneName %in% hk_genes,]
hk_df = hk_df %>% group_by(Sample,geneName) %>% arrange(start, .by_group=T)

Cs_id = unique(gene_key[gene_key$g_id %in% genes_common_alltrts,]$gene)



# test to find the best way to reformat the CLEAR results in order to use it for plotting with transcripts
# clear_res$circID = paste(clear_res$chrom,":",clear_res$start, "-",clear_res$end, sep='')
# clear_res$Sample = paste(clear_res$Sample,"_abund",sep='')
# circ.genename = unique(gene_key[gene_key$g_id==genes_common_alltrts[1],]$gene)
circ.subset = clear_res[clear_res$geneName==circ.genename,]
circ.subset.wide = circ.subset %>% group_by(circID) %>% select(circID,chrom,start,end,Sample,readNumber) %>% spread(Sample,readNumber,fill = 0)
colnames(circ.subset.wide) = c("circID","chrom","start","end","S14_abund",
                               "S15_abund","S16_abund","S18_abund",
                               "S19_abund","S20_abund","S21_abund",
                               "S22_abund","S23_abund","S24_abund",
                               "S25_abund","S26_abund","S27_abund",
                               "S28_abund","S29_abund")
no.testsamples=16
circ.subset.cts.matrix = as.matrix(circ.subset.wide[,5:no.testsamples])


paste(clear_res$Sample,"_abund",sep='')




hk_df1 = clear_res[clear_res$geneName==hk_genes[1],] 

plotTranscripts(genes_common_alltrts[1], bg, 
                samples=c('S14_abund', 'S15_abund'), 
                meas='FPKM', colorby='transcript')

overlayCircs(genes_common_alltrts[1], bg, samples=c('S14_abund'))


overlayCircsTest(gene=genes_common_alltrts[1], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
# overlayCircsTest(genes_common_alltrts[1], bg, samples=c('S14_abund', 'S15_abund'), gene_key, clear_res)





pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_MSTRG.11957_with.circ.overlay_w35h15.pdf", width =35 ,height = 15)
  overlayCircsTest(gene=genes_common_alltrts.reduced[2], bgobj=bg,
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'),
                   main=NULL, key=new_key, clear=clear_res)
dev.off()


pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[3],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[3], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()


pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[4],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[4], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()


pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[6],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[6], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()

# The corresponding transcript for the "MSTRG.6861" isoform does not match anything in the CLEAR output.
# "Csa11g034050" is not found in the CLEAR output
#
# Maybe there was a typo in the MSTRG id because it definitely matches up with that gene 
# I think it's supposed to be the TOPLESS2 gene: "Csa11g034060" (off by one number)

# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[7],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[7], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()



pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[8],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[8], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()


pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[9],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[9], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()

pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[10],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[10], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()

pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[11],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[11], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()

pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
               genes_common_alltrts.reduced[12],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
overlayCircsTest(gene=genes_common_alltrts.reduced[12], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=gene_key, clear=clear_res)
dev.off()

# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[13],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[13], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()


# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[14],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[14], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()

# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[15],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[15], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()
# 
# 
# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[16],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[16], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()


# pdf(file=paste("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/CommonGenes_",
#                genes_common_alltrts.reduced[17],"_with.circ.overlay.pdf",sep=''),width =18 ,height = 20)
# overlayCircsTest(gene=genes_common_alltrts.reduced[17], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=gene_key, clear=clear_res)
# dev.off()







# samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#           'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#           'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund')


## SUBSET & REFORMAT CIRCRNA DATA FROM CLEAR QUANT
# gene_keydf = key
# clear = clear_output
circ.gene0 = unique(gene_key[gene_key$g_id==genes_common_alltrts[1],]$gene) # match the stringtie g_id to the geneName
circ.sub0 = clear_res[clear_res$geneName==circ.gene0,]  # subset clear results by geneName

# count the number of unique circRNAs across all samples
# uniq.circs = length(unique(paste(circ.sub$chrom,":",circ.sub$start, "-",circ.sub$end, sep='')))
uniq.circs0 = length(unique(circ.sub0$circID))

# reformat the data and collapse rows into one per unique circRNA 
circ.wide0 = circ.sub0 %>% group_by(circID) %>% select(circID,chrom,start,end,Sample,readNumber) %>% spread(Sample,readNumber,fill = 0)
circ.cts.matrix0 = as.matrix(circ.wide0[,5:length(unique(circ.sub0$Sample))])
## do the same thing for the circ reads to establish colorscale for circ tracks
maxcol.circ0 = quantile(circ.cts.matrix0, 0.99)
colscale.circ0 = seq(0, maxcol.circ0, length.out=200) 

circIndex0 = which(colnames(circ.wide0)=='S14_abund')
circ_loop0 = unique(circ.wide0$circID)
cind0 = which(circ_loop0==circ_loop0[1]) 
circ.color0 = closestColor(circ.wide0[,circIndex0][cind0,],colscale.circ0)
circ.color1 = closestColor(as.numeric(circ.wide0[,circIndex0][cind0,]),colscale.circ0)
# 
# choices = rev(colorRampPalette(c("#009E73","#56B4E9","#FFFFFC"))(200))
# diffs = abs(x-colscale)
# return(choices[which.min(diffs)])

circ.gene2 = unique(gene_key[gene_key$g_id==MSTRGids.set4[5],]$gene) # match the stringtie g_id to the geneName
circ.sub3 = clear_df[clear_df$geneName==circ.gene2,]  # subset clear results by geneName
uniq.circs3 = length(unique(circ.sub3$circID)) # count the number of unique circRNAs across all samples
circ.wide3 = circ.sub3 %>% group_by(circID) %>% select(circID,chrom,start,end,Sample,readNumber) %>% spread(Sample,readNumber,fill = 0)
circ.cts.matrix3 = as.matrix(circ.wide3[,5:length(unique(circ.sub3$Sample))])

# write.csv(clear_res_complete, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/clear_res_complete.csv")
clear_df = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/clear_res_complete_df.csv")
clear_res = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/CLEAR/data_summaries/CLEAR_output_allsamples.merged.csv")

clear_df[is.na(clear_df$readNumber),]$readNumber = 0
clear_df$X.1=NULL
head(clear_df)
clear_df[is.na(clear_df$geneName),]

clear_df$new.geneName = NA
clear_df$new.geneName = clear_res[match(clear_df$circID, clear_res$circID),]$geneName
clear_df$geneName=clear_df$new.geneName
clear_df$new.geneName=NULL

Chr3:3385824-3387133

if(nrows(circ.sub2) == 1){
  maxcol.circ = quantile(circ.sub2$readNumber, 0.99)
  colscale.circ = seq(0, maxcol.circ, length.out=200) 
}else{
  # reformat the data and collapse rows into one per unique circRNA 
  circ.wide2 = circ.sub2 %>% group_by(circID) %>% select(circID,chrom,start,end,Sample,readNumber) %>% spread(Sample,readNumber,fill = 0)
  circ.cts.matrix2 = as.matrix(circ.wide2[,5:length(unique(circ.sub2$Sample))])
  ## do the same thing for the circ reads to establish colorscale for circ tracks
  maxcol.circ0 = quantile(circ.cts.matrix0, 0.99)
  colscale.circ0 = seq(0, maxcol.circ0, length.out=200) 
}

clear_res %>% complete(circID,Sample,fill=0)
  
  # group_by(circID) %>% select(circID,chrom,start,end,Sample,readNumber) %>% spread(Sample,readNumber,fill = 0)


# 
# > length(unique(clear_res$circID))
# [1] 3272




indexes(bg)$t2g$t_id[indexes(bg)$t2g$g_id=="MSTRG.11957"]




circ.geneX = unique(gene_key[gene_key$g_id=="MSTRG.11957",]$gene) # match the stringtie g_id to the geneName
circ.subX = clear[clear$geneName %in% circ.geneX,]

# circ.borderCol = "limegreen"
# 
# csub0 = circ.wide0[circ.wide0$circID==cind0,] ## creates a new df that contains the coordinates for that transcript id
# csub0 = csub0[order(csub0$start),] ## reorders the dataframe based on the transcript with the lowest start position
# 
# 
# cind = which(circ_loop==c) 
# circ.color = closestColor.circ(circ.wide[,circIndex][cind,],colscale.circ)
# csub = circ.wide0[circ.wide0$circID==cind0,] 
# 
## creates a new df that contains the coordinates for that transcript id
# csub = csub[order(csub$start),] ## reorders the dataframe based on the transcript with the lowest start position
#     
# for(cind in 1:dim(csub)[1]){
circ.borderCol = "gray20"
polygon(x=c(csub$start[cind], csub$start[cind], ## sets the top & bottom points / plotting coordinates for the leftmost corners of the rectangle/exon within that transcript
            csub$end[cind], csub$end[cind]), ## sets the top & bottom points / plotting coordinates for the rightmost corners of the rectangle/exon within that transcript
        y=c(cind-0.4,cind+0.4,cind+0.4,cind-0.4), ## sets the position within the plotting window for that transcript and the width of the rectangle
        col=circ.color, border=circ.borderCol)



for(c in circ_loop0){ ## plots each circRNA one-by-one 
  cind = which(circ_loop0==c) 
  # circ.color = closestColor.circ(circ.wide[,circIndex][cind,],colscale.circ)
  csub = circ.wide0[circ.wide0$circID==c,] ## creates a new df that contains the coordinates for that transcript id
  print(csub$start[c])
  
  #     csub = csub[order(csub$start),] ## reorders the dataframe based on the transcript with the lowest start position
  #     
  # for(cind in 1:dim(csub)[1]){
  # circ.borderCol = "gray20"
  # polygon(x=c(csub$start[cind], csub$start[cind], ## sets the top & bottom points / plotting coordinates for the leftmost corners of the rectangle/exon within that transcript
  #             csub$end[cind], csub$end[cind]), ## sets the top & bottom points / plotting coordinates for the rightmost corners of the rectangle/exon within that transcript
  #         y=c(cind-0.4,cind+0.4,cind+0.4,cind-0.4), ## sets the position within the plotting window for that transcript and the width of the rectangle
  #         col=circ.color, border=circ.borderCol) ## fill the rectangle with the color corresponding to the FPKM value and outline the rectangle in black
  # # }
}
################################################################################
# TO ADD CIRCRNA TO TRANSCRIPT PLOTS:
#
# First, check how many unique transcripts exist for that gene_id so that you get the circRNA positioned correctly
#
# in the background plotting code, 'ma' is first data structure that is created, which is:
# ma = IRanges::as.data.frame(structure(gown)$trans)\
#
# then for each plot (or gene_id), the ma IRanges dataframe is subset into a smaller dataframe called gtrans, which is: 
# gtrans = ma[ma$group_name %in% thetranscripts,]
#
# that dataframe is then subset into a smaller dataframe called 'transcript_loop' for each unique gene_id (or group_name) containing all transcripts with the same gene_id:
# transcript_loop = unique(gtrans$group_name)
#
#
polygon(x=c(6899458, 6899634,6900459,6900312), 
        y=c(1-0.2,1+0.2,1+0.2,1-0.2), 
        col="limegreen", border="grey60")




# Csa13g018960 = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/Csa13g018960.csv")
# test_df = Csa13g018960[Csa13g018960$Sample=="S14",]
# 
# test_gr = makeGRangesFromDataFrame(test_df,
#                          keep.extra.columns=T,
#                          ignore.strand=FALSE,
#                          seqinfo=NULL,
#                          seqnames.field="chrom",
#                          start.field="start",
#                          end.field="end",
#                          strand.field="strand",
#                          starts.in.df.are.0based=FALSE,
#                          na.rm=FALSE)
#



new_key = read.delim("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/GFFcomp.stringtie_merged.gtf.tmap")
t2g_key = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/transcript_gene_key.csv")
ATlip = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/lipidBS_genes_At.csv")
Cslip = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/CamRegBase_lipid_gene_search.csv")

unique(ATlip$identifier)
ATlip$At.gene = gsub("At","AT",ATlip$identifier)
ATlip$At.gene = gsub("at","AT",ATlip$At.gene)
ATlip$At.gene = gsub("g","G",ATlip$At.gene)

Cs_annot =  Cs_anno %>% separate(At_id, c("At_gene","At_iso"))
Atlip2Cs = Cs_annot[Cs_annot$At_gene %in% ATlip$At.gene,]
Cslip_genes = Atlip2Cs$Cs_gene


Cslip_from_both_databases = c(Cslip$Accession,Cslip_genes)
Cslip_from_both_databases_uniq = unique(Cslip_from_both_databases)
Cslip_from_both_databases_uniq_with_annotation = Cs_annot[Cs_annot$Cs_gene %in% Cslip_from_both_databases_uniq,]
write.csv(Cslip_from_both_databases_uniq_with_annotation, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/Heikes_Arabi_genes_and_CamRegBase_lipid_genes_search.csv")

# Cslip_from_both_databases_uniq_with_annotation_MSTRG = new_key[new_key$ref_gene_id %in% unique(Cslip_from_both_databases_uniq_with_annotation$Cs_gene),]$qry_gene_id
Cslip_from_both_databases_uniq_with_annotation$MSTRG.id = new_key[match(Cslip_from_both_databases_uniq_with_annotation$Cs_gene,new_key$ref_gene_id),]$qry_gene_id
# write.csv(Cslip_from_both_databases_uniq_with_annotation, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/Heikes_Arabi_genes_and_CamRegBase_lipid_genes_search.csv")
Cslip_desc = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/Heikes_Arabi_genes_and_CamRegBase_lipid_genes_search.csv")
CIRI_genes = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/circular_RNA/CIRI_genes_with_circRNA.csv")
CIRI_res = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/CIRI/CIRI2_output/CIRI2_output_files_merged_annotated_with_lncRNAs.csv")

pred_miRNA = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/predicted_miRNA_targets_gene.list.csv")
known_miRNA = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/known_miRNA_targets_gene.list.csv")


CIRI_res %>% group_by(day,trt,merge.gene.description) %>% summarise(n = sum(X.junction_reads),
                                                                    n.uniq.Csgenes = length(unique(merge.Cs_gene.id)))
length(unique(CIRI_res$merge.gene.description))

Cs_circs_inATlip = clear_res[clear_res$geneName %in% Cslip_genes,]
Cs_circs_inCslip = clear_res[clear_res$geneName %in% Cslip$Accession,] #from CamRegBase
CIRI_circs_inCslip = CIRI_res[CIRI_res$gene_id %in% Cslip$Accession,]$gene_id
CLEAR_circs_with_miRNA.pred = clear_res[clear_res$geneName %in% pred_miRNA$TargetGene_unique,]
CLEAR_circs_with_miRNA.knwn = clear_res[clear_res$geneName %in% known_miRNA$TargetGene_unique,]
CIRI_circs_with_miRNA.pred = CIRI_res[CIRI_res$gene_id %in% pred_miRNA$TargetGene_unique,]
CIRI_circs_with_miRNA.knwn = CIRI_res[CIRI_res$gene_id %in% known_miRNA$TargetGene_unique,]

  
  
  
CIRI_circs_inATlip = Atlip2Cs[Atlip2Cs$Cs_gene %in% CIRI_genes$gene_id,]
Aralip_inCIRI = Cslip[Cslip$Accession %in% CIRI_genes$gene_id,] #from CamRegBase
CIRI_circs_inCslip$At_gene = Cs_annot[match(CIRI_circs_inCslip$Accession,Cs_annot$Cs_gene),]$At_gene

write.csv(Cs_circs_inATlip, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/Cs_circs_inATlip.csv")
write.csv(Cs_circs_inCslip, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/Cs_circs_inCslip.csv")
write.csv(CIRI_circs_inCslip, file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/CIRI_circs_inAralip.csv")




MSTRGids_inATlip = new_key[new_key$ref_gene_id %in% Cs_circs_inATlip$geneName,]$qry_gene_id
MSTRGids_inCslip = new_key[new_key$ref_gene_id %in% Cs_circs_inCslip$geneName,]$qry_gene_id

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.from.AT.table.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids_inATlip)){
  overlayCircsTest(gene=MSTRGids_inATlip[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

MSTRGids = c(MSTRGids_inATlip,MSTRGids_inCslip)

MSTRGids = unique(MSTRGids) # 57 genes with circRNA

check = c("MSTRG.6","MSTRG.7","MSTRG.50","MSTRG.51")
MSTRGids %in% check


MSTRGids.set1 = MSTRGids[1:12]
MSTRGids.set2 = MSTRGids[13:25]
MSTRGids.set3 = MSTRGids[26:38]
MSTRGids.set4 = MSTRGids[39:48]
MSTRGids.set5 = MSTRGids[49:57]
MSTRGids.DElist = check

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set1_with.circ.overlay.plots_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set1)){
overlayCircsTest(gene=MSTRGids.set1[i], bgobj=bg, 
                 samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                           'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                          'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                          'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                 main=NULL, key=new_key, clear=clear_df)
}
dev.off()


pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set2_with.circ.overlay_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set2)){
  overlayCircsTest(gene=MSTRGids.set2[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set3_with.circ.overlay_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set3)){
  overlayCircsTest(gene=MSTRGids.set3[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set4_with.circ.overlay_test_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set4)){
  overlayCircsTest(gene=MSTRGids.set4[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set5_with.circ.overlay_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set5)){
  overlayCircsTest(gene=MSTRGids.set5[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/CommonGenes_DE.genes_with.circ.overlay_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.DElist)){
  overlayCircsTest(gene=MSTRGids.DElist[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()


pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/LipidGenes.set1_with.circ.overlay.plots_v2.pdf",width =20 ,height = 18)
for(i in 1:length(MSTRGids.set1)){
  overlayCircsTest(gene=MSTRGids.set1[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()

# Hk_genes_df = Cs_anno[Cs_anno$Cs_gene %in% hk_genes,]
# write.csv(Hk_genes_df,file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/Hk_genes_df.csv")
hk_genes_df = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/Hk_genes_df.csv")

hk_genes_mstrg = unique(new_key[new_key$ref_gene_id %in% hk_genes,]$qry_gene_id)
hk_genes

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/HouseKeepingGenes_with_circRNA.pdf",width =20 ,height = 18)
for(i in 1:length(hk_genes_mstrg)){
  overlayCircsTest(gene=hk_genes_mstrg[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()



trtsp_genes = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/treatment.specific.circs.csv")

cons.5L = trtsp_genes$cons.5L
cons.7L = trtsp_genes$cons.7L
cons.5D = trtsp_genes$cons.5D
cons.7D = trtsp_genes$cons.7D

cons.5L_mstrg = unique(new_key[new_key$ref_gene_id %in% cons.5L,]$qry_gene_id)
cons.7L_mstrg = unique(new_key[new_key$ref_gene_id %in% cons.7L,]$qry_gene_id)
cons.5D_mstrg = unique(new_key[new_key$ref_gene_id %in% cons.5D,]$qry_gene_id)
cons.7D_mstrg = unique(new_key[new_key$ref_gene_id %in% cons.7D,]$qry_gene_id)



pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/cons.5L_mstrg.pdf",width =20 ,height = 18)
for(i in 1:length(cons.5L_mstrg)){
  overlayCircsTest(gene=cons.5L_mstrg[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()



pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/cons.7L_mstrg.pdf",width =20 ,height = 18)
for(i in 1:length(cons.7L_mstrg)){
  overlayCircsTest(gene=cons.7L_mstrg[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()



pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/cons.5D_mstrg.pdf",width =20 ,height = 18)
for(i in 1:length(cons.5D_mstrg)){
  overlayCircsTest(gene=cons.5D_mstrg[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()



pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/linearRNA_with_circRNA_plots/cons.7D_mstrg.pdf",width =20 ,height = 18)
for(i in 1:length(cons.7D_mstrg)){
  overlayCircsTest(gene=cons.7D_mstrg[i], bgobj=bg, 
                   samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                             'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                             'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                             'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
                   main=NULL, key=new_key, clear=clear_df)
}
dev.off()




common_noncons = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/circRNA_to_overlay/common.noncons.csv")
common_noncons_mstrg = unique(new_key[new_key$ref_gene_id %in% common_noncons$all.noncons,]$qry_gene_id)



syntelogs = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/Cs_reference_genome/all_syntelogs.csv")



lin_res = read.csv("Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/stringtie_assembly/abundance/DE_texpression_results_adj.day.rep_texpr_matrix_annotated_with_tclass_FC.csv")

lin_res %>% group_by(t_class) %>% summarise(sum.5L = sum(FPKM.light.5),
                                            sum.7L = sum(FPKM.light.7),
                                            sum.5D = sum(FPKM.dark.5),
                                            sum.7D = sum(FPKM.dark.7),
                                            sum.dark = sum(c(FPKM.dark.5,FPKM.dark.7)),
                                            sum.light = sum(c(FPKM.light.5,FPKM.light.7)))

lin_res$FPKM.light.5 = lin_res$FPKM.S14_abund+lin_res$FPKM.S15_abund+lin_res$FPKM.S16_abund+lin_res$FPKM.S17_abund
lin_res$FPKM.light.7 = lin_res$FPKM.S22_abund+lin_res$FPKM.S23_abund+lin_res$FPKM.S24_abund+lin_res$FPKM.S25_abund
lin_res$FPKM.dark.5 = lin_res$FPKM.S18_abund+lin_res$FPKM.S19_abund+lin_res$FPKM.S20_abund+lin_res$FPKM.S21_abund
lin_res$FPKM.dark.7 = lin_res$FPKM.S26_abund+lin_res$FPKM.S27_abund+lin_res$FPKM.S28_abund+lin_res$FPKM.S29_abund


# clear_res_complete %>% select(geneName,readNumber)

# ma0 = IRanges::as.data.frame(structure(bg)$trans) ## creates IRanges data frame containing all transcripts
# thetranscripts0 = indexes(bg)$t2g$t_id[indexes(bg)$t2g$g_id=="MSTRG.37203"]
# g_id0 = texpr(bg, 'all')$gene_id ## pulls all the gene_id's from the texpr data
# smalldat0 = texpr(bg, 'FPKM')[which(g_id0 == "MSTRG.37203"),] ## pulls the texpr FPKM values for that gene of interest 
# t_id0 = texpr(bg, 'all')$t_id[which(g_id0 == "MSTRG.37203")] ## creates a vector of the t_id's for that gene 
# 
# if(numtx == 1){ ## if only one transcript exists for that gene_id... 
#   snames = names(smalldat) ## sets snames to NULL 
#   smalldat = matrix(smalldat, nrow=1) ## creates a matrix with all FPKM values collapsed into a single row
#   colnames(smalldat) = snames ## sets colnames for that new matrix
# }







# overlayCircsTest(gene=MSTRGids.set1[1], bgobj=bg, 
#                  samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
#                            'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
#                            'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
#                            'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'), 
#                  main=NULL, key=new_key, clear=clear_res)

# t2g_key[t2g_key$t_id==497,]
# new_key[new_key$qry_gene_id=='MSTRG.11957',]

# CROSS-REF BETWEEN GENE NAMES:
# 'qry_gene_id' = ballgown gene_id
# 'ref_gene_id' = CLEAR geneName



clusterTranscripts(gene='MSTRG.11957', gown=bg, k=2, method='kmeans')


lipid_genes = c("Csa05g058150","Csa16g038860","Csa03g023080","Csa07g046400","Csa11g033900")
lipid_genes_mstrg = unique(new_key[new_key$ref_gene_id %in% lipid_genes,]$qry_gene_id)
clusterTranscripts(gene=lipid_genes_mstrg[1], gown=bg, k=2,method='kmeans')
plotLatentTranscripts(gene=lipid_genes_mstrg[1], gown=bg, k=2, method='kmeans', returncluster=FALSE)

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/lipidBS_Genes.pdf",width =20 ,height = 18)
plotTranscripts(gene=lipid_genes_mstrg[1], gown=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund','S18_abund', 'S19_abund', 'S20_abund', 'S21_abund','S22_abund', 'S23_abund', 'S24_abund', 'S25_abund','S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'))
plotTranscripts(gene=lipid_genes_mstrg[2], gown=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund','S18_abund', 'S19_abund', 'S20_abund', 'S21_abund','S22_abund', 'S23_abund', 'S24_abund', 'S25_abund','S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'))
plotTranscripts(gene=lipid_genes_mstrg[3], gown=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund','S18_abund', 'S19_abund', 'S20_abund', 'S21_abund','S22_abund', 'S23_abund', 'S24_abund', 'S25_abund','S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'))
plotTranscripts(gene=lipid_genes_mstrg[4], gown=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund','S18_abund', 'S19_abund', 'S20_abund', 'S21_abund','S22_abund', 'S23_abund', 'S24_abund', 'S25_abund','S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'))
plotTranscripts(gene=lipid_genes_mstrg[5], gown=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund','S18_abund', 'S19_abund', 'S20_abund', 'S21_abund','S22_abund', 'S23_abund', 'S24_abund', 'S25_abund','S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'))
dev.off()

pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/lipidBS_Genes_clustered_kmeans2.pdf",width =20 ,height = 18)
plotLatentTranscripts(gene=lipid_genes_mstrg[1], gown=bg, k=2, method='kmeans', returncluster=FALSE)
plotLatentTranscripts(gene=lipid_genes_mstrg[2], gown=bg, k=2, method='kmeans', returncluster=FALSE)
# plotLatentTranscripts(gene=lipid_genes_mstrg[3], gown=bg, k=2, method='kmeans', returncluster=FALSE)
# plotLatentTranscripts(gene=lipid_genes_mstrg[4], gown=bg, k=2, method='kmeans', returncluster=FALSE)
# plotLatentTranscripts(gene=lipid_genes_mstrg[5], gown=bg, k=2, method='kmeans', returncluster=FALSE)
dev.off()


pdf(file="Q:/Shared drives/Sederoff Lab/Camelina circRNA Sequencing/linear_RNA/results/plots/DNALigaseI_with_circ_overlay",width =25 ,height = 18)
overlayCircsTest(gene="MSTRG.36359", bgobj=bg, samples=c('S14_abund', 'S15_abund', 'S16_abund', 'S17_abund',
                                                         'S22_abund', 'S23_abund', 'S24_abund', 'S25_abund',
                                                         'S18_abund', 'S19_abund', 'S20_abund', 'S21_abund',
                                                         'S26_abund', 'S27_abund', 'S28_abund', 'S29_abund'),main=NULL, key=new_key, clear=clear_df)
dev.off()



fpkm = texpr(bg,meas="FPKM")
# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(as.factor(paste(pData(bg)$trt,pData(bg)$day,sep='.')))+1,las=2,ylab='log2(FPKM+1)')
# ,col=c('#FFF2CC','#BDD7EE','#C5E0B4','#4472C4'))
boxplot(fpkm,col=c('#FFF2CC','#FFF2CC','#FFF2CC','#FFF2CC',
                   '#BDD7EE','#BDD7EE','#BDD7EE','#BDD7EE',
                   '#C5E0B4','#C5E0B4','#C5E0B4','#C5E0B4',
                   '#4472C4','#4472C4','#4472C4','#4472C4'),
        las=2,ylab='log2(FPKM+1)',
        xlab=rep(c('1','2','3','4'),4))
