##########################################################################################################################
#  Overview
##########################################################################################################################
# We load quantified total proteome and phosphoproteome data generated from MaxQuant data
# Data sets are annotated using current EntrezGene and GeneSymbol identifiers
# Finally, for both total proteome and phosphoproteome data we test if the averaged log2FC of the two replicates is significantly different from 0 using the R package LIMMA
#
# Note 1: Replicates 2 and 3 from raw data have been renamed to replicate 1 and 2 in the manuscript, respectively. There is no replicate 1 with data. 
# Note 2: the final phosphoproteome data has more than one line per phosphosite due to different values in column Multiplicity (i.e. number of phosphosites). For the differential epxression analyse, only the single row containing values for each phosphosite was used and then merged back to all sites, replicating data for different Multiplicity values. Therefore, all subsequent visualiations should also average results by site_id_unique.

library(data.table)
library(limma)

all_prot_diff_expr_tables = list()

##########################################################################################################################
#  PROTEOM 
##########################################################################################################################

proteom = fread("../../Data/Raw/Proteome/PZ Phos-Prot rep 2 and 3_whole.txt", sep="\t")
proteom$species = rep("C.trachomatis L2", nrow(proteom))
proteom$species[grep("Homo sapiens", proteom$"Fasta headers")] = "Homo sapiens"
nuc_val_cnt = apply(proteom[,c("Nuc_rep2",  "Nuc_rep3"),with=F], 1, function(x) sum(!is.na(x)))
total_val_cnt = apply(proteom[,c("Total_rep2",  "Total_rep3"),with=F], 1, function(x) sum(!is.na(x)))
proteom$Nuc_avg = rowSums(proteom[,c("Nuc_rep2",  "Nuc_rep3"),with=F],na.rm=T )
proteom$Nuc_avg[nuc_val_cnt==0] <- NA
proteom$Total_avg = rowSums(proteom[,c("Total_rep2",  "Total_rep3"),with=F],na.rm=T )
proteom$Total_avg[total_val_cnt==0] <- NA
#plot(proteom$Nuc_avg, proteom$Total_avg, col=as.numeric(as.factor(proteom$species)))
proteom$valid_cnt = nuc_val_cnt + total_val_cnt
proteom$uid = 1:nrow(proteom)

table(proteom$valid_cnt, proteom$species)
# why this distribution ??

setkey(phosphoproteom, "Proteins", "Positions within proteins")

library(org.Hs.eg.db)
# > org.Hs.eg.db
# OrgDb object:
# | DBSCHEMAVERSION: 2.1
# | Db type: OrgDb
# | Supporting package: AnnotationDbi
# | DBSCHEMA: HUMAN_DB
# | ORGANISM: Homo sapiens
# | SPECIES: Human
# | EGSOURCEDATE: 2015-Sep27
# | EGSOURCENAME: Entrez Gene
# | EGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
# | CENTRALID: EG
# | TAXID: 9606
# | GOSOURCENAME: Gene Ontology
# | GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
# | GOSOURCEDATE: 20150919
# | GOEGSOURCEDATE: 2015-Sep27
# | GOEGSOURCENAME: Entrez Gene
# | GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
# | KEGGSOURCENAME: KEGG GENOME
# | KEGGSOURCEURL: ftp://ftp.genome.jp/pub/kegg/genomes
# | KEGGSOURCEDATE: 2011-Mar15
# | GPSOURCENAME: UCSC Genome Bioinformatics (Homo sapiens)
# | GPSOURCEURL: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19
# | GPSOURCEDATE: 2010-Mar22
# | ENSOURCEDATE: 2015-Jul16
# | ENSOURCENAME: Ensembl
# | ENSOURCEURL: ftp://ftp.ensembl.org/pub/current_fasta
# | UPSOURCENAME: Uniprot
# | UPSOURCEURL: http://www.UniProt.org/
# | UPSOURCEDATE: Fri Oct  9 19:54:03 2015

all_Hs_entrezIDs = keys(org.Hs.eg.db)

sum(proteom$"Gene name"=="") # there are 566 proteins without gene name
table(proteom$species, proteom$"Gene name"=="") # most of those proteins without symbol are C.trachomatis proteins


# first try to get current Entrez IDs from database
prot_split = strsplit(proteom$'Protein IDs',";")
prot_split_fixed = lapply(prot_split, function(s) unique(sapply(strsplit(s,"-"), function(x) x[1]))) 
prot_cnts = sapply(prot_split_fixed, length)
expanded_prot_ids = data.frame(uid=rep(proteom$uid, times=prot_cnts), prot_id = unlist(prot_split_fixed), stringsAsFactors = F )

entrez_ids = select(org.Hs.eg.db, keys=expanded_prot_ids$prot_id, keytype = "UNIPROT", columns=c("ENTREZID"))
expanded_prot_ids_annotated = merge(expanded_prot_ids, entrez_ids, by.x="prot_id", by.y="UNIPROT", all.x=T, sort=F)
ll2 = tapply(expanded_prot_ids_annotated$ENTREZID, expanded_prot_ids_annotated$uid, function(x) length(unique(x[!is.na(x)])))
# there is quite a high number (~700) of proteins that show not match, but on closer inspection use unreviewed entries of known proteins/genes
# let's try harder - first idea would be to take the gene name from the fasta headers, but that column has truncated cells
# Just do a q&d hack and export all protein accessions, upload to the uniprot ID mapping webtool
# -> This gives:     8,664 out of 24,217 identifiers from UniProtKB AC/ID were successfully mapped to 7,553 Entrez Gene (GeneID) IDs.
# 7553 EntrezIds compares quite favorably to the 7464 (sum(ll2)) Ids we could obtain via org.Hs.eg.db
# So lets just go with what we got from org.Hs.eg.db

fasta_header_split_1 = strsplit(proteom$"Fasta headers", ";>")
fasta_header_split_2 = lapply(fasta_header_split_1, function(s) sapply(strsplit(s,"\\|"), function(x) x[3]) )
fasta_header_split_3 = lapply(fasta_header_split_2, function(s) unique(gsub("^GN=","",unlist(lapply(strsplit(s," "), function(x) x[grep("^GN=", x, perl=T)])))) )
ll3 = sapply(fasta_header_split_3, length)
table(ll3)
# we can fill in gene symbols for most of the proteins that had no entrez gene id. However, there are also 340 proteins with >1 symbol, approximately
# the same order of magnitude as the 450 from org.Hs.eg.db. Since Fasta headers also have been truncated, this lower number might just mean missed symbols
proteom$symbol_from_fasta = sapply(fasta_header_split_3, function(x) paste(x,collapse=","))
proteom$gene_symbol_curated = ifelse(proteom$"Gene name"=="", proteom$symbol_from_fasta, proteom$"Gene name")

# aggregate entrezID obtained via uniprot accessions from org.Hs.eg.db and add to proteom table
agg_entrezID = unlist(tapply(expanded_prot_ids_annotated$ENTREZID, expanded_prot_ids_annotated$uid, simplify=F ,
                      function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
proteom$entrezID_from_db = NA
proteom$entrezID_from_db[as.integer(names(agg_entrezID))] = as.character(as.vector(agg_entrezID))

# now use the provided gene symbols to obtain the missing entrez IDs via org.Hs.eg.db
gs_split = strsplit(proteom$gene_symbol_curated,";")
gs_cnts = sapply(gs_split, length)
expanded_gs = data.frame(uid=rep(proteom$uid, times=gs_cnts), gs = unlist(gs_split), stringsAsFactors = F )

entrez_ids_from_gs = select(org.Hs.eg.db, keys=expanded_gs$gs, keytype = "SYMBOL", columns=c("ENTREZID"))
expanded_gs_ids_annotated = merge(expanded_gs, unique(entrez_ids_from_gs), by.x="gs", by.y="SYMBOL", all.x=T, sort=F)

agg_entrezID = unlist(tapply(expanded_gs_ids_annotated$ENTREZID, expanded_gs_ids_annotated$uid, simplify=F ,
                             function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
proteom$entrezID_from_gs = NA
proteom$entrezID_from_gs[as.integer(names(agg_entrezID))] = as.character(as.vector(agg_entrezID))

# generate the final entrez IDs and count the number of proteins belonging to each of them as well as uniqueness of protein -> gene mapping
proteom$entrezID_curated = ifelse(is.na(proteom$entrezID_from_db), proteom$entrezID_from_gs, proteom$entrezID_from_db)
proteom$multiple_gene_ids = F
proteom$multiple_gene_ids[grep(",",proteom$entrezID_curated)] = T
proteom$gene_id_cnt = 0
rr = tapply(proteom$entrezID_curated, proteom$entrezID_curated, length)
proteom$gene_id_replicate_cnt = unlist(Map(function(x) nrow(proteom[entrezID_curated==x]), proteom$entrezID_curated))

# finally, get the symbol from the final entrez ID to make symbols and entrezID consistent
# now use the provided gene symbols to obtain the missing entrez IDs via org.Hs.eg.db
eID_split = strsplit(proteom$entrezID_curated,",")
eID_cnts = sapply(eID_split, length)
expanded_eID = data.frame(uid=rep(proteom$uid, times=eID_cnts), eID = unlist(eID_split), stringsAsFactors = F )

gs_from_eID = select(org.Hs.eg.db, keys=expanded_eID$eID, keytype = "ENTREZID", columns=c("SYMBOL"))
expanded_eID_annotated = merge(expanded_eID, unique(gs_from_eID), by.x="eID", by.y="ENTREZID", all.x=T, sort=F)

agg_gs = unlist(tapply(expanded_eID_annotated$SYMBOL, expanded_eID_annotated$uid, simplify=F ,
                             function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
#proteom$gene_symbol_final = NA
#proteom$gene_symbol_final[as.integer(names(agg_gs))] = as.character(as.vector(agg_gs))
setkey(proteom, "uid")
proteom[ i = as.integer(names(agg_gs)), gene_symbol_final := agg_gs]

# clean up
proteom[,symbol_from_fasta:=NULL]
proteom[,gene_symbol_curated:=NULL]
proteom[,"Gene names":=NULL]
proteom[,entrezID_from_db:=NULL]
proteom[,entrezID_from_gs:=NULL]
setnames(proteom, "Gene name","Gene_name_orig")
proteom_human = proteom[species=="Homo sapiens"]
proteom_ctr = proteom[species == "C.trachomatis L2"]

#### DPE

protein_expression_data = as.matrix(proteom_human[,c("Nuc_rep2","Nuc_rep3","Total_rep2","Total_rep3"), with=F])
rownames(protein_expression_data) = proteom_human$"entrezID_curated"
protein_expression_data = protein_expression_data[!is.na(rownames(protein_expression_data)),]

### Nuclear fraction
edata = avereps(protein_expression_data[,c(1,2)])
edata = edata[apply(edata,1,function(x) sum(is.na(x)))<2,]
design = c(1,1)
fit <- lmFit(edata,design)
fit <- eBayes(fit)
res = topTable(fit, adjust="BH", number=nrow(edata))

all_prot_diff_expr_tables[["proteom_nuclear"]] = res

colnames(res) = paste("proteom_nuclear", colnames(res), sep="_")
res$ID = rownames(res)

proteom_human = merge(proteom_human, res[,c(1,4,5,7)], by.x="entrezID_curated", by.y="ID", all.x=T, sort=F)

### Total fraction

edata = avereps(protein_expression_data[,c(3,4)])
edata = edata[apply(edata,1,function(x) sum(is.na(x)))<2,]
design = c(1,1)
fit <- lmFit(edata,design)
fit <- eBayes(fit)
res = topTable(fit, adjust="BH", number=nrow(edata))

all_prot_diff_expr_tables[["proteom_total"]] = res

colnames(res) = paste("proteom_total", colnames(res), sep="_")
res$ID = rownames(res)

proteom_human = merge(proteom_human, res[,c(1,4,5,7)], by.x="entrezID_curated", by.y="ID", all.x=T, sort=F)

##########################################################################################################################
# PHOSPHOPROTEOM
##########################################################################################################################
# This one does not contain the charge column, so unique identification of rows is impossible
# phosphoproteom = fread("./PZ Phos-Prot rep 2 and 3_phospho.txt", sep="\t")
# phosphoproteom$row_index = 1:nrow(phosphoproteom)
# We will load a previous annotated table with the same rows and more columns instead
load("../../Data/Processed/Proteome/Phosphoproteome_infection_annotated_with_phopsphosite_plus.Rdata")
phospho_proteom = data_merged
rm(data_merged)

phospho_proteom$species = rep("C.trachomatis L2", nrow(phospho_proteom))
phospho_proteom$species[grep("Homo sapiens", phospho_proteom$"Fasta headers")] = "Homo sapiens"

#setkey(phospho_proteom, "Proteins", "Positions within proteins")

sum(phospho_proteom$"Gene name"=="") # there are 288 proteins without gene name
table(phospho_proteom$species, phospho_proteom$"Gene name"=="") # most of those proteins without symbol are C.trachomatis proteins

# first try to get current Entrez IDs from database
prot_split = strsplit(phospho_proteom$'Proteins',";")
prot_split_fixed = lapply(prot_split, function(s) unique(sapply(strsplit(s,"-"), function(x) x[1]))) 
prot_cnts = sapply(prot_split_fixed, length)
expanded_prot_ids = data.frame(uid=rep(1:nrow(phospho_proteom), times=prot_cnts), prot_id = unlist(prot_split_fixed), stringsAsFactors = F )

entrez_ids = select(org.Hs.eg.db, keys=expanded_prot_ids$prot_id, keytype = "UNIPROT", columns=c("ENTREZID"))
expanded_prot_ids_annotated = merge(expanded_prot_ids, entrez_ids, by.x="prot_id", by.y="UNIPROT", all.x=T, sort=F)
ll2 = tapply(expanded_prot_ids_annotated$ENTREZID, expanded_prot_ids_annotated$uid, function(x) length(unique(x[!is.na(x)])))
# about 700 proteins where we couldn't get a entrez ID via Uniprot accessions

fasta_header_split_1 = strsplit(phospho_proteom$"Fasta headers", ";>")
fasta_header_split_2 = lapply(fasta_header_split_1, function(s) sapply(strsplit(s,"\\|"), function(x) x[3]) )
fasta_header_split_3 = lapply(fasta_header_split_2, function(s) unique(gsub("^GN=","",unlist(lapply(strsplit(s," "), function(x) x[grep("^GN=", x, perl=T)])))) )
ll3 = sapply(fasta_header_split_3, length)
table(ll3)
# we can fill in gene symbols for most of the proteins that had no entrez gene id. However, there are also 340 proteins with >1 symbol, approximately
# the same order of magnitude as the 450 from org.Hs.eg.db. Since Fasta headers also have been truncated, this lower number might just mean missed symbols
phospho_proteom$symbol_from_fasta = sapply(fasta_header_split_3, function(x) paste(x,collapse=","))
phospho_proteom$gene_symbol_curated = ifelse(phospho_proteom$"Gene name"=="", phospho_proteom$symbol_from_fasta, phospho_proteom$"Gene name")
# this brings us down to 21 proteins without gene symbol

# aggregate entrezID obtained via uniprot accessions from org.Hs.eg.db and add to phospho_proteom table
agg_entrezID = unlist(tapply(expanded_prot_ids_annotated$ENTREZID, expanded_prot_ids_annotated$uid, simplify=F ,
                             function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
phospho_proteom$entrezID_from_db = NA
phospho_proteom$entrezID_from_db[as.integer(names(agg_entrezID))] = as.character(as.vector(agg_entrezID))

# now use the provided gene symbols to obtain the missing entrez IDs via org.Hs.eg.db
gs_split = strsplit(phospho_proteom$gene_symbol_curated,";")
gs_cnts = sapply(gs_split, length)
expanded_gs = data.frame(uid=rep(1:nrow(phospho_proteom), times=gs_cnts), gs = unlist(gs_split), stringsAsFactors = F )

entrez_ids_from_gs = select(org.Hs.eg.db, keys=expanded_gs$gs, keytype = "SYMBOL", columns=c("ENTREZID"))
expanded_gs_ids_annotated = merge(expanded_gs, unique(entrez_ids_from_gs), by.x="gs", by.y="SYMBOL", all.x=T, sort=F)

agg_entrezID = unlist(tapply(expanded_gs_ids_annotated$ENTREZID, expanded_gs_ids_annotated$uid, simplify=F ,
                             function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
phospho_proteom$entrezID_from_gs = NA
phospho_proteom$entrezID_from_gs[as.integer(names(agg_entrezID))] = as.character(as.vector(agg_entrezID))

# generate the final entrez IDs and count the number of proteins belonging to each of them as well as uniqueness of protein -> gene mapping
phospho_proteom$entrezID_curated = ifelse(is.na(phospho_proteom$entrezID_from_db), phospho_proteom$entrezID_from_gs, phospho_proteom$entrezID_from_db)
phospho_proteom$multiple_gene_ids = F
phospho_proteom$multiple_gene_ids[grep(",",phospho_proteom$entrezID_curated)] = T
phospho_proteom$gene_id_cnt = 0
rr = table(phospho_proteom$entrezID_curated)/3
phospho_proteom$gene_id_replicate_cnt = rr[phospho_proteom$entrezID_curated]
#xx = unique(phospho_proteom[,c("entrezID_curated", "gene_id_replicate_cnt"), with=F])
#cumsum(table(xx$gene_id_replicate_cnt))
#barplot(table(xx$gene_id_replicate_cnt))

# finally, get the symbol from the final entrez ID to make symbols and entrezID consistent
# now use the provided gene symbols to obtain the missing entrez IDs via org.Hs.eg.db
eID_split = strsplit(phospho_proteom$entrezID_curated,",")
eID_cnts = sapply(eID_split, length)
expanded_eID = data.frame(uid=rep(1:nrow(phospho_proteom), times=eID_cnts), eID = unlist(eID_split), stringsAsFactors = F )

gs_from_eID = select(org.Hs.eg.db, keys=expanded_eID$eID, keytype = "ENTREZID", columns=c("SYMBOL"))
expanded_eID_annotated = merge(expanded_eID, unique(gs_from_eID), by.x="eID", by.y="ENTREZID", all.x=T, sort=F)

agg_gs = unlist(tapply(expanded_eID_annotated$SYMBOL, expanded_eID_annotated$uid, simplify=F ,
                       function(x) {v=unique(x[!is.na(x)]); if ( length(v)==0 ) {NA} else {paste(sort(v),collapse=",")} }  ) )
phospho_proteom$gene_symbol_final = NA
phospho_proteom$gene_symbol_final[as.integer(names(agg_gs))] = as.character(as.vector(agg_gs))

nuc_val_cnt = apply(phospho_proteom[,c("Nuc_rep2_phospho",  "Nuc_rep3_phospho"),with=F], 1, function(x) sum(!is.na(x)))
total_val_cnt = apply(phospho_proteom[,c("Total_rep2_phospho",  "Total_rep3_phospho"),with=F], 1, function(x) sum(!is.na(x)))
phospho_proteom$Nuc_avg = rowSums(phospho_proteom[,c("Nuc_rep2_phospho",  "Nuc_rep3_phospho"),with=F],na.rm=T )
phospho_proteom$Nuc_avg[nuc_val_cnt==0] <- NA
phospho_proteom$Total_avg = rowSums(phospho_proteom[,c("Total_rep2_phospho",  "Total_rep3_phospho"),with=F],na.rm=T )
phospho_proteom$Total_avg[total_val_cnt==0] <- NA
#plot(phospho_proteom$Nuc_avg, phospho_proteom$Total_avg, col=as.numeric(as.factor(phospho_proteom$species)))
phospho_proteom$valid_cnt = nuc_val_cnt + total_val_cnt
#barplot(table(phospho_proteom$valid_cnt/3))

phospho_proteom$site_id_unique = paste(phospho_proteom$entrezID_curated, phospho_proteom$"Sequence window", phospho_proteom$"Position", sep="_")

table(phospho_proteom$valid_cnt, phospho_proteom$species)
# why this distribution ??


# clean up
phospho_proteom[,symbol_from_fasta:=NULL]
phospho_proteom[,gene_symbol_curated:=NULL]
phospho_proteom[,"Gene names":=NULL]
phospho_proteom[,entrezID_from_db:=NULL]
phospho_proteom[,entrezID_from_gs:=NULL]
setnames(phospho_proteom, "Gene name","Gene_name_orig")
phospho_proteom_human = phospho_proteom[species=="Homo sapiens"]
phospho_proteom_ctr = phospho_proteom[species == "C.trachomatis L2"]



#### DPE
protein_expression_data = as.matrix(phospho_proteom_human[,c("Nuc_rep2_phospho","Nuc_rep3_phospho","Total_rep2_phospho","Total_rep3_phospho"), with=F])
rownames(protein_expression_data) = phospho_proteom_human$site_id_unique
protein_expression_data = protein_expression_data[!is.na(rownames(protein_expression_data)),]

### Nuclear fraction
edata = avereps(protein_expression_data[,c(1,2)])
edata = edata[apply(edata,1,function(x) sum(is.na(x)))<2,]
design = c(1,1)
fit <- lmFit(edata,design)
fit <- eBayes(fit)
res = topTable(fit, adjust="BH", number=nrow(edata))

all_prot_diff_expr_tables[["phospho_proteom_nuclear"]] = res

colnames(res) = paste("phospho_proteom_nuclear", colnames(res), sep="_")
res$ID = rownames(res)

phospho_proteom_human = merge(phospho_proteom_human, res[,c(1,4,5,7)], by.x="site_id_unique", by.y="ID", all.x=T, sort=F)

### Total fraction

edata = avereps(protein_expression_data[,c(3,4)])
design = c(1,1)
fit <- lmFit(edata,design)
fit <- eBayes(fit)
res = topTable(fit, adjust="BH", number=nrow(edata))

all_prot_diff_expr_tables[["phospho_proteom_total"]] = res

colnames(res) = paste("phospho_proteom_total", colnames(res), sep="_")
res$ID = rownames(res)

phospho_proteom_human = merge(phospho_proteom_human, res[,c(1,4,5,7)], by.x="site_id_unique", by.y="ID", all.x=T, sort=F)

write.table(phospho_proteom_human, file="../../Results/Proteome/Phosphoproteome_DiffExpression.txt", sep="\t", row.names=F, quote=F)
write.table(proteom_human, file="../../Results/Proteome/Proteome_DiffExpression.txt", sep="\t", quote=F, row.names=F)

