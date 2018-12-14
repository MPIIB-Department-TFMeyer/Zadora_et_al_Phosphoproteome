library(data.table)

setwd("/home/hilmar/Documents/Data/Protein/ModifiableSites/")

phosphosites = fread("gzip -dc Phosphorylation_site_dataset.gz", sep="\t", header=T)

# Columns that (should) exist in the table:
# "PROTEIN"   <- Protein ID (identifies ortholog proteins)
# "ACC_ID"    <- SwissProt/TrEMBL protein accession
# "GENE"      <- Gene name
# "HU_CHR_LOC" <- chromosomal location on human genome
# "MOD_TYPE"  <- Modification type
# "MOD_RSD"  <- Modified Residue
# "SITE_GRP_ID"  <- groups equivalent residue sites across species (which might have differen positions in AS seq)
# "ORGANISM"
# "MW_kD"  <- Protein weight
# "DOMAIN" <- Protein domain (if in a known domain)
# "SITE_+/-7_AA" <- SITE ID 
# "LT_LIT"  <- "The number of associated literature records derived using low-throughput experimental techniques. An LTP result may be more reliable than an MS2 result."
# "MS_LIT" <- "The number of published articles using proteomic MS experiments to locate residues that are modified."
# "MS_CST" <- "The number of associated shotgun proteomic experiments performed at Cell Signaling Technology (CST) in which the indicated PTM was observed."
# "CST_CAT#"  <- Catalog IDs for antibodies for that site
setnames(phosphosites, "SITE_+/-7_AA","SITE_ID")
#setkey(phosphosites, PROTEIN)
setkey(phosphosites, PROTEIN, SITE_ID)

# We want to map the data from phosphosites data base table to inhouse (human) data using mainly Uniprot accession and SITE ID. 
# We are mainly interested in a) if the site has been reported before, b) in which model system, c) with LT or MS methods (number)
# Adding d) a possibly existing antibody exists
# There might be sites where phosphorylation / PTM has only been described in one species but not tested or found yet in another
# non-PTM-sites for a given species are not listed in the tables, so we can't find them using exact string comparison with the site ID

human_ps=phosphosites[ORGANISM=="human"]
human_ps$LT_LIT[is.na(human_ps$LT_LIT)] <- 0
human_ps$MS_LIT[is.na(human_ps$MS_LIT)] <- 0
human_ps$MS_CST[is.na(human_ps$MS_CST)] <- 0

# aggregate information from other species
other_species_ps=phosphosites[ORGANISM %in% c("mouse","rat")] 
other_species_ps = other_species_ps[!(is.na(LT_LIT) & is.na(MS_LIT) & is.na(MS_CST))]
#table(phosphosites$ORGANISM)

tmp_lt_lit = dcast.data.table(other_species_ps,  PROTEIN + SITE_GRP_ID ~ ORGANISM, value.var="LT_LIT", fun.aggregate=function(x) paste(x,collapse=","))
tmp_lt_lit$mouse[tmp_lt_lit$mouse %in% c("","NA")] <- 0
tmp_lt_lit$rat[tmp_lt_lit$rat %in% c("","NA")] <- 0
setnames(tmp_lt_lit, c("mouse","rat"), paste("LT_LIT", c("mouse","rat"), sep="_"))

tmp_ms_lit = dcast.data.table(other_species_ps,  PROTEIN + SITE_GRP_ID ~ ORGANISM, value.var="MS_LIT", fun.aggregate=function(x) paste(x,collapse=","))
tmp_ms_lit$mouse[tmp_ms_lit$mouse %in% c("","NA")] <- 0
tmp_ms_lit$rat[tmp_ms_lit$rat %in% c("","NA")] <- 0
setnames(tmp_ms_lit, c("mouse","rat"), paste("MS_LIT", c("mouse","rat"), sep="_"))

tmp_ms_cst = dcast.data.table(other_species_ps,  PROTEIN + SITE_GRP_ID ~ ORGANISM, value.var="MS_CST", fun.aggregate=function(x) paste(x,collapse=","))
tmp_ms_cst$mouse[tmp_ms_cst$mouse %in% c("","NA")] <- 0
tmp_ms_cst$rat[tmp_ms_cst$rat %in% c("","NA")] <- 0
setnames(tmp_ms_cst, c("mouse","rat"), paste("MS_CST", c("mouse","rat"), sep="_"))

setkey(tmp_lt_lit, PROTEIN, SITE_GRP_ID)
setkey(tmp_ms_lit, PROTEIN, SITE_GRP_ID)
setkey(tmp_ms_cst, PROTEIN, SITE_GRP_ID)

other_species_combined = merge(tmp_lt_lit, tmp_ms_lit, all=T)
other_species_combined = merge(other_species_combined, tmp_ms_cst, all=T)

setkey(human_ps, PROTEIN, SITE_GRP_ID)
phosphosites_final = merge(human_ps, other_species_combined, all.x=T)

## disease associated sites 
disease_sites = fread("gzip -dc Disease-associated_sites.gz", sep="\t",  header=T)
ds_sites_short = disease_sites[,c("PROTEIN","SITE_GRP_ID","DISEASE","ALTERATION","PMIDs","MOD_RSD","ORGANISM"), with=FALSE]
ds_molten = melt(ds_sites_short, id.vars=c("PROTEIN","SITE_GRP_ID"), measure_vars = c("DISEASE","ALTERATION","PMIDs","MOD_RSD","ORGANISM"))
tmp_ds_disease = dcast.data.table(ds_molten,  PROTEIN + SITE_GRP_ID ~ variable, value.var="value", fun.aggregate=function(x) paste(unique(x),collapse=";"))

setnames(tmp_ds_disease,c("DISEASE","ALTERATION","PMIDs","MOD_RSD","ORGANISM"), paste("DISEASE",c("DISEASE","ALTERATION","PMIDs","MOD_RSD","ORGANISM"),sep="_"))
setkey(tmp_ds_disease, PROTEIN, SITE_GRP_ID)
phosphosites_final = merge(phosphosites_final, tmp_ds_disease, all.x=T)

## regulatory sites 
reg_sites = fread("gzip -dc Regulatory_sites.gz", sep="\t",  header=T)
reg_sites_short = reg_sites[,c("PROTEIN","SITE_GRP_ID","ON_FUNCTION","ON_PROCESS","ON_PROT_INTERACT","ON_OTHER_INTERACT","NOTES","PMIDs","MOD_RSD","ORGANISM"), with=FALSE]
rs_molten = melt(reg_sites_short, id.vars=c("PROTEIN","SITE_GRP_ID"), measure_vars = c("ON_FUNCTION","ON_PROCESS","ON_PROT_INTERACT","ON_OTHER_INTERACT","NOTES","PMIDs","MOD_RSD","ORGANISM"))
tmp_rs_disease = dcast.data.table(rs_molten,  PROTEIN + SITE_GRP_ID ~ variable, value.var="value", fun.aggregate=function(x) paste(unique(x),collapse=";"))

setnames(tmp_rs_disease,c("PMIDs","MOD_RSD","ORGANISM","NOTES"), paste("REGULATORY",c("PMIDs","MOD_RSD","ORGANISM","NOTES"),sep="_"))
setkey(tmp_rs_disease, PROTEIN, SITE_GRP_ID)
phosphosites_final = merge(phosphosites_final, tmp_rs_disease, all.x=T)

# write data set

phosphosites_final$prot_residue_id = paste(phosphosites_final$ACC_ID,phosphosites_final$MOD_RSD, sep="_")
setorder(phosphosites_final, prot_residue_id)
setcolorder(phosphosites_final, c(ncol(phosphosites_final), 1:(ncol(phosphosites_final)-1)) )
save(phosphosites_final, file="Phopshosites_diseases_and_regulatory_data_table.Rdata")

write.table(phosphosites_final, file="Phopshosites_diseases_and_regulatory.txt",sep="\t",dec=".",row.names=F, quote=F)
