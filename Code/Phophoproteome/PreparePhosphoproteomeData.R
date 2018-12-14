library(data.table)
phosphosites_final = fread("../../Data/External/Phosphosite/Phopshosites_diseases_and_regulatory.txt", sep="\t")
setorder(phosphosites_final, prot_residue_id)
setcolorder(phosphosites_final, c(ncol(phosphosites_final), 1:(ncol(phosphosites_final)-1)) )

full_data_file = "../../Data/Raw/Proteome/Perseus_rep2-3_20150330.txt.gz"
data = fread(paste("gzip -dc",full_data_file),sep="\t",header=T)
tmp = data[,100:ncol(data),with=F]
data = data[,!duplicated(colnames(data)), with=F]
data$phosphorylation_detected_cnt = apply(data[,c("Nuc_rep2_phospho", "Nuc_rep3_phospho", "Total_rep2_phospho", "Total_rep3_phospho"),with=F], 1, function(x) sum(!is.nan(x)) )
data$Gene_names = apply(tmp, 1, function(x) paste(unique(x[which((!is.na(x) & x!=""))]), collapse=";") )

setorder(data, "Protein", "Sequence window","Modification window")

# we usually find several protein identifiers per identified site. Only one or few will match Phopshosite plus data base IDs
# We therefore assume that the group of possible proteins for that site can be represented by the matching proteins
# from Phosphosite plus (which is usually the reviewed protein variant for a gene from SwissProt). 
# so here we first expand data so that each row contains a single protein accession
prot_and_pos = data[,c("Proteins","Positions within proteins","PhosphoSitePlus window","Unique identifier"),with=F]
uid = prot_and_pos$"Unique identifier"
prot_split = strsplit(data$Proteins,";")
pos_split = strsplit(data$"Positions within proteins",";")
prot_cnts = sapply(prot_split, length)
pos_cnts = sapply(pos_split, length)
if(any(which(prot_cnts!=pos_cnts)) ) stop("Conflicting protein ID and position numbers")
expanded_prot_ids = data.table(uid=rep(uid, times=prot_cnts), prot_id = unlist(prot_split), pos = unlist(pos_split), ppsite_window = rep(prot_and_pos$"PhosphoSitePlus window", times=prot_cnts) )

# unfortunately, there might also be more than one phosphosite IDs - will have to split them, too

ppsite_ids_split = lapply(strsplit(expanded_prot_ids$ppsite_window,";"), function(x) if(length(x)==0) {""} else {x} )
#ppsite_ids_split = strsplit(expanded_prot_ids$ppsite_window,";")
ppsite_ids_len = sapply(ppsite_ids_split, length)

expanded_prots_ids_and_site_windows = data.table(uid = rep(expanded_prot_ids$uid, times=ppsite_ids_len), 
                                                 prot_id = rep(expanded_prot_ids$prot_id, times=ppsite_ids_len ), 
                                                 pos = rep(expanded_prot_ids$pos, times=ppsite_ids_len ), 
                                                 ppsite_window = unlist(ppsite_ids_split) )

setkey(expanded_prots_ids_and_site_windows, "prot_id","ppsite_window", "pos")

phosphosites_final$pos = substr(phosphosites_final$MOD_RSD,2,nchar(phosphosites_final$MOD_RSD))
setkey(phosphosites_final, "ACC_ID","SITE_ID", "pos")

# merge protein site IDs with phosphosite plus data base
pp = expanded_prots_ids_and_site_windows[phosphosites_final, allow.cartesian=TRUE]
mpp = pp[uid %in% names(which(table(pp$uid)>1))]
setorder(mpp, "uid","prot_id","pos")
# some UIDs match to several entries in phosphosite. Aggregate them all on one line
pp_ts = melt(pp, id.vars = "uid")
pp_sw = dcast.data.table(pp_ts, uid ~ variable, fun.aggregate=function(x) paste(unique(x[which(!is.na(x))]), collapse=";"), value.var="value")

# now finally merge the original data table
setkey(pp_sw, "uid")
setkey(data, "Unique identifier")

data_merged = pp_sw[data]
setcolorder(data_merged, c(colnames(data)[!colnames(data)=="Unique identifier"], colnames(pp_sw)))
setorder(data_merged, "Leading proteins", "pos", "Charge")

save(data_merged, file= "../../Data/Processed/Proteome/Phosphoproteome_infection_annotated_with_phopsphosite_plus.Rdata")
