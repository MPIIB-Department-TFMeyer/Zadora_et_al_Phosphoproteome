library(readxl)
library(data.table)

# Documentation of iGPS result tables from 
# http://igps.biocuckoo.org/faq.php
# 1. Q: From your example sequence "P53667", the iGPS predicts severl results below:
# (1) S179 AGC/AKT SHGKRGLSVSIDPPH 2.011 1.21 String
# (2) S298 AGC/AKT KPVLRSCSIDRSPGA 1.64   1.21 String
# (3) S323 AGC/PKC KDLGRSESLRVVCRP 3.588 1.98 String
# How to interpret the results? For example, since the scores of (3) > (1) > (2), whether the predictions mean that the correct 
# probability is (3) > (1) > (2)?
# A: This is the mostly asked questions from users. There are two principles for interpreting the results. First, for the 
# same protein kinase (PK) group, higher score means higher probability that a phosphorylation sites can be modified by the 
# PK. Thus, since the score of (1) is greater than (2) for AGC/AKT, the (1) has higher probability to be a real hit than (2). 
# Second, the scores from different PK groups can not be compared, because different training data sets and procedures were 
# performed, and different thresholds were chosen. In this regard, comparison of scores of (3) to (1) and (2) is meaningless. 

##########################################################################################################################
# Conclusion: scores within each PK groups can be weighted by score. Different PK groups cannot be weighted, that's why we
# could only treat them as equivalent. On the other hand, if we would treat them as equivalent, e.g for a gene with 9/10 predictions being from PK
# group A and 1/10 from group B but average scores in B > A, that might introduce another bias, since one could argue that in this case group A should 
# get a higher weight than group B. 
#
# So, let's keep it simple and ignore score and just set the weight for all PKs to 1/N, with N being the number of predictions 
# for a given site. 

input_file = "../../Data/Processed/iGPS.xlsx"
predictions_total = read_xlsx(input_file, sheet = 1)
predictions_nuclear = read_xlsx(input_file, sheet = 2)

pt_dt = as.data.table(unique(predictions_total))
pt_nc = as.data.table(unique(predictions_nuclear))

pt_dt$site_id = paste(pt_dt$"# ID", "_", pt_dt$'Gene Name')
pt_nc$site_id = paste(pt_nc$"# ID", "_", pt_nc$'Gene Name')
# these are unique, i.e. no site belongs to more than one gene. Note that gene 42620 is actually "SEPT7"-

setkey(pt_dt, "site_id")
setkey(pt_nc, "site_id")

#pt_dt[,c("cnt", "score_sum"):= list(.N , sum(Score)), by=site_id]
pt_nc[,c("cnt", "rel"):= list(.N , 1/.N), by=site_id]
pt_dt[,c("cnt", "rel"):= list(.N , 1/.N), by=site_id]

weighted_sum_per_PK_dt = pt_dt[,.(weighted_sum=sum(rel),cnt=.N), by="Kinase Name"]
weighted_sum_per_PK_nc = pt_nc[,.(weighted_sum=sum(rel),cnt=.N), by="Kinase Name"]
setorderv(weighted_sum_per_PK_dt, "weighted_sum", order=-1)
setorderv(weighted_sum_per_PK_nc, "weighted_sum", order=-1)

write.table(weighted_sum_per_PK_dt, file="../../Results/KinomeTree/Weighted_relative_contributions_PK_total.txt", sep="\t", row.names=F, quote=F)
write.table(weighted_sum_per_PK_nc, file="../../Results/KinomeTree/Weighted_relative_contributions_PK_nuclear.txt", sep="\t", row.names=F, quote=F)
