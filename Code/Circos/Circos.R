library(circlize)

cairo_pdf(file="../../Results/Proteome/Circos_plot.pdf", width=20, height=20)

inp_data = read.table("../../Data/Processed/Final_list_v2_with_TF_for_Circos.csv", sep="\t", header=T, stringsAsFactors = F)
inp_data$site_ID = paste(inp_data$Gene.name, inp_data$Phosphosite)
inp_data$log2ratio_nuclear = -apply(inp_data[,c("Nuc_rep2_phospho", "Nuc_rep3_phospho")], 1, mean, na.rm=T)
inp_data$log2ratio_total = -apply(inp_data[,c("Total_rep2_phospho", "Total_rep3_phospho")], 1, mean, na.rm=T)
inp_data$log2ratio_total_limited = ifelse(abs(inp_data$log2ratio_total) > 2.5, sign(inp_data$log2ratio_total) * 2.5, inp_data$log2ratio_total )
inp_data$log2ratio_nuclear_limited = ifelse(abs(inp_data$log2ratio_nuclear) > 2.5, sign(inp_data$log2ratio_nuclear) * 2.5, inp_data$log2ratio_nuclear )
inp_data$nuclear_cnt = apply(!is.na(inp_data[,c("Nuc_rep2_phospho", "Nuc_rep3_phospho")]), 1, sum)
inp_data$total_cnt = apply(!is.na(inp_data[,c("Total_rep2_phospho", "Total_rep3_phospho")]), 1, sum)

# Collation in R3.2 is a mess - it ignores special characters completely when sorting. So here use a character as separator. 
inp_data = inp_data[order(paste(inp_data$Kinase.Family, inp_data$site_ID, sep="Z")), ]

all_groups = unique(inp_data$Kinase.Family)
sectors_ordered = factor(inp_data$Kinase.Family, levels=all_groups)
entries_per_group = tapply(inp_data$Kinase.Family, inp_data$Kinase.Family, length)[all_groups]

inp_data$order = unlist(lapply(entries_per_group, function(n) seq(1:n)))

par(mar = c(2, 2, 2, 2), lwd = 0.1, cex = 0.7)

circos.par("track.height" = 0.97, gap.degree = 0, canvas.xlim=c(-1.5,1.5), canvas.ylim=c(-1.5,1.5))

xlims = cbind(rep(0.5, length(entries_per_group)), entries_per_group + .5)
circos.initialize(factors = sectors_ordered, x = inp_data$order, xlim=xlims)

min_y = min(min(inp_data$log2ratio_nuclear_limited, na.rm=T),min(inp_data$log2ratio_total_limited, na.rm=T))
max_y = max(max(inp_data$log2ratio_nuclear_limited, na.rm=T),max(inp_data$log2ratio_total_limited, na.rm=T))
circos.trackPlotRegion(factors = inp_data$Kinase.Family, 
                       y = inp_data$log2ratio_nuclear_limited, ylim=c(min_y,max_y) )


cols = rainbow(length(all_groups))
names(cols) = all_groups

for (g in all_groups) {
  circos.updatePlotRegion(g, 1, bg.col = cols[g]) 
  sector_data = subset(inp_data, Kinase.Family==g)
  circos.axis(labels=paste(" ", sector_data$site_ID, " ", sep=""), labels.cex=1.5, major.at=sector_data$order, major.tick.percentage = 0.001, minor.ticks = 0, labels.facing="clockwise", direction="outside")
  circos.axis(labels=F, labels.cex=1.5, major.at=sector_data$order, major.tick.percentage = 0.02, minor.ticks = 0, labels.facing="clockwise", direction="inside")

    # Nuclear
  sector_data_valid = subset(sector_data, !is.na(log2ratio_nuclear_limited))
  if(nrow(sector_data_valid) > 0 ) {
    symbol = ifelse(sector_data_valid$nuclear_cnt==2, 16, 17)
    circos.lines(sector_data_valid$order, sector_data_valid$log2ratio_nuclear_limited, lwd=1.7)
    circos.points(sector_data_valid$order, sector_data_valid$log2ratio_nuclear_limited, pch=symbol, cex=2.3, col="black")
  }
  
  sector_data_invalid = subset(sector_data, is.na(log2ratio_nuclear_limited))
  if(nrow(sector_data_invalid)>0) {
    circos.points(sector_data_invalid$order, rep(0, nrow(sector_data_invalid)), pch=22, cex=2, bg="black")
  }

  # Total
  sector_data_valid2 = subset(sector_data, !is.na(log2ratio_total_limited))
  if( nrow(sector_data_valid2) > 0 ) {
    symbol = ifelse(sector_data_valid2$nuclear_cnt==2, 16, 17)
    circos.lines(sector_data_valid2$order, sector_data_valid2$log2ratio_total_limited, lwd=1.7, col="grey60")
    circos.points(sector_data_valid2$order, sector_data_valid2$log2ratio_total_limited, pch=symbol, cex=2.3, col="grey60")
  }
  sector_data_invalid2 = subset(sector_data, is.na(log2ratio_total_limited))
  if(nrow(sector_data_invalid2)>0) {
    circos.points(sector_data_invalid2$order, rep(0, nrow(sector_data_invalid2)), pch=25, cex=2, bg="grey50")
  }
  
  max_x = max(sector_data$order)
  x_100 = seq(0.5,max_x+0.5, max_x/100)
  circos.lines(x_100, rep(0.5, length(x_100)), lty=2, lwd=1.7)
  circos.lines(x_100, rep(-0.5, length(x_100)), lty=2, lwd=1.7)
  circos.lines(x_100, rep(0, length(x_100)), lty=1, lwd=1)
  
  circos.text(1, max_y-0.1, g, facing = "clockwise", niceFacing = T, cex = 1.5, adj=c(1,0.5))
  
  if(g=="NEK") {
    circos.yaxis("left", at=c(seq(-2, 2, by=1)), labels.cex = 1.5, tick.length=1)
  }
}

circos.clear()
dev.off()