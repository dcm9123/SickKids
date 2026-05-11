# Daniel Castaneda Mogollon, PhD
# May 8th, 2026

#install.packages("eulerr")
library(eulerr)


# NS1 vs S2 KOs
path = "/Users/danielcm/Desktop/SickKids/PICRUSt2.6/"
setwd(path)

fit = euler(c(
    "NS1" = 1240,
    "S2" = 466,
    "NS1&S2" = 3315
), shape = "ellipse", )

eulerr_options(fontsize = 38)

plot(fit,
    labels = list(fontsize = 6, fontface = "bold"),
    fills = list(fill = c("NS1" = "darkorange1", "S2" = "darkturquoise", "NS1&S2" = "brown"), alpha = 0.35),
    quantities = list(fontsize = 6, fontface = "bold"),
    edges = list(col = "black", lwd = 2),
    main = list(fontsize = 10, fontface = "bold", label = "KOs"))

#pdf("NS1_S2_KOs_Venn.pdf", width = 1.85, height = 1.85)
dev.off()





# NS1 vs NS6 enriched pathways
path = "/Users/danielcm/Desktop/SickKids/Maaslin2.6/"
setwd(path)

fit2 = euler(c(
    "NS1" = 11,
    "NS6" = 11,
    "NS1&NS6" = 16),
    shape = "ellipse")

png("NS1_NS6_enriched_pathways_when_compared_to_S2_S5_Venn.png", width = 40, height = 40, res = 800, units = "cm")


plot(fit2,
    labels = list(fontsize = 60, fontface = "bold"),
    fills = list(fill = c("NS1" = "red", "NS6" = "yellow", "NS1&NS6" = "orange"), alpha = 0.35),
    quantities = list(
        type = c("counts", "percent"),
        fontsize = 60, fontface = "bold"),
    edges = list(col = "black", lwd = 8))#,
    #main = list(fontsize = 48, fontface = "bold", label = "Enriched pathways in NS1 and NS6\n at w9w10"))
dev.off()


# S2 vs S5 enriched pathways
fit3 = euler(c(
    "S2" = 16,
    "S5" = 59,
    "S2&S5" = 7),
    shape = "ellipse")

png("S2_S5_enriched_pathways_when_compared_to_NS1_NS6_Venn.png", width = 40, height = 40, res = 800, units = "cm")

plot(fit3,
labels = list(fontsize = 60, fontface = "bold"),
fills = list(fill = c("S2" = "blue", "S5" = "green", "S2&S5" = "purple"), alpha = 0.35),
quantities = list(
    type = c("counts", "percent"),
    fontsize = 52, fontface = "bold"),
edges = list(col = "black", lwd = 8))#,
#main = list(fontsize = 48, fontface = "bold", label = "Enriched pathways in S2 and S5\n at w9w10"))
dev.off()
