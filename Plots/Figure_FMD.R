library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Polychrome)
library(ggsci)
library(ggthemes)
library(viridis)
library(RColorBrewer)
data <- read.csv("./FMD_data/FMD_Cumbria_Subset.csv")
data <- data[data$Infection < 100, ] #filter data with infection time less than 100
#cluster
chain_len<-2000
DPMM_statistics<-function(assignments, samples, centroids, chain_len){
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  means<-colMeans(samples[(chain_len/2):chain_len,])
  centroids_means<-colMeans(centroids[(chain_len/2):chain_len,])
  modes <- sapply(assignments[(chain_len/2):chain_len,], Mode)
  ci_95 <- apply(samples[(chain_len/2):chain_len, ], 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  res <- list("means"=means,"modes"=modes,"ci_95"=ci_95, "centroids"=centroids_means)
  return(res)
}
assignments<-read.csv("FMD_data/assignments_singlephi.csv")[,-1]
samples<-read.csv("FMD_data/samples_singlephi.csv")[,-1]
centroids<-read.csv("FMD_data/centroids_singlephi.csv")[,-1]
res<-DPMM_statistics(assignments, samples, centroids, chain_len)
data$DPMM_cluster<-res$modes
data$DPMM_cluster<-as.numeric(factor(data$DPMM_cluster))
#UK map Cumbria
uk_map <- ne_states(country = "United Kingdom", returnclass = "sf")
uk_map_osgb <- st_transform(uk_map, crs = 27700)
bounding_box <- st_bbox(c(xmin = 310000, xmax = 390000, ymin = 510000, ymax = 570000), crs = 27700)
cropped_map <- st_crop(uk_map_osgb, bounding_box)
points_osgb <- data.frame(
  easting = data$X * 1000,
  northing = data$Y * 1000,
  cull_status = ifelse(data$Cull == 0, "Not Culled", "Culled"),  # 0 = "Not Culled", 非0 = "Culled"
  cluster = data$DPMM_cluster
)
points_sf <- st_as_sf(points_osgb, coords = c("easting", "northing"), crs = 27700)
ggplot() +
  geom_sf(data = cropped_map, fill = "gray90", color = "black") +  # 地图
  geom_sf(data = points_sf, aes(shape = cull_status, color = cull_status), size = 2, alpha = 0.8) +  # 同时使用 shape 和 color
  scale_color_manual(values = c("Not Culled" = "blue", "Culled" = "red")) +  # 论文推荐配色
  scale_shape_manual(values = c("Not Culled" = 16, "Culled" = 17)) +  # 圆点 vs. 三角形
  labs(x = "Longitude", y = "Latitude", shape = "Cull Status", color = "Cull Status") +  # 统一图例名称
  guides(color = guide_legend(override.aes = list(shape = c(16, 17)))) +  # 让图例合
  theme(panel.background = element_rect(fill = "white", color = "black"),
        legend.position = c(0.8, 0.8),
        plot.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))  
ggsave("./Plots/Cumbria.pdf", width = 8, height = 6, device = "pdf")

#Plot Clustered 
points_sf$cluster <- as.factor(points_sf$cluster)  
colors_27 <- colorRampPalette(brewer.pal(12, "Paired"))(27)  
colors_27 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
               "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
               "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
               "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
               "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
               "#FF69B4", "#00CED1")
ggplot() +
  geom_sf(data = cropped_map, fill = "gray90", color = "black") +  
  geom_sf(data = points_sf, aes(color = factor(cluster)), size = 2, alpha = 0.9) + 
  #scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)+
  scale_color_manual(values = colors_27)+
  labs(x = "Longitude", y = "Latitude") +  # 统一图例名称
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "none")

ggsave("./Plots/Cumbria_clustered.pdf", width = 8, height = 6, device = "pdf")

