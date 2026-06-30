Satellite Remote Sensing, Wildfire Analytics & Genomics Portfolio

This repository contains four interconnected data science projects built in R, spanning spatial analysis, environmental monitoring, and computational biology.

---

📁 The Four Projects

---

## Project 1: Satellite Image Fetching & Forest Health Mapping (NDVI)

What it does:  
Programmatically connects to a live space API to find and download crisp, cloud-free satellite images of forest landscapes before and after wildfire events.

How it works:  
It extracts specific light wavelengths (Near-Infrared and Red) from the Sentinel-2 satellite network. These spectral bands are used to calculate NDVI (Normalized Difference Vegetation Index), a standardized scientific metric that measures vegetation health and density.

Skills shown:  
Connecting to geospatial data APIs (rstac), working with spatial coordinates, and applying remote sensing mathematics to quantify ecosystem health.

---

## Project 2: Wildfire Burn Scar & Severity Dashboard

What it does:  
Automates change detection to isolate and visualize wildfire burn severity and spatial extent.

How it works:  
It subtracts pre-fire NDVI values from post-fire NDVI values to quantify vegetation loss. The resulting spatial output highlights burn severity gradients and produces publication-quality comparative visualizations.

Skills shown:  
Change detection modeling, geospatial raster analysis, and advanced data visualization using ggplot2.

---

## Project 3: Interactive Web GIS Map

What it does:  
Transforms static geospatial outputs into an interactive web-based mapping application.

How it works:  
Satellite raster data is downsampled for performance, layered onto a dynamic basemap, and combined with administrative boundary shapefiles for interactive exploration.

Skills shown:  
Interactive mapping (leaflet), vector boundary processing (sf), and web-optimized geospatial data engineering.

---

## Project 4: Transcriptomic Analysis of Breast Cancer (Basal vs Normal)

What it does:  
Performs a full RNA-seq-style differential gene expression analysis comparing Basal breast cancer tissue to normal breast tissue using publicly available GEO data.

How it works:  
Gene expression data is retrieved from the NCBI GEO database and processed in R. Differential expression analysis is performed using linear modeling (limma). Significant genes are identified based on adjusted p-values and fold-change thresholds. Functional enrichment analysis (GO and KEGG) is then used to interpret the biological pathways associated with tumor progression. Results are visualized using PCA, volcano plots, heatmaps, and enrichment dotplots in a publication-style multi-panel figure.

Skills shown:  
Bioinformatics data retrieval (GEOquery), differential expression analysis (limma), gene annotation mapping, pathway enrichment analysis (clusterProfiler), and publication-grade visualization.

---

🛠️ Tools & Packages Used

Language: R  

Geospatial & Remote Sensing:
- terra  
- sf  
- tidyterra  
- leaflet  

Bioinformatics & Genomics:
- GEOquery  
- limma  
- clusterProfiler  
- org.Hs.eg.db  
- EnhancedVolcano  
- pheatmap  

Visualization:
- ggplot2  
- patchwork  
- enrichplot  

Data Streams:
- Microsoft Planetary Computer API (rstac)  
- NCBI Gene Expression Omnibus (GEO)
## 🛠️ Tools & Packages Used
* **Language:** R
* **Spatial Data Engines:** `terra` (for satellite grids) and `sf` (for boundary lines)
* **Graphics & Mapping:** `ggplot2`, `tidyterra`, and `leaflet`
* **Data Streams:** Microsoft Planetary Computer API (`rstac`)
