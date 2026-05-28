# Satellite Remote Sensing & Wildfire Analytics Portfolio

This repository contains three interconnected data science projects built in R. Together, they show how to use satellite data to measure the environmental damage caused by a major wildfire (The Cameron Peak Fire in Colorado).

---

## 📁 The Three Projects

### Project 1: Satellite Image Fetching & Forest Health Mapping (NDVI)
* **What it does:** Programmatically connects to a live space API to find and download crisp, cloud-free satellite images of the forest before and after the fire.
* **How it works:** It extracts specific light wavelengths (Near-Infrared and Red) from the Sentinel-2 satellite network. It then uses these bands to calculate **NDVI**, which is a scientific score that measures how green and healthy a forest canopy is.
* **Skills shown:** Connecting to data APIs (`rstac`), working with spatial map coordinates, and writing math equations to analyze satellite imagery.

### Project 2: Wildfire Burn Scar & Severity Dashboard
* **What it does:** Automates change detection to isolate the exact boundary and footprint of the wildfire's damage.
* **How it works:** It takes the pre-fire health map and subtracts the post-fire health map. The resulting math isolates exactly where the trees were destroyed. It organizes these maps into a beautiful, publication-ready side-by-side timeline graphic.
* **Skills shown:** Automated change-detection modeling, advanced data visualization (`ggplot2`), and fixing layout bugs like overlapping text labels.

### Project 3: Interactive Web GIS Map
* **What it does:** Transforms static data into an interactive, web-based map that anyone can click around and explore in a web browser.
* **How it works:** It shrinks the satellite data so it runs smoothly on a webpage, overlays it onto a live global satellite background, and automatically downloads and draws official government county borders (like Larimer County) right on top of the map.
* **Skills shown:** Interactive web mapping (`leaflet`), processing vector boundary files (`USAboundaries`), and data optimization for web applications.

---

## 🛠️ Tools & Packages Used
* **Language:** R
* **Spatial Data Engines:** `terra` (for satellite grids) and `sf` (for boundary lines)
* **Graphics & Mapping:** `ggplot2`, `tidyterra`, and `leaflet`
* **Data Streams:** Microsoft Planetary Computer API (`rstac`)
