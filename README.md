# Visualizing Spatial Patterns of Air Pollution in the Eastern United States  
### Impacts on Human Health  
**Date:** March 27 – April 10, 2024  

---

## 📌 Overview  
This project investigates the spatial distribution of key air pollutants across five Eastern U.S. states: **New York, Pennsylvania, New Jersey, Delaware, and Maryland**. Using a combination of **spatial interpolation** (IDW and Kriging) and **spatial clustering** (SKATER algorithm), we analyze the concentrations of:

- **Nitrogen Dioxide (NO₂)**  
- **Ground-level Ozone (O₃)**  
- **Fine Particulate Matter (PM2.5)**  

Our study reveals that pollutant concentrations across the region remain *within EPA safety standards* overall. However, clustering identifies notable hotspots, including elevated NO₂ in New Jersey, high O₃ levels in Northern Maryland, and rising PM2.5 levels in Pennsylvania.

---

## 🎯 Objectives  
- Estimate pollution levels at unmonitored locations using **IDW** and **Kriging**.  
- Detect spatial patterns and hotspots using **SKATER hierarchical clustering**.  
- Compare pollutant concentrations with **EPA standards** to assess potential health risks.  
- Visualize results through maps, variograms, and interpolated surfaces.

---

## 🗺️ Study Area  
The analysis focuses on the contiguous states of:  
**NY, PA, NJ, DE, MD**  
These states provide diverse climate, population density, and geographic characteristics critical for examining pollution distribution.

---

## 📂 Data Sources  
- **EPA Air Quality System (AQS) 2021 Annual Summary**  
  - NO₂  
  - O₃  
  - PM2.5  
- **State boundary shapefiles**  
  - 2021 U.S. Census Bureau TIGER/Line  
- Additional tools: QGIS, R (gstat, sp, tidyverse), GeoDa.

---

## 🧪 Methods  

### **1. Spatial Interpolation**
#### 🔹 Inverse Distance Weighting (IDW)
- Distance weight **k** optimized via **Leave-One-Out Cross-Validation (LOOCV)**.
- Selected *k* values:  
  - NO₂: **1**  
  - O₃: **1**  
  - PM2.5: **2**  

#### 🔹 Ordinary Kriging
- Variograms fitted for each pollutant (Spherical/Wave/Linear).
- LOOCV used to select models with lowest RMSE.  
- Assumptions checked: normality, stationarity, isotropy.

---

### **2. Spatial Clustering**  
Performed in **GeoDa** using the **SKATER** algorithm with Queen contiguity weights.

Cluster counts:
- **NO₂: 4 clusters**  
- **O₃: 10 clusters**  
- **PM2.5: 6 clusters**  

Minimum cluster size ensured to avoid single-element clusters.  
Results exported and visualized using R’s ggplot.

---

## 📊 Results Summary

### **NO₂**
- All interpolated values fall well **below EPA’s 53 ppb annual standard**.
- RMSE:
  - IDW: **3.82**
  - Kriging: **3.14**
- **Hotspot:** New Jersey shows relatively elevated concentrations despite safe absolute levels.

### **O₃**
- Estimated O₃ levels **0.036–0.043 ppm**, well below the **0.070 ppm EPA limit**.
- RMSE:
  - IDW: **0.0026**
  - Kriging: **0.002**
- **Hotspot:** Northern Maryland identified as the highest O₃ concentration zone.

### **PM2.5**
- All values remain below the **12.0 μg/m³ annual standard**.
- RMSE:
  - IDW: **1.10**
  - Kriging: **0.979**
- **Concern:** Pennsylvania approaching the EPA threshold.

---


---

## 🛠️ Tools & Technologies  
- **R** (gstat, sp, raster, tidyverse)  
- **Python** (optional visualization)  
- **GeoDa** (clustering)  
- **QGIS** (geospatial preprocessing)  
- **ggplot2** (visualizations)  

---

## 📈 Key Findings  
- Pollution levels overall remain within **EPA safe limits**.  
- Spatial clustering identifies **localized hotspots** not visible from summary statistics alone.  
- Results emphasize the importance of **continuous monitoring** and **proactive air-quality management**.

---

## 📜 License  
This project is available under the **MIT License** unless otherwise specified.

