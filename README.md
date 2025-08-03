# Gradient-Guided Search for Assured Contingency Landing Management

This repository provides a C/C++ implementation of a **gradient-guided search algorithm** for real-time emergency landing path planning in dense and constrained airspace. The framework enables air traffic deconfliction and minimizes operational disruptions during emergency mitigation by combining spatial risk modeling, 3D trajectory search, and heuristic optimization.

Developed at the **Autonomous Aerospace Systems Laboratory (A2Sys)**, Kevin T. Crofton Department of Aerospace and Ocean Engineering, Virginia Tech.

---

## 🔑 Key Features

- **Gradient-Guided 4D Search**  
  Enhances real-time emergency landing path planning with gradient-based expansion in 4D state space.

- **Air Traffic Deconfliction**  
  Incorporates traffic-aware risk to avoid congested regions during emergency descent.

- **3D Airspace Risk Modeling**  
  Considers dense arrival/departure corridors, urban air corridors, and restricted zones.

- **ADS-B Driven Risk Estimation**  
  Uses historical Automatic Dependent Surveillance–Broadcast (ADS-B) data for probabilistic traffic density risk.

- **Proximity-Based 3D Heatmaps**  
  Constructs spatial risk zones around critical urban air corridors and no-fly areas via computational geometry.

- **Cumulative Time-Exposure Metric**  
  Quantifies airspace risk based on path duration through high-density volumes.

- **Look-Ahead Heuristic in Congested Airspace**  
  Improves search efficiency by forecasting future expansion penalties in traffic-heavy zones.

- **Landing Site Selection for Operational Continuity**  
  Ranks candidate sites based on minimal disruption to active traffic corridors.

## 📊 Datasets Provided

The repository includes curated datasets for modeling urban airspace, risk-aware descent, and air traffic deconfliction in the Washington, D.C. area. These datasets support simulation, benchmarking, and deployment of the contingency landing planner.

### 1️⃣ Gliding Aircraft Performance

- **Cessna 182 Gliding Parameter Set and Flight Path Angle Look-Up Table**  
  Precomputed table for estimating glide descent angles as a function of steady wind and aircraft heading.
  **Source:** H. E. Tekaslan and E. M. Atkins, “Vehicle-to-Vehicle Approach to Assured Aircraft Emergency Road Landings,” _AIAA Journal of Guidance, Control, and Dynamics_, vol. 48, no. 8, pp. 1800–1817, 2025.

### 2️⃣ Ground Risk Data

- **Washington D.C. Area Census Dataset**  
  Provided as shapefiles, includes polygonal population zones used to compute ground risk exposure along descent paths using spatial R-tree queries.

### 3️⃣ Airport and Runway Data

- **Washington D.C. Airport and Runway Geometry**  
  Binary file containing runway endpoints, headings, and elevations for major airports in the region. Used to initialize candidate landing sites.

### 4️⃣ Helicopter Route Data

- **Urban Helicopter Route Geometry and Heatmaps**  
  Binary files containing low-altitude helicopter route shapes and associated airspace density heatmaps. Modeled using proximity-based risk for mixed-use corridors.

### 5️⃣ No-Fly Zone Data

- **Restricted and Prohibited Area Geometry and Heatmaps**  
  Binary files representing permanent and temporary no-fly zones in the Washington D.C. airspace. Includes 3D volumetric heatmaps for risk-aware planning.

### 6️⃣ Commercial Airport Corridor Heatmaps

- **Ronald Reagan Washington National Airport (DCA)**  
  Heatmaps capturing historical departure and arrival flows based on ADS-B data. Used to evaluate air traffic deconfliction and trajectory feasibility.

> All datasets are stored in the `data/` or `lib/<lib-name>/data`  directory. Heatmap binary files are preprocessed and ready to be utilized by the planner modules.


---

## 🔧 Build Instructions

### 1. Dependencies

This project uses standard C/C++ libraries and a few geospatial dependencies. On macOS, install via Homebrew:

```bash
brew install proj spatialindex shapelib
```
### 2. Build and Run
Configuration file `aclm.cfg` can be used for fundamental settings. To build and run the software:

```bash
cd <directory-to-GGS_ACLM> make && clear && ./test
```

### 3. Post-process
Results are written in `results/` folder.
