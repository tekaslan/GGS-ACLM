# Gradient-Guided Search for Assured Contingency Landing Management

This repository provides a C/C++ implementation of a **Gradient-guided Search for Assured Contingency Landing Planning** in dense and constrained airspace. The framework enables path planning with air traffic deconfliction and minimizes operational disruptions during emergency mitigation by combining spatial risk modeling and 4D trajectory search considering degraded fixed-wing aircraft performance.

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

The repository includes curated datasets for Cessna 182 gliding aircraft performance, as well as urban airspace and population density in the Washington, D.C. area. These datasets support simulation and benchmarking of emergency landing scenarios in congested environments.


### 1️⃣ Gliding Aircraft Performance

- **Cessna 182 Gliding Parameter Set and Flight Path Angle Look-Up Table**  
  Precomputed table for estimating glide descent angles as a function of steady wind and aircraft heading.

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
The provided test case simulates joint airspace and ground risk-aware emergency landing planning for an engine-out Cessna 182 over Washington, D.C. To build and run the test case:

```bash
cd <directory-to-GGS_ACLM> make && clear && ./test
```
### 3. Post-process
Results are written in `results/` folder.

### 4. Custom Case
To run a different case, users need to input emergency state coordinates and heading into the configuration file `aclm.cfg`.\
Make sure it is a reachable case defined within the modeled airspace environment (e.g., Washington D.C. by default).\
Goal and touchdown states may be left blank as they are set by the landing site selection module.\
Note: The default Cessna 182 dataset supports maximum 8 m/s (~ 15.6 knots) of wind speed.

### 5. Software Execution Exit Flags
0  : Search solver has found a solution to the best landing site.\
1  : Search solver has found a solution to an alternate landing site.\
2  : Search solver has not converged within the maximum allowed state expansions. Dubins solver has returned a fallback solution. Increasing MAX_ITER in `aclm.cfg` may help to find a search-based solution.\
-1 : Emergency state has been found unreachable.\
-2 : Emergency state has been initialized nearby a prohibited area.\
-3 : Search open-list (priority queue) is empty. Search-based solution has not been found.\

## 📚 Data Sources

References for datasets:

- **Gliding Cessna 182 Model and Flight Path Angle Look-Up Table**  
  H. E. Tekaslan and E. M. Atkins, “Vehicle-to-Vehicle Approach to Assured Aircraft Emergency Road Landings,” _AIAA Journal of Guidance, Control, and Dynamics_, vol. 48, no. 8, pp. 1800–1817, 2025.

- **Washington D.C. Census Data (Shapefiles)**  
  U.S. Census Bureau, 2020 TIGER/Line Shapefiles\
  🔗 https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html

- **Airport and Runway Data (Washington D.C.)**  
  Derived from FAA Airport GIS and FAA NASR datasets.\
  🔗 SkyVector, “Skyvector aeronautical charts,” https://skyvector.com, 2025, Accessed on March 2025.

- **Helicopter Route Data**  
  Constructed from publicly available FAA VFR Sectionals and regional helicopter corridors.\
  🔗 SkyVector, “Skyvector aeronautical charts,” https://skyvector.com, 2025, Accessed on March 2025.

- **No-Fly Zone Geometry**  
  Constructed from publicly available FAA VFR Sectionals and regional helicopter corridors.\
  🔗 SkyVector, “Skyvector aeronautical charts,” https://skyvector.com, 2025, Accessed on March 2025.

- **ADS-B Data**  
  Historical ADS-B data collected via OpenSky Network.\
M. Schäfer, M. Strohmeier, V. Lenders, I. Martinovic, and M. Wilhelm, “Bringing up opensky: A large-scale ads-b sensor network for research,” in _IPSN-14 Proceedings of the 13th International Symposium on Information Processing in Sensor Networks_, 2014, pp. 83–94\
  🔗 https://opensky-network.org

---

If you use this software in your research or projects, please cite:

Tekaslan, H. E., & Atkins, E. M. *Aircraft Contingency Landing Planning in Dense Airspace*. 2025, arXiv preprint arXiv:XXXX.XXXXX. https://arxiv.org/abs/XXXX.XXXXX

```bibtex
@misc{ggs-aclm,
  title     = {Aircraft Contingency Landing Planning in Dense Airspace},
  author    = {H. Emre Tekaslan and Ella M. Atkins},
  year      = {2025},
  eprint    = {},
  archivePrefix = {arXiv},
  primaryClass  = {},
  url       = {}
}
```

## 👤 Author

**H. Emre Tekaslan**  
Autonomous Aerospace Systems Laboratory (A2Sys)  
Kevin T. Crofton Department of Aerospace and Ocean Engineering  
Virginia Tech, Blacksburg, VA, USA

- 📧 Email: [tekaslan@vt.edu](mailto:tekaslan@vt.edu)  
- 🔗 Google Scholar: [https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en](https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en)  
- 🔗 LinkedIn: [https://www.linkedin.com/in/tekaslan/](https://www.linkedin.com/in/tekaslan/)
