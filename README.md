# Residential Wood Combustion Study

![Air Quality Modeling](https://img.shields.io/badge/domain-air_quality-blue)
![Python](https://img.shields.io/badge/language-python-green)

This repository contains tools for analyzing sensitivity of air quality modeling outputs from the SMOKE emissions modeling system and WRF meteorological modeling for this study analyzing the air quality and public health impacts of residential wood combustion (RWC) in the contiguous United States. Research was conducted by Kyan Shlipak and the Climate Change Research Group under P.I. Daniel Ethan Horton at the Northwestern University Department of Earth, Environmental, and Planetary Sciences.


## Project Structure
├── SMOKE_sensitivity_analyses/ # Main analysis scripts for SMOKE emissions

├── WRF_analyses/ # WRF meteorological model analysis

├── requirements.txt # Python dependencies

├── .gitignore # Git ignore rules

└── .gitattributes # Git configuration



## Key Features

- **SMOKE Emissions Analysis**
  - Sensitivity testing of emission inputs
  - Scenario comparisons (baseline vs. modified) for emissions
  - Emissions processing scripts to prepare for SMOKE

- **WRF Meteorology Analysis**
  - Meterological validation scripts
  - Chemistry validation scripts
  - Air pollution analyses
  - Population impact analyses
  - Health impact analyses
  - Distributional impact analyses
  - Further emissions visualizations and correlations

## Installation

1. Clone repository:
   ```bash
   git clone [repository_url]
   cd smoke-wrf-sensitivity
   ```
2. Install requirements:
  ```bash
  pip install -r requirements.txt
  ```
3. Run any analysis scripts.

## Contact
Researcher & corresponding author: Kyan Shlipak, kyanshlipak2026@u.northwestern.edu

## License
Distributed under the MIT License. See LICENSE for more information.
