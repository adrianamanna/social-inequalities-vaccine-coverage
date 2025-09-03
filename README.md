## Installation

## Install folder
```
git clone https://github.com/yourusername social-inequalities-vaccine-coverage.git
cd social-inequalities-vaccine-coverage
```

## Install requirements
```
pip install -r requirements.txt

```
## Project structure
```
├── code/
│   ├── get_epi_rates.py          # Calculate epidemiological rates
│   ├── get_prep_epi.py           # Epidemic simulations 
│   ├── make_figure_1&3.py        # Generate figures 1 & 3
│   ├── make_figure_2&4.py        # Generate figures 2 & 4 (NPIs)
│   ├── make_figure_5.py          # Generate figure 5 (true matrices)
│   └── utils/
│       ├── EpiModels.py          # Main SEIRD model class
│       ├── data_tools.py         # Data loading and processing
│       ├── fig_tools.py          # Visualization utilities
│       ├── init_EpiModel.py      # Model initialization
│       ├── init_EpiParams.py     # Parameter initialization
│       ├── init_MatrixPop.py     # Population and contact matrices
│       ├── matrix_tools_sy.py    # Synthetic matrix generation
│       ├── tools.py              # YAML utilities
│       └── LABS.py               # Constants and labels
├── configs/                      # YAML configuration files
├── res/                          # Results storage (pickle files)
├── figs/                         # Generated figures
└── data/                         # Input data files
```

# Usage

## 1. Epidemic simulations
### Generate epidemiological rates at the end of the epidemic f
python get_epi_rates.py

### Run and process the dynamics of the epidemic over time
python get_prep_epi.py

## 2. Generate figures
### Figures 1 & 3
python make_figure_1&3.py

### Figures 2 & 4: NPI analysis
python make_figure_2&4.py

### Figures 5: True contact matrices from Hungary
python make_figure_5.py


# Configuration
- config_epi.yaml: Core epidemiological parameters
- config_matrix.yaml: Contact matrix and population settings
- config_params_torun.yaml: Scenario definitions and model types
- parameters.yaml: Detailed scenario definitions and model parameters
- IFR_age.yaml: Age-stratified infection fatality rates


# Model Types
- Gab Model: Generalized contact matrices with age × socioeconomic status stratification
- Cij Model: Traditional age-stratified contact matrix model