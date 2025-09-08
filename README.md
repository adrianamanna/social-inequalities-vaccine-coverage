# Social inequalities in vaccine coverage and their effects on epidemic spreading
This project analyzes social inequalities in vaccine coverage using epidemiological modeling with different contact matrix approaches.

## Installation

### Clone Repository
```bash
git clone https://github.com/yourusername/social-inequalities-vaccine-coverage.git
cd social-inequalities-vaccine-coverage
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

## Project Structure

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

## Usage

### Configuration Setup

In `/configs/config_params_to_run.yaml`, set:
- `model_type`: Choose `'Cij'` or `'Gab'` to select the model type
- `params_to_run`: Choose from:
  - `'params1&3'` for figures 1 and 3
  - `'params2&4_NPIs'` for figures 2 and 4 
  - `'params5'` for figure 5


In `/configs/config_epi.yaml`, set `N_iter` equal to the number of simulations to run. 


### Model Requirements by Figure

#### Figures 1 & 3
Run simulations for **both models** (age-stratified and generalized):
```bash
python get_epi_rates.py
python get_prep_epi.py
```

#### Figures 2 & 4  
Run simulations for **generalized model only**:
```bash
python get_epi_rates.py
python get_prep_epi.py
```
*Note: Requires baseline scenario results from figures 1 & 3*

#### Figure 5
Run simulations for **both models** (age-stratified and generalized):
```bash
python get_epi_rates.py
python get_prep_epi.py
```

### Generate Figures
```bash
# Figures 1 & 3: Basic analysis
python make_figure_1&3.py

# Figures 2 & 4: Non-pharmaceutical interventions (NPI) analysis
python make_figure_2&4.py

# Figure 5: True contact matrices from Hungary
python make_figure_5.py
```

## Configuration Files

| File | Description |
|------|-------------|
| `config_epi.yaml` | Epidemiological parameters |
| `config_matrix.yaml` | Contact matrix and population settings |
| `config_params_to_run.yaml` | Scenario definitions and model types |
| `parameters.yaml` | Detailed scenario definitions and model parameters |
| `IFR_age.yaml` | Age-stratified infection fatality rates |

## Model Types

- **G_ab Model**: Generalized contact matrices with age × socioeconomic status stratification
- **C_ij Model**: Traditional age-stratified contact matrix model

## Output

- **Results**: Stored as pickle files in `/res/` directory
- **Figures**: Generated plots saved in `/figs/` directory