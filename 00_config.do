*=============================================================================
* 00_config.do
* Sets the folder paths used by every script in this repository.
* Run this first. Edit ROOT only (the ROOT line below), then run:  do 00_config.do
*
* Expected layout under ROOT (the folder where you saved this repository):
*   ROOT/NHANES_Harmonized_data.dta  analysis dataset (included in this repository)
*   ROOT/results/                    logs and output files are written here
*=============================================================================

global ROOT "/path/to/your/repository/folder"

global DATA    "$ROOT/NHANES_Harmonized_data.dta"
global RESULTS "$ROOT/results"

capture mkdir "$RESULTS"
cd "$RESULTS"

di "Config loaded. Working directory: $RESULTS"
