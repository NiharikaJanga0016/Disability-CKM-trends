*=============================================================================
* 00_config.do
* Sets the folder paths used by every script in this repository.
* Run this first. Edit ROOT only (the ROOT line below), then run:  do 00_config.do
*
* Expected folder layout under ROOT:
*   ROOT/raw/         (optional) NHANES .xpt files downloaded from CDC
*   ROOT/harmonized/  nhanes_harmonized_primary_v2.dta (included in this repository)
*   ROOT/results/     logs and output files are written here
*=============================================================================

global ROOT "/path/to/your/NHANES/folder"

global RAW     "$ROOT/raw"
global DATA    "$ROOT/harmonized/nhanes_harmonized_primary_v2.dta"
global RESULTS "$ROOT/results"

capture mkdir "$ROOT/harmonized"
capture mkdir "$RESULTS"
cd "$RESULTS"

di "Config loaded. Working directory: $RESULTS"
