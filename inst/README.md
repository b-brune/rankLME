# Description of the supplemental code 

## Main Folder

* `simulation_setup.Rmd`: provides a basic description of how the model can be simulated
* `simulation_helpers.R`: helper functions loaded by `simulation_setup.Rmd`
* `results_leverage_with_hat.Rdata` simulation results from `simulation_setup.Rmd`

## Folder `data_examples`

* `sleepstudy_example.R`: provides code to reproduce the sleepstudy example in the paper
* `pv_examples.R`: provides code to reproduce the two applications to photovoltaik data in the paper


## Folder `simulation_study`

* `registry_setup_simulation.R`: Sets up the simulation registry in the `batchtools` framework. Running this code will carry out the simulation studies reported in the paper. Running all experiments on 120 cores took about 8-12 hours.
* `setup_simulation_evaluation.R`: provides code preprocess the simulation results 
* `simulation_evaluation_plots.R`: Code that produces the plots in the paper and supplement
