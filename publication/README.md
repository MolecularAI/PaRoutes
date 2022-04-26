# PaRoutes and AiZynthFinder

This folder contains a few scripts and input files to use PaRoutes with AiZynthFinder. 

The scripts should be sufficient to reproduce the results presented in the PaRoutes publication.

## Running AiZynthFinder

If you have installed AiZynthFinder in your environment, you can run the route predictions using one of the six provided input files.
One per search algorithm and reference set. For instance

    aizynthcli --smiles ../data/n1-targets.txt --config aizynthfinder_config_retro_n1.yml --output aizynthfinder_output_retro_n1.hdf5

To run AiZynthFinder with Retro* on the **n1** reference set.

## Analyzing the output

To use the scripts in the `analysis` folder, the routes need to be extracted from the AiZynthFinder output file. This can be done using for instance


    python extract_trees.py aizynthfinder_output_retro_n1.hdf5

And then you can run the analysis scripts using for instance

    python ../analysis/route_quality.py --routes aizynthfinder_output_retro_n1_trees.json.gz --references ../data/n1-routes.json --output route_quality_retro_n1.csv

for route quality analysis, and

    python ../analysis/route_clusters.py --routes aizynthfinder_output_retro_n1_trees.json.gz --model ../data/chembl_10k_route_distance_model.ckpt --nclusters 0 --min_density 2 --output cluster_analysis_retro_n1.json

for route clustering analysis.

Please consult the README file in the `analysis` folder for further information.


## Files included

* *aizynth_n1_stock.txt* and *aizynth_n5_stock.txt* - the stock molecules in InChI-key format
* *aizynthfinder_config_X_Y.yml* - the AiZynthFinder input configuration for method X on set Y
* *extract_trees.py* - script to extract trees from the AiZynthFinder output
* *retrostar_value_model.pickle* - the weights for the Retro* value model in AiZynthFinder format
* *torch2np_vm_model.py* - a Jupyter notebook that can be used to create *retrostar_value_model.pickle*

