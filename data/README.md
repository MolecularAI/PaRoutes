# PaRoutes dataset and models

PaRoutes provides two sets for benchmarking: **n1** and **n5**.

### n1-files

* `n1-targets.txt` - SMILE strings of the target molecules
* `n1-stock.txt` - SMILE strings of the stock molecules
* `n1-routes.json` - the reference routes, used by the `route_quality.py` script to compute accuracies

### n5-files

* `n5-targets.txt` - SMILE strings of the target molecules
* `n5-stock.txt` - SMILE strings of the stock molecules
* `n5-routes.json` - the reference routes, used by the `route_quality.py` script to compute accuracies

## Files for modelling

To allow training of one-step retrosynthesis models, PaRoutes provides:

* `uspto_raw_template_library.csv` - curated subset of USPTO dataset that was used to extract the reference routes
* `uspto_rxn_n1_raw_template_library.csv` - curated subset of USPTO with reactions not found in the **n1** routes
* `uspto_rxn_n5_raw_template_library.csv` - curated subset of USPTO with reactions not found in the **n5** routes

There are also the following files with trained template-based one-step models:

* `uspto_rxn_n1_keras_model.hdf5` and `uspto_rxn_n1_unique_templates.hdf5` - Keras model and associated templates trained on data not found in the **n1** routes
* `uspto_rxn_n5_keras_model.hdf5` and `uspto_rxn_n5_unique_templates.hdf5` - Keras model and associated templates trained on data not found in the **n5** routes

## Other files

* `chembl_10k_route_distance_model.ckpt` - trained LSTM model to compute distances between synthetic routes
* `150k_routes.json.gz` - all ~150k routes extracted from the USPTO dataset