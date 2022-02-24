# Setting up PaRoutes

This folder provides the script used to extract the PaRoutes datasets. It is provided for reference. For most users, the data provided in the `../data`
folder should be sufficient.

Roughly, these commands were used to extract routes from the USPTO dataset provided in the `../data` folder.

    python extract_uspto_data.py --filename ../data/uspto_raw_template_library.csv
    
    python extract_routes.py --max-workers 32
    
    python analyze_routes.py

Up to this point the procedure is the same for the **n1** and **n5** benchmark sets. The output is a pickled dictionary of extracted, re-formatted and analyzed routes.

Next we perform a number of additional sets to find non-overlapping routes that
are diverse. This will depend on how many route that are extracted from each patent.
    
    python find_non_overlaps.py --output non_overlapping_routes.pickle --routes-per-patent 1
    
    python select_routes.py --filename non_overlapping_routes.pickle --model ../data/chembl_10k_route_distance_model.ckpt --output selected_routes.pickle

    python extract_training_data.py --filename ../data/uspto_raw_template_library.csv --routes ref_routes.json  --output_tmpl uspto_tmpl_raw_template_library.csv --output_rxn uspto_rxn_raw_template_library.csv

For more information on each script, please read the documentation at the top
of each script.