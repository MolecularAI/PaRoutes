# Analysis of route predictions

Provided that you have run a route prediction on either the **n1** or **n5** set,
this folder provide scripts and instructions for how to analyze those predictions.

In order to analyze the route predictions, one need to convert them into a suitable format. This is explained below. 

## Route format

The routes or reaction tree are input as dictionaries that can be loaded from e.g. a JSON file. Multiple-routes that needs be read by the analysis scripts are supplied as a list of dictionaries.

The structure of the dictionary is checked by the code, and if the validation fails, an exception will be raised.

The input structure should be a the type `MoleculeNode` which is a dictionary with the following fields

* **smiles** - the SMILES string of the molecule represented by the node
* **type** - should be the string "mol"
* **in_stock** - an *optional* bolean to indicate if a molecule is in stock. 
* **children** - an *optional* list, containing at most one dictionary of type `ReactionNode`

**in_stock** is necessary for accurate route scoring and analysis, but is not required by the calculations.

The `ReactionNode` type is a dictionary with the following fields

* **type** - should be the string "reaction"
* **children** - a list of dictionaries of type `MoleculeNode`

It is easy to realize that this is a recursive definition. All extra fields in the dictionaries are ignored.

And example dictionary is shown below


    {
        "smiles": "CCCCOc1ccc(CC(=O)N(C)O)cc1",
        "type": "mol",
        "in_stock": false
        "children": [
            {
                "type": "reaction",
                "children": [
                    {
                        "smiles": "CCCCOc1ccc(CC(=O)Cl)cc1",
                        "type": "mol",
                        "in_stock": true
                    },
                    {
                        "smiles": "CNO",
                        "type": "mol",
                        "in_stock": true
                    }
                ]
            }
        ]
    }


## Analysis scripts

This folder contains two scripts for analyses:

* `route_quality.py` - computes mainly the top-n accuracies but also a range of other metrics
* `route_clusters.py` - perform pairwise distance calculations and cluster the routes

### Route quality analysis

Example usage:

    python route_quality.py --routes ~/output_routes.json --references ../data/n1-routes.json --output ~/route_analyses.csv

where `~/output_routes.json` is the correctly formated route predictions for the **n1** set.

The script will print out the key metrics but the created CSV file can be further analysed.

The script will do some analysis of the reference route, and then it will do analysis for the routes up to rank 1, 5 and 10 (by default). The columns `best-1`, `best-5`, `best-10` holds the minimum distance between the reference route and the predicted routes up to rank 1, 5 and 10. To compute, e.g. top-1 accuracy one has to check how many of the `best-1` is equal to 0. See this python snippet:

    import pandas as pd
    data = pd.read_csv("route_analyses.csv")
    print("top-1: ", (data["best-1"]==0).mean())

For more info on the columns in the CSV-file have a look at the documentation at the top of the Python script.

### Route diversity analysis

Example usage:

    python route_clusters.py --routes ~/output_routes.json --model ../data/chembl_10k_route_distance_model.ckpt --nclusters 0 --min_density 2 --output ~/cluster_analyses.json

The output is a JSON file, which holds a list of dictionaries, one for each target. This JSON file can easily be loaded as a pandas Dataframe. For instance, 
to compute the mean number of clusters this small snippet can be used

    import pandas as pd
    data = pd.read_json("cluster_analyses.json")
    nclusters = data.apply(
        lambda row: max(row.cluster_labels)+1 if row.cluster_labels else 0, axis=1
    )
    print("Average number of clusters: ", nclusters.mean())

For more info on the content of the JSON-file, have a look at the documentation at the top of the Python script.