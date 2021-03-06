<img src="assets/logo.png" width="100">

PaRoutes is a framework for benchmarking multi-step retrosynthesis methods,
i.e. route predictions. 

It provides:

* A curated reaction dataset for building one-step retrosynthesis models
* Two sets of 10,000 routes
* Two sets of stock molecules to use as stop-criterion for the search
* Scripts to compute route quality and route diversity metrics

## Prerequisites

Before you begin, ensure you have met the following requirements:

* Linux, Windows or macOS platforms are supported - as long as the dependencies are supported on these platforms.

* You have installed [anaconda](https://www.anaconda.com/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) with python 3.7 - 3.9

The tool has been developed on a Linux platform.

## Installation

First clone the repository using Git.

Then execute the following commands in the root of the repository 

    conda env create -f env.yml
    conda activate paroutes-env
    python data/download_data.py

Now all the dependencies and datasets are setup.


## Usage

### Performing route predictions

PaRoutes provide a list of targets and stock molecules in SMILES format
for two sets **n1** and **n5**.

For **n1** you find in the `data/` folder of the repository
* `n1-targets.txt` - the target molecules
* `n1-stock.txt` - the stock molecules

For **n5** you find in the `data/` folder of the repository
* `n5-targets.txt` - the target molecules
* `n5-stock.txt` - the stock molecules

For more information on the files in the `data/` folder, please read the README file for that folder.

### Analysing predictions

The predicted route exported by your software need to be converted to a format
that can be read by the analysis tool. This format is outlined in the `analysis\README.md`

The following command for analysis assumes:

1. The current directory is the root of the `paroutes` repo
2. Your route predictions for the **n1** targets in a JSON format is located at `~/output_routes.json`

Then you can type

    python analysis/route_quality.py --routes ~/output_routes.json --references data/n1-routes.json --output ~/route_analyses.csv


to calculate the route quality metrics. It will print out how many of the targets were solved and the top-1, top-5 and top-10 accuracies (by default). For further details have a look in the `data/README.md` file.

To perform clustering on the same dataset, you can type

    python analysis/route_clusters.py --routes ~/output_routes.json --model data/chembl_10k_route_distance_model.ckpt --min_density 2 --output ~/cluster_analyses.json

The script will print out the average number of clusters formed for each target. For further details have a look in the `data/README.md` file. 


## Contributing

We welcome contributions, in the form of issues or pull requests.

If you have a question or want to report a bug, please submit an issue.


To contribute with code to the project, follow these steps:

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Make your changes and commit them: `git commit -m '<commit_message>'`
4. Push to the remote branch: `git push`
5. Create the pull request.

Please use ``black`` package for formatting, and follow ``pep8`` style guide.


## Contributors

* [@SGenheden](https://www.github.com/SGenheden)
* [@EBjerrum](https://www.github.com/EBjerrum)

Yasmine Nahal is acknowledged for the creation of the PaRoutes logo.

The contributors have limited time for support questions, but please do not hesitate to submit an issue (see above).

## License

The software is licensed under the Apache 2.0 license (see LICENSE file), and is free and provided as-is.

## References