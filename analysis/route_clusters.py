""" 
Module containing CLI tool calculate route distances and do clustering 

Example:
    python route_clusters.py --routes output_trees.json.gz --model chembl_10k_route_distance_model.ckpt --min_density 2 --output cluster_analysis.json

The output JSON-file will contain a list of dictionaries, one for each target. The dictionary will
have these files:
    * distance_matrix - the pair-wise route distances
    * distances_time - the time to compute the distances
    * cluster_labels - the cluster label for each target
    * cluster_time - the time to compute the clusters

The JSON-file can easily be read into a pandas Dataframe:

    import pandas as pd
    data = pd.read_json("cluster_analysis.json")

"""
from __future__ import annotations
import json
import gzip
import argparse
import time
import math
from typing import List, Dict, Any, Optional, Sequence

from tqdm import tqdm

import route_distances.lstm.defaults as defaults
from route_distances.route_distances import route_distances_calculator
from route_distances.clustering import ClusteringHelper
from route_distances.utils.type_utils import RouteDistancesCalculator


def _get_args(optional_args: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        "Tool to calculate pairwise distances between routes and cluster them"
    )
    parser.add_argument(
        "--routes", nargs="+", required=True, help="the routes in JSON format"
    )
    parser.add_argument(
        "--fp_size",
        type=int,
        default=defaults.FP_SIZE,
        help="the fingerprint size for the LSTM model",
    )
    parser.add_argument(
        "--lstm_size", type=int, default=defaults.LSTM_SIZE, help="the LSTM layer size"
    )
    parser.add_argument(
        "--model",
        required=True,
        help="specify 'ted' to cluster with TED or a filename to an LSTM model",
    )
    parser.add_argument(
        "--nclusters",
        type=int,
        default=0,
        help="number of clusters. provide 0 to opimize.",
    )
    parser.add_argument(
        "--min_density",
        type=int,
        default=None,
        help="the if optimizing, specifies the minimum cluster density",
    )
    parser.add_argument(
        "--output", default="route_clustering.json", help="the output file"
    )
    return parser.parse_args(optional_args)


def _calc_distances(
    routes: List[Dict[str, Any]], calculator: RouteDistancesCalculator
) -> Dict[str, Any]:
    if len(routes) == 1:
        return {"distance_matrix": [[0.0]], "distances_time": 0}

    time0 = time.perf_counter_ns()
    distances = calculator(routes)
    dict_ = {
        "distance_matrix": distances.tolist(),
        "distances_time": (time.perf_counter_ns() - time0) * 1e-9,
    }
    return dict_


def _do_clustering(
    nroutes: int,
    distances: List[List[float]],
    nclusters: int,
    min_density: int = None,
) -> Dict[str, Any]:
    if distances == [[0.0]] or nroutes < 3:
        return {"cluster_labels": [], "cluster_time": 0}

    if min_density is None:
        max_clusters = min(nroutes, 10)
    else:
        max_clusters = int(math.ceil(nroutes / min_density))

    time0 = time.perf_counter_ns()
    labels = ClusteringHelper.cluster(
        distances, nclusters, max_clusters=max_clusters
    ).tolist()
    cluster_time = (time.perf_counter_ns() - time0) * 1e-9
    return {"cluster_labels": labels, "cluster_time": cluster_time}


def main(optional_args: Optional[Sequence[str]] = None) -> None:
    """ Entry-point for CLI tool """
    args = _get_args(optional_args)

    routes_list = []
    for filename in args.routes:
        if filename.endswith(".gz"):
            with gzip.open(filename, "rt", encoding="UTF-8") as fileobj:
                routes = json.load(fileobj)
        else:
            with open(filename) as fileobj:
                routes = json.load(fileobj)
        routes_list.extend(routes)

    if args.model == "ted":
        calculator = route_distances_calculator("ted", content="both")
    else:
        calculator = route_distances_calculator(
            "lstm",
            model_path=args.model,
            fp_size=args.fp_size,
            lstm_size=args.lstm_size,
        )

    dist_data = []
    for routes in tqdm(routes_list):
        dist_data.append(_calc_distances(routes=routes, calculator=calculator))

    if args.nclusters is not None:
        nroutes_sum = 0
        for routes, data in zip(tqdm(routes_list), dist_data):
            cluster_data = _do_clustering(
                nroutes=len(routes),
                distances=data["distance_matrix"],
                nclusters=args.nclusters,
                min_density=args.min_density,
            )
            if cluster_data["cluster_labels"]:
                nroutes_sum += max(cluster_data["cluster_labels"]) + 1
            data.update(cluster_data)

    print(f"Average number of clusters: {nroutes_sum/len(dist_data):.2f}")

    with open(args.output, "w") as fileobj:
        json.dump(dist_data, fileobj, indent=4)


if __name__ == "__main__":
    main()
