"""
Select a diverse set of routes that are not deeper than a maximum threshold.
Diversity calculation re-implemented from diversipy package.

Example:
    python select_routes.py --model ../data/chembl_10k_route_distance_model.ckpt


The input and output is a pickled list of an internal route structure format that
is common to many of the setup scripts. 

The script will output the selected routes in JSON-format, the selected targets in
a text file with SMILES strings and the stock compounds in a textfile with SMILES strings
or InChI-keys.
"""
import argparse
import pickle
import json
from typing import List, Dict, Any, Set

import numpy as np
import pandas as pd

from route_distances.route_distances import route_distances_calculator
from route_distances.clustering import ClusteringHelper
from route_distances.utils.routes import extract_leaves


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool extract non-overlapping routes")
    parser.add_argument(
        "--filename",
        default="non_overlapping_routes.pickle",
        help="the path to an output file generated by the 'find_non_overlaps.py' script",
    )
    parser.add_argument(
        "--output",
        default="selected_routes.pickle",
        help="the output filename",
    )
    parser.add_argument(
        "--model", required=True, help="the path to a LSTM route distance model"
    )
    parser.add_argument(
        "--stock", default="stock.txt", help="the filename of the stock file"
    )
    parser.add_argument(
        "--targets", default="targets.txt", help="the filename of the targets file"
    )
    parser.add_argument(
        "--ref-routes",
        default="ref_routes.json",
        help="the filename of the file with the reference routes",
    )
    parser.add_argument(
        "--size", type=int, default=10000, help="the number of routes to select"
    )
    parser.add_argument(
        "--max-reaction", type=int, default=10, help="the maximum depth of a route"
    )
    parser.add_argument(
        "--stock_kind",
        choices=["smi", "inchi"],
        default="smi",
        help="the format of the stock, either SMILES strings or InChI keys",
    )
    return parser.parse_args()


def _select_routes_greedy_maxmin(distances: np.ndarray, size: int) -> List[int]:
    aggregated_dist_criteria = distances[0, :]
    previous_index = np.argmax(aggregated_dist_criteria)
    selected_indices = [previous_index]
    while len(selected_indices) < size:
        sel_distances = distances[previous_index, :]
        aggregated_dist_criteria = np.minimum(
            aggregated_dist_criteria, sel_distances.ravel()
        )
        previous_index = np.argmax(aggregated_dist_criteria)
        selected_indices.append(previous_index)
    return selected_indices


def main() -> None:
    args = _get_args()

    with open(args.filename, "rb") as fileobj:
        routes = pickle.load(fileobj)
    print(f"Read {len(routes)}  routes in total")

    routes = [route for route in routes if route["nreactions"] <= args.max_reaction]
    print(f"Keeping {len(routes)} with at most {args.max_reaction} reactions")

    tree_dicts = [route["rt"] for route in routes]
    calculator = route_distances_calculator(model="lstm", model_path=args.model)
    distances = calculator(tree_dicts)
    indices = _select_routes_greedy_maxmin(distances, args.size)

    selected_routes = [routes[idx] for idx in indices]
    with open(args.output, "wb") as fileobj:
        pickle.dump(selected_routes, fileobj)

    count_data = pd.DataFrame(
        {
            "nleaves": [route["nleaves"] for route in selected_routes],
            "nmols": [route["nmols"] for route in selected_routes],
            "nreactions": [route["nreactions"] for route in selected_routes],
            "llr": [route["llr"] for route in selected_routes],
        }
    )
    print(count_data.describe(include="all"))
    branched_routes = [route["llr"] != route["nreactions"] for route in selected_routes]
    print(f"Number of branched routes: {sum(branched_routes)}")

    leaves = set()
    for tree in selected_routes:
        if args.stock_kind == "inchi":
            leaves = leaves.union(tree["leaves"])
        else:
            leaves = leaves.union(extract_leaves(tree["rt"]))
    with open(args.stock, "w") as fileobj:
        fileobj.write("\n".join(leaves))

    with open(args.targets, "w") as fileobj:
        fileobj.write("\n".join([route["rt"]["smiles"] for route in selected_routes]))

    with open(args.ref_routes, "w") as fileobj:
        json.dump([route["rt"] for route in selected_routes], fileobj, indent=4)


if __name__ == "__main__":
    main()
