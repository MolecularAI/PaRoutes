""" 
Module containing CLI tool calculate route quality metrics 

Example:
    python route_quality.py --routes output_trees.json.gz --references n1-routes.json --output route_analyses.csv

The output CSV file will have the following columns:
    * ref lrr - the longer linear route (LRR) of the reference route
    * ref nleaves - the number of leaves in the reference route
    * solved target - if at least one predicted route has all starting material in stock
    * max llr-1 - the maximumum LRR of the routes ranked at 1
    * min lrr-1 - the minimum LRR of the routes ranked at 1
    * min nleaves-1 - the minimum number of leaves in the routes ranked at 1
    * max leaves overlap-1 - the maximum leaf overlap of the routes ranked at 1
    * mean solved-1 - the average number of solved routes of the routes ranked at 1
    * best-1 - the minimum route distance between the reference at the routes ranked at 1

the same quanities are computed for the other rankes provided by the `--ks`` command argument,
that is 1, 5 and 10, by default.
"""
import argparse
import gzip
import json
from typing import Dict, Any, Set, List, Tuple, Optional, Sequence

from tqdm import tqdm
import numpy as np
import pandas as pd

from route_distances.ted.reactiontree import ReactionTreeWrapper
from route_distances.utils.routes import (
    calc_llr,
    extract_leaves,
    is_solved,
    route_scorer,
    route_ranks,
)


def _analyze_routes(
    routes: List[Dict[str, Any]],
    reference: Dict[str, Any],
    ks: List[int],
) -> Dict[str, float]:
    max_k = max(ks)
    routes, scores = route_scorer(routes)
    ranks = route_ranks(scores)

    ref_route = ReactionTreeWrapper(reference, content="both")
    ref_lrr = calc_llr(reference)
    ref_leaves = extract_leaves(reference)
    nref_leaves = len(ref_leaves)

    llr = []
    solved = []
    leaves_overlap = []
    nleafs = []
    true_rank = None
    distances = []
    solved_target = False
    for rank, tree_dict in zip(ranks, routes):
        this_is_solved = is_solved(tree_dict)
        if this_is_solved:
            solved_target = True
        if rank > max_k and (true_rank or not this_is_solved):
            break

        solved.append(int(this_is_solved))
        llr.append(calc_llr(tree_dict))
        this_leaves = extract_leaves(tree_dict)
        nleafs.append(len(this_leaves))
        leaves_overlap.append(
            len(ref_leaves.intersection(this_leaves))
            / ((nref_leaves + len(this_leaves)) / 2)
        )

        try:
            other_route = ReactionTreeWrapper(tree_dict, content="both")
        except ValueError:
            distances.append(1e6)
        else:
            dist = ref_route.distance_to(other_route)
            distances.append(dist)
            if dist < 1e-6:
                true_rank = rank

    stats = {}
    stats["ref lrr"] = ref_lrr
    stats["ref nleaves"] = len(ref_leaves)
    stats["solved target"] = solved_target
    for this_k in ks:
        if this_k + 1 in ranks:
            end = ranks.index(this_k + 1)
        else:
            end = len(ranks)
        stats[f"max llr-{this_k}"] = max(llr[:end])
        stats[f"min llr-{this_k}"] = min(llr[:end])
        stats[f"min nleaves-{this_k}"] = min(nleafs[:end])
        stats[f"max leaves overlap-{this_k}"] = max(leaves_overlap[:end])
        stats[f"mean solved-{this_k}"] = np.mean(solved[:end])
        stats[f"best-{this_k}"] = min(distances[:end])
    stats["true_rank"] = true_rank
    return stats


def _get_args(optional_args: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool to analyse route predictions")
    parser.add_argument(
        "--routes", nargs="+", required=True, help="the routes in JSON format"
    )
    parser.add_argument("--references", required=True, help="the reference routes")
    parser.add_argument(
        "--ks",
        type=int,
        nargs="+",
        default=[1, 5, 10],
        help="the ranks at which to calculate metrics",
    )
    parser.add_argument(
        "--output",
        default="route_analyses.csv",
        help="the filename of the output of the analysis",
    )
    return parser.parse_args(optional_args)


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

    with open(args.references, "r") as fileobj:
        references = json.load(fileobj)

    stats = []
    for routes, reference in zip(routes_list, tqdm(references)):
        stats.append(_analyze_routes(routes, reference, args.ks))
        if len(stats) == 10:
            break
    stats = pd.DataFrame(stats)

    print(f"Number of solved targets: {stats['solved target'].sum()}")
    for rank in args.ks:
        topn = (stats[f"best-{rank}"] == 0).mean()
        print(f"top-{rank}: {topn:.2f}")

    stats.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()