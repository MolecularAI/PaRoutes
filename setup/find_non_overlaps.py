"""
Overlap analysis of the extract routes. Is used to create a set
of non-overlapping routes. Two routes are non-overlapping if
the set of intermediates in one route is not among the set of leaves
in the other route, and vice versa.

Example:

  python find_non_overlaps.py --output non_overlapping_routes.pickle

"""
import argparse
import pickle
import random
from typing import List, Dict, Any

from tqdm import tqdm


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool extract non-overlapping routes")
    parser.add_argument(
        "--filename",
        default="loaded_routes.pickle",
        help="the output from the analyze_route script",
    )
    parser.add_argument(
        "--output", required=True, help="the extracted non-overlapping routes"
    )
    parser.add_argument("--seed", type=int, default=1689, help="random seed")
    parser.add_argument(
        "--routes-per-patent",
        type=int,
        default=1,
        help="how many routes to consider per patent",
    )
    return parser.parse_args()


def _find_non_overlappings(
    routes: List[Dict[str, Any]],
    show_progress: bool,
    max_routes: int = None,
):
    taken = set()
    non_overlapping_routes = []
    routes2 = tqdm(routes) if show_progress else routes
    for route in routes2:
        if route["id"] in taken:
            continue
        taken.add(route["id"])

        found_overlap = False
        for route2 in routes:
            if route2["id"] in taken:
                continue

            intersect1 = route["leaves"].intersection(route2["intermediates"])
            intersect2 = route2["leaves"].intersection(route["intermediates"])
            same_root = route["root"] == route2["root"]
            root1_in_other = route["root"] in route2["intermediates"]
            root2_in_other = route2["root"] in route["intermediates"]
            if (
                same_root
                or root1_in_other
                or root2_in_other
                or (intersect1 or intersect2)
            ):
                found_overlap = True
                taken.add(route2["id"])
                break

        if not found_overlap:
            non_overlapping_routes.append(route)

    # Keeping only one route per root
    routes_by_root = {route["root"]: route for route in non_overlapping_routes}
    non_overlapping_routes = list(routes_by_root.values())
    if len(non_overlapping_routes) > 0:
        if max_routes:
            return random.sample(
                non_overlapping_routes, k=min(max_routes, len(non_overlapping_routes))
            )
        else:
            return non_overlapping_routes
    return random.sample(routes, k=1)


def main() -> None:
    args = _get_args()

    random.seed(args.seed)

    with open(args.filename, "rb") as fileobj:
        routes_by_patent = pickle.load(fileobj)

    non_overlapping_routes_within_patents = []
    print(f"Extracting {args.routes_per_patent} routes per patent", flush=True)
    for routes in tqdm(routes_by_patent.values()):
        non_overlapping_routes_within_patents.extend(
            _find_non_overlappings(
                routes, show_progress=False, max_routes=args.routes_per_patent
            )
        )
    print(
        f"Number of non-overlapping routes within patents: {len(non_overlapping_routes_within_patents)}",
        flush=True,
    )

    non_overlapping_routes = _find_non_overlappings(
        non_overlapping_routes_within_patents, show_progress=True
    )
    print(f"Number of non-overlapping routes: {len(non_overlapping_routes)}")

    with open(args.output, "wb") as fileobj:
        pickle.dump(non_overlapping_routes, fileobj)


if __name__ == "__main__":
    main()
