"""
Analyze and re-format the extracted routes

The script will do the following for each extracted route:

    * Re-format it to comply with analysis scripts
    * Extract InChI keys of leaves and intermediates
    * Extract statistics on molecules, reactions and route length (LLR)

It will remove all routes with just one leaf.

The output will be a pickled dictionary where each key is a patent ID
and the value is a list of loaded and analysed routes.

Example:

  python analyze_routes.py

"""
import argparse
import pickle
from typing import Tuple, List, Dict, Any, Set

import pandas as pd
from rdkit import Chem
from tqdm import tqdm


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool to analyze and re-format extracted routes")
    parser.add_argument(
        "--filename",
        default="uspto_routes.pickle",
        help="the output from the extract_routes script",
    )
    parser.add_argument(
        "--output",
        default="loaded_routes.pickle",
        help="the loaded and analyzed routes",
    )
    return parser.parse_args()


def _calc_depth(tree_dict: Dict[str, Any], depth: int = 0) -> int:
    children = tree_dict.get("children", [])
    if children:
        return max(_calc_depth(child, depth + 1) for child in children)
    return depth


def _extract_nodes(
    tree_dict: Dict[str, Any],
    reactions: List[Dict[str, Any]],
    leaves: Set[str],
    intermediates: Set[str],
):
    if tree_dict["type"] == "mol":
        children = tree_dict.get("children", [])
        if children:
            intermediates.add(tree_dict["smiles"])
            for child in children:
                _extract_nodes(child, reactions, leaves, intermediates)
        else:
            leaves.add(tree_dict["smiles"])

    else:
        reactions.append(tree_dict["metadata"])
        for child in tree_dict["children"]:
            _extract_nodes(child, reactions, leaves, intermediates)


def _extracted_route_to_tree_dict(route: Dict[str, Any]) -> Dict[str, Any]:
    tree_dict = {"smiles": route["smiles"], "type": "mol"}

    children = route.get("child", [])
    tree_dict["in_stock"] = not bool(children)

    if not children:
        return tree_dict

    tree_dict["children"] = [
        {
            "type": "reaction",
            "smiles": "",
            "metadata": route["record_data"]["data"],
            "children": [_extracted_route_to_tree_dict(child) for child in children],
        }
    ]
    return tree_dict


def _load_route(tree_dict: Dict[str, Any], id_str: str) -> Dict[str, Any]:
    rt = _extracted_route_to_tree_dict(tree_dict["tree"])

    reactions = []
    leaves = set()
    intermediates = set()
    _extract_nodes(rt, reactions, leaves, intermediates)

    leaf_inchis = set(Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in leaves)
    inter_inchis = set(
        Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in intermediates
    )
    return {
        "rt": rt,
        "root": Chem.MolToInchiKey(Chem.MolFromSmiles(rt["smiles"])),
        "leaves": leaf_inchis,
        "intermediates": inter_inchis,
        "id": id_str,
        "nreactions": len(reactions),
        "nleaves": len(leaf_inchis),
        "nmols": len(leaf_inchis) + len(inter_inchis),
        "llr": _calc_depth(rt) // 2,
    }


def main() -> None:
    args = _get_args()

    with open(args.filename, "rb") as fileobj:
        routes_by_patent = pickle.load(fileobj)

    sum_routes = sum(len(routes) for routes in routes_by_patent)
    print(f"Total number of reaction groups: {len(routes_by_patent)}")
    print(f"Total number of routes: {sum_routes}")
    print(f"Average number of routes: {sum_routes/len(routes_by_patent):.2f}")
    print(
        f"Maximun number of routes for a group: {min(len(routes) for routes in routes_by_patent)}"
    )
    print(
        f"Minimum number of routes for a group: {max(len(routes) for routes in routes_by_patent)}"
    )
    print(
        f"Number of groups having at least one routes: {sum(bool(routes) for routes in routes_by_patent)}\n"
    )

    print("Loading and transforming extracted routes")
    loaded_routes_by_patent = {}
    for patent_id, routes in tqdm(routes_by_patent.items()):
        if not routes:
            continue
        loaded_routes = []
        for route_idx, route in enumerate(routes):
            loaded_route = _load_route(route, f"{patent_id}@{route_idx}")
            if loaded_route["nleaves"] > 1:
                loaded_routes.append(loaded_route)
        if loaded_routes:
            loaded_routes_by_patent[patent_id] = loaded_routes

    count_data = pd.DataFrame(
        {
            "nleaves": [
                route["nleaves"]
                for routes in loaded_routes_by_patent.values()
                for route in routes
            ],
            "nmols": [
                route["nmols"]
                for routes in loaded_routes_by_patent.values()
                for route in routes
            ],
            "nreactions": [
                route["nreactions"]
                for routes in loaded_routes_by_patent.values()
                for route in routes
            ],
            "llr": [
                route["llr"]
                for routes in loaded_routes_by_patent.values()
                for route in routes
            ],
        }
    )
    print(count_data.describe(include="all"))
    branched_routes = [
        route["llr"] != route["nreactions"]
        for routes in loaded_routes_by_patent.values()
        for route in routes
    ]
    print(f"Number of branched routes: {sum(branched_routes)}")

    with open(args.output, "wb") as fileobj:
        pickle.dump(loaded_routes_by_patent, fileobj)


if __name__ == "__main__":
    main()
