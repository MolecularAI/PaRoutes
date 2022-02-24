"""
Subset the reaction data for modelling to only include templates that
are included in the reference routes, but not the reactions included
in the reference routes.

Example:

  python extract_training_data.py --filename ../data/uspto_raw_template_library.csv --routes ref_routes.json  --output_tmpl uspto_tmpl_raw_template_library.csv --output_rxn uspto_rxn_raw_template_library.csv

"""
import argparse
import json
from typing import Dict, Any, Set

import pandas as pd


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool extract training data")
    parser.add_argument("--filename", required=True, help="the raw USPTO dataset")
    parser.add_argument(
        "--routes",
        required=True,
        help="the extracted reference routes created by the select_routes script",
    )
    parser.add_argument(
        "--output_tmpl",
        required=True,
        help="the selected USPTO dataset with templates only in the reference routes",
    )
    parser.add_argument(
        "--output_rxn",
        required=True,
        help="the selected USPTO dataset for any template",
    )
    return parser.parse_args()


def _traverse_route(
    route: Dict[str, Any], template_hashes: Set[str], reaction_hashes: Set[str]
) -> None:
    if route["type"] == "reaction":
        template_hashes.add(route["metadata"]["template_hash"])
        reaction_hashes.add(route["metadata"]["reaction_hash"])
    for child in route.get("children", []):
        _traverse_route(child, template_hashes, reaction_hashes)


def main() -> None:
    args = _get_args()

    template_data = pd.read_csv(args.filename)
    noriginal_reactions = len(template_data)
    noriginal_templates = len(set(template_data.template_hash))
    print(f"Number of original reactions: {noriginal_reactions}")
    print(f"Number of original unique templates: {noriginal_templates}")

    with open(args.routes, "r") as fileobj:
        routes = json.load(fileobj)

    template_hashes = set()
    reaction_hashes = set()
    for route in routes:
        _traverse_route(route, template_hashes, reaction_hashes)
    print(
        f"Number of unique reactions in routes: {len(reaction_hashes)} ({len(reaction_hashes)/noriginal_reactions*100:.2f}%)"
    )
    print(
        f"Number of unique templates in routes: {len(template_hashes)} ({len(template_hashes)/noriginal_templates*100:.2f}%)"
    )

    sel1 = ~template_data["reaction_hash"].isin(reaction_hashes)
    sel2 = template_data["template_hash"].isin(template_hashes)
    template_data_sel = template_data[sel1 & sel2]
    print(
        f"Selecting {len(template_data_sel)} ({len(template_data_sel)/noriginal_reactions*100:.2f}%) rows for training"
    )
    nmin = template_data_sel.groupby("template_hash").count()["ID"].min()
    print(f"The minimum number of reactions per template: {nmin}")

    assert len(set(template_data_sel["template_hash"])) == len(template_hashes)
    template_data_sel.to_csv(args.output_tmpl, index=False)

    template_data_sel = template_data[sel1]
    print(
        f"Selecting {len(template_data_sel)} ({len(template_data_sel)/noriginal_reactions*100:.2f}%) rows for training"
    )
    template_data_sel.to_csv(args.output_rxn, index=False)


if __name__ == "__main__":
    main()
