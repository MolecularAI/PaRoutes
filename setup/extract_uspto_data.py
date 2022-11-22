"""
Extract and re-format the USPTO template library data to be
parsable by the route extraction script. The reactions are
grouped by patent.

The script will remove all rows with template occurrence less than a limit + 1,
which leaves enough data to train a model on reactions from the reference routes.
A random sample of three row for each unique template is then left out from template library.
All selected reaction from USPTO that correspond to the rest of the template library are then
processed to create a grouping by patent.

Update in v2.0 - The template library example reactions are extended with all selected reactions with the same reaction hash

Example:

  python extract_uspto_data.py --template_library ../data/uspto_template_library.csv --all_reactions ../data/uspto_selected_reactions_all.csv

"""
import argparse
import pickle
from collections import defaultdict
from typing import Tuple, Dict

import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm import tqdm


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool to extract USPTO data for route extraction")
    parser.add_argument(
        "--template_library",
        required=True,
        help="the path to the USPTO template library",
    )
    parser.add_argument(
        "--all_reactions", required=True, help="the path to the raw USPTO dataset"
    )
    parser.add_argument(
        "--output",
        default="uspto_data.pickle",
        help="extracted USPTO reaction grouped by patent",
    )
    parser.add_argument(
        "--min-occurrence",
        type=int,
        default=3,
        help="the minimum template occurance to keep for training",
    )
    parser.add_argument("--seed", type=int, default=1689, help="a random seed")
    return parser.parse_args()


def _extract_one_row(row: pd.Series) -> Tuple[str, Dict[str, str]]:
    patent_id = row["id"].split(";")[0]
    data = {
        "smiles": ">>".join(row.rsmi_processed.split(">")[::2]),
        "rsmi": row.rsmi_processed,
        "reaction_hash": row.PseudoHash,
        "ID": row["id"],
        "RingBreaker": row.RingBreaker,
    }
    return patent_id, data


def _split_template_group(group_data, train_indices, routes_indices, seed):
    train_arr, route_arr = train_test_split(
        group_data.index, train_size=3, random_state=seed, shuffle=True
    )
    train_indices.extend(train_arr)
    routes_indices.extend(route_arr)


def main() -> None:
    args = _get_args()

    template_data = pd.read_csv(
        args.template_library, sep="\t", usecols=["reaction_hash", "template_hash"]
    )
    print(f"Number of reactions in template library: {len(template_data)}")

    template_group = template_data.groupby("template_hash")
    template_group = template_group.size().sort_values(ascending=False)
    min_index = template_group[template_group >= args.min_occurrence + 1].index
    template_data = template_data[
        template_data["template_hash"].isin(min_index)
    ].reset_index()
    print(
        f"Number of reactions with template occurrence >= {args.min_occurrence+1}: {len(template_data)}"
    )

    train_indices = []
    routes_indices = []
    template_group = template_data.groupby("template_hash").apply(
        _split_template_group,
        train_indices=train_indices,
        routes_indices=routes_indices,
        seed=args.seed,
    )
    print(f"Number of templates left out: {len(train_indices)}")
    print(f"Number of templates kept for route extraction: {len(routes_indices)}")

    reaction_data = pd.read_csv(
        args.all_reactions,
        sep="\t",
        usecols=["id", "classification", "rsmi_processed", "PseudoHash", "RingBreaker"],
    )
    reaction_data = reaction_data[
        reaction_data["PseudoHash"].isin(
            template_data.loc[routes_indices, "reaction_hash"]
        )
    ]
    print(f"Number of original reactions for route extraction: {len(reaction_data)}")

    uspto_data = defaultdict(list)
    for _, row in tqdm(reaction_data.iterrows(), total=len(reaction_data)):
        key, data = _extract_one_row(row)
        uspto_data[key].append({"data": data})

    print(f"Number of patents: {len(uspto_data)}")
    with open(args.output, "wb") as fileobj:
        pickle.dump(uspto_data, fileobj)


if __name__ == "__main__":
    main()
