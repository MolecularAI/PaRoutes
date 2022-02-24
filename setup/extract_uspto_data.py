"""
Extract and re-format the USPTO template library data to be
parsable by the route extraction script. The reactions are
grouped by patent.

The script will remove all rows with template occurrence less than a limit + 1,
which leaves enough data to train a model on reactions from the reference routes.
A random sample of three reactions for each unique template is then left out from the
reaction list and saved to a separate CSV.

Example:

  python extract_uspto_data.py --filename ../data/uspto_raw_template_library.csv

"""
import argparse
import pickle
from collections import defaultdict
from typing import Tuple, Dict

import pandas as pd
from tqdm import tqdm


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Tool to extract USPTO data for route extraction")
    parser.add_argument(
        "--filename", required=True, help="the path to the raw USPTO dataset"
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
    patent_id = row.ID.split(";")[0]
    # Keep only atom-mapped reactants
    reactant_smi = ".".join(smi for smi in row.reactants.split(".") if ":" in smi)

    data = {
        "smiles": reactant_smi + ">>" + row.products,
        "reaction_hash": row.reaction_hash,
        "ID": row.ID,
        "retro_template": row.retro_template,
        "template_hash": row.template_hash,
    }
    return patent_id, data


def main() -> None:
    args = _get_args()

    template_data = pd.read_csv(args.filename)
    print(f"Number of original reactions: {len(template_data)}")

    template_group = template_data.groupby("template_hash")
    template_group = template_group.size().sort_values(ascending=False)
    min_index = template_group[template_group >= args.min_occurrence + 1].index
    template_data = template_data[template_data["template_hash"].isin(min_index)]
    print(
        f"Number of reactions with template occurrence >= {args.min_occurrence+1}: {len(template_data)}"
    )

    template_group = template_data.groupby("template_hash")
    leave_out_data = template_group.sample(n=3, random_state=args.seed)
    leave_out_data.to_csv(args.filename.replace(".csv", "_left_out.csv"), index=False)
    print(f"Number of reactions left out: {len(leave_out_data)}")

    template_data = template_data[
        ~template_data["reaction_hash"].isin(leave_out_data["reaction_hash"])
    ]
    print(f"Number of original reactions for route extraction: {len(template_data)}")

    uspto_data = defaultdict(list)
    for _, row in tqdm(template_data.iterrows(), total=len(template_data)):
        key, data = _extract_one_row(row)
        uspto_data[key].append({"data": data})

    print(f"Number of patents: {len(uspto_data)}")
    with open(args.output, "wb") as fileobj:
        pickle.dump(uspto_data, fileobj)


if __name__ == "__main__":
    main()
