"""
Extract routes from grouped reaction data using a method of Mo et al. (2021)

The extracted routes will be saved to a pickled dictionary where the key
is the patent ID and the value is a list of routes.

Example:

  python extract_routes.py --max-workers 2

"""
import argparse
import pickle
from concurrent.futures import ProcessPoolExecutor
from typing import Tuple, List, Dict, Any

from external.pathway_extraction import extract_one_patent


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        "Tool to extract routes from grouped reaction data"
    )
    parser.add_argument(
        "--filename",
        default="uspto_data.pickle",
        help="the output from the extract_uspto_data script",
    )
    parser.add_argument(
        "--output", default="uspto_routes.pickle", help="the extracted routes"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="the maximum number of parallel workers",
    )
    return parser.parse_args()


def _extract_pathways_wrapper(args: Tuple[str, List[Dict[str, str]]]) -> Dict[str, Any]:
    patent_id, records = args
    return extract_one_patent(records, patent_id)


def main() -> None:
    args = _get_args()

    with open(args.filename, "rb") as f:
        data = pickle.load(f)

    data_list = list(data.items())
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        extracted_pathways = [
            output for output in executor.map(_extract_pathways_wrapper, data_list)
        ]

    pathway_dict = {}
    for pathway in extracted_pathways:
        pathway_dict[pathway["patentID"]] = pathway["trees"]

    with open(args.output, "wb") as fileobj:
        pickle.dump(pathway_dict, fileobj)


if __name__ == "__main__":
    main()
