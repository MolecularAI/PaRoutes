""" Module with script to download public data
"""
import os
import sys
from pathlib import Path

import requests
import tqdm

FILES_TO_DOWNLOAD = {}
FILES_TO_DOWNLOAD["V1"] = [
    {
        "filename": "chembl_10k_route_distance_model.ckpt",
        "url": "https://zenodo.org/record/4925903/files/chembl_10k_route_distance_model.ckpt?download=1",
    },
    {
        "filename": "n1-routes.json",
        "url": "https://zenodo.org/record/6275421/files/n1-routes.json?download=1",
    },
    {
        "filename": "n1-targets.txt",
        "url": "https://zenodo.org/record/6275421/files/n1-targets.txt?download=1",
    },
    {
        "filename": "n1-stock.txt",
        "url": "https://zenodo.org/record/6275421/files/n1-stock.txt?download=1",
    },
    {
        "filename": "n5-routes.json",
        "url": "https://zenodo.org/record/6275421/files/n5-routes.json?download=1",
    },
    {
        "filename": "n5-targets.txt",
        "url": "https://zenodo.org/record/6275421/files/n5-targets.txt?download=1",
    },
    {
        "filename": "n5-stock.txt",
        "url": "https://zenodo.org/record/6275421/files/n5-stock.txt?download=1",
    },
    {
        "filename": "uspto_raw_template_library.csv",
        "url": "https://zenodo.org/record/6275421/files/uspto_raw_template_library.csv?download=1",
    },
    {
        "filename": "uspto_rxn_n1_raw_template_library.csv",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n1_raw_template_library.csv?download=1",
    },
    {
        "filename": "uspto_rxn_n5_raw_template_library.csv",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n5_raw_template_library.csv?download=1",
    },
    {
        "filename": "uspto_rxn_n1_keras_model.hdf5",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n1_keras_model.hdf5?download=1",
    },
    {
        "filename": "uspto_rxn_n1_unique_templates.hdf5",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n1_unique_templates.hdf5?download=1",
    },
    {
        "filename": "uspto_rxn_n5_keras_model.hdf5",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n5_keras_model.hdf5?download=1",
    },
    {
        "filename": "uspto_rxn_n5_unique_templates.hdf5",
        "url": "https://zenodo.org/record/6275421/files/uspto_rxn_n5_unique_templates.hdf5?download=1",
    },
    {
        "filename": "150k_routes.json.gz",
        "url": "https://zenodo.org/record/6275421/files/150k_routes.json.gz?download=1",
    },
]


FILES_TO_DOWNLOAD["V2"] = [
    {
        "filename": "chembl_10k_route_distance_model.ckpt",
        "url": "https://zenodo.org/record/4925903/files/chembl_10k_route_distance_model.ckpt?download=1",
    },
    {
        "filename": "n1-routes.json",
        "url": "https://zenodo.org/record/7341155/files/ref_routes_n1.json?download=1",
    },
    {
        "filename": "n1-targets.txt",
        "url": "https://zenodo.org/record/7341155/files/targets_n1.txt?download=1",
    },
    {
        "filename": "n1-stock.txt",
        "url": "https://zenodo.org/record/7341155/files/stock_n1.txt?download=1",
    },
    {
        "filename": "n5-routes.json",
        "url": "https://zenodo.org/record/7341155/files/ref_routes_n5.json?download=1",
    },
    {
        "filename": "n5-targets.txt",
        "url": "https://zenodo.org/record/7341155/files/targets_n5.txt?download=1",
    },
    {
        "filename": "n5-stock.txt",
        "url": "https://zenodo.org/record/7341155/files/stock_n5.txt?download=1",
    },
    {
        "filename": "uspto_template_library.csv",
        "url": "https://zenodo.org/record/7341155/files/uspto_template_library.csv?download=1",
    },
    {
        "filename": "selected_reactions_all.csv",
        "url": "https://zenodo.org/record/7341155/files/selected_reactions_all.csv?download=1",
    },
    {
        "filename": "uspto_keras_model.hdf5",
        "url": "https://zenodo.org/record/7341155/files/uspto_keras_model.hdf5?download=1",
    },
    {
        "filename": "uspto_unique_templates.csv.gz",
        "url": "https://zenodo.org/record/7341155/files/uspto_unique_templates.csv.gz?download=1",
    },
    {
        "filename": "all_routes.json.gz",
        "url": "https://zenodo.org/record/7341155/files/all_loaded_routes.json.gz?download=1",
    },
]


def _download_file(url: str, filename: str) -> None:
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        pbar = tqdm.tqdm(
            total=total_size, desc=os.path.basename(filename), unit="B", unit_scale=True
        )
        with open(filename, "wb") as fileobj:
            for chunk in response.iter_content(chunk_size=1024):
                fileobj.write(chunk)
                pbar.update(len(chunk))
        pbar.close()


def main() -> None:
    """Entry-point for CLI"""
    version = sys.argv[1] if len(sys.argv) > 1 else "V2"
    base_path = Path(__file__).parent
    for filespec in FILES_TO_DOWNLOAD[version]:
        try:
            _download_file(
                filespec["url"], os.path.join(base_path, filespec["filename"])
            )
        except requests.HTTPError as err:
            print(f"Download failed with message {str(err)}")
            sys.exit(1)


if __name__ == "__main__":
    main()
