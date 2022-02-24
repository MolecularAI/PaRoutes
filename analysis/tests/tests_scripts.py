import os
import json

import pandas as pd
import numpy as np

from ..route_clusters import main as route_clusters_main
from ..route_quality import main as route_quality_main


def test_clustering(tmpdir, shared_datadir, capsys):
    routes_filename = str(shared_datadir / "routes.json")
    output_filename = str(tmpdir / "output.json")

    route_clusters_main(
        [
            "--routes",
            routes_filename,
            "--model",
            "ted",
            "--min_density",
            "2",
            "--output",
            output_filename,
        ]
    )

    assert os.path.exists(output_filename)

    with open(output_filename, "r") as fileobj:
        results = json.load(fileobj)

    assert len(results) == 2
    expected_keys = [
        "distance_matrix",
        "distances_time",
        "cluster_labels",
        "cluster_time",
    ]
    assert all(key in results[0] for key in expected_keys)
    assert results[0]["cluster_labels"] == []
    assert results[1]["cluster_labels"] == [0, 0, 0, 0, 0, 0, 1]

    output = capsys.readouterr()
    assert "Average number of clusters: 1.00" in output.out


def test_quality_analysis(tmpdir, shared_datadir, capsys):
    routes_filename = str(shared_datadir / "routes.json")
    ref_filename = str(shared_datadir / "ref-routes.json")
    output_filename = str(tmpdir / "output.json")

    route_quality_main(
        [
            "--routes",
            routes_filename,
            "--references",
            ref_filename,
            "--output",
            output_filename,
        ]
    )

    assert os.path.exists(output_filename)

    results = pd.read_csv(output_filename)
    assert results["ref lrr"].tolist() == [3, 4]
    assert results["ref nleaves"].tolist() == [4, 4]
    assert results["solved target"].tolist() == [True, True]

    assert results["max llr-1"].tolist() == [3, 2]
    assert results["min llr-1"].tolist() == [3, 2]
    assert results["min nleaves-1"].tolist() == [4, 3]
    assert np.round(results["max leaves overlap-1"], 3).tolist() == [1, 0.857]
    assert results["mean solved-1"].tolist() == [1.0, 1.0]
    assert np.round(results["best-1"], 3).tolist() == [0.0, 7.045]

    assert results["max llr-5"].tolist() == [3, 4]
    assert results["min llr-5"].tolist() == [3, 2]
    assert results["min nleaves-5"].tolist() == [4, 2]
    assert np.round(results["max leaves overlap-5"], 3).tolist() == [1, 1.0]
    assert results["mean solved-5"].tolist() == [1.0, 1.0]
    assert np.round(results["best-5"], 3).tolist() == [0.0, 0.0]

    output = capsys.readouterr()
    assert "Number of solved targets: 2" in output.out
    assert "top-1: 0.50" in output.out
    assert "top-5: 1.00" in output.out
    assert "top-10: 1.00" in output.out
