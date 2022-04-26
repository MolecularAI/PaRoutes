import json
import gzip
import sys

import pandas as pd

for filename in sys.argv[1:]:
    df = pd.read_hdf(filename, "table")
    with gzip.open(
        filename.replace(".hdf5", "_trees.json.gz"), "wt", encoding="UTF-8"
    ) as fileobj:
        json.dump(df.trees.values.tolist(), fileobj)
