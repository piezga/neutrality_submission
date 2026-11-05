"""
This script takes the raw data from https://datadryad.org/dataset/doi:10.15146/5xcp-0d46
and converts it to a .csv file with only spatial coordinates + name
"""

from pathlib import Path
import rdata
import pandas as pd

for file in Path("data/raw").glob("*.rdata"):
    data = rdata.read_rda(file)

    for name, obj in data.items():
        df = pd.DataFrame(obj)
        df = df[df["status"] == "A"]                    # alive stems only
        df = df[["gx", "gy", "sp", "dbh"]].copy()
        df.columns = ["x", "y", "name", "dbh"]
        df = df.dropna(subset=["x","y"])
        out = file.with_name(f"{file.stem}.csv")
        df.to_csv(out, index=False)

        print(f"{file} -> {out}")

