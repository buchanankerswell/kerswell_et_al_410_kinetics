#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_arguments() -> Namespace:
    """"""
    parser = ArgumentParser(description="Write depth profile summary markdown table.")
    parser.add_argument("--in-dir", type=str, help="Output directory")

    return parser.parse_args()


#######################################################
## .1. Main                                      !!! ##
#######################################################
def main():
    """"""
    args = parse_arguments()

    in_dir = Path(args.in_dir)
    in_path = in_dir / "depth-profile-data.csv"
    out_path = in_dir / "depth-profile-summary-table.md"

    print("    --------------------------------------------------")
    print("    Writing depth profile summary md table")
    print("    --------------------------------------------------")

    if out_path.exists():
        print(f" -- Found table: {out_path.name}!")
        return

    print(f" -> {out_path.name}")

    df = pd.read_csv(in_path)
    df_100 = cast(pd.DataFrame, df[df["timestep"] == 100])
    df_100 = df_100.drop_duplicates(subset=["model_id", "timestep"], keep="first")

    df_100[["model_id", "Z"]] = df_100["model_id"].str.extract(r"^(plume|slab)_(\d+\.\d+e[+-]?\d+)$").astype(str)
    df_100["Z"] = df_100["Z"].apply(lambda x: f"`{x}`")
    df_100["displacement"] = (df_100["displacement"] / 1e3).round(0)
    df_100["width"] = (df_100["width"] / 1e3).round(0)
    df_100["max_velocity"] = df_100["max_velocity"].round(2)
    df_100["max_reaction_rate"] = np.log10(df_100["max_reaction_rate"]).round(2)

    cols = ["model_id", "Z", "displacement", "width", "max_velocity", "max_reaction_rate"]
    df_100 = df_100[cols]

    df_100 = df_100.rename(  # pyright: ignore
        columns={
            "model_id": "Model",
            "Z": "$Z$",
            "displacement": "410 Displacement",
            "width": "410 Width",
            "max_velocity": "Max $\\vec{u}_y$",
            "max_reaction_rate": "Log$_{10}$ Max $\\dot{X}$",
        }
    )

    markdown_table: str = df_100.to_markdown(index=False) or ""

    with open(out_path, "w") as f:
        f.write(markdown_table)


if __name__ == "__main__":
    main()
