import pathlib
import pandas as pd
import argparse
import re
import plotly.express as px 


def get_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser()
  parser.add_argument("path", help="Path to the csv")

  return parser



def main():
    parser = get_parser()
    args = parser.parse_args()
    path = args.path
    df = pd.read_csv(path)

    _df = df.set_index(["application", "task", "machine", "conf"])["duration_w_stress"].unstack()
    _df = _df.divide(_df["real"], axis=0)
    _df = _df.reset_index()
    _df = _df.melt(id_vars=["application", "task", "machine"], value_vars=_df.columns, var_name="conf", value_name="duration_w_stress_ratio")
    print(_df)
    
    


if __name__ =="__main__":
    main()