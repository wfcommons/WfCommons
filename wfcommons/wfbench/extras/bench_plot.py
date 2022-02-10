import pathlib
import pandas as pd
import argparse
import re
import plotly.express as px 


this_dir = pathlib.Path(__file__).resolve().parent

def convert_time(text: str) -> float:
  return sum([float(p) * 60**i for i, p in enumerate(str(text).split(":")[::-1])])

def get_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser()
  parser.add_argument("path", help="Path to the csv")
  parser.add_argument("-m", "--machine", help="Machine used")
  return parser

colors = {
    "1_9": "#e69f00",
    "2_8": "#56b4e9",
    "3_7": "#009e73",
    "4_6": "#f0e442",
    "5_5": "#0072b2",
    "6_4": "#d55e00",
    "7_3": "#cc79a7",
    "8_2": "#000000",
    "9_1": "#aa3377",
    "real":"#332288"
}

symbols = ["circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "star", "hexagon", "pentagon"]
def main():
  parser = get_parser()
  args = parser.parse_args()
  path = pathlib.Path(args.path)
  machine = args.machine
  files = []

  for file in path.glob("*.csv"):
    task, _ = re.match(r"^(.+)_(^.+)?.*?_.+?.*?.*?$", str(file)).groups()
    task = pathlib.Path(task).stem
    df = pd.read_csv(str(file), index_col=0)
    df["task"] = task
    files.append(df)
  
  df_all = pd.concat(files) 
  df_all["tool"] = df_all["type"]
  df_all.loc[df_all["tool"] != "real", "tool"] = "sysbench"

  df_all["time"] = df_all["time"].apply(convert_time)
  df_all["label"] = df_all["tool"] #+ "_" + df_all["server"] 

  fig = px.strip(
    df_all,
    x="task", y="time",
    color="type",
    width=1500, 
    height=750,
    color_discrete_map=colors,
    title=machine,
    category_orders={
      "task": sorted(df_all["task"].unique())
    }
    # symbol="tool",
    # symbol_sequence=symbols
  ).update_traces(
    marker={
      "size": 15,
      "line": {
        "width": 2,
        "color": "DarkSlateGrey"
      },
      # "symbol": symbols
    },
    jitter=1,
  )
  

  fig.update_layout(
    legend=dict(
      font_size=20
    ),
    font_size = 30,
    yaxis_title="Time (s)",
    xaxis_title="Task"
  )

  fig.update_xaxes(
    tickangle = 45
  )
  savedir = this_dir.joinpath(f"new_test/{machine}")
  fig.write_image(savedir.joinpath(f"{machine}_time_plot.png"))

if __name__ == '__main__':
  main()
