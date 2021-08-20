import pandas as pd
import pathlib
import argparse
import re

this_dir = pathlib.Path(__file__).resolve().parent

def find_mean(df: pd.Series) -> pd.DataFrame:
    return df.mean()

def convert_time(text: str) -> float:
  return sum([float(p) * 60**i for i, p in enumerate(str(text).split(":")[::-1])])
    
def percent_error(real_mean: float, bm_mean: float) -> float:
    return ((bm_mean - real_mean)/real_mean)*100

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to the csv")
    parser.add_argument("-r", "--ratio", help="Correct CPU percentage for task")
    parser.add_argument("-t", "--task", help="Task")

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    path = pathlib.Path(args.path)
    ratio = float(args.ratio)
    cpu_thread = int(ratio*10) 
    mem_thread = 10 - cpu_thread
    label = f'{cpu_thread}_{mem_thread}'  

    lines = []
    for file in path.glob("*.csv"):
        df = pd.read_csv(file, index_col=0)
        df["time"] = df["time"].apply(convert_time)
        real = df[df["type"] == "real"]
        real = real.set_index("type")
        bm = df[~(df["type"] == "real")]
        bm = bm.set_index("type")


        real_mean = round(float(find_mean(real["time"])), 3)
        bm_thread = bm.loc[label]
        bm_mean = round(float(find_mean(bm_thread)), 3)

        error = round(percent_error(real_mean, bm_mean), 3)

        savedir = file.parent.joinpath(f"error")
        savedir.mkdir(exist_ok=True, parents=True)

        _, machine = re.match(r"^.+_(^.+)?.*?_(.+).*?.*?$", str(file.stem), re.DOTALL).groups()
        lines.append(f'{machine} {cpu_thread} {mem_thread} {real_mean} {bm_mean} {error} \n')
    
    with savedir.joinpath(f"{args.task}_error.txt").open("w+") as fp:
        fp.writelines(lines)


if __name__ == "__main__":
    main()