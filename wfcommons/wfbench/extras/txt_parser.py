import argparse
import pathlib
import pandas as pd
import re

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=pathlib.Path, help="Path to the txt files")
    parser.add_argument("-m", "--machine", help="Machine used for the experiment")
    parser.add_argument("-w", "--workflow", help="Workflow used for the experiment")

    return parser

def dirt_parser(path, name, machine, workflow):
    rows = []
    for file in path.glob("*.txt"):
        real, num = re.match(r"^(real)?.*?(\d_\d)?$", file.stem, re.DOTALL).groups()
        _type = real if real else num
        for match in re.finditer(r"\b([^\s]+)elapsed", file.read_text()):
            rows.append([_type, match.group(1)])
            

    df = pd.DataFrame(rows, columns=["type", "time"])

    savedir = path.joinpath("csv")
    savedir.mkdir(parents=True, exist_ok=True)

    savedir.joinpath(f"{name}_{workflow}_{machine}.csv").write_text(df.to_csv())

def chameleon_real_parser(path, workflow):
    rows = []
    for file in path.glob("*.txt"):
        print(f"Checking file:{str(file.stem)}")
        fullname = str(file.stem).split("-")
        name = fullname[0]
        if name == "madd":
            name = "mAdd"
        elif name == "mbgmodel":
            name = "mBgModel"
        elif name == "mviewer":
            name = "mViewer"
        elif name == "mproject":
            name = "mProject"
        machine = fullname[-1]
        
        with file.open("r") as fp:
            for line in fp.readlines():
                if workflow == "montage":
                    if "real" in line:
                        _, min, sec, ms = re.match(r".*?(([0-9]+)m([0-9]+).([0-9]+)).*?$", line, re.DOTALL).groups()
                        time = f"{min}:{sec}.{ms}"
                        rows.append([name, machine, "real", time])      
                elif workflow == "genome":
                    if "user" in line:
                        _, min, sec, ms = re.match(r".*?(([0-9]+)m([0-9]+).([0-9]+)).*?$", line, re.DOTALL).groups()
                        time = f"{min}:{sec}.{ms}"
                        rows.append([name, machine, "real", time])  

    df = pd.DataFrame(rows, columns=["name","machine", "type", "time"])

    savedir = path.joinpath("csv")
    savedir.mkdir(parents=True, exist_ok=True)

    savedir.joinpath(f"all_{workflow}_real.csv").write_text(df.to_csv())
 
def genome_parser(path, name, machine, workflow):
    rows = []
    for file in path.glob("*.txt"):
        real, num = re.match(r"^(real)?.*?(\d_\d)?$", file.stem, re.DOTALL).groups()
        _type = real if real else num
        for match in re.finditer(r"\b([^\s]+)user", file.read_text()):
            rows.append([_type, match.group(1)])
            

        df = pd.DataFrame(rows, columns=["type", "time"])

        savedir = path.joinpath("csv")
        savedir.mkdir(parents=True, exist_ok=True)

        savedir.joinpath(f"{name}_{workflow}_{machine}.csv").write_text(df.to_csv())

def main():
    parser = get_parser()
    args = parser.parse_args()
    path: pathlib.Path = args.path
    name = path.stem 
    machine = args.machine
    workflow = args.workflow
    
    # dirt_parser(path, name, machine, workflow)
    # chameleon_real_parser(path, workflow)
    genome_parser(path, name, machine, workflow)    
        

if __name__ == "__main__":
    main()