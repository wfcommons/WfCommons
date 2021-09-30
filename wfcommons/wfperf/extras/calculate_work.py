from importlib.resources import read_text
import math
import pathlib
import time
import warnings 
import pandas as pd
import numpy as np
import subprocess
import plotly.express as px 
import plotly.graph_objs as go 
from typing import Callable, Dict, List, Optional, Tuple, Union

this_dir = pathlib.Path(__file__).resolve().parent
DRAW_FIT = False 
DO_TEST = False



tasks = {
    "mAdd":         (0.6, 1560989, 13.432413101196289), # max-prime found with one thread
    "mBgModel":     (0.7, 2513944, 15.435254573822021),
    "mProject":     (0.7, 5800800, 49.057353496551514),
    "mViewer":      (0.6, 11307825,107.97233176231384),
    "frequency":    (0.6, 29348298,415.41742062568665),
    "individuals":  (0.7, 6200800, 54.11822032928467) ,
    "mutation":     (0.7, 13300800,158.35056829452515)
}

def max_prime_average(path: pathlib.Path):
    df = pd.read_csv(path, index_col=0)
    df["duration"] = df["duration"].astype(float)
    df = df.drop(columns=['sample'])
    df = df.groupby(["max_prime"]).mean().reset_index()
    return df          
   
def find_max_prime_for_thread(work: Dict):  
    emp = fit_mx_distr()
    rows = []
    for task, (cpu, max_prime, duration) in tasks.items():
        rows.append([task, duration, cpu*10, max_prime, True])
        work = cpu * 10 * duration
        for i in range(1, 10):
            cpu_time = work / i
            rows.append([task, cpu_time, i, round(emp(cpu_time)), False])
    
    return pd.DataFrame(rows, columns=["task", "duration", "cpu", "max_prime", "real"])

def compress_legend(fig):
    group1_base, group2_base  = fig.data[0].task.split(",")
    lines_marker_name = []
    for i, trace in enumerate(fig.data):
        part1,part2 = trace.task.split(',')
        if part1 == group1_base:
            lines_marker_name.append(
                {"line": trace.line.to_plotly_json(), 
                 "marker": trace.marker.to_plotly_json(), 
                 "mode": trace.mode, 
                 "task": part2.lstrip(" ")
                }
            )
        if part2 != group2_base:
            trace['task'] = ''
            trace['showlegend']=False
        else:
            trace['task'] = part1
    
    ## Add the line/markers for the 2nd group
    for lmn in lines_marker_name:
        lmn["line"]["color"] = "green"
        lmn["marker"]["color"] = "green"
        fig.add_trace(go.Scatter(y=[None], **lmn))
    fig.update_layout(legend_title_text='', 
                        legend_itemclick=False,
                        legend_itemdoubleclick= False)

def draw(results: pd.DataFrame):
    results = results[~results["real"]]
    # results["duration"] = np.log2(results["duration"])
    # results["max_prime"] = np.log2(results["max_prime"])


    results["cpu"] = results["cpu"].apply(lambda x: f"{x*10}")

    wrap = 2
    fig: go.Figure
    fig = px.scatter(
        results, 
        y="max_prime", x="duration",
        color="cpu",
        facet_col="task", facet_col_wrap=wrap,
        facet_col_spacing=0.06,
        width=75*15, height=50*15,
        template="simple_white",
        labels={
            "duration": "Time (s)",
            "max_prime": "Max Prime",
            "cpu": "     CPU  "
        },
        category_orders={
            "task": sorted(results["task"].unique())
        }
    )   

    names = set(results["task"].unique())
    for i, task in enumerate(sorted(names)):
        row = math.ceil(len(names) / wrap) - (i // wrap)
        col = i%wrap+1
        fig.add_vline(
            x=tasks[task][2], 
            row=row, 
            col=col,
            line_width=3,
            line_color="green",
            line_dash="dot", 
            annotation_text="  Real", 
            annotation_position="top right"
        )

    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    for axis in fig.layout:
        if type(fig.layout[axis]) == go.layout.YAxis:
            fig.layout[axis].title.text = ''
    fig.update_layout(
        # keep the original annotations and add a list of new annotations:
        annotations = list(fig.layout.annotations) + 
        [go.layout.Annotation(
                x=-0.1,
                y=0.5,
                font=dict(
                    size=14
                ),
                showarrow=False,
                text="Max Prime",
                textangle=-90,
                xref="paper",
                yref="paper"
            )
        ]
    )

  
    fig.update_xaxes(matches=None)
    fig.update_yaxes(matches=None)
    
    fig.update_xaxes(showticklabels=True)
    fig.update_yaxes(showticklabels=True)
    fig.write_image("max_prime.png")



def fit_mx_distr(do_test: bool = False, do_draw: bool = False) -> Callable[[Union[np.ndarray, float]], Union[np.ndarray, float]]:
    csv = this_dir.joinpath("max_prime.csv")
    prime_df = max_prime_average(csv)
    prime_df = prime_df[prime_df["duration"] > 10]
    prime_df["duration"] = prime_df["duration"]
    # prime_df["max_prime"] = np.log2(prime_df["max_prime"])

    p = np.polyfit(
        x=prime_df["duration"].values, 
        y=prime_df["max_prime"].values,
        deg=3
    )
    est_max_prime = np.poly1d(p)

    if do_test:
        with this_dir.joinpath('original_max.txt').open('w+') as fp:
            for task, (duration) in tasks.items():
                start = time.time()
                proc = subprocess.Popen(
                    [
                        "sysbench", 
                        "cpu",
                        f"--cpu-max-prime={est_max_prime(duration)}", 
                        "run"
                    ], 
                    stdout=fp, 
                    stderr=fp
                )
                proc.wait()
                actual_duration = time.time() - start 
                print(f"{task} ({duration}s) -> sysbench (~{actual_duration}s)")
    
    if do_draw:
        fig = px.line(prime_df, x="duration", y="max_prime")
        xs = np.linspace(prime_df["duration"].min(), prime_df["duration"].max(), 1000)
        ys = est_max_prime(xs)
        fig.add_trace(
            go.Scatter(
                x=xs, y=ys, 
                mode="lines" 
            )
        )

        fig.write_image(this_dir.joinpath("max_prime_fit.png"))

    return np.poly1d(p)

def main():
    results = find_max_prime_for_thread(tasks)
    results.to_csv(this_dir.joinpath("thread_description.csv"))
    draw(results)
    fit = fit_mx_distr(do_test=False, do_draw=True)


if __name__ == "__main__":
    main() 