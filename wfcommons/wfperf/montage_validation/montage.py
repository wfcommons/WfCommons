from wfcommons.wfperf.montage_validation.montage_perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import MontageRecipe
import pathlib

this_dir = pathlib.Path(__file__).resolve().parent

def total_tasks():
    num_tasks = {'mProject': 64, 
                 'mDiffFit': 2016, 
                 'mConcatFit': 1, 
                 'mBgModel': 1, 
                 'mBackground': 64, 
                 'mImgtbl': 1, 
                 'mAdd': 1, 
                 'mViewer': 1}
    total = 0
    for (k, v) in num_tasks.items():
        if k == 'mViewer':
            num_tasks[k] = v*3 + 1 
        else: 
            num_tasks[k] = v*3

    for values in num_tasks.values():
        total += values

    return total

def main():
     
    num_tasks = total_tasks()
    tasks = {'mProject': (12800000, 7, 120), 
             'mDiffFit': (24900000 , 7, 1), 
             'mConcatFit': (24900000 , 7, 5), 
             'mBgModel': (1910000, 7, 120), 
             'mBackground': (24900000 , 7, 1), 
             'mImgtbl': (24900000 , 7, 2),
             'mAdd': (1050000, 6, 120),
             'mViewer': (7400000, 6, 120)}

    
    bench = WorkflowBenchmark(MontageRecipe, num_tasks)
    bench.create("/home/tgcoleman/tests/Montage", tasks, verbose=True)
    # return
    # locked = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt.lock")
    # coresfile = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt")
    

    
    
if __name__ == "__main__":
    main()