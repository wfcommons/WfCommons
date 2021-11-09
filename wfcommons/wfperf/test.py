import os
from Pegasus.api import *


def which(file):
    for path in os.environ['PATH'].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return os.path.join(path, file)
    return None


wf = Workflow('Blast-Benchmark', infer_dependencies=True)
tc = TransformationCatalog()
rc = ReplicaCatalog()

task_output_files = {}

transformation_path = which('/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
if transformation_path is None:
    raise RuntimeError('Unable to find /home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
transformation = Transformation('split_fasta', site='local',
                                pfn='/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py',
                                is_stageable=True)
transformation.add_env(PATH='/usr/bin:/bin:.')
transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')
tc.add_transformations(transformation)

transformation_path = which('/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
if transformation_path is None:
    raise RuntimeError('Unable to find /home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
transformation = Transformation('blastall', site='local',
                                pfn=',/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py',
                                is_stageable=True)
transformation.add_env(PATH='/usr/bin:/bin:.')
transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')
tc.add_transformations(transformation)

transformation_path = which('/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
if transformation_path is None:
    raise RuntimeError('Unable to find /home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
transformation = Transformation('cat_blast', site='local',
                                pfn=',/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py',
                                is_stageable=True)
transformation.add_env(PATH='/usr/bin:/bin:.')
transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')
tc.add_transformations(transformation)

transformation_path = which('/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
if transformation_path is None:
    raise RuntimeError('Unable to find /home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py')
transformation = Transformation('cat', site='local',
                                pfn='/home/cc/test/wfcommons/wfcommons/wfperf/wfperf_benchmark.py',
                                is_stageable=True)
transformation.add_env(PATH='/usr/bin:/bin:.')
transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')
tc.add_transformations(transformation)

job_1 = Job('split_fasta', _id='split_fasta_00000001')
out_file_1 = File('split_fasta_00000001_output.txt')
task_output_files['job_1'] = out_file_1
job_1.add_outputs(out_file_1, stage_out=False, register_replica=False)
job_1.add_args('split', 'split_fasta_00000001', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=split_fasta_00000001_output.txt')
wf.add_jobs(job_1)

job_2 = Job('blastall', _id='blastall_00000002')
out_file_2 = File('blastall_00000002_output.txt')
task_output_files['job_2'] = out_file_2
job_2.add_outputs(out_file_2, stage_out=False, register_replica=False)
job_2.add_args('blastall', 'blastall_00000002', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000002_output.txt')
wf.add_jobs(job_2)

job_3 = Job('cat_blast', _id='cat_blast_00000042')
out_file_3 = File('cat_blast_00000042_output.txt')
task_output_files['job_3'] = out_file_3
job_3.add_outputs(out_file_3, stage_out=True, register_replica=True)
job_3.add_args('cat', 'cat_blast_00000042', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=cat_blast_00000042_output.txt')
wf.add_jobs(job_3)

if 'job_2' in task_output_files:
  job_3.add_inputs(task_output_files['job_2'])
wf.add_dependency(job_3, parents=[job_2])

job_4 = Job('cat', _id='cat_00000043')
out_file_4 = File('cat_00000043_output.txt')
task_output_files['job_4'] = out_file_4
job_4.add_outputs(out_file_4, stage_out=True, register_replica=True)
job_4.add_args('cat', 'cat_00000043', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=cat_00000043_output.txt')
wf.add_jobs(job_4)

if 'job_2' in task_output_files:
  job_4.add_inputs(task_output_files['job_2'])
wf.add_dependency(job_4, parents=[job_2])

if 'job_1' in task_output_files:
  job_2.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_2, parents=[job_1])

job_5 = Job('blastall', _id='blastall_00000003')
out_file_5 = File('blastall_00000003_output.txt')
task_output_files['job_5'] = out_file_5
job_5.add_outputs(out_file_5, stage_out=False, register_replica=False)
job_5.add_args('blastall', 'blastall_00000003', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000003_output.txt')
wf.add_jobs(job_5)

if 'job_5' in task_output_files:
  job_3.add_inputs(task_output_files['job_5'])
wf.add_dependency(job_3, parents=[job_5])

if 'job_5' in task_output_files:
  job_4.add_inputs(task_output_files['job_5'])
wf.add_dependency(job_4, parents=[job_5])

if 'job_1' in task_output_files:
  job_5.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_5, parents=[job_1])

job_6 = Job('blastall', _id='blastall_00000004')
out_file_6 = File('blastall_00000004_output.txt')
task_output_files['job_6'] = out_file_6
job_6.add_outputs(out_file_6, stage_out=False, register_replica=False)
job_6.add_args('blastall', 'blastall_00000004', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000004_output.txt')
wf.add_jobs(job_6)

if 'job_6' in task_output_files:
  job_3.add_inputs(task_output_files['job_6'])
wf.add_dependency(job_3, parents=[job_6])

if 'job_6' in task_output_files:
  job_4.add_inputs(task_output_files['job_6'])
wf.add_dependency(job_4, parents=[job_6])

if 'job_1' in task_output_files:
  job_6.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_6, parents=[job_1])

job_7 = Job('blastall', _id='blastall_00000005')
out_file_7 = File('blastall_00000005_output.txt')
task_output_files['job_7'] = out_file_7
job_7.add_outputs(out_file_7, stage_out=False, register_replica=False)
job_7.add_args('blastall', 'blastall_00000005', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000005_output.txt')
wf.add_jobs(job_7)

if 'job_7' in task_output_files:
  job_3.add_inputs(task_output_files['job_7'])
wf.add_dependency(job_3, parents=[job_7])

if 'job_7' in task_output_files:
  job_4.add_inputs(task_output_files['job_7'])
wf.add_dependency(job_4, parents=[job_7])

if 'job_1' in task_output_files:
  job_7.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_7, parents=[job_1])

job_8 = Job('blastall', _id='blastall_00000006')
out_file_8 = File('blastall_00000006_output.txt')
task_output_files['job_8'] = out_file_8
job_8.add_outputs(out_file_8, stage_out=False, register_replica=False)
job_8.add_args('blastall', 'blastall_00000006', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000006_output.txt')
wf.add_jobs(job_8)

if 'job_8' in task_output_files:
  job_3.add_inputs(task_output_files['job_8'])
wf.add_dependency(job_3, parents=[job_8])

if 'job_8' in task_output_files:
  job_4.add_inputs(task_output_files['job_8'])
wf.add_dependency(job_4, parents=[job_8])

if 'job_1' in task_output_files:
  job_8.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_8, parents=[job_1])

job_9 = Job('blastall', _id='blastall_00000007')
out_file_9 = File('blastall_00000007_output.txt')
task_output_files['job_9'] = out_file_9
job_9.add_outputs(out_file_9, stage_out=False, register_replica=False)
job_9.add_args('blastall', 'blastall_00000007', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000007_output.txt')
wf.add_jobs(job_9)

if 'job_9' in task_output_files:
  job_3.add_inputs(task_output_files['job_9'])
wf.add_dependency(job_3, parents=[job_9])

if 'job_9' in task_output_files:
  job_4.add_inputs(task_output_files['job_9'])
wf.add_dependency(job_4, parents=[job_9])

if 'job_1' in task_output_files:
  job_9.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_9, parents=[job_1])

job_10 = Job('blastall', _id='blastall_00000008')
out_file_10 = File('blastall_00000008_output.txt')
task_output_files['job_10'] = out_file_10
job_10.add_outputs(out_file_10, stage_out=False, register_replica=False)
job_10.add_args('blastall', 'blastall_00000008', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000008_output.txt')
wf.add_jobs(job_10)

if 'job_10' in task_output_files:
  job_3.add_inputs(task_output_files['job_10'])
wf.add_dependency(job_3, parents=[job_10])

if 'job_10' in task_output_files:
  job_4.add_inputs(task_output_files['job_10'])
wf.add_dependency(job_4, parents=[job_10])

if 'job_1' in task_output_files:
  job_10.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_10, parents=[job_1])

job_11 = Job('blastall', _id='blastall_00000009')
out_file_11 = File('blastall_00000009_output.txt')
task_output_files['job_11'] = out_file_11
job_11.add_outputs(out_file_11, stage_out=False, register_replica=False)
job_11.add_args('blastall', 'blastall_00000009', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000009_output.txt')
wf.add_jobs(job_11)

if 'job_11' in task_output_files:
  job_3.add_inputs(task_output_files['job_11'])
wf.add_dependency(job_3, parents=[job_11])

if 'job_11' in task_output_files:
  job_4.add_inputs(task_output_files['job_11'])
wf.add_dependency(job_4, parents=[job_11])

if 'job_1' in task_output_files:
  job_11.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_11, parents=[job_1])

job_12 = Job('blastall', _id='blastall_00000010')
out_file_12 = File('blastall_00000010_output.txt')
task_output_files['job_12'] = out_file_12
job_12.add_outputs(out_file_12, stage_out=False, register_replica=False)
job_12.add_args('blastall', 'blastall_00000010', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000010_output.txt')
wf.add_jobs(job_12)

if 'job_12' in task_output_files:
  job_3.add_inputs(task_output_files['job_12'])
wf.add_dependency(job_3, parents=[job_12])

if 'job_12' in task_output_files:
  job_4.add_inputs(task_output_files['job_12'])
wf.add_dependency(job_4, parents=[job_12])

if 'job_1' in task_output_files:
  job_12.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_12, parents=[job_1])

job_13 = Job('blastall', _id='blastall_00000011')
out_file_13 = File('blastall_00000011_output.txt')
task_output_files['job_13'] = out_file_13
job_13.add_outputs(out_file_13, stage_out=False, register_replica=False)
job_13.add_args('blastall', 'blastall_00000011', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000011_output.txt')
wf.add_jobs(job_13)

if 'job_13' in task_output_files:
  job_3.add_inputs(task_output_files['job_13'])
wf.add_dependency(job_3, parents=[job_13])

if 'job_13' in task_output_files:
  job_4.add_inputs(task_output_files['job_13'])
wf.add_dependency(job_4, parents=[job_13])

if 'job_1' in task_output_files:
  job_13.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_13, parents=[job_1])

job_14 = Job('blastall', _id='blastall_00000012')
out_file_14 = File('blastall_00000012_output.txt')
task_output_files['job_14'] = out_file_14
job_14.add_outputs(out_file_14, stage_out=False, register_replica=False)
job_14.add_args('blastall', 'blastall_00000012', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000012_output.txt')
wf.add_jobs(job_14)

if 'job_14' in task_output_files:
  job_3.add_inputs(task_output_files['job_14'])
wf.add_dependency(job_3, parents=[job_14])

if 'job_14' in task_output_files:
  job_4.add_inputs(task_output_files['job_14'])
wf.add_dependency(job_4, parents=[job_14])

if 'job_1' in task_output_files:
  job_14.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_14, parents=[job_1])

job_15 = Job('blastall', _id='blastall_00000013')
out_file_15 = File('blastall_00000013_output.txt')
task_output_files['job_15'] = out_file_15
job_15.add_outputs(out_file_15, stage_out=False, register_replica=False)
job_15.add_args('blastall', 'blastall_00000013', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000013_output.txt')
wf.add_jobs(job_15)

if 'job_15' in task_output_files:
  job_3.add_inputs(task_output_files['job_15'])
wf.add_dependency(job_3, parents=[job_15])

if 'job_15' in task_output_files:
  job_4.add_inputs(task_output_files['job_15'])
wf.add_dependency(job_4, parents=[job_15])

if 'job_1' in task_output_files:
  job_15.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_15, parents=[job_1])

job_16 = Job('blastall', _id='blastall_00000014')
out_file_16 = File('blastall_00000014_output.txt')
task_output_files['job_16'] = out_file_16
job_16.add_outputs(out_file_16, stage_out=False, register_replica=False)
job_16.add_args('blastall', 'blastall_00000014', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000014_output.txt')
wf.add_jobs(job_16)

if 'job_16' in task_output_files:
  job_3.add_inputs(task_output_files['job_16'])
wf.add_dependency(job_3, parents=[job_16])

if 'job_16' in task_output_files:
  job_4.add_inputs(task_output_files['job_16'])
wf.add_dependency(job_4, parents=[job_16])

if 'job_1' in task_output_files:
  job_16.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_16, parents=[job_1])

job_17 = Job('blastall', _id='blastall_00000015')
out_file_17 = File('blastall_00000015_output.txt')
task_output_files['job_17'] = out_file_17
job_17.add_outputs(out_file_17, stage_out=False, register_replica=False)
job_17.add_args('blastall', 'blastall_00000015', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000015_output.txt')
wf.add_jobs(job_17)

if 'job_17' in task_output_files:
  job_3.add_inputs(task_output_files['job_17'])
wf.add_dependency(job_3, parents=[job_17])

if 'job_17' in task_output_files:
  job_4.add_inputs(task_output_files['job_17'])
wf.add_dependency(job_4, parents=[job_17])

if 'job_1' in task_output_files:
  job_17.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_17, parents=[job_1])

job_18 = Job('blastall', _id='blastall_00000016')
out_file_18 = File('blastall_00000016_output.txt')
task_output_files['job_18'] = out_file_18
job_18.add_outputs(out_file_18, stage_out=False, register_replica=False)
job_18.add_args('blastall', 'blastall_00000016', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000016_output.txt')
wf.add_jobs(job_18)

if 'job_18' in task_output_files:
  job_3.add_inputs(task_output_files['job_18'])
wf.add_dependency(job_3, parents=[job_18])

if 'job_18' in task_output_files:
  job_4.add_inputs(task_output_files['job_18'])
wf.add_dependency(job_4, parents=[job_18])

if 'job_1' in task_output_files:
  job_18.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_18, parents=[job_1])

job_19 = Job('blastall', _id='blastall_00000017')
out_file_19 = File('blastall_00000017_output.txt')
task_output_files['job_19'] = out_file_19
job_19.add_outputs(out_file_19, stage_out=False, register_replica=False)
job_19.add_args('blastall', 'blastall_00000017', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000017_output.txt')
wf.add_jobs(job_19)

if 'job_19' in task_output_files:
  job_3.add_inputs(task_output_files['job_19'])
wf.add_dependency(job_3, parents=[job_19])

if 'job_19' in task_output_files:
  job_4.add_inputs(task_output_files['job_19'])
wf.add_dependency(job_4, parents=[job_19])

if 'job_1' in task_output_files:
  job_19.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_19, parents=[job_1])

job_20 = Job('blastall', _id='blastall_00000018')
out_file_20 = File('blastall_00000018_output.txt')
task_output_files['job_20'] = out_file_20
job_20.add_outputs(out_file_20, stage_out=False, register_replica=False)
job_20.add_args('blastall', 'blastall_00000018', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000018_output.txt')
wf.add_jobs(job_20)

if 'job_20' in task_output_files:
  job_3.add_inputs(task_output_files['job_20'])
wf.add_dependency(job_3, parents=[job_20])

if 'job_20' in task_output_files:
  job_4.add_inputs(task_output_files['job_20'])
wf.add_dependency(job_4, parents=[job_20])

if 'job_1' in task_output_files:
  job_20.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_20, parents=[job_1])

job_21 = Job('blastall', _id='blastall_00000019')
out_file_21 = File('blastall_00000019_output.txt')
task_output_files['job_21'] = out_file_21
job_21.add_outputs(out_file_21, stage_out=False, register_replica=False)
job_21.add_args('blastall', 'blastall_00000019', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000019_output.txt')
wf.add_jobs(job_21)

if 'job_21' in task_output_files:
  job_3.add_inputs(task_output_files['job_21'])
wf.add_dependency(job_3, parents=[job_21])

if 'job_21' in task_output_files:
  job_4.add_inputs(task_output_files['job_21'])
wf.add_dependency(job_4, parents=[job_21])

if 'job_1' in task_output_files:
  job_21.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_21, parents=[job_1])

job_22 = Job('blastall', _id='blastall_00000020')
out_file_22 = File('blastall_00000020_output.txt')
task_output_files['job_22'] = out_file_22
job_22.add_outputs(out_file_22, stage_out=False, register_replica=False)
job_22.add_args('blastall', 'blastall_00000020', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000020_output.txt')
wf.add_jobs(job_22)

if 'job_22' in task_output_files:
  job_3.add_inputs(task_output_files['job_22'])
wf.add_dependency(job_3, parents=[job_22])

if 'job_22' in task_output_files:
  job_4.add_inputs(task_output_files['job_22'])
wf.add_dependency(job_4, parents=[job_22])

if 'job_1' in task_output_files:
  job_22.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_22, parents=[job_1])

job_23 = Job('blastall', _id='blastall_00000021')
out_file_23 = File('blastall_00000021_output.txt')
task_output_files['job_23'] = out_file_23
job_23.add_outputs(out_file_23, stage_out=False, register_replica=False)
job_23.add_args('blastall', 'blastall_00000021', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000021_output.txt')
wf.add_jobs(job_23)

if 'job_23' in task_output_files:
  job_3.add_inputs(task_output_files['job_23'])
wf.add_dependency(job_3, parents=[job_23])

if 'job_23' in task_output_files:
  job_4.add_inputs(task_output_files['job_23'])
wf.add_dependency(job_4, parents=[job_23])

if 'job_1' in task_output_files:
  job_23.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_23, parents=[job_1])

job_24 = Job('blastall', _id='blastall_00000022')
out_file_24 = File('blastall_00000022_output.txt')
task_output_files['job_24'] = out_file_24
job_24.add_outputs(out_file_24, stage_out=False, register_replica=False)
job_24.add_args('blastall', 'blastall_00000022', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000022_output.txt')
wf.add_jobs(job_24)

if 'job_24' in task_output_files:
  job_3.add_inputs(task_output_files['job_24'])
wf.add_dependency(job_3, parents=[job_24])

if 'job_24' in task_output_files:
  job_4.add_inputs(task_output_files['job_24'])
wf.add_dependency(job_4, parents=[job_24])

if 'job_1' in task_output_files:
  job_24.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_24, parents=[job_1])

job_25 = Job('blastall', _id='blastall_00000023')
out_file_25 = File('blastall_00000023_output.txt')
task_output_files['job_25'] = out_file_25
job_25.add_outputs(out_file_25, stage_out=False, register_replica=False)
job_25.add_args('blastall', 'blastall_00000023', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000023_output.txt')
wf.add_jobs(job_25)

if 'job_25' in task_output_files:
  job_3.add_inputs(task_output_files['job_25'])
wf.add_dependency(job_3, parents=[job_25])

if 'job_25' in task_output_files:
  job_4.add_inputs(task_output_files['job_25'])
wf.add_dependency(job_4, parents=[job_25])

if 'job_1' in task_output_files:
  job_25.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_25, parents=[job_1])

job_26 = Job('blastall', _id='blastall_00000024')
out_file_26 = File('blastall_00000024_output.txt')
task_output_files['job_26'] = out_file_26
job_26.add_outputs(out_file_26, stage_out=False, register_replica=False)
job_26.add_args('blastall', 'blastall_00000024', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000024_output.txt')
wf.add_jobs(job_26)

if 'job_26' in task_output_files:
  job_3.add_inputs(task_output_files['job_26'])
wf.add_dependency(job_3, parents=[job_26])

if 'job_26' in task_output_files:
  job_4.add_inputs(task_output_files['job_26'])
wf.add_dependency(job_4, parents=[job_26])

if 'job_1' in task_output_files:
  job_26.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_26, parents=[job_1])

job_27 = Job('blastall', _id='blastall_00000025')
out_file_27 = File('blastall_00000025_output.txt')
task_output_files['job_27'] = out_file_27
job_27.add_outputs(out_file_27, stage_out=False, register_replica=False)
job_27.add_args('blastall', 'blastall_00000025', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000025_output.txt')
wf.add_jobs(job_27)

if 'job_27' in task_output_files:
  job_3.add_inputs(task_output_files['job_27'])
wf.add_dependency(job_3, parents=[job_27])

if 'job_27' in task_output_files:
  job_4.add_inputs(task_output_files['job_27'])
wf.add_dependency(job_4, parents=[job_27])

if 'job_1' in task_output_files:
  job_27.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_27, parents=[job_1])

job_28 = Job('blastall', _id='blastall_00000026')
out_file_28 = File('blastall_00000026_output.txt')
task_output_files['job_28'] = out_file_28
job_28.add_outputs(out_file_28, stage_out=False, register_replica=False)
job_28.add_args('blastall', 'blastall_00000026', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000026_output.txt')
wf.add_jobs(job_28)

if 'job_28' in task_output_files:
  job_3.add_inputs(task_output_files['job_28'])
wf.add_dependency(job_3, parents=[job_28])

if 'job_28' in task_output_files:
  job_4.add_inputs(task_output_files['job_28'])
wf.add_dependency(job_4, parents=[job_28])

if 'job_1' in task_output_files:
  job_28.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_28, parents=[job_1])

job_29 = Job('blastall', _id='blastall_00000027')
out_file_29 = File('blastall_00000027_output.txt')
task_output_files['job_29'] = out_file_29
job_29.add_outputs(out_file_29, stage_out=False, register_replica=False)
job_29.add_args('blastall', 'blastall_00000027', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000027_output.txt')
wf.add_jobs(job_29)

if 'job_29' in task_output_files:
  job_3.add_inputs(task_output_files['job_29'])
wf.add_dependency(job_3, parents=[job_29])

if 'job_29' in task_output_files:
  job_4.add_inputs(task_output_files['job_29'])
wf.add_dependency(job_4, parents=[job_29])

if 'job_1' in task_output_files:
  job_29.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_29, parents=[job_1])

job_30 = Job('blastall', _id='blastall_00000028')
out_file_30 = File('blastall_00000028_output.txt')
task_output_files['job_30'] = out_file_30
job_30.add_outputs(out_file_30, stage_out=False, register_replica=False)
job_30.add_args('blastall', 'blastall_00000028', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000028_output.txt')
wf.add_jobs(job_30)

if 'job_30' in task_output_files:
  job_3.add_inputs(task_output_files['job_30'])
wf.add_dependency(job_3, parents=[job_30])

if 'job_30' in task_output_files:
  job_4.add_inputs(task_output_files['job_30'])
wf.add_dependency(job_4, parents=[job_30])

if 'job_1' in task_output_files:
  job_30.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_30, parents=[job_1])

job_31 = Job('blastall', _id='blastall_00000029')
out_file_31 = File('blastall_00000029_output.txt')
task_output_files['job_31'] = out_file_31
job_31.add_outputs(out_file_31, stage_out=False, register_replica=False)
job_31.add_args('blastall', 'blastall_00000029', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000029_output.txt')
wf.add_jobs(job_31)

if 'job_31' in task_output_files:
  job_3.add_inputs(task_output_files['job_31'])
wf.add_dependency(job_3, parents=[job_31])

if 'job_31' in task_output_files:
  job_4.add_inputs(task_output_files['job_31'])
wf.add_dependency(job_4, parents=[job_31])

if 'job_1' in task_output_files:
  job_31.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_31, parents=[job_1])

job_32 = Job('blastall', _id='blastall_00000030')
out_file_32 = File('blastall_00000030_output.txt')
task_output_files['job_32'] = out_file_32
job_32.add_outputs(out_file_32, stage_out=False, register_replica=False)
job_32.add_args('blastall', 'blastall_00000030', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000030_output.txt')
wf.add_jobs(job_32)

if 'job_32' in task_output_files:
  job_3.add_inputs(task_output_files['job_32'])
wf.add_dependency(job_3, parents=[job_32])

if 'job_32' in task_output_files:
  job_4.add_inputs(task_output_files['job_32'])
wf.add_dependency(job_4, parents=[job_32])

if 'job_1' in task_output_files:
  job_32.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_32, parents=[job_1])

job_33 = Job('blastall', _id='blastall_00000031')
out_file_33 = File('blastall_00000031_output.txt')
task_output_files['job_33'] = out_file_33
job_33.add_outputs(out_file_33, stage_out=False, register_replica=False)
job_33.add_args('blastall', 'blastall_00000031', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000031_output.txt')
wf.add_jobs(job_33)

if 'job_33' in task_output_files:
  job_3.add_inputs(task_output_files['job_33'])
wf.add_dependency(job_3, parents=[job_33])

if 'job_33' in task_output_files:
  job_4.add_inputs(task_output_files['job_33'])
wf.add_dependency(job_4, parents=[job_33])

if 'job_1' in task_output_files:
  job_33.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_33, parents=[job_1])

job_34 = Job('blastall', _id='blastall_00000032')
out_file_34 = File('blastall_00000032_output.txt')
task_output_files['job_34'] = out_file_34
job_34.add_outputs(out_file_34, stage_out=False, register_replica=False)
job_34.add_args('blastall', 'blastall_00000032', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000032_output.txt')
wf.add_jobs(job_34)

if 'job_34' in task_output_files:
  job_3.add_inputs(task_output_files['job_34'])
wf.add_dependency(job_3, parents=[job_34])

if 'job_34' in task_output_files:
  job_4.add_inputs(task_output_files['job_34'])
wf.add_dependency(job_4, parents=[job_34])

if 'job_1' in task_output_files:
  job_34.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_34, parents=[job_1])

job_35 = Job('blastall', _id='blastall_00000033')
out_file_35 = File('blastall_00000033_output.txt')
task_output_files['job_35'] = out_file_35
job_35.add_outputs(out_file_35, stage_out=False, register_replica=False)
job_35.add_args('blastall', 'blastall_00000033', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000033_output.txt')
wf.add_jobs(job_35)

if 'job_35' in task_output_files:
  job_3.add_inputs(task_output_files['job_35'])
wf.add_dependency(job_3, parents=[job_35])

if 'job_35' in task_output_files:
  job_4.add_inputs(task_output_files['job_35'])
wf.add_dependency(job_4, parents=[job_35])

if 'job_1' in task_output_files:
  job_35.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_35, parents=[job_1])

job_36 = Job('blastall', _id='blastall_00000034')
out_file_36 = File('blastall_00000034_output.txt')
task_output_files['job_36'] = out_file_36
job_36.add_outputs(out_file_36, stage_out=False, register_replica=False)
job_36.add_args('blastall', 'blastall_00000034', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000034_output.txt')
wf.add_jobs(job_36)

if 'job_36' in task_output_files:
  job_3.add_inputs(task_output_files['job_36'])
wf.add_dependency(job_3, parents=[job_36])

if 'job_36' in task_output_files:
  job_4.add_inputs(task_output_files['job_36'])
wf.add_dependency(job_4, parents=[job_36])

if 'job_1' in task_output_files:
  job_36.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_36, parents=[job_1])

job_37 = Job('blastall', _id='blastall_00000035')
out_file_37 = File('blastall_00000035_output.txt')
task_output_files['job_37'] = out_file_37
job_37.add_outputs(out_file_37, stage_out=False, register_replica=False)
job_37.add_args('blastall', 'blastall_00000035', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000035_output.txt')
wf.add_jobs(job_37)

if 'job_37' in task_output_files:
  job_3.add_inputs(task_output_files['job_37'])
wf.add_dependency(job_3, parents=[job_37])

if 'job_37' in task_output_files:
  job_4.add_inputs(task_output_files['job_37'])
wf.add_dependency(job_4, parents=[job_37])

if 'job_1' in task_output_files:
  job_37.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_37, parents=[job_1])

job_38 = Job('blastall', _id='blastall_00000036')
out_file_38 = File('blastall_00000036_output.txt')
task_output_files['job_38'] = out_file_38
job_38.add_outputs(out_file_38, stage_out=False, register_replica=False)
job_38.add_args('blastall', 'blastall_00000036', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000036_output.txt')
wf.add_jobs(job_38)

if 'job_38' in task_output_files:
  job_3.add_inputs(task_output_files['job_38'])
wf.add_dependency(job_3, parents=[job_38])

if 'job_38' in task_output_files:
  job_4.add_inputs(task_output_files['job_38'])
wf.add_dependency(job_4, parents=[job_38])

if 'job_1' in task_output_files:
  job_38.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_38, parents=[job_1])

job_39 = Job('blastall', _id='blastall_00000037')
out_file_39 = File('blastall_00000037_output.txt')
task_output_files['job_39'] = out_file_39
job_39.add_outputs(out_file_39, stage_out=False, register_replica=False)
job_39.add_args('blastall', 'blastall_00000037', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000037_output.txt')
wf.add_jobs(job_39)

if 'job_39' in task_output_files:
  job_3.add_inputs(task_output_files['job_39'])
wf.add_dependency(job_3, parents=[job_39])

if 'job_39' in task_output_files:
  job_4.add_inputs(task_output_files['job_39'])
wf.add_dependency(job_4, parents=[job_39])

if 'job_1' in task_output_files:
  job_39.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_39, parents=[job_1])

job_40 = Job('blastall', _id='blastall_00000038')
out_file_40 = File('blastall_00000038_output.txt')
task_output_files['job_40'] = out_file_40
job_40.add_outputs(out_file_40, stage_out=False, register_replica=False)
job_40.add_args('blastall', 'blastall_00000038', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000038_output.txt')
wf.add_jobs(job_40)

if 'job_40' in task_output_files:
  job_3.add_inputs(task_output_files['job_40'])
wf.add_dependency(job_3, parents=[job_40])

if 'job_40' in task_output_files:
  job_4.add_inputs(task_output_files['job_40'])
wf.add_dependency(job_4, parents=[job_40])

if 'job_1' in task_output_files:
  job_40.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_40, parents=[job_1])

job_41 = Job('blastall', _id='blastall_00000039')
out_file_41 = File('blastall_00000039_output.txt')
task_output_files['job_41'] = out_file_41
job_41.add_outputs(out_file_41, stage_out=False, register_replica=False)
job_41.add_args('blastall', 'blastall_00000039', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000039_output.txt')
wf.add_jobs(job_41)

if 'job_41' in task_output_files:
  job_3.add_inputs(task_output_files['job_41'])
wf.add_dependency(job_3, parents=[job_41])

if 'job_41' in task_output_files:
  job_4.add_inputs(task_output_files['job_41'])
wf.add_dependency(job_4, parents=[job_41])

if 'job_1' in task_output_files:
  job_41.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_41, parents=[job_1])

job_42 = Job('blastall', _id='blastall_00000040')
out_file_42 = File('blastall_00000040_output.txt')
task_output_files['job_42'] = out_file_42
job_42.add_outputs(out_file_42, stage_out=False, register_replica=False)
job_42.add_args('blastall', 'blastall_00000040', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000040_output.txt')
wf.add_jobs(job_42)

if 'job_42' in task_output_files:
  job_3.add_inputs(task_output_files['job_42'])
wf.add_dependency(job_3, parents=[job_42])

if 'job_42' in task_output_files:
  job_4.add_inputs(task_output_files['job_42'])
wf.add_dependency(job_4, parents=[job_42])

if 'job_1' in task_output_files:
  job_42.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_42, parents=[job_1])

job_43 = Job('blastall', _id='blastall_00000041')
out_file_43 = File('blastall_00000041_output.txt')
task_output_files['job_43'] = out_file_43
job_43.add_outputs(out_file_43, stage_out=False, register_replica=False)
job_43.add_args('blastall', 'blastall_00000041', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000041_output.txt')
wf.add_jobs(job_43)

if 'job_43' in task_output_files:
  job_3.add_inputs(task_output_files['job_43'])
wf.add_dependency(job_3, parents=[job_43])

if 'job_43' in task_output_files:
  job_4.add_inputs(task_output_files['job_43'])
wf.add_dependency(job_4, parents=[job_43])

if 'job_1' in task_output_files:
  job_43.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_43, parents=[job_1])

job_44 = Job('blastall', _id='blastall_00000044')
out_file_44 = File('blastall_00000044_output.txt')
task_output_files['job_44'] = out_file_44
job_44.add_outputs(out_file_44, stage_out=False, register_replica=False)
job_44.add_args('blastall', 'blastall_00000044', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000044_output.txt')
wf.add_jobs(job_44)

if 'job_44' in task_output_files:
  job_3.add_inputs(task_output_files['job_44'])
wf.add_dependency(job_3, parents=[job_44])

if 'job_44' in task_output_files:
  job_4.add_inputs(task_output_files['job_44'])
wf.add_dependency(job_4, parents=[job_44])

if 'job_1' in task_output_files:
  job_44.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_44, parents=[job_1])

job_45 = Job('blastall', _id='blastall_00000045')
out_file_45 = File('blastall_00000045_output.txt')
task_output_files['job_45'] = out_file_45
job_45.add_outputs(out_file_45, stage_out=False, register_replica=False)
job_45.add_args('blastall', 'blastall_00000045', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000045_output.txt')
wf.add_jobs(job_45)

if 'job_45' in task_output_files:
  job_3.add_inputs(task_output_files['job_45'])
wf.add_dependency(job_3, parents=[job_45])

if 'job_45' in task_output_files:
  job_4.add_inputs(task_output_files['job_45'])
wf.add_dependency(job_4, parents=[job_45])

if 'job_1' in task_output_files:
  job_45.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_45, parents=[job_1])

job_46 = Job('blastall', _id='blastall_00000046')
out_file_46 = File('blastall_00000046_output.txt')
task_output_files['job_46'] = out_file_46
job_46.add_outputs(out_file_46, stage_out=False, register_replica=False)
job_46.add_args('blastall', 'blastall_00000046', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000046_output.txt')
wf.add_jobs(job_46)

if 'job_46' in task_output_files:
  job_3.add_inputs(task_output_files['job_46'])
wf.add_dependency(job_3, parents=[job_46])

if 'job_46' in task_output_files:
  job_4.add_inputs(task_output_files['job_46'])
wf.add_dependency(job_4, parents=[job_46])

if 'job_1' in task_output_files:
  job_46.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_46, parents=[job_1])

job_47 = Job('blastall', _id='blastall_00000047')
out_file_47 = File('blastall_00000047_output.txt')
task_output_files['job_47'] = out_file_47
job_47.add_outputs(out_file_47, stage_out=False, register_replica=False)
job_47.add_args('blastall', 'blastall_00000047', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000047_output.txt')
wf.add_jobs(job_47)

if 'job_47' in task_output_files:
  job_3.add_inputs(task_output_files['job_47'])
wf.add_dependency(job_3, parents=[job_47])

if 'job_47' in task_output_files:
  job_4.add_inputs(task_output_files['job_47'])
wf.add_dependency(job_4, parents=[job_47])

if 'job_1' in task_output_files:
  job_47.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_47, parents=[job_1])

job_48 = Job('blastall', _id='blastall_00000048')
out_file_48 = File('blastall_00000048_output.txt')
task_output_files['job_48'] = out_file_48
job_48.add_outputs(out_file_48, stage_out=False, register_replica=False)
job_48.add_args('blastall', 'blastall_00000048', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000048_output.txt')
wf.add_jobs(job_48)

if 'job_48' in task_output_files:
  job_3.add_inputs(task_output_files['job_48'])
wf.add_dependency(job_3, parents=[job_48])

if 'job_48' in task_output_files:
  job_4.add_inputs(task_output_files['job_48'])
wf.add_dependency(job_4, parents=[job_48])

if 'job_1' in task_output_files:
  job_48.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_48, parents=[job_1])

job_49 = Job('blastall', _id='blastall_00000049')
out_file_49 = File('blastall_00000049_output.txt')
task_output_files['job_49'] = out_file_49
job_49.add_outputs(out_file_49, stage_out=False, register_replica=False)
job_49.add_args('blastall', 'blastall_00000049', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000049_output.txt')
wf.add_jobs(job_49)

if 'job_49' in task_output_files:
  job_3.add_inputs(task_output_files['job_49'])
wf.add_dependency(job_3, parents=[job_49])

if 'job_49' in task_output_files:
  job_4.add_inputs(task_output_files['job_49'])
wf.add_dependency(job_4, parents=[job_49])

if 'job_1' in task_output_files:
  job_49.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_49, parents=[job_1])

job_50 = Job('blastall', _id='blastall_00000050')
out_file_50 = File('blastall_00000050_output.txt')
task_output_files['job_50'] = out_file_50
job_50.add_outputs(out_file_50, stage_out=False, register_replica=False)
job_50.add_args('blastall', 'blastall_00000050', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000050_output.txt')
wf.add_jobs(job_50)

if 'job_50' in task_output_files:
  job_3.add_inputs(task_output_files['job_50'])
wf.add_dependency(job_3, parents=[job_50])

if 'job_50' in task_output_files:
  job_4.add_inputs(task_output_files['job_50'])
wf.add_dependency(job_4, parents=[job_50])

if 'job_1' in task_output_files:
  job_50.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_50, parents=[job_1])

job_51 = Job('blastall', _id='blastall_00000051')
out_file_51 = File('blastall_00000051_output.txt')
task_output_files['job_51'] = out_file_51
job_51.add_outputs(out_file_51, stage_out=False, register_replica=False)
job_51.add_args('blastall', 'blastall_00000051', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000051_output.txt')
wf.add_jobs(job_51)

if 'job_51' in task_output_files:
  job_3.add_inputs(task_output_files['job_51'])
wf.add_dependency(job_3, parents=[job_51])

if 'job_51' in task_output_files:
  job_4.add_inputs(task_output_files['job_51'])
wf.add_dependency(job_4, parents=[job_51])

if 'job_1' in task_output_files:
  job_51.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_51, parents=[job_1])

job_52 = Job('blastall', _id='blastall_00000052')
out_file_52 = File('blastall_00000052_output.txt')
task_output_files['job_52'] = out_file_52
job_52.add_outputs(out_file_52, stage_out=False, register_replica=False)
job_52.add_args('blastall', 'blastall_00000052', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000052_output.txt')
wf.add_jobs(job_52)

if 'job_52' in task_output_files:
  job_3.add_inputs(task_output_files['job_52'])
wf.add_dependency(job_3, parents=[job_52])

if 'job_52' in task_output_files:
  job_4.add_inputs(task_output_files['job_52'])
wf.add_dependency(job_4, parents=[job_52])

if 'job_1' in task_output_files:
  job_52.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_52, parents=[job_1])

job_53 = Job('blastall', _id='blastall_00000053')
out_file_53 = File('blastall_00000053_output.txt')
task_output_files['job_53'] = out_file_53
job_53.add_outputs(out_file_53, stage_out=False, register_replica=False)
job_53.add_args('blastall', 'blastall_00000053', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000053_output.txt')
wf.add_jobs(job_53)

if 'job_53' in task_output_files:
  job_3.add_inputs(task_output_files['job_53'])
wf.add_dependency(job_3, parents=[job_53])

if 'job_53' in task_output_files:
  job_4.add_inputs(task_output_files['job_53'])
wf.add_dependency(job_4, parents=[job_53])

if 'job_1' in task_output_files:
  job_53.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_53, parents=[job_1])

job_54 = Job('blastall', _id='blastall_00000054')
out_file_54 = File('blastall_00000054_output.txt')
task_output_files['job_54'] = out_file_54
job_54.add_outputs(out_file_54, stage_out=False, register_replica=False)
job_54.add_args('blastall', 'blastall_00000054', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000054_output.txt')
wf.add_jobs(job_54)

if 'job_54' in task_output_files:
  job_3.add_inputs(task_output_files['job_54'])
wf.add_dependency(job_3, parents=[job_54])

if 'job_54' in task_output_files:
  job_4.add_inputs(task_output_files['job_54'])
wf.add_dependency(job_4, parents=[job_54])

if 'job_1' in task_output_files:
  job_54.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_54, parents=[job_1])

job_55 = Job('blastall', _id='blastall_00000055')
out_file_55 = File('blastall_00000055_output.txt')
task_output_files['job_55'] = out_file_55
job_55.add_outputs(out_file_55, stage_out=False, register_replica=False)
job_55.add_args('blastall', 'blastall_00000055', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000055_output.txt')
wf.add_jobs(job_55)

if 'job_55' in task_output_files:
  job_3.add_inputs(task_output_files['job_55'])
wf.add_dependency(job_3, parents=[job_55])

if 'job_55' in task_output_files:
  job_4.add_inputs(task_output_files['job_55'])
wf.add_dependency(job_4, parents=[job_55])

if 'job_1' in task_output_files:
  job_55.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_55, parents=[job_1])

job_56 = Job('blastall', _id='blastall_00000056')
out_file_56 = File('blastall_00000056_output.txt')
task_output_files['job_56'] = out_file_56
job_56.add_outputs(out_file_56, stage_out=False, register_replica=False)
job_56.add_args('blastall', 'blastall_00000056', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000056_output.txt')
wf.add_jobs(job_56)

if 'job_56' in task_output_files:
  job_3.add_inputs(task_output_files['job_56'])
wf.add_dependency(job_3, parents=[job_56])

if 'job_56' in task_output_files:
  job_4.add_inputs(task_output_files['job_56'])
wf.add_dependency(job_4, parents=[job_56])

if 'job_1' in task_output_files:
  job_56.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_56, parents=[job_1])

job_57 = Job('blastall', _id='blastall_00000057')
out_file_57 = File('blastall_00000057_output.txt')
task_output_files['job_57'] = out_file_57
job_57.add_outputs(out_file_57, stage_out=False, register_replica=False)
job_57.add_args('blastall', 'blastall_00000057', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000057_output.txt')
wf.add_jobs(job_57)

if 'job_57' in task_output_files:
  job_3.add_inputs(task_output_files['job_57'])
wf.add_dependency(job_3, parents=[job_57])

if 'job_57' in task_output_files:
  job_4.add_inputs(task_output_files['job_57'])
wf.add_dependency(job_4, parents=[job_57])

if 'job_1' in task_output_files:
  job_57.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_57, parents=[job_1])

job_58 = Job('blastall', _id='blastall_00000058')
out_file_58 = File('blastall_00000058_output.txt')
task_output_files['job_58'] = out_file_58
job_58.add_outputs(out_file_58, stage_out=False, register_replica=False)
job_58.add_args('blastall', 'blastall_00000058', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000058_output.txt')
wf.add_jobs(job_58)

if 'job_58' in task_output_files:
  job_3.add_inputs(task_output_files['job_58'])
wf.add_dependency(job_3, parents=[job_58])

if 'job_58' in task_output_files:
  job_4.add_inputs(task_output_files['job_58'])
wf.add_dependency(job_4, parents=[job_58])

if 'job_1' in task_output_files:
  job_58.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_58, parents=[job_1])

job_59 = Job('blastall', _id='blastall_00000059')
out_file_59 = File('blastall_00000059_output.txt')
task_output_files['job_59'] = out_file_59
job_59.add_outputs(out_file_59, stage_out=False, register_replica=False)
job_59.add_args('blastall', 'blastall_00000059', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000059_output.txt')
wf.add_jobs(job_59)

if 'job_59' in task_output_files:
  job_3.add_inputs(task_output_files['job_59'])
wf.add_dependency(job_3, parents=[job_59])

if 'job_59' in task_output_files:
  job_4.add_inputs(task_output_files['job_59'])
wf.add_dependency(job_4, parents=[job_59])

if 'job_1' in task_output_files:
  job_59.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_59, parents=[job_1])

job_60 = Job('blastall', _id='blastall_00000060')
out_file_60 = File('blastall_00000060_output.txt')
task_output_files['job_60'] = out_file_60
job_60.add_outputs(out_file_60, stage_out=False, register_replica=False)
job_60.add_args('blastall', 'blastall_00000060', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000060_output.txt')
wf.add_jobs(job_60)

if 'job_60' in task_output_files:
  job_3.add_inputs(task_output_files['job_60'])
wf.add_dependency(job_3, parents=[job_60])

if 'job_60' in task_output_files:
  job_4.add_inputs(task_output_files['job_60'])
wf.add_dependency(job_4, parents=[job_60])

if 'job_1' in task_output_files:
  job_60.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_60, parents=[job_1])

job_61 = Job('blastall', _id='blastall_00000061')
out_file_61 = File('blastall_00000061_output.txt')
task_output_files['job_61'] = out_file_61
job_61.add_outputs(out_file_61, stage_out=False, register_replica=False)
job_61.add_args('blastall', 'blastall_00000061', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000061_output.txt')
wf.add_jobs(job_61)

if 'job_61' in task_output_files:
  job_3.add_inputs(task_output_files['job_61'])
wf.add_dependency(job_3, parents=[job_61])

if 'job_61' in task_output_files:
  job_4.add_inputs(task_output_files['job_61'])
wf.add_dependency(job_4, parents=[job_61])

if 'job_1' in task_output_files:
  job_61.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_61, parents=[job_1])

job_62 = Job('blastall', _id='blastall_00000062')
out_file_62 = File('blastall_00000062_output.txt')
task_output_files['job_62'] = out_file_62
job_62.add_outputs(out_file_62, stage_out=False, register_replica=False)
job_62.add_args('blastall', 'blastall_00000062', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000062_output.txt')
wf.add_jobs(job_62)

if 'job_62' in task_output_files:
  job_3.add_inputs(task_output_files['job_62'])
wf.add_dependency(job_3, parents=[job_62])

if 'job_62' in task_output_files:
  job_4.add_inputs(task_output_files['job_62'])
wf.add_dependency(job_4, parents=[job_62])

if 'job_1' in task_output_files:
  job_62.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_62, parents=[job_1])

job_63 = Job('blastall', _id='blastall_00000063')
out_file_63 = File('blastall_00000063_output.txt')
task_output_files['job_63'] = out_file_63
job_63.add_outputs(out_file_63, stage_out=False, register_replica=False)
job_63.add_args('blastall', 'blastall_00000063', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000063_output.txt')
wf.add_jobs(job_63)

if 'job_63' in task_output_files:
  job_3.add_inputs(task_output_files['job_63'])
wf.add_dependency(job_3, parents=[job_63])

if 'job_63' in task_output_files:
  job_4.add_inputs(task_output_files['job_63'])
wf.add_dependency(job_4, parents=[job_63])

if 'job_1' in task_output_files:
  job_63.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_63, parents=[job_1])

job_64 = Job('blastall', _id='blastall_00000064')
out_file_64 = File('blastall_00000064_output.txt')
task_output_files['job_64'] = out_file_64
job_64.add_outputs(out_file_64, stage_out=False, register_replica=False)
job_64.add_args('blastall', 'blastall_00000064', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000064_output.txt')
wf.add_jobs(job_64)

if 'job_64' in task_output_files:
  job_3.add_inputs(task_output_files['job_64'])
wf.add_dependency(job_3, parents=[job_64])

if 'job_64' in task_output_files:
  job_4.add_inputs(task_output_files['job_64'])
wf.add_dependency(job_4, parents=[job_64])

if 'job_1' in task_output_files:
  job_64.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_64, parents=[job_1])

job_65 = Job('blastall', _id='blastall_00000065')
out_file_65 = File('blastall_00000065_output.txt')
task_output_files['job_65'] = out_file_65
job_65.add_outputs(out_file_65, stage_out=False, register_replica=False)
job_65.add_args('blastall', 'blastall_00000065', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000065_output.txt')
wf.add_jobs(job_65)

if 'job_65' in task_output_files:
  job_3.add_inputs(task_output_files['job_65'])
wf.add_dependency(job_3, parents=[job_65])

if 'job_65' in task_output_files:
  job_4.add_inputs(task_output_files['job_65'])
wf.add_dependency(job_4, parents=[job_65])

if 'job_1' in task_output_files:
  job_65.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_65, parents=[job_1])

job_66 = Job('blastall', _id='blastall_00000066')
out_file_66 = File('blastall_00000066_output.txt')
task_output_files['job_66'] = out_file_66
job_66.add_outputs(out_file_66, stage_out=False, register_replica=False)
job_66.add_args('blastall', 'blastall_00000066', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000066_output.txt')
wf.add_jobs(job_66)

if 'job_66' in task_output_files:
  job_3.add_inputs(task_output_files['job_66'])
wf.add_dependency(job_3, parents=[job_66])

if 'job_66' in task_output_files:
  job_4.add_inputs(task_output_files['job_66'])
wf.add_dependency(job_4, parents=[job_66])

if 'job_1' in task_output_files:
  job_66.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_66, parents=[job_1])

job_67 = Job('blastall', _id='blastall_00000067')
out_file_67 = File('blastall_00000067_output.txt')
task_output_files['job_67'] = out_file_67
job_67.add_outputs(out_file_67, stage_out=False, register_replica=False)
job_67.add_args('blastall', 'blastall_00000067', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000067_output.txt')
wf.add_jobs(job_67)

if 'job_67' in task_output_files:
  job_3.add_inputs(task_output_files['job_67'])
wf.add_dependency(job_3, parents=[job_67])

if 'job_67' in task_output_files:
  job_4.add_inputs(task_output_files['job_67'])
wf.add_dependency(job_4, parents=[job_67])

if 'job_1' in task_output_files:
  job_67.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_67, parents=[job_1])

job_68 = Job('blastall', _id='blastall_00000068')
out_file_68 = File('blastall_00000068_output.txt')
task_output_files['job_68'] = out_file_68
job_68.add_outputs(out_file_68, stage_out=False, register_replica=False)
job_68.add_args('blastall', 'blastall_00000068', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000068_output.txt')
wf.add_jobs(job_68)

if 'job_68' in task_output_files:
  job_3.add_inputs(task_output_files['job_68'])
wf.add_dependency(job_3, parents=[job_68])

if 'job_68' in task_output_files:
  job_4.add_inputs(task_output_files['job_68'])
wf.add_dependency(job_4, parents=[job_68])

if 'job_1' in task_output_files:
  job_68.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_68, parents=[job_1])

job_69 = Job('blastall', _id='blastall_00000069')
out_file_69 = File('blastall_00000069_output.txt')
task_output_files['job_69'] = out_file_69
job_69.add_outputs(out_file_69, stage_out=False, register_replica=False)
job_69.add_args('blastall', 'blastall_00000069', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000069_output.txt')
wf.add_jobs(job_69)

if 'job_69' in task_output_files:
  job_3.add_inputs(task_output_files['job_69'])
wf.add_dependency(job_3, parents=[job_69])

if 'job_69' in task_output_files:
  job_4.add_inputs(task_output_files['job_69'])
wf.add_dependency(job_4, parents=[job_69])

if 'job_1' in task_output_files:
  job_69.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_69, parents=[job_1])

job_70 = Job('blastall', _id='blastall_00000070')
out_file_70 = File('blastall_00000070_output.txt')
task_output_files['job_70'] = out_file_70
job_70.add_outputs(out_file_70, stage_out=False, register_replica=False)
job_70.add_args('blastall', 'blastall_00000070', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000070_output.txt')
wf.add_jobs(job_70)

if 'job_70' in task_output_files:
  job_3.add_inputs(task_output_files['job_70'])
wf.add_dependency(job_3, parents=[job_70])

if 'job_70' in task_output_files:
  job_4.add_inputs(task_output_files['job_70'])
wf.add_dependency(job_4, parents=[job_70])

if 'job_1' in task_output_files:
  job_70.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_70, parents=[job_1])

job_71 = Job('blastall', _id='blastall_00000071')
out_file_71 = File('blastall_00000071_output.txt')
task_output_files['job_71'] = out_file_71
job_71.add_outputs(out_file_71, stage_out=False, register_replica=False)
job_71.add_args('blastall', 'blastall_00000071', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000071_output.txt')
wf.add_jobs(job_71)

if 'job_71' in task_output_files:
  job_3.add_inputs(task_output_files['job_71'])
wf.add_dependency(job_3, parents=[job_71])

if 'job_71' in task_output_files:
  job_4.add_inputs(task_output_files['job_71'])
wf.add_dependency(job_4, parents=[job_71])

if 'job_1' in task_output_files:
  job_71.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_71, parents=[job_1])

job_72 = Job('blastall', _id='blastall_00000072')
out_file_72 = File('blastall_00000072_output.txt')
task_output_files['job_72'] = out_file_72
job_72.add_outputs(out_file_72, stage_out=False, register_replica=False)
job_72.add_args('blastall', 'blastall_00000072', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000072_output.txt')
wf.add_jobs(job_72)

if 'job_72' in task_output_files:
  job_3.add_inputs(task_output_files['job_72'])
wf.add_dependency(job_3, parents=[job_72])

if 'job_72' in task_output_files:
  job_4.add_inputs(task_output_files['job_72'])
wf.add_dependency(job_4, parents=[job_72])

if 'job_1' in task_output_files:
  job_72.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_72, parents=[job_1])

job_73 = Job('blastall', _id='blastall_00000073')
out_file_73 = File('blastall_00000073_output.txt')
task_output_files['job_73'] = out_file_73
job_73.add_outputs(out_file_73, stage_out=False, register_replica=False)
job_73.add_args('blastall', 'blastall_00000073', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000073_output.txt')
wf.add_jobs(job_73)

if 'job_73' in task_output_files:
  job_3.add_inputs(task_output_files['job_73'])
wf.add_dependency(job_3, parents=[job_73])

if 'job_73' in task_output_files:
  job_4.add_inputs(task_output_files['job_73'])
wf.add_dependency(job_4, parents=[job_73])

if 'job_1' in task_output_files:
  job_73.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_73, parents=[job_1])

job_74 = Job('blastall', _id='blastall_00000074')
out_file_74 = File('blastall_00000074_output.txt')
task_output_files['job_74'] = out_file_74
job_74.add_outputs(out_file_74, stage_out=False, register_replica=False)
job_74.add_args('blastall', 'blastall_00000074', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000074_output.txt')
wf.add_jobs(job_74)

if 'job_74' in task_output_files:
  job_3.add_inputs(task_output_files['job_74'])
wf.add_dependency(job_3, parents=[job_74])

if 'job_74' in task_output_files:
  job_4.add_inputs(task_output_files['job_74'])
wf.add_dependency(job_4, parents=[job_74])

if 'job_1' in task_output_files:
  job_74.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_74, parents=[job_1])

job_75 = Job('blastall', _id='blastall_00000075')
out_file_75 = File('blastall_00000075_output.txt')
task_output_files['job_75'] = out_file_75
job_75.add_outputs(out_file_75, stage_out=False, register_replica=False)
job_75.add_args('blastall', 'blastall_00000075', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000075_output.txt')
wf.add_jobs(job_75)

if 'job_75' in task_output_files:
  job_3.add_inputs(task_output_files['job_75'])
wf.add_dependency(job_3, parents=[job_75])

if 'job_75' in task_output_files:
  job_4.add_inputs(task_output_files['job_75'])
wf.add_dependency(job_4, parents=[job_75])

if 'job_1' in task_output_files:
  job_75.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_75, parents=[job_1])

job_76 = Job('blastall', _id='blastall_00000076')
out_file_76 = File('blastall_00000076_output.txt')
task_output_files['job_76'] = out_file_76
job_76.add_outputs(out_file_76, stage_out=False, register_replica=False)
job_76.add_args('blastall', 'blastall_00000076', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000076_output.txt')
wf.add_jobs(job_76)

if 'job_76' in task_output_files:
  job_3.add_inputs(task_output_files['job_76'])
wf.add_dependency(job_3, parents=[job_76])

if 'job_76' in task_output_files:
  job_4.add_inputs(task_output_files['job_76'])
wf.add_dependency(job_4, parents=[job_76])

if 'job_1' in task_output_files:
  job_76.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_76, parents=[job_1])

job_77 = Job('blastall', _id='blastall_00000077')
out_file_77 = File('blastall_00000077_output.txt')
task_output_files['job_77'] = out_file_77
job_77.add_outputs(out_file_77, stage_out=False, register_replica=False)
job_77.add_args('blastall', 'blastall_00000077', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000077_output.txt')
wf.add_jobs(job_77)

if 'job_77' in task_output_files:
  job_3.add_inputs(task_output_files['job_77'])
wf.add_dependency(job_3, parents=[job_77])

if 'job_77' in task_output_files:
  job_4.add_inputs(task_output_files['job_77'])
wf.add_dependency(job_4, parents=[job_77])

if 'job_1' in task_output_files:
  job_77.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_77, parents=[job_1])

job_78 = Job('blastall', _id='blastall_00000078')
out_file_78 = File('blastall_00000078_output.txt')
task_output_files['job_78'] = out_file_78
job_78.add_outputs(out_file_78, stage_out=False, register_replica=False)
job_78.add_args('blastall', 'blastall_00000078', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000078_output.txt')
wf.add_jobs(job_78)

if 'job_78' in task_output_files:
  job_3.add_inputs(task_output_files['job_78'])
wf.add_dependency(job_3, parents=[job_78])

if 'job_78' in task_output_files:
  job_4.add_inputs(task_output_files['job_78'])
wf.add_dependency(job_4, parents=[job_78])

if 'job_1' in task_output_files:
  job_78.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_78, parents=[job_1])

job_79 = Job('blastall', _id='blastall_00000079')
out_file_79 = File('blastall_00000079_output.txt')
task_output_files['job_79'] = out_file_79
job_79.add_outputs(out_file_79, stage_out=False, register_replica=False)
job_79.add_args('blastall', 'blastall_00000079', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000079_output.txt')
wf.add_jobs(job_79)

if 'job_79' in task_output_files:
  job_3.add_inputs(task_output_files['job_79'])
wf.add_dependency(job_3, parents=[job_79])

if 'job_79' in task_output_files:
  job_4.add_inputs(task_output_files['job_79'])
wf.add_dependency(job_4, parents=[job_79])

if 'job_1' in task_output_files:
  job_79.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_79, parents=[job_1])

job_80 = Job('blastall', _id='blastall_00000080')
out_file_80 = File('blastall_00000080_output.txt')
task_output_files['job_80'] = out_file_80
job_80.add_outputs(out_file_80, stage_out=False, register_replica=False)
job_80.add_args('blastall', 'blastall_00000080', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000080_output.txt')
wf.add_jobs(job_80)

if 'job_80' in task_output_files:
  job_3.add_inputs(task_output_files['job_80'])
wf.add_dependency(job_3, parents=[job_80])

if 'job_80' in task_output_files:
  job_4.add_inputs(task_output_files['job_80'])
wf.add_dependency(job_4, parents=[job_80])

if 'job_1' in task_output_files:
  job_80.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_80, parents=[job_1])

job_81 = Job('blastall', _id='blastall_00000081')
out_file_81 = File('blastall_00000081_output.txt')
task_output_files['job_81'] = out_file_81
job_81.add_outputs(out_file_81, stage_out=False, register_replica=False)
job_81.add_args('blastall', 'blastall_00000081', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000081_output.txt')
wf.add_jobs(job_81)

if 'job_81' in task_output_files:
  job_3.add_inputs(task_output_files['job_81'])
wf.add_dependency(job_3, parents=[job_81])

if 'job_81' in task_output_files:
  job_4.add_inputs(task_output_files['job_81'])
wf.add_dependency(job_4, parents=[job_81])

if 'job_1' in task_output_files:
  job_81.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_81, parents=[job_1])

job_82 = Job('blastall', _id='blastall_00000082')
out_file_82 = File('blastall_00000082_output.txt')
task_output_files['job_82'] = out_file_82
job_82.add_outputs(out_file_82, stage_out=False, register_replica=False)
job_82.add_args('blastall', 'blastall_00000082', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000082_output.txt')
wf.add_jobs(job_82)

if 'job_82' in task_output_files:
  job_3.add_inputs(task_output_files['job_82'])
wf.add_dependency(job_3, parents=[job_82])

if 'job_82' in task_output_files:
  job_4.add_inputs(task_output_files['job_82'])
wf.add_dependency(job_4, parents=[job_82])

if 'job_1' in task_output_files:
  job_82.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_82, parents=[job_1])

job_83 = Job('blastall', _id='blastall_00000083')
out_file_83 = File('blastall_00000083_output.txt')
task_output_files['job_83'] = out_file_83
job_83.add_outputs(out_file_83, stage_out=False, register_replica=False)
job_83.add_args('blastall', 'blastall_00000083', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000083_output.txt')
wf.add_jobs(job_83)

if 'job_83' in task_output_files:
  job_3.add_inputs(task_output_files['job_83'])
wf.add_dependency(job_3, parents=[job_83])

if 'job_83' in task_output_files:
  job_4.add_inputs(task_output_files['job_83'])
wf.add_dependency(job_4, parents=[job_83])

if 'job_1' in task_output_files:
  job_83.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_83, parents=[job_1])

job_84 = Job('blastall', _id='blastall_00000084')
out_file_84 = File('blastall_00000084_output.txt')
task_output_files['job_84'] = out_file_84
job_84.add_outputs(out_file_84, stage_out=False, register_replica=False)
job_84.add_args('blastall', 'blastall_00000084', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000084_output.txt')
wf.add_jobs(job_84)

if 'job_84' in task_output_files:
  job_3.add_inputs(task_output_files['job_84'])
wf.add_dependency(job_3, parents=[job_84])

if 'job_84' in task_output_files:
  job_4.add_inputs(task_output_files['job_84'])
wf.add_dependency(job_4, parents=[job_84])

if 'job_1' in task_output_files:
  job_84.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_84, parents=[job_1])

job_85 = Job('blastall', _id='blastall_00000085')
out_file_85 = File('blastall_00000085_output.txt')
task_output_files['job_85'] = out_file_85
job_85.add_outputs(out_file_85, stage_out=False, register_replica=False)
job_85.add_args('blastall', 'blastall_00000085', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000085_output.txt')
wf.add_jobs(job_85)

if 'job_85' in task_output_files:
  job_3.add_inputs(task_output_files['job_85'])
wf.add_dependency(job_3, parents=[job_85])

if 'job_85' in task_output_files:
  job_4.add_inputs(task_output_files['job_85'])
wf.add_dependency(job_4, parents=[job_85])

if 'job_1' in task_output_files:
  job_85.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_85, parents=[job_1])

job_86 = Job('blastall', _id='blastall_00000086')
out_file_86 = File('blastall_00000086_output.txt')
task_output_files['job_86'] = out_file_86
job_86.add_outputs(out_file_86, stage_out=False, register_replica=False)
job_86.add_args('blastall', 'blastall_00000086', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000086_output.txt')
wf.add_jobs(job_86)

if 'job_86' in task_output_files:
  job_3.add_inputs(task_output_files['job_86'])
wf.add_dependency(job_3, parents=[job_86])

if 'job_86' in task_output_files:
  job_4.add_inputs(task_output_files['job_86'])
wf.add_dependency(job_4, parents=[job_86])

if 'job_1' in task_output_files:
  job_86.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_86, parents=[job_1])

job_87 = Job('blastall', _id='blastall_00000087')
out_file_87 = File('blastall_00000087_output.txt')
task_output_files['job_87'] = out_file_87
job_87.add_outputs(out_file_87, stage_out=False, register_replica=False)
job_87.add_args('blastall', 'blastall_00000087', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000087_output.txt')
wf.add_jobs(job_87)

if 'job_87' in task_output_files:
  job_3.add_inputs(task_output_files['job_87'])
wf.add_dependency(job_3, parents=[job_87])

if 'job_87' in task_output_files:
  job_4.add_inputs(task_output_files['job_87'])
wf.add_dependency(job_4, parents=[job_87])

if 'job_1' in task_output_files:
  job_87.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_87, parents=[job_1])

job_88 = Job('blastall', _id='blastall_00000088')
out_file_88 = File('blastall_00000088_output.txt')
task_output_files['job_88'] = out_file_88
job_88.add_outputs(out_file_88, stage_out=False, register_replica=False)
job_88.add_args('blastall', 'blastall_00000088', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000088_output.txt')
wf.add_jobs(job_88)

if 'job_88' in task_output_files:
  job_3.add_inputs(task_output_files['job_88'])
wf.add_dependency(job_3, parents=[job_88])

if 'job_88' in task_output_files:
  job_4.add_inputs(task_output_files['job_88'])
wf.add_dependency(job_4, parents=[job_88])

if 'job_1' in task_output_files:
  job_88.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_88, parents=[job_1])

job_89 = Job('blastall', _id='blastall_00000089')
out_file_89 = File('blastall_00000089_output.txt')
task_output_files['job_89'] = out_file_89
job_89.add_outputs(out_file_89, stage_out=False, register_replica=False)
job_89.add_args('blastall', 'blastall_00000089', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000089_output.txt')
wf.add_jobs(job_89)

if 'job_89' in task_output_files:
  job_3.add_inputs(task_output_files['job_89'])
wf.add_dependency(job_3, parents=[job_89])

if 'job_89' in task_output_files:
  job_4.add_inputs(task_output_files['job_89'])
wf.add_dependency(job_4, parents=[job_89])

if 'job_1' in task_output_files:
  job_89.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_89, parents=[job_1])

job_90 = Job('blastall', _id='blastall_00000090')
out_file_90 = File('blastall_00000090_output.txt')
task_output_files['job_90'] = out_file_90
job_90.add_outputs(out_file_90, stage_out=False, register_replica=False)
job_90.add_args('blastall', 'blastall_00000090', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000090_output.txt')
wf.add_jobs(job_90)

if 'job_90' in task_output_files:
  job_3.add_inputs(task_output_files['job_90'])
wf.add_dependency(job_3, parents=[job_90])

if 'job_90' in task_output_files:
  job_4.add_inputs(task_output_files['job_90'])
wf.add_dependency(job_4, parents=[job_90])

if 'job_1' in task_output_files:
  job_90.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_90, parents=[job_1])

job_91 = Job('blastall', _id='blastall_00000091')
out_file_91 = File('blastall_00000091_output.txt')
task_output_files['job_91'] = out_file_91
job_91.add_outputs(out_file_91, stage_out=False, register_replica=False)
job_91.add_args('blastall', 'blastall_00000091', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000091_output.txt')
wf.add_jobs(job_91)

if 'job_91' in task_output_files:
  job_3.add_inputs(task_output_files['job_91'])
wf.add_dependency(job_3, parents=[job_91])

if 'job_91' in task_output_files:
  job_4.add_inputs(task_output_files['job_91'])
wf.add_dependency(job_4, parents=[job_91])

if 'job_1' in task_output_files:
  job_91.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_91, parents=[job_1])

job_92 = Job('blastall', _id='blastall_00000092')
out_file_92 = File('blastall_00000092_output.txt')
task_output_files['job_92'] = out_file_92
job_92.add_outputs(out_file_92, stage_out=False, register_replica=False)
job_92.add_args('blastall', 'blastall_00000092', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000092_output.txt')
wf.add_jobs(job_92)

if 'job_92' in task_output_files:
  job_3.add_inputs(task_output_files['job_92'])
wf.add_dependency(job_3, parents=[job_92])

if 'job_92' in task_output_files:
  job_4.add_inputs(task_output_files['job_92'])
wf.add_dependency(job_4, parents=[job_92])

if 'job_1' in task_output_files:
  job_92.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_92, parents=[job_1])

job_93 = Job('blastall', _id='blastall_00000093')
out_file_93 = File('blastall_00000093_output.txt')
task_output_files['job_93'] = out_file_93
job_93.add_outputs(out_file_93, stage_out=False, register_replica=False)
job_93.add_args('blastall', 'blastall_00000093', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000093_output.txt')
wf.add_jobs(job_93)

if 'job_93' in task_output_files:
  job_3.add_inputs(task_output_files['job_93'])
wf.add_dependency(job_3, parents=[job_93])

if 'job_93' in task_output_files:
  job_4.add_inputs(task_output_files['job_93'])
wf.add_dependency(job_4, parents=[job_93])

if 'job_1' in task_output_files:
  job_93.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_93, parents=[job_1])

job_94 = Job('blastall', _id='blastall_00000094')
out_file_94 = File('blastall_00000094_output.txt')
task_output_files['job_94'] = out_file_94
job_94.add_outputs(out_file_94, stage_out=False, register_replica=False)
job_94.add_args('blastall', 'blastall_00000094', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000094_output.txt')
wf.add_jobs(job_94)

if 'job_94' in task_output_files:
  job_3.add_inputs(task_output_files['job_94'])
wf.add_dependency(job_3, parents=[job_94])

if 'job_94' in task_output_files:
  job_4.add_inputs(task_output_files['job_94'])
wf.add_dependency(job_4, parents=[job_94])

if 'job_1' in task_output_files:
  job_94.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_94, parents=[job_1])

job_95 = Job('blastall', _id='blastall_00000095')
out_file_95 = File('blastall_00000095_output.txt')
task_output_files['job_95'] = out_file_95
job_95.add_outputs(out_file_95, stage_out=False, register_replica=False)
job_95.add_args('blastall', 'blastall_00000095', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000095_output.txt')
wf.add_jobs(job_95)

if 'job_95' in task_output_files:
  job_3.add_inputs(task_output_files['job_95'])
wf.add_dependency(job_3, parents=[job_95])

if 'job_95' in task_output_files:
  job_4.add_inputs(task_output_files['job_95'])
wf.add_dependency(job_4, parents=[job_95])

if 'job_1' in task_output_files:
  job_95.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_95, parents=[job_1])

job_96 = Job('blastall', _id='blastall_00000096')
out_file_96 = File('blastall_00000096_output.txt')
task_output_files['job_96'] = out_file_96
job_96.add_outputs(out_file_96, stage_out=False, register_replica=False)
job_96.add_args('blastall', 'blastall_00000096', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000096_output.txt')
wf.add_jobs(job_96)

if 'job_96' in task_output_files:
  job_3.add_inputs(task_output_files['job_96'])
wf.add_dependency(job_3, parents=[job_96])

if 'job_96' in task_output_files:
  job_4.add_inputs(task_output_files['job_96'])
wf.add_dependency(job_4, parents=[job_96])

if 'job_1' in task_output_files:
  job_96.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_96, parents=[job_1])

job_97 = Job('blastall', _id='blastall_00000097')
out_file_97 = File('blastall_00000097_output.txt')
task_output_files['job_97'] = out_file_97
job_97.add_outputs(out_file_97, stage_out=False, register_replica=False)
job_97.add_args('blastall', 'blastall_00000097', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000097_output.txt')
wf.add_jobs(job_97)

if 'job_97' in task_output_files:
  job_3.add_inputs(task_output_files['job_97'])
wf.add_dependency(job_3, parents=[job_97])

if 'job_97' in task_output_files:
  job_4.add_inputs(task_output_files['job_97'])
wf.add_dependency(job_4, parents=[job_97])

if 'job_1' in task_output_files:
  job_97.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_97, parents=[job_1])

job_98 = Job('blastall', _id='blastall_00000098')
out_file_98 = File('blastall_00000098_output.txt')
task_output_files['job_98'] = out_file_98
job_98.add_outputs(out_file_98, stage_out=False, register_replica=False)
job_98.add_args('blastall', 'blastall_00000098', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000098_output.txt')
wf.add_jobs(job_98)

if 'job_98' in task_output_files:
  job_3.add_inputs(task_output_files['job_98'])
wf.add_dependency(job_3, parents=[job_98])

if 'job_98' in task_output_files:
  job_4.add_inputs(task_output_files['job_98'])
wf.add_dependency(job_4, parents=[job_98])

if 'job_1' in task_output_files:
  job_98.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_98, parents=[job_1])

job_99 = Job('blastall', _id='blastall_00000099')
out_file_99 = File('blastall_00000099_output.txt')
task_output_files['job_99'] = out_file_99
job_99.add_outputs(out_file_99, stage_out=False, register_replica=False)
job_99.add_args('blastall', 'blastall_00000099', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000099_output.txt')
wf.add_jobs(job_99)

if 'job_99' in task_output_files:
  job_3.add_inputs(task_output_files['job_99'])
wf.add_dependency(job_3, parents=[job_99])

if 'job_99' in task_output_files:
  job_4.add_inputs(task_output_files['job_99'])
wf.add_dependency(job_4, parents=[job_99])

if 'job_1' in task_output_files:
  job_99.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_99, parents=[job_1])

job_100 = Job('blastall', _id='blastall_00000100')
out_file_100 = File('blastall_00000100_output.txt')
task_output_files['job_100'] = out_file_100
job_100.add_outputs(out_file_100, stage_out=False, register_replica=False)
job_100.add_args('blastall', 'blastall_00000100', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000100_output.txt')
wf.add_jobs(job_100)

if 'job_100' in task_output_files:
  job_3.add_inputs(task_output_files['job_100'])
wf.add_dependency(job_3, parents=[job_100])

if 'job_100' in task_output_files:
  job_4.add_inputs(task_output_files['job_100'])
wf.add_dependency(job_4, parents=[job_100])

if 'job_1' in task_output_files:
  job_100.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_100, parents=[job_1])

job_101 = Job('blastall', _id='blastall_00000101')
out_file_101 = File('blastall_00000101_output.txt')
task_output_files['job_101'] = out_file_101
job_101.add_outputs(out_file_101, stage_out=False, register_replica=False)
job_101.add_args('blastall', 'blastall_00000101', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000101_output.txt')
wf.add_jobs(job_101)

if 'job_101' in task_output_files:
  job_3.add_inputs(task_output_files['job_101'])
wf.add_dependency(job_3, parents=[job_101])

if 'job_101' in task_output_files:
  job_4.add_inputs(task_output_files['job_101'])
wf.add_dependency(job_4, parents=[job_101])

if 'job_1' in task_output_files:
  job_101.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_101, parents=[job_1])

job_102 = Job('blastall', _id='blastall_00000102')
out_file_102 = File('blastall_00000102_output.txt')
task_output_files['job_102'] = out_file_102
job_102.add_outputs(out_file_102, stage_out=False, register_replica=False)
job_102.add_args('blastall', 'blastall_00000102', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000102_output.txt')
wf.add_jobs(job_102)

if 'job_102' in task_output_files:
  job_3.add_inputs(task_output_files['job_102'])
wf.add_dependency(job_3, parents=[job_102])

if 'job_102' in task_output_files:
  job_4.add_inputs(task_output_files['job_102'])
wf.add_dependency(job_4, parents=[job_102])

if 'job_1' in task_output_files:
  job_102.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_102, parents=[job_1])

job_103 = Job('blastall', _id='blastall_00000103')
out_file_103 = File('blastall_00000103_output.txt')
task_output_files['job_103'] = out_file_103
job_103.add_outputs(out_file_103, stage_out=False, register_replica=False)
job_103.add_args('blastall', 'blastall_00000103', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000103_output.txt')
wf.add_jobs(job_103)

if 'job_103' in task_output_files:
  job_3.add_inputs(task_output_files['job_103'])
wf.add_dependency(job_3, parents=[job_103])

if 'job_103' in task_output_files:
  job_4.add_inputs(task_output_files['job_103'])
wf.add_dependency(job_4, parents=[job_103])

if 'job_1' in task_output_files:
  job_103.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_103, parents=[job_1])

job_104 = Job('blastall', _id='blastall_00000104')
out_file_104 = File('blastall_00000104_output.txt')
task_output_files['job_104'] = out_file_104
job_104.add_outputs(out_file_104, stage_out=False, register_replica=False)
job_104.add_args('blastall', 'blastall_00000104', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000104_output.txt')
wf.add_jobs(job_104)

if 'job_104' in task_output_files:
  job_3.add_inputs(task_output_files['job_104'])
wf.add_dependency(job_3, parents=[job_104])

if 'job_104' in task_output_files:
  job_4.add_inputs(task_output_files['job_104'])
wf.add_dependency(job_4, parents=[job_104])

if 'job_1' in task_output_files:
  job_104.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_104, parents=[job_1])

job_105 = Job('blastall', _id='blastall_00000105')
out_file_105 = File('blastall_00000105_output.txt')
task_output_files['job_105'] = out_file_105
job_105.add_outputs(out_file_105, stage_out=False, register_replica=False)
job_105.add_args('blastall', 'blastall_00000105', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000105_output.txt')
wf.add_jobs(job_105)

if 'job_105' in task_output_files:
  job_3.add_inputs(task_output_files['job_105'])
wf.add_dependency(job_3, parents=[job_105])

if 'job_105' in task_output_files:
  job_4.add_inputs(task_output_files['job_105'])
wf.add_dependency(job_4, parents=[job_105])

if 'job_1' in task_output_files:
  job_105.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_105, parents=[job_1])

job_106 = Job('blastall', _id='blastall_00000106')
out_file_106 = File('blastall_00000106_output.txt')
task_output_files['job_106'] = out_file_106
job_106.add_outputs(out_file_106, stage_out=False, register_replica=False)
job_106.add_args('blastall', 'blastall_00000106', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000106_output.txt')
wf.add_jobs(job_106)

if 'job_106' in task_output_files:
  job_3.add_inputs(task_output_files['job_106'])
wf.add_dependency(job_3, parents=[job_106])

if 'job_106' in task_output_files:
  job_4.add_inputs(task_output_files['job_106'])
wf.add_dependency(job_4, parents=[job_106])

if 'job_1' in task_output_files:
  job_106.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_106, parents=[job_1])

job_107 = Job('blastall', _id='blastall_00000107')
out_file_107 = File('blastall_00000107_output.txt')
task_output_files['job_107'] = out_file_107
job_107.add_outputs(out_file_107, stage_out=False, register_replica=False)
job_107.add_args('blastall', 'blastall_00000107', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000107_output.txt')
wf.add_jobs(job_107)

if 'job_107' in task_output_files:
  job_3.add_inputs(task_output_files['job_107'])
wf.add_dependency(job_3, parents=[job_107])

if 'job_107' in task_output_files:
  job_4.add_inputs(task_output_files['job_107'])
wf.add_dependency(job_4, parents=[job_107])

if 'job_1' in task_output_files:
  job_107.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_107, parents=[job_1])

job_108 = Job('blastall', _id='blastall_00000108')
out_file_108 = File('blastall_00000108_output.txt')
task_output_files['job_108'] = out_file_108
job_108.add_outputs(out_file_108, stage_out=False, register_replica=False)
job_108.add_args('blastall', 'blastall_00000108', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000108_output.txt')
wf.add_jobs(job_108)

if 'job_108' in task_output_files:
  job_3.add_inputs(task_output_files['job_108'])
wf.add_dependency(job_3, parents=[job_108])

if 'job_108' in task_output_files:
  job_4.add_inputs(task_output_files['job_108'])
wf.add_dependency(job_4, parents=[job_108])

if 'job_1' in task_output_files:
  job_108.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_108, parents=[job_1])

job_109 = Job('blastall', _id='blastall_00000109')
out_file_109 = File('blastall_00000109_output.txt')
task_output_files['job_109'] = out_file_109
job_109.add_outputs(out_file_109, stage_out=False, register_replica=False)
job_109.add_args('blastall', 'blastall_00000109', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000109_output.txt')
wf.add_jobs(job_109)

if 'job_109' in task_output_files:
  job_3.add_inputs(task_output_files['job_109'])
wf.add_dependency(job_3, parents=[job_109])

if 'job_109' in task_output_files:
  job_4.add_inputs(task_output_files['job_109'])
wf.add_dependency(job_4, parents=[job_109])

if 'job_1' in task_output_files:
  job_109.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_109, parents=[job_1])

job_110 = Job('blastall', _id='blastall_00000110')
out_file_110 = File('blastall_00000110_output.txt')
task_output_files['job_110'] = out_file_110
job_110.add_outputs(out_file_110, stage_out=False, register_replica=False)
job_110.add_args('blastall', 'blastall_00000110', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000110_output.txt')
wf.add_jobs(job_110)

if 'job_110' in task_output_files:
  job_3.add_inputs(task_output_files['job_110'])
wf.add_dependency(job_3, parents=[job_110])

if 'job_110' in task_output_files:
  job_4.add_inputs(task_output_files['job_110'])
wf.add_dependency(job_4, parents=[job_110])

if 'job_1' in task_output_files:
  job_110.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_110, parents=[job_1])

job_111 = Job('blastall', _id='blastall_00000111')
out_file_111 = File('blastall_00000111_output.txt')
task_output_files['job_111'] = out_file_111
job_111.add_outputs(out_file_111, stage_out=False, register_replica=False)
job_111.add_args('blastall', 'blastall_00000111', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000111_output.txt')
wf.add_jobs(job_111)

if 'job_111' in task_output_files:
  job_3.add_inputs(task_output_files['job_111'])
wf.add_dependency(job_3, parents=[job_111])

if 'job_111' in task_output_files:
  job_4.add_inputs(task_output_files['job_111'])
wf.add_dependency(job_4, parents=[job_111])

if 'job_1' in task_output_files:
  job_111.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_111, parents=[job_1])

job_112 = Job('blastall', _id='blastall_00000112')
out_file_112 = File('blastall_00000112_output.txt')
task_output_files['job_112'] = out_file_112
job_112.add_outputs(out_file_112, stage_out=False, register_replica=False)
job_112.add_args('blastall', 'blastall_00000112', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000112_output.txt')
wf.add_jobs(job_112)

if 'job_112' in task_output_files:
  job_3.add_inputs(task_output_files['job_112'])
wf.add_dependency(job_3, parents=[job_112])

if 'job_112' in task_output_files:
  job_4.add_inputs(task_output_files['job_112'])
wf.add_dependency(job_4, parents=[job_112])

if 'job_1' in task_output_files:
  job_112.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_112, parents=[job_1])

job_113 = Job('blastall', _id='blastall_00000113')
out_file_113 = File('blastall_00000113_output.txt')
task_output_files['job_113'] = out_file_113
job_113.add_outputs(out_file_113, stage_out=False, register_replica=False)
job_113.add_args('blastall', 'blastall_00000113', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000113_output.txt')
wf.add_jobs(job_113)

if 'job_113' in task_output_files:
  job_3.add_inputs(task_output_files['job_113'])
wf.add_dependency(job_3, parents=[job_113])

if 'job_113' in task_output_files:
  job_4.add_inputs(task_output_files['job_113'])
wf.add_dependency(job_4, parents=[job_113])

if 'job_1' in task_output_files:
  job_113.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_113, parents=[job_1])

job_114 = Job('blastall', _id='blastall_00000114')
out_file_114 = File('blastall_00000114_output.txt')
task_output_files['job_114'] = out_file_114
job_114.add_outputs(out_file_114, stage_out=False, register_replica=False)
job_114.add_args('blastall', 'blastall_00000114', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000114_output.txt')
wf.add_jobs(job_114)

if 'job_114' in task_output_files:
  job_3.add_inputs(task_output_files['job_114'])
wf.add_dependency(job_3, parents=[job_114])

if 'job_114' in task_output_files:
  job_4.add_inputs(task_output_files['job_114'])
wf.add_dependency(job_4, parents=[job_114])

if 'job_1' in task_output_files:
  job_114.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_114, parents=[job_1])

job_115 = Job('blastall', _id='blastall_00000115')
out_file_115 = File('blastall_00000115_output.txt')
task_output_files['job_115'] = out_file_115
job_115.add_outputs(out_file_115, stage_out=False, register_replica=False)
job_115.add_args('blastall', 'blastall_00000115', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000115_output.txt')
wf.add_jobs(job_115)

if 'job_115' in task_output_files:
  job_3.add_inputs(task_output_files['job_115'])
wf.add_dependency(job_3, parents=[job_115])

if 'job_115' in task_output_files:
  job_4.add_inputs(task_output_files['job_115'])
wf.add_dependency(job_4, parents=[job_115])

if 'job_1' in task_output_files:
  job_115.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_115, parents=[job_1])

job_116 = Job('blastall', _id='blastall_00000116')
out_file_116 = File('blastall_00000116_output.txt')
task_output_files['job_116'] = out_file_116
job_116.add_outputs(out_file_116, stage_out=False, register_replica=False)
job_116.add_args('blastall', 'blastall_00000116', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000116_output.txt')
wf.add_jobs(job_116)

if 'job_116' in task_output_files:
  job_3.add_inputs(task_output_files['job_116'])
wf.add_dependency(job_3, parents=[job_116])

if 'job_116' in task_output_files:
  job_4.add_inputs(task_output_files['job_116'])
wf.add_dependency(job_4, parents=[job_116])

if 'job_1' in task_output_files:
  job_116.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_116, parents=[job_1])

job_117 = Job('blastall', _id='blastall_00000117')
out_file_117 = File('blastall_00000117_output.txt')
task_output_files['job_117'] = out_file_117
job_117.add_outputs(out_file_117, stage_out=False, register_replica=False)
job_117.add_args('blastall', 'blastall_00000117', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000117_output.txt')
wf.add_jobs(job_117)

if 'job_117' in task_output_files:
  job_3.add_inputs(task_output_files['job_117'])
wf.add_dependency(job_3, parents=[job_117])

if 'job_117' in task_output_files:
  job_4.add_inputs(task_output_files['job_117'])
wf.add_dependency(job_4, parents=[job_117])

if 'job_1' in task_output_files:
  job_117.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_117, parents=[job_1])

job_118 = Job('blastall', _id='blastall_00000118')
out_file_118 = File('blastall_00000118_output.txt')
task_output_files['job_118'] = out_file_118
job_118.add_outputs(out_file_118, stage_out=False, register_replica=False)
job_118.add_args('blastall', 'blastall_00000118', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000118_output.txt')
wf.add_jobs(job_118)

if 'job_118' in task_output_files:
  job_3.add_inputs(task_output_files['job_118'])
wf.add_dependency(job_3, parents=[job_118])

if 'job_118' in task_output_files:
  job_4.add_inputs(task_output_files['job_118'])
wf.add_dependency(job_4, parents=[job_118])

if 'job_1' in task_output_files:
  job_118.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_118, parents=[job_1])

job_119 = Job('blastall', _id='blastall_00000119')
out_file_119 = File('blastall_00000119_output.txt')
task_output_files['job_119'] = out_file_119
job_119.add_outputs(out_file_119, stage_out=False, register_replica=False)
job_119.add_args('blastall', 'blastall_00000119', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000119_output.txt')
wf.add_jobs(job_119)

if 'job_119' in task_output_files:
  job_3.add_inputs(task_output_files['job_119'])
wf.add_dependency(job_3, parents=[job_119])

if 'job_119' in task_output_files:
  job_4.add_inputs(task_output_files['job_119'])
wf.add_dependency(job_4, parents=[job_119])

if 'job_1' in task_output_files:
  job_119.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_119, parents=[job_1])

job_120 = Job('blastall', _id='blastall_00000120')
out_file_120 = File('blastall_00000120_output.txt')
task_output_files['job_120'] = out_file_120
job_120.add_outputs(out_file_120, stage_out=False, register_replica=False)
job_120.add_args('blastall', 'blastall_00000120', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000120_output.txt')
wf.add_jobs(job_120)

if 'job_120' in task_output_files:
  job_3.add_inputs(task_output_files['job_120'])
wf.add_dependency(job_3, parents=[job_120])

if 'job_120' in task_output_files:
  job_4.add_inputs(task_output_files['job_120'])
wf.add_dependency(job_4, parents=[job_120])

if 'job_1' in task_output_files:
  job_120.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_120, parents=[job_1])

job_121 = Job('blastall', _id='blastall_00000121')
out_file_121 = File('blastall_00000121_output.txt')
task_output_files['job_121'] = out_file_121
job_121.add_outputs(out_file_121, stage_out=False, register_replica=False)
job_121.add_args('blastall', 'blastall_00000121', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000121_output.txt')
wf.add_jobs(job_121)

if 'job_121' in task_output_files:
  job_3.add_inputs(task_output_files['job_121'])
wf.add_dependency(job_3, parents=[job_121])

if 'job_121' in task_output_files:
  job_4.add_inputs(task_output_files['job_121'])
wf.add_dependency(job_4, parents=[job_121])

if 'job_1' in task_output_files:
  job_121.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_121, parents=[job_1])

job_122 = Job('blastall', _id='blastall_00000122')
out_file_122 = File('blastall_00000122_output.txt')
task_output_files['job_122'] = out_file_122
job_122.add_outputs(out_file_122, stage_out=False, register_replica=False)
job_122.add_args('blastall', 'blastall_00000122', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000122_output.txt')
wf.add_jobs(job_122)

if 'job_122' in task_output_files:
  job_3.add_inputs(task_output_files['job_122'])
wf.add_dependency(job_3, parents=[job_122])

if 'job_122' in task_output_files:
  job_4.add_inputs(task_output_files['job_122'])
wf.add_dependency(job_4, parents=[job_122])

if 'job_1' in task_output_files:
  job_122.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_122, parents=[job_1])

job_123 = Job('blastall', _id='blastall_00000123')
out_file_123 = File('blastall_00000123_output.txt')
task_output_files['job_123'] = out_file_123
job_123.add_outputs(out_file_123, stage_out=False, register_replica=False)
job_123.add_args('blastall', 'blastall_00000123', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000123_output.txt')
wf.add_jobs(job_123)

if 'job_123' in task_output_files:
  job_3.add_inputs(task_output_files['job_123'])
wf.add_dependency(job_3, parents=[job_123])

if 'job_123' in task_output_files:
  job_4.add_inputs(task_output_files['job_123'])
wf.add_dependency(job_4, parents=[job_123])

if 'job_1' in task_output_files:
  job_123.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_123, parents=[job_1])

job_124 = Job('blastall', _id='blastall_00000124')
out_file_124 = File('blastall_00000124_output.txt')
task_output_files['job_124'] = out_file_124
job_124.add_outputs(out_file_124, stage_out=False, register_replica=False)
job_124.add_args('blastall', 'blastall_00000124', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000124_output.txt')
wf.add_jobs(job_124)

if 'job_124' in task_output_files:
  job_3.add_inputs(task_output_files['job_124'])
wf.add_dependency(job_3, parents=[job_124])

if 'job_124' in task_output_files:
  job_4.add_inputs(task_output_files['job_124'])
wf.add_dependency(job_4, parents=[job_124])

if 'job_1' in task_output_files:
  job_124.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_124, parents=[job_1])

job_125 = Job('blastall', _id='blastall_00000125')
out_file_125 = File('blastall_00000125_output.txt')
task_output_files['job_125'] = out_file_125
job_125.add_outputs(out_file_125, stage_out=False, register_replica=False)
job_125.add_args('blastall', 'blastall_00000125', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000125_output.txt')
wf.add_jobs(job_125)

if 'job_125' in task_output_files:
  job_3.add_inputs(task_output_files['job_125'])
wf.add_dependency(job_3, parents=[job_125])

if 'job_125' in task_output_files:
  job_4.add_inputs(task_output_files['job_125'])
wf.add_dependency(job_4, parents=[job_125])

if 'job_1' in task_output_files:
  job_125.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_125, parents=[job_1])

job_126 = Job('blastall', _id='blastall_00000126')
out_file_126 = File('blastall_00000126_output.txt')
task_output_files['job_126'] = out_file_126
job_126.add_outputs(out_file_126, stage_out=False, register_replica=False)
job_126.add_args('blastall', 'blastall_00000126', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000126_output.txt')
wf.add_jobs(job_126)

if 'job_126' in task_output_files:
  job_3.add_inputs(task_output_files['job_126'])
wf.add_dependency(job_3, parents=[job_126])

if 'job_126' in task_output_files:
  job_4.add_inputs(task_output_files['job_126'])
wf.add_dependency(job_4, parents=[job_126])

if 'job_1' in task_output_files:
  job_126.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_126, parents=[job_1])

job_127 = Job('blastall', _id='blastall_00000127')
out_file_127 = File('blastall_00000127_output.txt')
task_output_files['job_127'] = out_file_127
job_127.add_outputs(out_file_127, stage_out=False, register_replica=False)
job_127.add_args('blastall', 'blastall_00000127', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000127_output.txt')
wf.add_jobs(job_127)

if 'job_127' in task_output_files:
  job_3.add_inputs(task_output_files['job_127'])
wf.add_dependency(job_3, parents=[job_127])

if 'job_127' in task_output_files:
  job_4.add_inputs(task_output_files['job_127'])
wf.add_dependency(job_4, parents=[job_127])

if 'job_1' in task_output_files:
  job_127.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_127, parents=[job_1])

job_128 = Job('blastall', _id='blastall_00000128')
out_file_128 = File('blastall_00000128_output.txt')
task_output_files['job_128'] = out_file_128
job_128.add_outputs(out_file_128, stage_out=False, register_replica=False)
job_128.add_args('blastall', 'blastall_00000128', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000128_output.txt')
wf.add_jobs(job_128)

if 'job_128' in task_output_files:
  job_3.add_inputs(task_output_files['job_128'])
wf.add_dependency(job_3, parents=[job_128])

if 'job_128' in task_output_files:
  job_4.add_inputs(task_output_files['job_128'])
wf.add_dependency(job_4, parents=[job_128])

if 'job_1' in task_output_files:
  job_128.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_128, parents=[job_1])

job_129 = Job('blastall', _id='blastall_00000129')
out_file_129 = File('blastall_00000129_output.txt')
task_output_files['job_129'] = out_file_129
job_129.add_outputs(out_file_129, stage_out=False, register_replica=False)
job_129.add_args('blastall', 'blastall_00000129', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000129_output.txt')
wf.add_jobs(job_129)

if 'job_129' in task_output_files:
  job_3.add_inputs(task_output_files['job_129'])
wf.add_dependency(job_3, parents=[job_129])

if 'job_129' in task_output_files:
  job_4.add_inputs(task_output_files['job_129'])
wf.add_dependency(job_4, parents=[job_129])

if 'job_1' in task_output_files:
  job_129.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_129, parents=[job_1])

job_130 = Job('blastall', _id='blastall_00000130')
out_file_130 = File('blastall_00000130_output.txt')
task_output_files['job_130'] = out_file_130
job_130.add_outputs(out_file_130, stage_out=False, register_replica=False)
job_130.add_args('blastall', 'blastall_00000130', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000130_output.txt')
wf.add_jobs(job_130)

if 'job_130' in task_output_files:
  job_3.add_inputs(task_output_files['job_130'])
wf.add_dependency(job_3, parents=[job_130])

if 'job_130' in task_output_files:
  job_4.add_inputs(task_output_files['job_130'])
wf.add_dependency(job_4, parents=[job_130])

if 'job_1' in task_output_files:
  job_130.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_130, parents=[job_1])

job_131 = Job('blastall', _id='blastall_00000131')
out_file_131 = File('blastall_00000131_output.txt')
task_output_files['job_131'] = out_file_131
job_131.add_outputs(out_file_131, stage_out=False, register_replica=False)
job_131.add_args('blastall', 'blastall_00000131', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000131_output.txt')
wf.add_jobs(job_131)

if 'job_131' in task_output_files:
  job_3.add_inputs(task_output_files['job_131'])
wf.add_dependency(job_3, parents=[job_131])

if 'job_131' in task_output_files:
  job_4.add_inputs(task_output_files['job_131'])
wf.add_dependency(job_4, parents=[job_131])

if 'job_1' in task_output_files:
  job_131.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_131, parents=[job_1])

job_132 = Job('blastall', _id='blastall_00000132')
out_file_132 = File('blastall_00000132_output.txt')
task_output_files['job_132'] = out_file_132
job_132.add_outputs(out_file_132, stage_out=False, register_replica=False)
job_132.add_args('blastall', 'blastall_00000132', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000132_output.txt')
wf.add_jobs(job_132)

if 'job_132' in task_output_files:
  job_3.add_inputs(task_output_files['job_132'])
wf.add_dependency(job_3, parents=[job_132])

if 'job_132' in task_output_files:
  job_4.add_inputs(task_output_files['job_132'])
wf.add_dependency(job_4, parents=[job_132])

if 'job_1' in task_output_files:
  job_132.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_132, parents=[job_1])

job_133 = Job('blastall', _id='blastall_00000133')
out_file_133 = File('blastall_00000133_output.txt')
task_output_files['job_133'] = out_file_133
job_133.add_outputs(out_file_133, stage_out=False, register_replica=False)
job_133.add_args('blastall', 'blastall_00000133', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000133_output.txt')
wf.add_jobs(job_133)

if 'job_133' in task_output_files:
  job_3.add_inputs(task_output_files['job_133'])
wf.add_dependency(job_3, parents=[job_133])

if 'job_133' in task_output_files:
  job_4.add_inputs(task_output_files['job_133'])
wf.add_dependency(job_4, parents=[job_133])

if 'job_1' in task_output_files:
  job_133.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_133, parents=[job_1])

job_134 = Job('blastall', _id='blastall_00000134')
out_file_134 = File('blastall_00000134_output.txt')
task_output_files['job_134'] = out_file_134
job_134.add_outputs(out_file_134, stage_out=False, register_replica=False)
job_134.add_args('blastall', 'blastall_00000134', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000134_output.txt')
wf.add_jobs(job_134)

if 'job_134' in task_output_files:
  job_3.add_inputs(task_output_files['job_134'])
wf.add_dependency(job_3, parents=[job_134])

if 'job_134' in task_output_files:
  job_4.add_inputs(task_output_files['job_134'])
wf.add_dependency(job_4, parents=[job_134])

if 'job_1' in task_output_files:
  job_134.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_134, parents=[job_1])

job_135 = Job('blastall', _id='blastall_00000135')
out_file_135 = File('blastall_00000135_output.txt')
task_output_files['job_135'] = out_file_135
job_135.add_outputs(out_file_135, stage_out=False, register_replica=False)
job_135.add_args('blastall', 'blastall_00000135', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000135_output.txt')
wf.add_jobs(job_135)

if 'job_135' in task_output_files:
  job_3.add_inputs(task_output_files['job_135'])
wf.add_dependency(job_3, parents=[job_135])

if 'job_135' in task_output_files:
  job_4.add_inputs(task_output_files['job_135'])
wf.add_dependency(job_4, parents=[job_135])

if 'job_1' in task_output_files:
  job_135.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_135, parents=[job_1])

job_136 = Job('blastall', _id='blastall_00000136')
out_file_136 = File('blastall_00000136_output.txt')
task_output_files['job_136'] = out_file_136
job_136.add_outputs(out_file_136, stage_out=False, register_replica=False)
job_136.add_args('blastall', 'blastall_00000136', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000136_output.txt')
wf.add_jobs(job_136)

if 'job_136' in task_output_files:
  job_3.add_inputs(task_output_files['job_136'])
wf.add_dependency(job_3, parents=[job_136])

if 'job_136' in task_output_files:
  job_4.add_inputs(task_output_files['job_136'])
wf.add_dependency(job_4, parents=[job_136])

if 'job_1' in task_output_files:
  job_136.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_136, parents=[job_1])

job_137 = Job('blastall', _id='blastall_00000137')
out_file_137 = File('blastall_00000137_output.txt')
task_output_files['job_137'] = out_file_137
job_137.add_outputs(out_file_137, stage_out=False, register_replica=False)
job_137.add_args('blastall', 'blastall_00000137', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000137_output.txt')
wf.add_jobs(job_137)

if 'job_137' in task_output_files:
  job_3.add_inputs(task_output_files['job_137'])
wf.add_dependency(job_3, parents=[job_137])

if 'job_137' in task_output_files:
  job_4.add_inputs(task_output_files['job_137'])
wf.add_dependency(job_4, parents=[job_137])

if 'job_1' in task_output_files:
  job_137.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_137, parents=[job_1])

job_138 = Job('blastall', _id='blastall_00000138')
out_file_138 = File('blastall_00000138_output.txt')
task_output_files['job_138'] = out_file_138
job_138.add_outputs(out_file_138, stage_out=False, register_replica=False)
job_138.add_args('blastall', 'blastall_00000138', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000138_output.txt')
wf.add_jobs(job_138)

if 'job_138' in task_output_files:
  job_3.add_inputs(task_output_files['job_138'])
wf.add_dependency(job_3, parents=[job_138])

if 'job_138' in task_output_files:
  job_4.add_inputs(task_output_files['job_138'])
wf.add_dependency(job_4, parents=[job_138])

if 'job_1' in task_output_files:
  job_138.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_138, parents=[job_1])

job_139 = Job('blastall', _id='blastall_00000139')
out_file_139 = File('blastall_00000139_output.txt')
task_output_files['job_139'] = out_file_139
job_139.add_outputs(out_file_139, stage_out=False, register_replica=False)
job_139.add_args('blastall', 'blastall_00000139', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000139_output.txt')
wf.add_jobs(job_139)

if 'job_139' in task_output_files:
  job_3.add_inputs(task_output_files['job_139'])
wf.add_dependency(job_3, parents=[job_139])

if 'job_139' in task_output_files:
  job_4.add_inputs(task_output_files['job_139'])
wf.add_dependency(job_4, parents=[job_139])

if 'job_1' in task_output_files:
  job_139.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_139, parents=[job_1])

job_140 = Job('blastall', _id='blastall_00000140')
out_file_140 = File('blastall_00000140_output.txt')
task_output_files['job_140'] = out_file_140
job_140.add_outputs(out_file_140, stage_out=False, register_replica=False)
job_140.add_args('blastall', 'blastall_00000140', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000140_output.txt')
wf.add_jobs(job_140)

if 'job_140' in task_output_files:
  job_3.add_inputs(task_output_files['job_140'])
wf.add_dependency(job_3, parents=[job_140])

if 'job_140' in task_output_files:
  job_4.add_inputs(task_output_files['job_140'])
wf.add_dependency(job_4, parents=[job_140])

if 'job_1' in task_output_files:
  job_140.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_140, parents=[job_1])

job_141 = Job('blastall', _id='blastall_00000141')
out_file_141 = File('blastall_00000141_output.txt')
task_output_files['job_141'] = out_file_141
job_141.add_outputs(out_file_141, stage_out=False, register_replica=False)
job_141.add_args('blastall', 'blastall_00000141', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000141_output.txt')
wf.add_jobs(job_141)

if 'job_141' in task_output_files:
  job_3.add_inputs(task_output_files['job_141'])
wf.add_dependency(job_3, parents=[job_141])

if 'job_141' in task_output_files:
  job_4.add_inputs(task_output_files['job_141'])
wf.add_dependency(job_4, parents=[job_141])

if 'job_1' in task_output_files:
  job_141.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_141, parents=[job_1])

job_142 = Job('blastall', _id='blastall_00000142')
out_file_142 = File('blastall_00000142_output.txt')
task_output_files['job_142'] = out_file_142
job_142.add_outputs(out_file_142, stage_out=False, register_replica=False)
job_142.add_args('blastall', 'blastall_00000142', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000142_output.txt')
wf.add_jobs(job_142)

if 'job_142' in task_output_files:
  job_3.add_inputs(task_output_files['job_142'])
wf.add_dependency(job_3, parents=[job_142])

if 'job_142' in task_output_files:
  job_4.add_inputs(task_output_files['job_142'])
wf.add_dependency(job_4, parents=[job_142])

if 'job_1' in task_output_files:
  job_142.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_142, parents=[job_1])

job_143 = Job('blastall', _id='blastall_00000143')
out_file_143 = File('blastall_00000143_output.txt')
task_output_files['job_143'] = out_file_143
job_143.add_outputs(out_file_143, stage_out=False, register_replica=False)
job_143.add_args('blastall', 'blastall_00000143', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000143_output.txt')
wf.add_jobs(job_143)

if 'job_143' in task_output_files:
  job_3.add_inputs(task_output_files['job_143'])
wf.add_dependency(job_3, parents=[job_143])

if 'job_143' in task_output_files:
  job_4.add_inputs(task_output_files['job_143'])
wf.add_dependency(job_4, parents=[job_143])

if 'job_1' in task_output_files:
  job_143.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_143, parents=[job_1])

job_144 = Job('blastall', _id='blastall_00000144')
out_file_144 = File('blastall_00000144_output.txt')
task_output_files['job_144'] = out_file_144
job_144.add_outputs(out_file_144, stage_out=False, register_replica=False)
job_144.add_args('blastall', 'blastall_00000144', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000144_output.txt')
wf.add_jobs(job_144)

if 'job_144' in task_output_files:
  job_3.add_inputs(task_output_files['job_144'])
wf.add_dependency(job_3, parents=[job_144])

if 'job_144' in task_output_files:
  job_4.add_inputs(task_output_files['job_144'])
wf.add_dependency(job_4, parents=[job_144])

if 'job_1' in task_output_files:
  job_144.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_144, parents=[job_1])

job_145 = Job('blastall', _id='blastall_00000145')
out_file_145 = File('blastall_00000145_output.txt')
task_output_files['job_145'] = out_file_145
job_145.add_outputs(out_file_145, stage_out=False, register_replica=False)
job_145.add_args('blastall', 'blastall_00000145', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000145_output.txt')
wf.add_jobs(job_145)

if 'job_145' in task_output_files:
  job_3.add_inputs(task_output_files['job_145'])
wf.add_dependency(job_3, parents=[job_145])

if 'job_145' in task_output_files:
  job_4.add_inputs(task_output_files['job_145'])
wf.add_dependency(job_4, parents=[job_145])

if 'job_1' in task_output_files:
  job_145.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_145, parents=[job_1])

job_146 = Job('blastall', _id='blastall_00000146')
out_file_146 = File('blastall_00000146_output.txt')
task_output_files['job_146'] = out_file_146
job_146.add_outputs(out_file_146, stage_out=False, register_replica=False)
job_146.add_args('blastall', 'blastall_00000146', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000146_output.txt')
wf.add_jobs(job_146)

if 'job_146' in task_output_files:
  job_3.add_inputs(task_output_files['job_146'])
wf.add_dependency(job_3, parents=[job_146])

if 'job_146' in task_output_files:
  job_4.add_inputs(task_output_files['job_146'])
wf.add_dependency(job_4, parents=[job_146])

if 'job_1' in task_output_files:
  job_146.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_146, parents=[job_1])

job_147 = Job('blastall', _id='blastall_00000147')
out_file_147 = File('blastall_00000147_output.txt')
task_output_files['job_147'] = out_file_147
job_147.add_outputs(out_file_147, stage_out=False, register_replica=False)
job_147.add_args('blastall', 'blastall_00000147', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000147_output.txt')
wf.add_jobs(job_147)

if 'job_147' in task_output_files:
  job_3.add_inputs(task_output_files['job_147'])
wf.add_dependency(job_3, parents=[job_147])

if 'job_147' in task_output_files:
  job_4.add_inputs(task_output_files['job_147'])
wf.add_dependency(job_4, parents=[job_147])

if 'job_1' in task_output_files:
  job_147.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_147, parents=[job_1])

job_148 = Job('blastall', _id='blastall_00000148')
out_file_148 = File('blastall_00000148_output.txt')
task_output_files['job_148'] = out_file_148
job_148.add_outputs(out_file_148, stage_out=False, register_replica=False)
job_148.add_args('blastall', 'blastall_00000148', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000148_output.txt')
wf.add_jobs(job_148)

if 'job_148' in task_output_files:
  job_3.add_inputs(task_output_files['job_148'])
wf.add_dependency(job_3, parents=[job_148])

if 'job_148' in task_output_files:
  job_4.add_inputs(task_output_files['job_148'])
wf.add_dependency(job_4, parents=[job_148])

if 'job_1' in task_output_files:
  job_148.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_148, parents=[job_1])

job_149 = Job('blastall', _id='blastall_00000149')
out_file_149 = File('blastall_00000149_output.txt')
task_output_files['job_149'] = out_file_149
job_149.add_outputs(out_file_149, stage_out=False, register_replica=False)
job_149.add_args('blastall', 'blastall_00000149', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000149_output.txt')
wf.add_jobs(job_149)

if 'job_149' in task_output_files:
  job_3.add_inputs(task_output_files['job_149'])
wf.add_dependency(job_3, parents=[job_149])

if 'job_149' in task_output_files:
  job_4.add_inputs(task_output_files['job_149'])
wf.add_dependency(job_4, parents=[job_149])

if 'job_1' in task_output_files:
  job_149.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_149, parents=[job_1])

job_150 = Job('blastall', _id='blastall_00000150')
out_file_150 = File('blastall_00000150_output.txt')
task_output_files['job_150'] = out_file_150
job_150.add_outputs(out_file_150, stage_out=False, register_replica=False)
job_150.add_args('blastall', 'blastall_00000150', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000150_output.txt')
wf.add_jobs(job_150)

if 'job_150' in task_output_files:
  job_3.add_inputs(task_output_files['job_150'])
wf.add_dependency(job_3, parents=[job_150])

if 'job_150' in task_output_files:
  job_4.add_inputs(task_output_files['job_150'])
wf.add_dependency(job_4, parents=[job_150])

if 'job_1' in task_output_files:
  job_150.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_150, parents=[job_1])

job_151 = Job('blastall', _id='blastall_00000151')
out_file_151 = File('blastall_00000151_output.txt')
task_output_files['job_151'] = out_file_151
job_151.add_outputs(out_file_151, stage_out=False, register_replica=False)
job_151.add_args('blastall', 'blastall_00000151', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000151_output.txt')
wf.add_jobs(job_151)

if 'job_151' in task_output_files:
  job_3.add_inputs(task_output_files['job_151'])
wf.add_dependency(job_3, parents=[job_151])

if 'job_151' in task_output_files:
  job_4.add_inputs(task_output_files['job_151'])
wf.add_dependency(job_4, parents=[job_151])

if 'job_1' in task_output_files:
  job_151.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_151, parents=[job_1])

job_152 = Job('blastall', _id='blastall_00000152')
out_file_152 = File('blastall_00000152_output.txt')
task_output_files['job_152'] = out_file_152
job_152.add_outputs(out_file_152, stage_out=False, register_replica=False)
job_152.add_args('blastall', 'blastall_00000152', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000152_output.txt')
wf.add_jobs(job_152)

if 'job_152' in task_output_files:
  job_3.add_inputs(task_output_files['job_152'])
wf.add_dependency(job_3, parents=[job_152])

if 'job_152' in task_output_files:
  job_4.add_inputs(task_output_files['job_152'])
wf.add_dependency(job_4, parents=[job_152])

if 'job_1' in task_output_files:
  job_152.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_152, parents=[job_1])

job_153 = Job('blastall', _id='blastall_00000153')
out_file_153 = File('blastall_00000153_output.txt')
task_output_files['job_153'] = out_file_153
job_153.add_outputs(out_file_153, stage_out=False, register_replica=False)
job_153.add_args('blastall', 'blastall_00000153', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000153_output.txt')
wf.add_jobs(job_153)

if 'job_153' in task_output_files:
  job_3.add_inputs(task_output_files['job_153'])
wf.add_dependency(job_3, parents=[job_153])

if 'job_153' in task_output_files:
  job_4.add_inputs(task_output_files['job_153'])
wf.add_dependency(job_4, parents=[job_153])

if 'job_1' in task_output_files:
  job_153.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_153, parents=[job_1])

job_154 = Job('blastall', _id='blastall_00000154')
out_file_154 = File('blastall_00000154_output.txt')
task_output_files['job_154'] = out_file_154
job_154.add_outputs(out_file_154, stage_out=False, register_replica=False)
job_154.add_args('blastall', 'blastall_00000154', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000154_output.txt')
wf.add_jobs(job_154)

if 'job_154' in task_output_files:
  job_3.add_inputs(task_output_files['job_154'])
wf.add_dependency(job_3, parents=[job_154])

if 'job_154' in task_output_files:
  job_4.add_inputs(task_output_files['job_154'])
wf.add_dependency(job_4, parents=[job_154])

if 'job_1' in task_output_files:
  job_154.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_154, parents=[job_1])

job_155 = Job('blastall', _id='blastall_00000155')
out_file_155 = File('blastall_00000155_output.txt')
task_output_files['job_155'] = out_file_155
job_155.add_outputs(out_file_155, stage_out=False, register_replica=False)
job_155.add_args('blastall', 'blastall_00000155', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000155_output.txt')
wf.add_jobs(job_155)

if 'job_155' in task_output_files:
  job_3.add_inputs(task_output_files['job_155'])
wf.add_dependency(job_3, parents=[job_155])

if 'job_155' in task_output_files:
  job_4.add_inputs(task_output_files['job_155'])
wf.add_dependency(job_4, parents=[job_155])

if 'job_1' in task_output_files:
  job_155.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_155, parents=[job_1])

job_156 = Job('blastall', _id='blastall_00000156')
out_file_156 = File('blastall_00000156_output.txt')
task_output_files['job_156'] = out_file_156
job_156.add_outputs(out_file_156, stage_out=False, register_replica=False)
job_156.add_args('blastall', 'blastall_00000156', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000156_output.txt')
wf.add_jobs(job_156)

if 'job_156' in task_output_files:
  job_3.add_inputs(task_output_files['job_156'])
wf.add_dependency(job_3, parents=[job_156])

if 'job_156' in task_output_files:
  job_4.add_inputs(task_output_files['job_156'])
wf.add_dependency(job_4, parents=[job_156])

if 'job_1' in task_output_files:
  job_156.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_156, parents=[job_1])

job_157 = Job('blastall', _id='blastall_00000157')
out_file_157 = File('blastall_00000157_output.txt')
task_output_files['job_157'] = out_file_157
job_157.add_outputs(out_file_157, stage_out=False, register_replica=False)
job_157.add_args('blastall', 'blastall_00000157', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000157_output.txt')
wf.add_jobs(job_157)

if 'job_157' in task_output_files:
  job_3.add_inputs(task_output_files['job_157'])
wf.add_dependency(job_3, parents=[job_157])

if 'job_157' in task_output_files:
  job_4.add_inputs(task_output_files['job_157'])
wf.add_dependency(job_4, parents=[job_157])

if 'job_1' in task_output_files:
  job_157.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_157, parents=[job_1])

job_158 = Job('blastall', _id='blastall_00000158')
out_file_158 = File('blastall_00000158_output.txt')
task_output_files['job_158'] = out_file_158
job_158.add_outputs(out_file_158, stage_out=False, register_replica=False)
job_158.add_args('blastall', 'blastall_00000158', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000158_output.txt')
wf.add_jobs(job_158)

if 'job_158' in task_output_files:
  job_3.add_inputs(task_output_files['job_158'])
wf.add_dependency(job_3, parents=[job_158])

if 'job_158' in task_output_files:
  job_4.add_inputs(task_output_files['job_158'])
wf.add_dependency(job_4, parents=[job_158])

if 'job_1' in task_output_files:
  job_158.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_158, parents=[job_1])

job_159 = Job('blastall', _id='blastall_00000159')
out_file_159 = File('blastall_00000159_output.txt')
task_output_files['job_159'] = out_file_159
job_159.add_outputs(out_file_159, stage_out=False, register_replica=False)
job_159.add_args('blastall', 'blastall_00000159', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000159_output.txt')
wf.add_jobs(job_159)

if 'job_159' in task_output_files:
  job_3.add_inputs(task_output_files['job_159'])
wf.add_dependency(job_3, parents=[job_159])

if 'job_159' in task_output_files:
  job_4.add_inputs(task_output_files['job_159'])
wf.add_dependency(job_4, parents=[job_159])

if 'job_1' in task_output_files:
  job_159.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_159, parents=[job_1])

job_160 = Job('blastall', _id='blastall_00000160')
out_file_160 = File('blastall_00000160_output.txt')
task_output_files['job_160'] = out_file_160
job_160.add_outputs(out_file_160, stage_out=False, register_replica=False)
job_160.add_args('blastall', 'blastall_00000160', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000160_output.txt')
wf.add_jobs(job_160)

if 'job_160' in task_output_files:
  job_3.add_inputs(task_output_files['job_160'])
wf.add_dependency(job_3, parents=[job_160])

if 'job_160' in task_output_files:
  job_4.add_inputs(task_output_files['job_160'])
wf.add_dependency(job_4, parents=[job_160])

if 'job_1' in task_output_files:
  job_160.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_160, parents=[job_1])

job_161 = Job('blastall', _id='blastall_00000161')
out_file_161 = File('blastall_00000161_output.txt')
task_output_files['job_161'] = out_file_161
job_161.add_outputs(out_file_161, stage_out=False, register_replica=False)
job_161.add_args('blastall', 'blastall_00000161', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000161_output.txt')
wf.add_jobs(job_161)

if 'job_161' in task_output_files:
  job_3.add_inputs(task_output_files['job_161'])
wf.add_dependency(job_3, parents=[job_161])

if 'job_161' in task_output_files:
  job_4.add_inputs(task_output_files['job_161'])
wf.add_dependency(job_4, parents=[job_161])

if 'job_1' in task_output_files:
  job_161.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_161, parents=[job_1])

job_162 = Job('blastall', _id='blastall_00000162')
out_file_162 = File('blastall_00000162_output.txt')
task_output_files['job_162'] = out_file_162
job_162.add_outputs(out_file_162, stage_out=False, register_replica=False)
job_162.add_args('blastall', 'blastall_00000162', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000162_output.txt')
wf.add_jobs(job_162)

if 'job_162' in task_output_files:
  job_3.add_inputs(task_output_files['job_162'])
wf.add_dependency(job_3, parents=[job_162])

if 'job_162' in task_output_files:
  job_4.add_inputs(task_output_files['job_162'])
wf.add_dependency(job_4, parents=[job_162])

if 'job_1' in task_output_files:
  job_162.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_162, parents=[job_1])

job_163 = Job('blastall', _id='blastall_00000163')
out_file_163 = File('blastall_00000163_output.txt')
task_output_files['job_163'] = out_file_163
job_163.add_outputs(out_file_163, stage_out=False, register_replica=False)
job_163.add_args('blastall', 'blastall_00000163', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000163_output.txt')
wf.add_jobs(job_163)

if 'job_163' in task_output_files:
  job_3.add_inputs(task_output_files['job_163'])
wf.add_dependency(job_3, parents=[job_163])

if 'job_163' in task_output_files:
  job_4.add_inputs(task_output_files['job_163'])
wf.add_dependency(job_4, parents=[job_163])

if 'job_1' in task_output_files:
  job_163.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_163, parents=[job_1])

job_164 = Job('blastall', _id='blastall_00000164')
out_file_164 = File('blastall_00000164_output.txt')
task_output_files['job_164'] = out_file_164
job_164.add_outputs(out_file_164, stage_out=False, register_replica=False)
job_164.add_args('blastall', 'blastall_00000164', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000164_output.txt')
wf.add_jobs(job_164)

if 'job_164' in task_output_files:
  job_3.add_inputs(task_output_files['job_164'])
wf.add_dependency(job_3, parents=[job_164])

if 'job_164' in task_output_files:
  job_4.add_inputs(task_output_files['job_164'])
wf.add_dependency(job_4, parents=[job_164])

if 'job_1' in task_output_files:
  job_164.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_164, parents=[job_1])

job_165 = Job('blastall', _id='blastall_00000165')
out_file_165 = File('blastall_00000165_output.txt')
task_output_files['job_165'] = out_file_165
job_165.add_outputs(out_file_165, stage_out=False, register_replica=False)
job_165.add_args('blastall', 'blastall_00000165', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000165_output.txt')
wf.add_jobs(job_165)

if 'job_165' in task_output_files:
  job_3.add_inputs(task_output_files['job_165'])
wf.add_dependency(job_3, parents=[job_165])

if 'job_165' in task_output_files:
  job_4.add_inputs(task_output_files['job_165'])
wf.add_dependency(job_4, parents=[job_165])

if 'job_1' in task_output_files:
  job_165.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_165, parents=[job_1])

job_166 = Job('blastall', _id='blastall_00000166')
out_file_166 = File('blastall_00000166_output.txt')
task_output_files['job_166'] = out_file_166
job_166.add_outputs(out_file_166, stage_out=False, register_replica=False)
job_166.add_args('blastall', 'blastall_00000166', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000166_output.txt')
wf.add_jobs(job_166)

if 'job_166' in task_output_files:
  job_3.add_inputs(task_output_files['job_166'])
wf.add_dependency(job_3, parents=[job_166])

if 'job_166' in task_output_files:
  job_4.add_inputs(task_output_files['job_166'])
wf.add_dependency(job_4, parents=[job_166])

if 'job_1' in task_output_files:
  job_166.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_166, parents=[job_1])

job_167 = Job('blastall', _id='blastall_00000167')
out_file_167 = File('blastall_00000167_output.txt')
task_output_files['job_167'] = out_file_167
job_167.add_outputs(out_file_167, stage_out=False, register_replica=False)
job_167.add_args('blastall', 'blastall_00000167', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000167_output.txt')
wf.add_jobs(job_167)

if 'job_167' in task_output_files:
  job_3.add_inputs(task_output_files['job_167'])
wf.add_dependency(job_3, parents=[job_167])

if 'job_167' in task_output_files:
  job_4.add_inputs(task_output_files['job_167'])
wf.add_dependency(job_4, parents=[job_167])

if 'job_1' in task_output_files:
  job_167.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_167, parents=[job_1])

job_168 = Job('blastall', _id='blastall_00000168')
out_file_168 = File('blastall_00000168_output.txt')
task_output_files['job_168'] = out_file_168
job_168.add_outputs(out_file_168, stage_out=False, register_replica=False)
job_168.add_args('blastall', 'blastall_00000168', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000168_output.txt')
wf.add_jobs(job_168)

if 'job_168' in task_output_files:
  job_3.add_inputs(task_output_files['job_168'])
wf.add_dependency(job_3, parents=[job_168])

if 'job_168' in task_output_files:
  job_4.add_inputs(task_output_files['job_168'])
wf.add_dependency(job_4, parents=[job_168])

if 'job_1' in task_output_files:
  job_168.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_168, parents=[job_1])

job_169 = Job('blastall', _id='blastall_00000169')
out_file_169 = File('blastall_00000169_output.txt')
task_output_files['job_169'] = out_file_169
job_169.add_outputs(out_file_169, stage_out=False, register_replica=False)
job_169.add_args('blastall', 'blastall_00000169', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000169_output.txt')
wf.add_jobs(job_169)

if 'job_169' in task_output_files:
  job_3.add_inputs(task_output_files['job_169'])
wf.add_dependency(job_3, parents=[job_169])

if 'job_169' in task_output_files:
  job_4.add_inputs(task_output_files['job_169'])
wf.add_dependency(job_4, parents=[job_169])

if 'job_1' in task_output_files:
  job_169.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_169, parents=[job_1])

job_170 = Job('blastall', _id='blastall_00000170')
out_file_170 = File('blastall_00000170_output.txt')
task_output_files['job_170'] = out_file_170
job_170.add_outputs(out_file_170, stage_out=False, register_replica=False)
job_170.add_args('blastall', 'blastall_00000170', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000170_output.txt')
wf.add_jobs(job_170)

if 'job_170' in task_output_files:
  job_3.add_inputs(task_output_files['job_170'])
wf.add_dependency(job_3, parents=[job_170])

if 'job_170' in task_output_files:
  job_4.add_inputs(task_output_files['job_170'])
wf.add_dependency(job_4, parents=[job_170])

if 'job_1' in task_output_files:
  job_170.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_170, parents=[job_1])

job_171 = Job('blastall', _id='blastall_00000171')
out_file_171 = File('blastall_00000171_output.txt')
task_output_files['job_171'] = out_file_171
job_171.add_outputs(out_file_171, stage_out=False, register_replica=False)
job_171.add_args('blastall', 'blastall_00000171', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000171_output.txt')
wf.add_jobs(job_171)

if 'job_171' in task_output_files:
  job_3.add_inputs(task_output_files['job_171'])
wf.add_dependency(job_3, parents=[job_171])

if 'job_171' in task_output_files:
  job_4.add_inputs(task_output_files['job_171'])
wf.add_dependency(job_4, parents=[job_171])

if 'job_1' in task_output_files:
  job_171.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_171, parents=[job_1])

job_172 = Job('blastall', _id='blastall_00000172')
out_file_172 = File('blastall_00000172_output.txt')
task_output_files['job_172'] = out_file_172
job_172.add_outputs(out_file_172, stage_out=False, register_replica=False)
job_172.add_args('blastall', 'blastall_00000172', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000172_output.txt')
wf.add_jobs(job_172)

if 'job_172' in task_output_files:
  job_3.add_inputs(task_output_files['job_172'])
wf.add_dependency(job_3, parents=[job_172])

if 'job_172' in task_output_files:
  job_4.add_inputs(task_output_files['job_172'])
wf.add_dependency(job_4, parents=[job_172])

if 'job_1' in task_output_files:
  job_172.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_172, parents=[job_1])

job_173 = Job('blastall', _id='blastall_00000173')
out_file_173 = File('blastall_00000173_output.txt')
task_output_files['job_173'] = out_file_173
job_173.add_outputs(out_file_173, stage_out=False, register_replica=False)
job_173.add_args('blastall', 'blastall_00000173', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000173_output.txt')
wf.add_jobs(job_173)

if 'job_173' in task_output_files:
  job_3.add_inputs(task_output_files['job_173'])
wf.add_dependency(job_3, parents=[job_173])

if 'job_173' in task_output_files:
  job_4.add_inputs(task_output_files['job_173'])
wf.add_dependency(job_4, parents=[job_173])

if 'job_1' in task_output_files:
  job_173.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_173, parents=[job_1])

job_174 = Job('blastall', _id='blastall_00000174')
out_file_174 = File('blastall_00000174_output.txt')
task_output_files['job_174'] = out_file_174
job_174.add_outputs(out_file_174, stage_out=False, register_replica=False)
job_174.add_args('blastall', 'blastall_00000174', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000174_output.txt')
wf.add_jobs(job_174)

if 'job_174' in task_output_files:
  job_3.add_inputs(task_output_files['job_174'])
wf.add_dependency(job_3, parents=[job_174])

if 'job_174' in task_output_files:
  job_4.add_inputs(task_output_files['job_174'])
wf.add_dependency(job_4, parents=[job_174])

if 'job_1' in task_output_files:
  job_174.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_174, parents=[job_1])

job_175 = Job('blastall', _id='blastall_00000175')
out_file_175 = File('blastall_00000175_output.txt')
task_output_files['job_175'] = out_file_175
job_175.add_outputs(out_file_175, stage_out=False, register_replica=False)
job_175.add_args('blastall', 'blastall_00000175', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000175_output.txt')
wf.add_jobs(job_175)

if 'job_175' in task_output_files:
  job_3.add_inputs(task_output_files['job_175'])
wf.add_dependency(job_3, parents=[job_175])

if 'job_175' in task_output_files:
  job_4.add_inputs(task_output_files['job_175'])
wf.add_dependency(job_4, parents=[job_175])

if 'job_1' in task_output_files:
  job_175.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_175, parents=[job_1])

job_176 = Job('blastall', _id='blastall_00000176')
out_file_176 = File('blastall_00000176_output.txt')
task_output_files['job_176'] = out_file_176
job_176.add_outputs(out_file_176, stage_out=False, register_replica=False)
job_176.add_args('blastall', 'blastall_00000176', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000176_output.txt')
wf.add_jobs(job_176)

if 'job_176' in task_output_files:
  job_3.add_inputs(task_output_files['job_176'])
wf.add_dependency(job_3, parents=[job_176])

if 'job_176' in task_output_files:
  job_4.add_inputs(task_output_files['job_176'])
wf.add_dependency(job_4, parents=[job_176])

if 'job_1' in task_output_files:
  job_176.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_176, parents=[job_1])

job_177 = Job('blastall', _id='blastall_00000177')
out_file_177 = File('blastall_00000177_output.txt')
task_output_files['job_177'] = out_file_177
job_177.add_outputs(out_file_177, stage_out=False, register_replica=False)
job_177.add_args('blastall', 'blastall_00000177', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000177_output.txt')
wf.add_jobs(job_177)

if 'job_177' in task_output_files:
  job_3.add_inputs(task_output_files['job_177'])
wf.add_dependency(job_3, parents=[job_177])

if 'job_177' in task_output_files:
  job_4.add_inputs(task_output_files['job_177'])
wf.add_dependency(job_4, parents=[job_177])

if 'job_1' in task_output_files:
  job_177.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_177, parents=[job_1])

job_178 = Job('blastall', _id='blastall_00000178')
out_file_178 = File('blastall_00000178_output.txt')
task_output_files['job_178'] = out_file_178
job_178.add_outputs(out_file_178, stage_out=False, register_replica=False)
job_178.add_args('blastall', 'blastall_00000178', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000178_output.txt')
wf.add_jobs(job_178)

if 'job_178' in task_output_files:
  job_3.add_inputs(task_output_files['job_178'])
wf.add_dependency(job_3, parents=[job_178])

if 'job_178' in task_output_files:
  job_4.add_inputs(task_output_files['job_178'])
wf.add_dependency(job_4, parents=[job_178])

if 'job_1' in task_output_files:
  job_178.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_178, parents=[job_1])

job_179 = Job('blastall', _id='blastall_00000179')
out_file_179 = File('blastall_00000179_output.txt')
task_output_files['job_179'] = out_file_179
job_179.add_outputs(out_file_179, stage_out=False, register_replica=False)
job_179.add_args('blastall', 'blastall_00000179', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000179_output.txt')
wf.add_jobs(job_179)

if 'job_179' in task_output_files:
  job_3.add_inputs(task_output_files['job_179'])
wf.add_dependency(job_3, parents=[job_179])

if 'job_179' in task_output_files:
  job_4.add_inputs(task_output_files['job_179'])
wf.add_dependency(job_4, parents=[job_179])

if 'job_1' in task_output_files:
  job_179.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_179, parents=[job_1])

job_180 = Job('blastall', _id='blastall_00000180')
out_file_180 = File('blastall_00000180_output.txt')
task_output_files['job_180'] = out_file_180
job_180.add_outputs(out_file_180, stage_out=False, register_replica=False)
job_180.add_args('blastall', 'blastall_00000180', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000180_output.txt')
wf.add_jobs(job_180)

if 'job_180' in task_output_files:
  job_3.add_inputs(task_output_files['job_180'])
wf.add_dependency(job_3, parents=[job_180])

if 'job_180' in task_output_files:
  job_4.add_inputs(task_output_files['job_180'])
wf.add_dependency(job_4, parents=[job_180])

if 'job_1' in task_output_files:
  job_180.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_180, parents=[job_1])

job_181 = Job('blastall', _id='blastall_00000181')
out_file_181 = File('blastall_00000181_output.txt')
task_output_files['job_181'] = out_file_181
job_181.add_outputs(out_file_181, stage_out=False, register_replica=False)
job_181.add_args('blastall', 'blastall_00000181', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000181_output.txt')
wf.add_jobs(job_181)

if 'job_181' in task_output_files:
  job_3.add_inputs(task_output_files['job_181'])
wf.add_dependency(job_3, parents=[job_181])

if 'job_181' in task_output_files:
  job_4.add_inputs(task_output_files['job_181'])
wf.add_dependency(job_4, parents=[job_181])

if 'job_1' in task_output_files:
  job_181.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_181, parents=[job_1])

job_182 = Job('blastall', _id='blastall_00000182')
out_file_182 = File('blastall_00000182_output.txt')
task_output_files['job_182'] = out_file_182
job_182.add_outputs(out_file_182, stage_out=False, register_replica=False)
job_182.add_args('blastall', 'blastall_00000182', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000182_output.txt')
wf.add_jobs(job_182)

if 'job_182' in task_output_files:
  job_3.add_inputs(task_output_files['job_182'])
wf.add_dependency(job_3, parents=[job_182])

if 'job_182' in task_output_files:
  job_4.add_inputs(task_output_files['job_182'])
wf.add_dependency(job_4, parents=[job_182])

if 'job_1' in task_output_files:
  job_182.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_182, parents=[job_1])

job_183 = Job('blastall', _id='blastall_00000183')
out_file_183 = File('blastall_00000183_output.txt')
task_output_files['job_183'] = out_file_183
job_183.add_outputs(out_file_183, stage_out=False, register_replica=False)
job_183.add_args('blastall', 'blastall_00000183', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000183_output.txt')
wf.add_jobs(job_183)

if 'job_183' in task_output_files:
  job_3.add_inputs(task_output_files['job_183'])
wf.add_dependency(job_3, parents=[job_183])

if 'job_183' in task_output_files:
  job_4.add_inputs(task_output_files['job_183'])
wf.add_dependency(job_4, parents=[job_183])

if 'job_1' in task_output_files:
  job_183.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_183, parents=[job_1])

job_184 = Job('blastall', _id='blastall_00000184')
out_file_184 = File('blastall_00000184_output.txt')
task_output_files['job_184'] = out_file_184
job_184.add_outputs(out_file_184, stage_out=False, register_replica=False)
job_184.add_args('blastall', 'blastall_00000184', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000184_output.txt')
wf.add_jobs(job_184)

if 'job_184' in task_output_files:
  job_3.add_inputs(task_output_files['job_184'])
wf.add_dependency(job_3, parents=[job_184])

if 'job_184' in task_output_files:
  job_4.add_inputs(task_output_files['job_184'])
wf.add_dependency(job_4, parents=[job_184])

if 'job_1' in task_output_files:
  job_184.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_184, parents=[job_1])

job_185 = Job('blastall', _id='blastall_00000185')
out_file_185 = File('blastall_00000185_output.txt')
task_output_files['job_185'] = out_file_185
job_185.add_outputs(out_file_185, stage_out=False, register_replica=False)
job_185.add_args('blastall', 'blastall_00000185', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000185_output.txt')
wf.add_jobs(job_185)

if 'job_185' in task_output_files:
  job_3.add_inputs(task_output_files['job_185'])
wf.add_dependency(job_3, parents=[job_185])

if 'job_185' in task_output_files:
  job_4.add_inputs(task_output_files['job_185'])
wf.add_dependency(job_4, parents=[job_185])

if 'job_1' in task_output_files:
  job_185.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_185, parents=[job_1])

job_186 = Job('blastall', _id='blastall_00000186')
out_file_186 = File('blastall_00000186_output.txt')
task_output_files['job_186'] = out_file_186
job_186.add_outputs(out_file_186, stage_out=False, register_replica=False)
job_186.add_args('blastall', 'blastall_00000186', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000186_output.txt')
wf.add_jobs(job_186)

if 'job_186' in task_output_files:
  job_3.add_inputs(task_output_files['job_186'])
wf.add_dependency(job_3, parents=[job_186])

if 'job_186' in task_output_files:
  job_4.add_inputs(task_output_files['job_186'])
wf.add_dependency(job_4, parents=[job_186])

if 'job_1' in task_output_files:
  job_186.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_186, parents=[job_1])

job_187 = Job('blastall', _id='blastall_00000187')
out_file_187 = File('blastall_00000187_output.txt')
task_output_files['job_187'] = out_file_187
job_187.add_outputs(out_file_187, stage_out=False, register_replica=False)
job_187.add_args('blastall', 'blastall_00000187', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000187_output.txt')
wf.add_jobs(job_187)

if 'job_187' in task_output_files:
  job_3.add_inputs(task_output_files['job_187'])
wf.add_dependency(job_3, parents=[job_187])

if 'job_187' in task_output_files:
  job_4.add_inputs(task_output_files['job_187'])
wf.add_dependency(job_4, parents=[job_187])

if 'job_1' in task_output_files:
  job_187.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_187, parents=[job_1])

job_188 = Job('blastall', _id='blastall_00000188')
out_file_188 = File('blastall_00000188_output.txt')
task_output_files['job_188'] = out_file_188
job_188.add_outputs(out_file_188, stage_out=False, register_replica=False)
job_188.add_args('blastall', 'blastall_00000188', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000188_output.txt')
wf.add_jobs(job_188)

if 'job_188' in task_output_files:
  job_3.add_inputs(task_output_files['job_188'])
wf.add_dependency(job_3, parents=[job_188])

if 'job_188' in task_output_files:
  job_4.add_inputs(task_output_files['job_188'])
wf.add_dependency(job_4, parents=[job_188])

if 'job_1' in task_output_files:
  job_188.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_188, parents=[job_1])

job_189 = Job('blastall', _id='blastall_00000189')
out_file_189 = File('blastall_00000189_output.txt')
task_output_files['job_189'] = out_file_189
job_189.add_outputs(out_file_189, stage_out=False, register_replica=False)
job_189.add_args('blastall', 'blastall_00000189', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000189_output.txt')
wf.add_jobs(job_189)

if 'job_189' in task_output_files:
  job_3.add_inputs(task_output_files['job_189'])
wf.add_dependency(job_3, parents=[job_189])

if 'job_189' in task_output_files:
  job_4.add_inputs(task_output_files['job_189'])
wf.add_dependency(job_4, parents=[job_189])

if 'job_1' in task_output_files:
  job_189.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_189, parents=[job_1])

job_190 = Job('blastall', _id='blastall_00000190')
out_file_190 = File('blastall_00000190_output.txt')
task_output_files['job_190'] = out_file_190
job_190.add_outputs(out_file_190, stage_out=False, register_replica=False)
job_190.add_args('blastall', 'blastall_00000190', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000190_output.txt')
wf.add_jobs(job_190)

if 'job_190' in task_output_files:
  job_3.add_inputs(task_output_files['job_190'])
wf.add_dependency(job_3, parents=[job_190])

if 'job_190' in task_output_files:
  job_4.add_inputs(task_output_files['job_190'])
wf.add_dependency(job_4, parents=[job_190])

if 'job_1' in task_output_files:
  job_190.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_190, parents=[job_1])

job_191 = Job('blastall', _id='blastall_00000191')
out_file_191 = File('blastall_00000191_output.txt')
task_output_files['job_191'] = out_file_191
job_191.add_outputs(out_file_191, stage_out=False, register_replica=False)
job_191.add_args('blastall', 'blastall_00000191', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000191_output.txt')
wf.add_jobs(job_191)

if 'job_191' in task_output_files:
  job_3.add_inputs(task_output_files['job_191'])
wf.add_dependency(job_3, parents=[job_191])

if 'job_191' in task_output_files:
  job_4.add_inputs(task_output_files['job_191'])
wf.add_dependency(job_4, parents=[job_191])

if 'job_1' in task_output_files:
  job_191.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_191, parents=[job_1])

job_192 = Job('blastall', _id='blastall_00000192')
out_file_192 = File('blastall_00000192_output.txt')
task_output_files['job_192'] = out_file_192
job_192.add_outputs(out_file_192, stage_out=False, register_replica=False)
job_192.add_args('blastall', 'blastall_00000192', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000192_output.txt')
wf.add_jobs(job_192)

if 'job_192' in task_output_files:
  job_3.add_inputs(task_output_files['job_192'])
wf.add_dependency(job_3, parents=[job_192])

if 'job_192' in task_output_files:
  job_4.add_inputs(task_output_files['job_192'])
wf.add_dependency(job_4, parents=[job_192])

if 'job_1' in task_output_files:
  job_192.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_192, parents=[job_1])

job_193 = Job('blastall', _id='blastall_00000193')
out_file_193 = File('blastall_00000193_output.txt')
task_output_files['job_193'] = out_file_193
job_193.add_outputs(out_file_193, stage_out=False, register_replica=False)
job_193.add_args('blastall', 'blastall_00000193', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000193_output.txt')
wf.add_jobs(job_193)

if 'job_193' in task_output_files:
  job_3.add_inputs(task_output_files['job_193'])
wf.add_dependency(job_3, parents=[job_193])

if 'job_193' in task_output_files:
  job_4.add_inputs(task_output_files['job_193'])
wf.add_dependency(job_4, parents=[job_193])

if 'job_1' in task_output_files:
  job_193.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_193, parents=[job_1])

job_194 = Job('blastall', _id='blastall_00000194')
out_file_194 = File('blastall_00000194_output.txt')
task_output_files['job_194'] = out_file_194
job_194.add_outputs(out_file_194, stage_out=False, register_replica=False)
job_194.add_args('blastall', 'blastall_00000194', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000194_output.txt')
wf.add_jobs(job_194)

if 'job_194' in task_output_files:
  job_3.add_inputs(task_output_files['job_194'])
wf.add_dependency(job_3, parents=[job_194])

if 'job_194' in task_output_files:
  job_4.add_inputs(task_output_files['job_194'])
wf.add_dependency(job_4, parents=[job_194])

if 'job_1' in task_output_files:
  job_194.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_194, parents=[job_1])

job_195 = Job('blastall', _id='blastall_00000195')
out_file_195 = File('blastall_00000195_output.txt')
task_output_files['job_195'] = out_file_195
job_195.add_outputs(out_file_195, stage_out=False, register_replica=False)
job_195.add_args('blastall', 'blastall_00000195', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000195_output.txt')
wf.add_jobs(job_195)

if 'job_195' in task_output_files:
  job_3.add_inputs(task_output_files['job_195'])
wf.add_dependency(job_3, parents=[job_195])

if 'job_195' in task_output_files:
  job_4.add_inputs(task_output_files['job_195'])
wf.add_dependency(job_4, parents=[job_195])

if 'job_1' in task_output_files:
  job_195.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_195, parents=[job_1])

job_196 = Job('blastall', _id='blastall_00000196')
out_file_196 = File('blastall_00000196_output.txt')
task_output_files['job_196'] = out_file_196
job_196.add_outputs(out_file_196, stage_out=False, register_replica=False)
job_196.add_args('blastall', 'blastall_00000196', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000196_output.txt')
wf.add_jobs(job_196)

if 'job_196' in task_output_files:
  job_3.add_inputs(task_output_files['job_196'])
wf.add_dependency(job_3, parents=[job_196])

if 'job_196' in task_output_files:
  job_4.add_inputs(task_output_files['job_196'])
wf.add_dependency(job_4, parents=[job_196])

if 'job_1' in task_output_files:
  job_196.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_196, parents=[job_1])

job_197 = Job('blastall', _id='blastall_00000197')
out_file_197 = File('blastall_00000197_output.txt')
task_output_files['job_197'] = out_file_197
job_197.add_outputs(out_file_197, stage_out=False, register_replica=False)
job_197.add_args('blastall', 'blastall_00000197', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000197_output.txt')
wf.add_jobs(job_197)

if 'job_197' in task_output_files:
  job_3.add_inputs(task_output_files['job_197'])
wf.add_dependency(job_3, parents=[job_197])

if 'job_197' in task_output_files:
  job_4.add_inputs(task_output_files['job_197'])
wf.add_dependency(job_4, parents=[job_197])

if 'job_1' in task_output_files:
  job_197.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_197, parents=[job_1])

job_198 = Job('blastall', _id='blastall_00000198')
out_file_198 = File('blastall_00000198_output.txt')
task_output_files['job_198'] = out_file_198
job_198.add_outputs(out_file_198, stage_out=False, register_replica=False)
job_198.add_args('blastall', 'blastall_00000198', '--path-lock=/tmp/cores.txt.lock', '--path-cores=/tmp/cores.txt', '--percent-cpu=0.5', '--cpu-work=100', '--data', '--file-size=50', '--out=blastall_00000198_output.txt')
wf.add_jobs(job_198)

if 'job_198' in task_output_files:
  job_3.add_inputs(task_output_files['job_198'])
wf.add_dependency(job_3, parents=[job_198])

if 'job_198' in task_output_files:
  job_4.add_inputs(task_output_files['job_198'])
wf.add_dependency(job_4, parents=[job_198])

if 'job_1' in task_output_files:
  job_198.add_inputs(task_output_files['job_1'])
wf.add_dependency(job_198, parents=[job_1])

in_file_199 = File('sys_input_0.txt')
rc.add_replica('local', 'sys_input_0.txt', 'file://' + os.getcwd() + '/data/sys_input_0.txt')
job_1.add_inputs(in_file_199)
print('Using input data: ' + os.getcwd() + '/data/sys_input_0.txt')

wf.add_replica_catalog(rc)
wf.add_transformation_catalog(tc)
wf.write('Blast-Benchmark-benchmark-workflow.yml')
