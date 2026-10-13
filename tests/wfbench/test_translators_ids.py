#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2026 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Regression tests (no Docker, no workflow systems) for creating workflow benchmarks and translating
them, with the task IDs, task names, and file IDs allowed by WfFormat 1.6.
"""

import ast
import json
import pathlib
import re
import shutil
import subprocess
import sys
import types

import pytest
import yaml

from wfcommons import BlastRecipe
from wfcommons.common import File, Task, Workflow
from wfcommons.wfbench import WorkflowBenchmark
from wfcommons.wfbench import AirflowTranslator, BashTranslator, CWLTranslator, DaskTranslator, \
    MakeflowTranslator, NextflowTranslator, ParslTranslator, PegasusTranslator, PyCompssTranslator, \
    SnakemakeTranslator, StreamflowTranslator, SwiftTTranslator, TaskVineTranslator
from wfcommons.wfbench.translator.abstract_translator import Translator
from wfcommons.wfinstances import SchemaValidator

ALL_TRANSLATORS = [AirflowTranslator, BashTranslator, CWLTranslator, DaskTranslator, MakeflowTranslator,
                   NextflowTranslator, ParslTranslator, PegasusTranslator, PyCompssTranslator,
                   SnakemakeTranslator, StreamflowTranslator, SwiftTTranslator, TaskVineTranslator]

# task IDs that are valid in WfFormat 1.6 but not identifiers in the target languages
TRICKY_IDS = ["#1", "-x", "a:b", "a/b", "a.b", "a-b", "m#1", "class", "2"]
TRICKY_NAMES = ["split-fasta", "proc.step", "proc.step", "it's a name", "merge all", "proc.step", "m#1", "2", "x:y"]
# the first task is the parent, and the last task is the child, of all other tasks
TRICKY_EDGES = [(0, i) for i in range(1, len(TRICKY_IDS) - 1)] + \
               [(i, len(TRICKY_IDS) - 1) for i in range(1, len(TRICKY_IDS) - 1)]


@pytest.fixture(autouse=True)
def _no_wfbench_binary(monkeypatch):
    # translations copy the wfbench executable, which may not be installed
    monkeypatch.setattr(Translator, "_copy_binary_files",
                        lambda self, output_folder: output_folder.joinpath("bin").mkdir(exist_ok=True))


def _build_workflow(ids, names, edges, with_files=True, file_ids=None, runtimes=None) -> Workflow:
    """Build a workflow in which each task produces one file, used by all its children."""
    workflow = Workflow(name="tricky")
    file_ids = file_ids or ["workflow_input"] + [f"f{i}" for i in range(len(ids))]
    for i, (task_id, name) in enumerate(zip(ids, names)):
        input_files, output_files = [], []
        if with_files:
            parents = [parent for parent, child in edges if child == i]
            input_files = [File(file_ids[parent + 1], 100 + parent) for parent in parents] or [File(file_ids[0], 1000)]
            output_files = [File(file_ids[i + 1], 100 + i)]
        workflow.add_task(Task(name=name, task_id=task_id, runtime=runtimes[i] if runtimes else 10.0 + i,
                               input_files=input_files, output_files=output_files,
                               bytes_read=100, bytes_written=200, avg_cpu=50.0, memory=1024))
    for parent, child in edges:
        workflow.add_dependency(ids[parent], ids[child])
    return workflow


def _create_benchmark(tmp_path: pathlib.Path, workflow: Workflow, **kwargs) -> (WorkflowBenchmark, dict):
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=len(workflow.tasks))
    kwargs.setdefault("cpu_work", 100)
    json_path = benchmark.create_benchmark_from_synthetic_workflow(tmp_path / "benchmark", workflow, **kwargs)
    with open(json_path) as f:
        data = json.load(f)
    # the generated benchmark is a valid WfFormat instance
    SchemaValidator().validate_instance(data)
    return benchmark, data


def _arguments(data: dict) -> dict:
    return {t["id"]: t["command"]["arguments"] for t in data["workflow"]["execution"]["tasks"]}


class TestWorkflowBenchmark:

    @pytest.mark.unit
    def test_argument_order(self, tmp_path) -> None:
        workflow = _build_workflow(["t1", "t2"], ["a", "b"], [(0, 1)])
        _, data = _create_benchmark(tmp_path, workflow, time=30, num_chunks=7)
        for args in _arguments(data).values():
            assert "--num-chunks 7" in args
            assert "--time 30" in args

        workflow = _build_workflow(["t1", "t2"], ["a", "b"], [(0, 1)])
        _, data = _create_benchmark(tmp_path / "no_time", workflow)
        for args in _arguments(data).values():
            assert "--num-chunks 10" in args
            assert not any(a.startswith("--time") or "None" in a for a in args)

    @pytest.mark.unit
    def test_file_renaming(self, tmp_path) -> None:
        # file IDs that are substrings of task IDs and of other file IDs
        ids = ["s:1", "p/1", "p/2", "m#1"]
        edges = [(0, 1), (0, 2), (1, 3), (2, 3)]
        workflow = _build_workflow(ids, ["s", "p", "p", "m"], edges, file_ids=["a", "1", "11", "a1", "1a"])
        _, data = _create_benchmark(tmp_path, workflow)
        arguments = _arguments(data)
        tasks = {t["id"]: t for t in data["workflow"]["specification"]["tasks"]}

        for task_id, args in arguments.items():
            # the task name argument is not modified
            assert args[0] == f"--name {task_id}"
            # file arguments match the (renamed) files of the task
            output_files = ast.literal_eval(next(a for a in args if a.startswith("--output-files"))[15:])
            input_files = ast.literal_eval(next(a for a in args if a.startswith("--input-files"))[14:])
            assert list(output_files) == tasks[task_id]["outputFiles"]
            assert input_files == tasks[task_id]["inputFiles"]
            # file names are safe to use as paths
            for file_id in tasks[task_id]["outputFiles"] + tasks[task_id]["inputFiles"]:
                assert re.fullmatch(r"[0-9A-Za-z_.-]+", file_id)

        # dependencies are still expressed by files
        for parent, child in edges:
            assert set(tasks[ids[parent]]["outputFiles"]) & set(tasks[ids[child]]["inputFiles"])

    @pytest.mark.unit
    def test_zero_runtimes(self, tmp_path) -> None:
        workflow = _build_workflow(["t1", "t2", "t3"], ["a", "b", "b"], [(0, 1), (0, 2)], runtimes=[0, 0, 0])
        _, data = _create_benchmark(tmp_path, workflow, percent_cpu={"a": 0.5, "b": 0.4})
        for args in _arguments(data).values():
            assert "--cpu-work 100" in args

    @pytest.mark.unit
    def test_execution_data_and_spec_only_workflow(self, tmp_path) -> None:
        workflow = _build_workflow(["t1", "t2"], ["a", "b"], [(0, 1)])
        workflow.include_execution = False
        _, data = _create_benchmark(tmp_path, workflow, mem=100)
        # the benchmark commands are in the execution section
        assert len(data["workflow"]["execution"]["tasks"]) == 2
        for task in data["workflow"]["execution"]["tasks"]:
            assert task["command"]["program"] == "wfbench"
            assert task["runtimeInSeconds"] == 0
            # measured data of the original workflow are not kept
            for key in ("readBytes", "writtenBytes", "avgCPU"):
                assert key not in task
            # requested resources are set by the benchmark
            assert task["coreCount"] >= 1
            assert task["memoryInBytes"] > 1024

    @pytest.mark.unit
    def test_values_per_task_name(self, tmp_path, monkeypatch) -> None:
        # recipe tasks have no category: values are given per task name
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=45)
        monkeypatch.setattr("builtins.input", lambda *args: "")
        input_file = tmp_path / "input.json"
        benchmark.generate_input_file(input_file)
        params = json.loads(input_file.read_text())
        assert set(params["cpu_work"]) == {"split_fasta", "blastall", "cat_blast", "cat"}

        json_path = benchmark.create_benchmark(tmp_path / "benchmark",
                                               percent_cpu=params["percent_cpu"],
                                               cpu_work=params["cpu_work"],
                                               gpu_work=params["gpu_work"])
        data = json.loads(json_path.read_text())
        SchemaValidator().validate_instance(data)
        for args in _arguments(data).values():
            assert "--percent-cpu 0.6" in args and "--cpu-work 1000" in args and "--gpu-work 100" in args


class TestTranslators:

    @pytest.mark.unit
    @pytest.mark.parametrize("translator_class", ALL_TRANSLATORS)
    def test_missing_command(self, translator_class) -> None:
        workflow = _build_workflow(["t1", "t2"], ["a", "b"], [(0, 1)])
        with pytest.raises(ValueError, match=r"'t1', 't2'.*WorkflowBenchmark"):
            translator_class(workflow)

    @pytest.mark.unit
    def test_unique_identifiers(self) -> None:
        names = {"k1": "valid_name", "k2": "a.b", "k3": "a-b", "k4": "a_b", "k5": "1abc", "k6": "class",
                 "k7": "", "k8": "valid_name"}
        identifiers = Translator._unique_identifiers(names, reserved={"class"})
        assert identifiers["k1"] == "valid_name"
        assert identifiers["k4"] == "a_b"
        assert identifiers["k2"] == "a_b_1" and identifiers["k3"] == "a_b_2"
        assert identifiers["k5"] == "_1abc"
        assert identifiers["k6"] == "class_1"
        assert identifiers["k7"] == "_"
        assert identifiers["k8"] == "valid_name_1"
        assert len(set(identifiers.values())) == len(names)
        # deterministic
        assert identifiers == Translator._unique_identifiers(names, reserved={"class"})

    @pytest.mark.unit
    def test_graph(self) -> None:
        workflow = BlastRecipe.from_num_tasks(45).build_workflow()
        workflow.generate_json()
        expected = workflow.workflow_json["workflow"]["specification"]["tasks"]
        translator = BashTranslator(workflow)
        assert list(translator.tasks) == [task["id"] for task in expected]
        assert translator.root_task_names == [task["id"] for task in expected if not task["parents"]]
        for task in expected:
            assert translator.task_parents[task["id"]] == task["parents"]
            assert translator.task_children[task["id"]] == task["children"]

    @pytest.mark.unit
    def test_safe_ids_are_kept(self, tmp_path) -> None:
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=45)
        benchmark.create_benchmark(tmp_path / "benchmark", cpu_work=10, data=1)
        task_ids = set(benchmark.workflow.tasks)
        assert all(task_id.isidentifier() for task_id in task_ids)

        ParslTranslator(benchmark.workflow).translate(tmp_path / "parsl")
        code = (tmp_path / "parsl/parsl_workflow.py").read_text()
        assert set(re.findall(r"^(\S+) = generic_shell_app", code, re.M)) == task_ids

        SnakemakeTranslator(benchmark.workflow).translate(tmp_path / "snakemake")
        rules = re.findall(r"^rule (\S+):", (tmp_path / "snakemake/workflow.smk").read_text(), re.M)
        assert set(rules) == task_ids | {"all_wfbench_tasks"}
        assert ".wfbench_done" not in (tmp_path / "snakemake/workflow.smk").read_text()

        NextflowTranslator(benchmark.workflow).translate(tmp_path / "nextflow")
        processes = re.findall(r"^process (\S+)\(\)", (tmp_path / "nextflow/modules/tasks_0.nf").read_text(), re.M)
        assert set(processes) == task_ids
        assert {p.name for p in (tmp_path / "nextflow/bin").glob("script_*.sh")} == \
               {f"script_{task_id}.sh" for task_id in task_ids}

    @pytest.mark.unit
    def test_python_translators(self, tmp_path) -> None:
        workflow = _build_workflow(TRICKY_IDS, TRICKY_NAMES, TRICKY_EDGES)
        benchmark, _ = _create_benchmark(tmp_path, workflow)

        outputs = {AirflowTranslator: "workflow.py", DaskTranslator: "dask_workflow.py",
                   ParslTranslator: "parsl_workflow.py", PegasusTranslator: "pegasus_workflow.py",
                   PyCompssTranslator: "pycompss_workflow.py", TaskVineTranslator: "taskvine_workflow.py"}
        for translator_class, file_name in outputs.items():
            output_folder = tmp_path / translator_class.__name__
            translator_class(benchmark.workflow).translate(output_folder)
            code = (output_folder / file_name).read_text()
            ast.parse(code)
            # task IDs starting with '#' or '-' are attached to the '--name' flag
            assert "--name=#1" in code and "--name=-x" in code

        airflow_ids = re.findall(r'task_id="([^"]+)"', (tmp_path / "AirflowTranslator/workflow.py").read_text())
        assert len(set(airflow_ids)) == len(TRICKY_IDS)
        assert all(re.fullmatch(r"[\w.-]+", task_id) for task_id in airflow_ids)

        pegasus_ids = re.findall(r"_id='([^']+)'", (tmp_path / "PegasusTranslator/pegasus_workflow.py").read_text())
        assert len(set(pegasus_ids)) == len(TRICKY_IDS)
        assert all(re.fullmatch(r"[-0-9a-zA-Z_]+", task_id) for task_id in pegasus_ids)

    @pytest.mark.unit
    def test_shell_and_file_based_translators(self, tmp_path) -> None:
        workflow = _build_workflow(TRICKY_IDS, TRICKY_NAMES, TRICKY_EDGES)
        benchmark, _ = _create_benchmark(tmp_path, workflow)

        BashTranslator(benchmark.workflow).translate(tmp_path / "bash")
        script = tmp_path / "bash/run_workflow.sh"
        assert "--name=#1" in script.read_text()
        if shutil.which("bash"):
            assert subprocess.run(["bash", "-n", str(script)]).returncode == 0

        SnakemakeTranslator(benchmark.workflow).translate(tmp_path / "snakemake")
        rules = re.findall(r"^rule (\S+):", (tmp_path / "snakemake/workflow.smk").read_text(), re.M)
        assert len(set(rules)) == len(TRICKY_IDS) + 1
        assert all(r.isidentifier() for r in rules)

        MakeflowTranslator(benchmark.workflow).translate(tmp_path / "makeflow")
        lines = (tmp_path / "makeflow/workflow.makeflow").read_text().splitlines()
        rules = [line for line in lines if line and not line.startswith(("#", "\t"))]
        assert len(rules) == len(TRICKY_IDS)
        assert all(line.count(":") == 1 and not line.startswith(" :") for line in rules)
        # '#' starts a comment anywhere in a Makeflow command, unless quoted
        commands = [line for line in lines if line.startswith("\t")]
        assert any(" '--name=#1' " in command for command in commands)
        assert any(" --name 'm#1' " in command for command in commands)

        NextflowTranslator(benchmark.workflow).translate(tmp_path / "nextflow")
        module = (tmp_path / "nextflow/modules/tasks_0.nf").read_text()
        processes = re.findall(r"^process (\S+)\(\)", module, re.M)
        assert len(set(processes)) == len(TRICKY_IDS)
        assert all(re.fullmatch(r"[A-Za-z_]\w*", p) for p in processes) and "class" not in processes
        # map entries whose key is not an identifier use the subscript notation
        assert all(re.fullmatch(r"[A-Za-z_]\w*", key) for key in re.findall(r"results\.(\S+) =", module))
        assert len(list((tmp_path / "nextflow/bin").glob("script_*.sh"))) == len(TRICKY_IDS)

        for translator_class in (CWLTranslator, StreamflowTranslator):
            output_folder = tmp_path / translator_class.__name__
            translator_class(benchmark.workflow).translate(output_folder)
            cwl = yaml.safe_load((output_folder / "main.cwl").read_text())
            assert all(re.fullmatch(r"[A-Za-z_][\w.-]*", str(step)) for step in cwl["steps"])
            sources = [s for step in cwl["steps"].values() if isinstance(step["in"], dict)
                       for s in (step["in"]["input_files"].get("source", [])
                                 if isinstance(step["in"].get("input_files"), dict) else [])]
            assert all(source.split("/")[0] in cwl["steps"] for source in sources)
        inputs = yaml.safe_load((tmp_path / "StreamflowTranslator/inputs.yml").read_text())
        assert list(inputs) == list(yaml.safe_load((tmp_path / "StreamflowTranslator/main.cwl").read_text())["inputs"])

        SwiftTTranslator(benchmark.workflow).translate(tmp_path / "swiftt")
        code = (tmp_path / "swiftt/workflow.swift").read_text()
        assert all(re.fullmatch(r"[A-Za-z_]\w*", v) for v in re.findall(r"^int (\S+)__out\[\];", code, re.M))

    @pytest.mark.unit
    def test_parsl_file_paths(self, tmp_path, monkeypatch) -> None:
        # a workflow input file ID with '/', and short file IDs that are part of other names
        ids = ["s:1", "p/1", "p/2", "m#1"]
        workflow = _build_workflow(ids, ["s", "p", "p", "m"], [(0, 1), (0, 2), (1, 3), (2, 3)],
                                   file_ids=["inputs/raw.dat", "1", "11", "a1", "1a"])
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=4)
        benchmark.workflow = workflow
        for task in workflow.tasks.values():
            outputs = {f.file_id: f.size for f in task.output_files}
            task.program = "wfbench"
            task.args = [f"--name {task.task_id}", "--cpu-work 10", f"--output-files {outputs}",
                         f"--input-files {[f.file_id for f in task.input_files]}"]
        ParslTranslator(workflow).translate(tmp_path / "parsl")
        assert (tmp_path / "parsl/data/inputs/raw.dat").exists()

        commands = self._run_parsl_workflow(tmp_path / "parsl", monkeypatch)
        assert len(commands) == 4
        folder = (tmp_path / "parsl").resolve()
        for task_id, command in zip(ids, commands):
            task = workflow.tasks[task_id]
            assert f"--name {task_id} " in command
            files = json.loads(re.search(r"--input-files '(\[.*?\])'", command).group(1))
            assert files == [str(folder / ("data" if task_id == "s:1" else "output") / f.file_id)
                             for f in task.input_files]
            files = json.loads(re.search(r"--output-files '(\{.*?\})'", command).group(1))
            assert list(files) == [str(folder / "output" / f.file_id) for f in task.output_files]

    @staticmethod
    def _run_parsl_workflow(folder: pathlib.Path, monkeypatch) -> list:
        """Run a generated Parsl workflow with a fake Parsl module that records the commands."""
        commands = []

        class FakeFile:
            def __init__(self, path):
                self.filepath = str(path)

        class FakeAppFuture:
            def __init__(self, outputs):
                self.outputs = outputs
                self.stdout = self.stderr = None

            def result(self):
                return 0

        def bash_app(function):
            def submit(cmd, inputs=None, outputs=None, stdout=None, stderr=None):
                # futures are replaced by their results (app futures) or files (data futures)
                inputs = [i.result() if isinstance(i, FakeAppFuture) else i for i in inputs or []]
                commands.append(function(cmd, inputs=inputs, outputs=outputs or []))
                return FakeAppFuture(outputs or [])
            return submit

        modules = {name: types.ModuleType(name) for name in
                   ["parsl", "parsl.app", "parsl.app.app", "parsl.config", "parsl.executors",
                    "parsl.executors.threads", "parsl.providers", "parsl.data_provider",
                    "parsl.data_provider.files"]}
        modules["parsl"].clear = modules["parsl"].load = lambda *args: None
        modules["parsl"].dfk = lambda: types.SimpleNamespace(cleanup=lambda: None)
        modules["parsl.app.app"].bash_app = bash_app
        modules["parsl.config"].Config = lambda **kwargs: None
        modules["parsl.executors"].HighThroughputExecutor = lambda **kwargs: None
        modules["parsl.executors.threads"].ThreadPoolExecutor = lambda **kwargs: None
        modules["parsl.providers"].LocalProvider = lambda **kwargs: None
        modules["parsl.data_provider.files"].File = FakeFile
        for name, module in modules.items():
            monkeypatch.setitem(sys.modules, name, module)

        monkeypatch.chdir(folder)
        code = (folder / "parsl_workflow.py").read_text()
        exec(compile(code, str(folder / "parsl_workflow.py"), "exec"), {"__name__": "__main__"})
        return commands

    @pytest.mark.unit
    def test_workflow_without_files(self, tmp_path) -> None:
        workflow = _build_workflow(["t1", "t2", "t3", "t4"], ["a", "b", "b", "c"],
                                   [(0, 1), (0, 2), (1, 3), (2, 3)], with_files=False)
        benchmark, _ = _create_benchmark(tmp_path, workflow)

        for translator_class in ALL_TRANSLATORS:
            translator_class(benchmark.workflow).translate(tmp_path / translator_class.__name__)

        # dependencies are expressed with marker files for systems that infer them from files
        smk = (tmp_path / "SnakemakeTranslator/workflow.smk").read_text()
        rule_t4 = smk[smk.index("rule t4:"):]
        assert '"data/t2.wfbench_done"' in rule_t4 and '"data/t3.wfbench_done"' in rule_t4
        makeflow = (tmp_path / "MakeflowTranslator/workflow.makeflow").read_text()
        assert "data/t4.wfbench_done : data/t2.wfbench_done data/t3.wfbench_done" in makeflow
        script = (tmp_path / "NextflowTranslator/bin/script_t4.sh").read_text()
        assert "t2.wfbench_done" in script and "t4.wfbench_done" in script

        # workflow without inputs
        cwl = yaml.safe_load((tmp_path / "CWLTranslator/main.cwl").read_text())
        assert cwl["inputs"] == []
        assert "input_files" not in cwl["steps"]["t1"]["in"]
        assert yaml.safe_load((tmp_path / "CWLTranslator/config.yml").read_text()) == {}

    @pytest.mark.unit
    def test_swift_t_multiple_outputs(self, tmp_path) -> None:
        workflow = _build_workflow(["t1", "t2", "t3"], ["a", "b", "b"], [(0, 1), (0, 2)])
        workflow.tasks["t1"].output_files.append(File("extra", 10))
        benchmark, _ = _create_benchmark(tmp_path, workflow)
        with pytest.raises(ValueError, match="multiple outputs"):
            SwiftTTranslator(benchmark.workflow)
