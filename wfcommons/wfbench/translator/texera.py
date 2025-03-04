import logging
from typing import Union, Optional
import pathlib
import json
import ast

from .abstract_translator import Translator
from ...common import Workflow, Task

class TexeraTranslator(Translator):
    """
    Translate a WfFormat workflow into a Texera workflow (JSON format),
    using Python UDF operators to invoke `wfbench` or other shell commands
    via subprocess.
    """

    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[logging.Logger] = None):
        super().__init__(workflow, logger)

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        1) 创建输出目录
        2) 可选：复制所需可执行文件 (wfbench等) 到 bin/
        3) 可选：生成 Root Tasks 的输入文件 data/
        4) 把 WfCommons 的 Tasks 转成 Texera Operators (Python UDF)
        5) 用 self._write_output_file(...) 把 JSON 写到 workflow-texera.json
        """
        output_folder.mkdir(parents=True, exist_ok=True)

        # 如果需要 wfbench 等可执行文件
        self._copy_binary_files(output_folder)
        # 如果需要在 data/ 下生成 root-task 的输入文件
        self._generate_input_files(output_folder)

        # 构造 Texera workflow JSON
        texera_workflow = self._build_texera_workflow(output_folder)

        # 写出到 workflow-texera.json
        json_file_path = output_folder.joinpath("workflow-texera.json")
        self._write_output_file(json.dumps(texera_workflow, indent=2), json_file_path)

        self.logger.info(f"Texera workflow has been created at {json_file_path}.")

    def _build_texera_workflow(self, output_folder: pathlib.Path) -> dict:
        """
        将 self.tasks + 依赖  => Texera operators + links
        """
        operators = []
        links = []

        # 遍历所有 task，创建一个 PythonUDF operator
        for task_id, task_obj in self.tasks.items():
            op = {
                "operatorID": task_id,
                "operatorType": "PythonUDF",  # 假设 Texera 用 "PythonUDF" 来表示 Python UDF Operator
                "inputPorts": [0],   # 举例: 默认一个输入端口
                "outputPorts": [0],  # 默认一个输出端口
                "operatorProperties": {
                    # 这里是最关键的：把 Python 脚本源码放进来
                    "code": self._generate_python_udf_code(task_obj, output_folder),
                }
            }
            operators.append(op)

            # 构造父 -> 子的 links
            parent_list = self.task_parents.get(task_id, [])
            for p in parent_list:
                link = {
                    "origin": {
                        "operatorID": p,
                        "portIndex": 0
                    },
                    "destination": {
                        "operatorID": task_id,
                        "portIndex": 0
                    }
                }
                links.append(link)

        workflow_json = {
            "operators": operators,
            "links": links
        }
        return workflow_json

    def _generate_python_udf_code(self, task: Task, output_folder: pathlib.Path) -> str:
        """
        这里写 Python UDF Operator 在 Texera 中执行时的脚本。
        目标：用 subprocess.run(...) 执行 bin/wfbench --xxx --yyy 之类的命令。
        由于 Texera 会在运行时把 code 字符串当成一个 Python 脚本执行，
        我们要在这个脚本里自己决定如何处理输入/输出, 以及如何调用命令。
        """

        # 1) 先把 Task 的 args、input_files、output_files 转成可执行的命令行
        #    可以参考 BashTranslator 里的处理方法, 这里只是简化示例:
        cmd_str = f"./bin/{task.program} "  # e.g. ./bin/wfbench

        # 遍历 task.args，把 --input-files / --output-files 里的路径改成 data/xxx
        # （如果你要和 BashTranslator 一样的处理）
        final_args = []
        for arg in task.args:
            if arg.startswith("--output-files"):
                flag, outfiles_str = arg.split(" ", 1)
                outfiles = ast.literal_eval(outfiles_str)
                # 把文件路径前加上 data/
                outfiles_dict = { f"data/{k}": v for (k,v) in outfiles.items() }
                # 转成 JSON 字符串, 并转义内部双引号
                outfiles_json = json.dumps(outfiles_dict).replace('"', '\\"')
                final_args.append(f'{flag} "{outfiles_json}"')
            elif arg.startswith("--input-files"):
                flag, infiles_str = arg.split(" ", 1)
                infiles = ast.literal_eval(infiles_str)
                # 每个文件前加 data/
                infiles_arr = [f"data/{f}" for f in infiles]
                infiles_json = json.dumps(infiles_arr).replace('"', '\\"')
                final_args.append(f'{flag} "{infiles_json}"')
            else:
                final_args.append(arg)

        cmd_str += " ".join(final_args)

        # 2) 构造 Python 脚本，使用 subprocess.run(cmd, shell=True)
        #    在 Texera 的 Python UDF 里，会把这个 "code" 当成一个大字符串脚本执行
        #    你还可以灵活加 process_tuple() / open() / stdout 等逻辑
        #    这里是最小可行示例:
        python_code = f'''
import subprocess

# Texera Python UDF 需要一个 process_tuple(input_tuple) 函数 (取决于版本，不同版本可能略有差异)
def process_tuple(tuple_):
    # 我们这里简单地执行子进程并返回原样
    cmd = "{cmd_str}"
    subprocess.run(cmd, shell=True, check=True)
    return tuple_

# 还需要一个 keep_open 或 process_init, process_close 等视 Texera 版本情况
def process_init():
    pass

def process_close():
    pass
'''
        return python_code
