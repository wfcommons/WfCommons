import pathlib
import shutil

from logging import Logger
from typing import Dict, Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class PyCompssTranslator(Translator):
    """
    A WfFormat parser for creating PyCOMPSs workflow applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)
        self.parsed_tasks = []
        self.task_counter = 1
        self.output_files_map = {}

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        self.script = "\n# workflow tasks\n"
        # PyCOMPSs translator
        self._pycompss_code()

        # Generates pycompss workflow file: template + script
        with open(this_dir.joinpath("templates/pycompss_template.py")) as fp:
            run_workflow_code = fp.read()
        run_workflow_code = run_workflow_code.replace("# Generated code goes here", self.script)
        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("pycompss_workflow.py"), "w") as fp:
            fp.write(run_workflow_code)
        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

    def _pycompss_code(self) -> None:
        # GENERATES PYCOMPSS TASKS (functions)
        all_pycompss_tasks_as_functions = []
        for task in self.tasks.values():
            # @task parameters
            if len(task.input_files) > 0 and len(task.output_files) > 0:
                self.script += f"@task(filePath=FILE_INOUT)\n"
            elif len(task.input_files) > 0:
                self.script += f"@task(filePath=FILE_IN)\n"
            elif len(task.output_files) > 0:
                self.script += f"@task(filePath=FILE_OUT)\n"
            else:
                self.script += f"@task\n"
            # function name
            # function_name = f"{task.name}_{task.task_id}"
            function_name = task.name
            # function parameters
            function_parameter_names = ""
            for i in range(len(task.input_files)):
                if len(task.input_files) == 1:
                    function_parameter_names += f"file{i}"
                else:
                    if i == 0:
                        function_parameter_names += f"file{i}"
                    else:
                        function_parameter_names += f", file{i}"
            self.script += f"def {function_name}"
            self.script += "("
            self.script += function_parameter_names
            self.script += "):\n"
            # function body
            self.script += f"\tpass\n\n"
            # PYCOMPSS TASKS METHOD CALL DEFINITION
            function_parameters = ""
            for i in range(len(task.input_files)):
                if len(task.input_files) == 1:
                    function_parameters += f"{task.input_files[i].file_id}"
                else:
                    if i == 0:
                        function_parameters += f"{task.input_files[i].file_id}"
                    else:
                        function_parameters += f", {task.input_files[i].file_id}"
            all_pycompss_tasks_as_functions.append(f"{function_name}({function_parameters})")

        # INVOKE PYCOMPSS TASKS (functions)
        self.script += f"\n\ndef main_program():\n"
        for func in all_pycompss_tasks_as_functions:
            self.script += f"\t{func}\n"

        # CALL TO MAIN METHOD
        self.script += f"\n\nif __name__ == \"__main__\":\n"
        self.script += f"\tmain_program()\n"