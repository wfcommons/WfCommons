# Adding Flowcept to a new translator

1. For some translators (e.g., Dask), you can instrument the template with:

    ```
    # FLOWCEPT_INIT
    # FLOWCEPT_END
    ```

    These `# tags` will be replaced by the translator, and you must know where to place them in your template file.
    
    For some translators (e.g., Parsl), these tags will be placed within the translator.
    For some translators (e.g., Pycompss), these tags will NOT be used, but you will use just `self._flowcept_init_python()` and `self._flowcept_stop_python()`.

2. Add the following code to the `Translator.translate` method

    ```
    # generate Flowcept code
    if self.workflow.workflow_id is not None:
        run_workflow_code = run_workflow_code.replace("# FLOWCEPT_INIT",
                                                    self._flowcept_init_python(self.workflow.workflow_id,
                                                                                self.workflow.name))
        run_workflow_code = run_workflow_code.replace("# FLOWCEPT_END", self._flowcept_stop_python())
    ```

    Note that workflow_id is assigned by flowcept. So, if workflow.workflow_id is not None, flowcept is going to be used.

3. Using: 

    # create a workflow benchmark object to generate specifications based on a recipe **with flowcept.**

    `benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=45, with_flowcept=True)`
