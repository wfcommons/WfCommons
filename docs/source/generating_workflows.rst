Generating Workflows
====================

Examples
--------

The following example generates a *Seismology* synthetic workflow trace
based on the number of pair of signals to estimate earthquake STFs (`num_pairs`). ::

    from workflowhub import WorkflowGenerator
    from workflowhub.generator import SeismologyRecipe

    # creating a Seismology workflow recipe based on the number
    # of pair of signals to estimate earthquake STFs
    recipe = SeismologyRecipe.from_num_pairs(num_pairs=10)

    # creating an instance of the workflow generator with the
    # Seismology workflow recipe
    generator = WorkflowGenerator(recipe)

    # generating a synthetic workflow trace of the Seismology workflow
    workflow = generator.build_workflow()

    # writing the synthetic workflow trace into a JSON file
    workflow.write_json('seismology-workflow.json')

