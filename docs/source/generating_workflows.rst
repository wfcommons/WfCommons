Generating Workflows
====================

Generating a *Seismology* synthetic workflow trace::

    from workflowhub import WorkflowGenerator
    from workflowhub.generator import SeismologyRecipe

    recipe = SeismologyRecipe.from_num_pairs(num_pairs=2)

    generator = WorkflowGenerator(recipe)
    workflow = generator.build_workflow()

    workflow.write_json('seismology-workflow.json')

