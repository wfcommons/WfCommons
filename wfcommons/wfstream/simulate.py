#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simulate a synthetic instance. Not built yet.

The last step of the new-workflow flow: hand a generated WfFormat instance to a
simulator and return its prediction for the size the user asked about. Kept as a
named seam so the flow in `pipeline.py` is complete and the gap is explicit.
"""

import pathlib


def simulate(instance: pathlib.Path, **kwargs) -> dict:
    """Simulate one synthetic instance; returns the simulator's result."""
    raise NotImplementedError(
        f"simulation is not built yet; the instance is ready at {instance}"
    )
