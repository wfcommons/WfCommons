#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import ndcctools.taskvine as vine


# Create a new manager
m = vine.Manager(9123)
print(f"listening on port {m.port}")


# helper function to report final status of a task
def process_result(t):
    if t:
        if t.successful():
            print(f"task {t.id} done: {t.command}")
        elif t.completed():
            print(f"task {t.id} completed with an execution error, exit code {t.exit_code}")
        else:
            print(f"task {t.id} failed with status {t.result}")


def wait_for_tasks_completion():
    print("waiting for tasks to complete...")
    while not m.empty():
        t = m.wait(2)
        if t:
            process_result(t)
    print("all tasks complete!")


wfbench = m.declare_file("bin/wfbech")
cpu_bench = m.declare_file("bin/cpu-benchmark")

# Generated code goes here