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

# Generated code goes here
