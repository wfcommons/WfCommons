#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


class TraceError(Exception):
	"""
		Represents the most basic error in the trace package of the WorkflowHub project
	"""
	pass


class TraceNotValid(TraceError):
	"""
		Error when the trace is not valided against the provided JSON schema
	"""
	pass
