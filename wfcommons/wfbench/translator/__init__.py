#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from .dask import DaskTranslator
from .nextflow import NextflowTranslator
from .parsl import ParslTranslator
from .pegasus import PegasusTranslator
from .swift_t import SwiftTTranslator
from .taskvine import TaskVineTranslator
from .cwl import CWLTranslator
from .bash import BashTranslator
from .texera import TexeraTranslator