#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from typing import Dict, Union, Optional
from logging import Logger

from ..utils import NoValue


class MachineSystem(NoValue):
    """Machine system type."""
    LINUX = 'linux'
    MACOS = 'macos'
    WINDOWS = 'windows'


class Machine:
    """Representation of one compute machine.

    :param name: Machine node name.
    :type name: str
    :param cpu: A dictionary containing information about the CPU specification.
                Must at least contains two fields: *count* (number of CPU cores)
                and speed (CPU speed of each core in MHz).

                .. code-block:: python

                    cpu = {
                        'coreCount': 48,
                        'speedInMHz': 1200,
                        'vendor': 'Vendor Name'
                    }

    :type cpu: Dict[str, Union[int, str]]
    :param system: Machine system (linux, macos, windows).
    :type system: MachineSystem
    :param architecture: Machine architecture (e.g., x86_64, ppc).
    :type architecture: str
    :param memory: Total machine's RAM memory in bytes.
    :type memory: int
    :param release: Machine release.
    :type release: str
    :param hashcode: MD5 Hashcode for the Machine.
    :type hashcode: str
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self,
                 name: str,
                 cpu: Optional[Dict[str, Union[int, str]]] = None,
                 system: Optional[MachineSystem] = None,
                 architecture: Optional[str] = None,
                 memory: Optional[int] = None,
                 release: Optional[str] = None,
                 hashcode: Optional[str] = None,
                 logger: Optional[Logger] = None
                 ) -> None:
        """A machine from a workflow."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.name: str = name
        self.system: MachineSystem = system
        self.architecture: str = architecture
        self.memory: int = memory
        self.release: str = release
        self.hashcode = hashcode

        # the CPU information is optional in WfFormat
        cpu = cpu if cpu else {}
        self.cpu_cores: Optional[int] = cpu.get('coreCount', None)
        self.cpu_speed: int = cpu.get('speedInMHz', 0)
        self.cpu_flops: int = self.cpu_cores * self.cpu_speed * 10 ** 6 if self.cpu_cores and self.cpu_speed else 0
        self.cpu_vendor: Optional[str] = cpu.get('vendor', None)

        self.logger.debug(f"created machine: {self.name} with {self.cpu_cores} cores and {self.cpu_flops} FLOPS.")

    def as_dict(self) -> Dict[str, Union[int, str]]:
        """A JSON representation of the machine.

        :return: A JSON object representation of the machine.
        :rtype: Dict[str, Union[int, str]]
        """
        machine = {"nodeName": self.name}
        if self.system:
            machine['system'] = self.system.value
        if self.architecture:
            machine['architecture'] = self.architecture
        if self.memory:
            machine['memoryInBytes'] = self.memory
        if self.release:
            machine['release'] = self.release
        cpu = {}
        if self.cpu_cores:
            cpu['coreCount'] = self.cpu_cores
        if self.cpu_speed:
            cpu['speedInMHz'] = self.cpu_speed
        if self.cpu_vendor:
            cpu['vendor'] = self.cpu_vendor
        if cpu:
            machine['cpu'] = cpu
        return machine
