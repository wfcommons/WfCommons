#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
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
                        'count': 48,
                        'speed': 1200
                    }

    :type cpu: Dict[str, Union[int, str]]
    :param system: Machine system (linux, macos, windows).
    :type system: MachineSystem
    :param architecture: Machine architecture (e.g., x86_64, ppc).
    :type architecture: str
    :param memory: Total machine's RAM memory in KB.
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
                 cpu: Dict[str, Union[int, str]],
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

        self.cpu_cores: int = cpu['count']
        self.cpu_speed: int = cpu['speed'] if 'speed' in cpu else 0
        self.cpu_flops: int = cpu['count'] * cpu['speed'] * 10 ^ 6 if 'speed' in cpu else 0
        self.cpu_vendor: str = cpu['vendor'] if 'vendor' in cpu else None

        self.logger.debug("created machine: {0} with {1} cores and {2} FLOPS.".format(
            self.name, self.cpu_cores, self.cpu_flops)
        )

    def as_dict(self) -> Dict:
        """A JSON representation of the machine.

        :return: A JSON object representation of the machine.
        :rtype: Dict
        """
        machine = {"nodeName": self.name}
        if self.system:
            machine['system'] = self.system.value
        if self.architecture:
            machine['architecture'] = self.architecture
        if self.memory:
            machine['memory'] = self.memory
        if self.release:
            machine['release'] = self.release
        if self.cpu_cores:
            machine['cpu'] = {'count': self.cpu_cores}
        if self.cpu_speed:
            machine['cpu']['speed'] = self.cpu_speed
        if self.cpu_vendor:
            machine['cpu']['vendor'] = self.cpu_vendor
        return machine
