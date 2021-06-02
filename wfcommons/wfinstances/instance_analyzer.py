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
import math
import numpy
import scipy.stats

from logging import Logger
from matplotlib import pyplot
from os import path
from typing import Any, Dict, List, Optional, Tuple

from .instance import Instance
from ..common.task import Task
from ..common.file import FileLink
from ..utils import best_fit_distribution, NoValue


class InstanceElement(NoValue):
    RUNTIME = ('runtime', 'Runtime (s)')
    INPUT = ('input', 'Input File Size (bytes)')
    OUTPUT = ('output', 'Input File Size (bytes)')


class InstanceAnalyzer:
    """Set of tools for analyzing collections of instances.

    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self, logger: Optional[Logger] = None) -> None:
        """Create an object of the instance analyzer."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.instances: List[Instance] = []
        self.tasks_summary: Dict[str, List:[Task]] = {}
        self.instances_summary: Dict[str, Dict[str, Any]] = {}

    def append_instance(self, instance: Instance) -> None:
        """Append a workflow instance object to the instance analyzer.

        .. code-block:: python

            instance = Instance(input_instance = 'instance.json', schema = 'schema.json')
            instance_analyzer = InstanceAnalyzer()
            instance_analyzer.append_instance(instance)

        :param instance: A workflow instance object.
        :type instance: Instance
        """
        if instance not in self.instances:
            self.instances.append(instance)
            self.logger.debug('Appended instance: {} ({} tasks)'.format(instance.name, len(instance.workflow.nodes)))

    def build_summary(self, tasks_list: List[str], include_raw_data: Optional[bool] = True) -> Dict[
        str, Dict[str, Any]]:
        """Analyzes appended instances and produce a summary of the analysis per task prefix.

        .. code-block:: python

            workflow_tasks = ['sG1IterDecon', 'wrapper_siftSTFByMisfit']
            instances_summary = instance_analyzer.build_summary(workflow_tasks, include_raw_data=False)

        :param tasks_list: List of workflow tasks prefix (e.g., mProject, sol2sanger, add_replace)
        :type tasks_list: List[str]
        :param include_raw_data: Whether to include the raw data in the instance summary.
        :type include_raw_data: bool

        :return: A summary of the analysis of instances in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Dict[str, Any]]
        """
        self.logger.debug('Building summary for {} instances'.format(len(self.instances)))

        tasks_list = sorted(list(tasks_list), key=len, reverse=True)  # had to sorted so it would get all cases
        # build tasks summary
        for instance in self.instances:
            self.logger.debug('Parsing instance: {} ({} tasks)'.format(instance.name, len(instance.workflow.nodes)))

            for node in instance.workflow.nodes.data():
                task: Task = node[1]['task']
                task_name: str = [j for j in tasks_list if task.name.startswith(j)][
                    0]  # it was eliminating bwa_index because bwa came before it
                if task_name not in self.tasks_summary:
                    self.tasks_summary[task_name] = []
                self.tasks_summary[task_name].append(task)

        # build instances summary
        for task_name in self.tasks_summary:
            runtime_list: List[float] = []
            inputs_dict: Dict[str, Any] = {}
            outputs_dict: Dict[str, Any] = {}

            for task in self.tasks_summary[task_name]:
                runtime_list.append(task.runtime)

                for file in task.files:
                    extension: str = path.splitext(file.name)[1] if '.' in file.name else file.name
                    if extension[1:].isnumeric():
                        extension = path.splitext(file.name.replace(extension, ''))[1]

                    if file.link == FileLink.INPUT:
                        _append_file_to_dict(extension, inputs_dict, file.size)
                    elif file.link == FileLink.OUTPUT:
                        _append_file_to_dict(extension, outputs_dict, file.size)

            _best_fit_distribution_for_file(inputs_dict, include_raw_data)
            _best_fit_distribution_for_file(outputs_dict, include_raw_data)

            self.instances_summary[task_name] = {
                'runtime': {
                    'min': min(runtime_list),
                    'max': max(runtime_list),
                    'distribution': _json_format_distribution_fit(best_fit_distribution(runtime_list))
                },
                'input': inputs_dict,
                'output': outputs_dict
            }
            if include_raw_data:
                self.instances_summary[task_name]['runtime']['data'] = runtime_list

        return self.instances_summary

    def generate_fit_plots(self, instance_element: InstanceElement, outfile_prefix: Optional[str] = None) -> None:
        """
        Produce fit plots as images for each entry of an instance element generated by the summary analysis. For
        entries in which there are no distribution (i.e., constant value), no plot will be generated.

        :param instance_element: Workflow element for which the fit plots will be generated.
        :type instance_element: InstanceElement
        :param outfile_prefix: Prefix to be attached to each generated plot file name (optional).
        :type outfile_prefix: str
        """
        self.logger.info('Generating fit plots ({}).'.format(instance_element.value[0]))
        outfile_prefix = outfile_prefix + '_' if outfile_prefix else ''

        for task_summary in self.instances_summary:
            outfile = outfile_prefix + task_summary.lower() + '-' + instance_element.value[0]
            el = self.instances_summary[task_summary][instance_element.value[0]]

            if instance_element == InstanceElement.RUNTIME:
                _generate_fit_plots(el, task_summary + ' (' + instance_element.value[0] + ')',
                                    xlabel=instance_element.value[1], outfile=outfile + '.png', logger=self.logger)
            else:
                for k in el:
                    ext = k if '.' not in k else k[1:]
                    _generate_fit_plots(el[k], task_summary + ' (' + instance_element.value[0] + '): ' + ext,
                                        xlabel=instance_element.value[1], outfile=outfile + '-' + ext + '.png',
                                        logger=self.logger)

    def generate_all_fit_plots(self, outfile_prefix: Optional[str] = None) -> None:
        """
        Produce fit plots as images for each entry of the summary analysis. For entries in which there are no
        distribution (i.e., constant value), no plot will be generated.

        :param outfile_prefix: Prefix to be attached to each generated plot file name (optional).
        :type outfile_prefix: str
        """
        for instance_element in InstanceElement:
            self.generate_fit_plots(instance_element, outfile_prefix)


def _append_file_to_dict(extension: str, dict_obj: Dict[str, Any], file_size: int) -> None:
    """Add a file size to a file type extension dictionary.

    :param extension: File type extension.
    :type extension: str
    :param dict_obj: Dictionary of file type extensions.
    :type dict_obj: Dict[str, Any]
    :param file_size: File size in bytes.
    :type file_size: int
    """
    if extension not in dict_obj:
        dict_obj[extension] = {'data': [], 'distribution': None}
    dict_obj[extension]['data'].append(file_size)


def _best_fit_distribution_for_file(dict_obj, include_raw_data) -> None:
    """Find the best fit distribution for a file.

    :param dict_obj: Dictionary of file type extensions.
    :type dict_obj: Dict[str, Any]
    :param include_raw_data:
    :type include_raw_data: bool
    """
    for ext in dict_obj:
        dict_obj[ext]['min'] = min(dict_obj[ext]['data'])
        dict_obj[ext]['max'] = max(dict_obj[ext]['data'])
        if dict_obj[ext]['min'] != dict_obj[ext]['max']:
            dict_obj[ext]['distribution'] = _json_format_distribution_fit(
                best_fit_distribution(dict_obj[ext]['data']))
        if not include_raw_data:
            del dict_obj[ext]['data']


def _json_format_distribution_fit(dist_tuple: Tuple) -> Dict[str, Any]:
    """Format the best fit distribution data into a dictionary

    :param dist_tuple: Tuple containing best fit distribution name and parameters.
    :type dist_tuple: Tuple

    :return:
    :rtype: Dict[str, Any]
    """
    formatted_entry = {'name': dist_tuple[0], 'params': []}
    for p in dist_tuple[1]:
        formatted_entry['params'].append(p)
    return formatted_entry


def _generate_fit_plots(el: Dict, title: str, xlabel: str, outfile: str, font_size: Optional[int] = None,
                        logger: Optional[Logger] = None) -> None:
    """Produce a fit plot as an image for an entry of an instance element generated by the summary analysis.

    :param el: Entry of an instance element generated by the summary analysis.
    :type el: Dict
    :param title: Plot title.
    :type title: str
    :param xlabel: X-axis label.
    :type xlabel: str
    :param outfile: Plot file name.
    :type outfile: str
    :param font_size: Size of the font.
    :type outfile: Optional[int]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """
    if not el['distribution']:
        return

    distribution = getattr(scipy.stats, el['distribution']['name'])
    params = el['distribution']['params']
    kwargs = params[:-2]
    loc = params[-2]
    scale = params[-1]

    raw_data = el['data']
    bins = math.ceil(len(raw_data) / 10)
    normalized = (raw_data - numpy.min(raw_data)) / (numpy.max(raw_data) - numpy.min(raw_data))
    y, x = numpy.histogram(normalized, bins=bins, density=True)

    if font_size:
        old_font_size = pyplot.rcParams["font.size"]
        pyplot.rcParams["font.size"] = str(font_size)

    pyplot.grid(True)
    pyplot.plot(x, distribution.cdf(x, *kwargs, loc=loc, scale=scale), 'r-', lw=2, alpha=0.6,
                label=distribution.name)
    pyplot.hist(normalized, bins=bins, density=True, cumulative=True)

    pyplot.title(title)
    pyplot.xlabel(xlabel)
    pyplot.ylabel('CDF')
    pyplot.legend(loc='best')
    pyplot.savefig(outfile)
    pyplot.clf()
    logger.info('Generated fit plot: {}'.format(outfile))

    if font_size:
        pyplot.rcParams["font.size"] = old_font_size
