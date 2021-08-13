#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import jsonschema
import logging
import os
import requests

from logging import Logger
from typing import Any, Dict, Optional


class SchemaValidator:
    """
    Validate JSON files against WfCommons schema. If schema file path is not
    provided, it will look for a local copy of the WfCommons schema, and if
    not available it will fetch the latest schema from the
    `WfCommons schema GitHub <https://github.com/wfcommons/workflow-schema>`_
    repository.

    :param schema_file: JSON schema file path.
    :type schema_file: str
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self, schema_file: Optional[str] = None, logger: Optional[Logger] = None) -> None:
        """Create an object of the schema validator class."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.schema = self._load_schema(schema_file)

    def validate_instance(self, data: Dict[str, Any]):
        """Perform syntax validation against the schema, and semantic validation.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        self._syntax_validation(data)
        self._semantic_validation(data)

    def _load_schema(self, schema_file: Optional[str] = None):
        """
        Load the schema file. If schema file path is not provided, it will look for
        a local copy of the WfCommons schema, and if not available it will fetch
        the latest schema from the GitHub repository.

        :param schema_file: JSON schema file path.
        :type schema_file: str

        :return: The JSON schema.
        :rtype: json
        """
        if schema_file:
            # schema file provided
            self.logger.info('Using schema file: {}'.format(schema_file))
            return json.loads(open(schema_file).read())

        # looking for local copy of schema file
        schema_path = os.getcwd() + '/wfcommons-schema.json'
        if os.path.exists(schema_path):
            self.logger.info('Using schema file: {}'.format(schema_path))
            return json.loads(open(schema_path).read())

        # fetching latest schema file from GitHub repository
        url = 'https://raw.githubusercontent.com/wfcommons/wfformat/master/wfcommons-schema.json'
        response = requests.get(url)
        schema = json.loads(response.content)
        with open(schema_path, 'w') as outfile:
            json.dump(schema, outfile)
        self.logger.info(
            'Using latest schema file from GitHub repository (saved local copy into {}).'.format(schema_path))
        return schema

    def _syntax_validation(self, data: Dict[str, Any]):
        """Validate the JSON workflow execution instance against the schema.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        v = jsonschema.Draft4Validator(self.schema)
        has_error = False
        for error in sorted(v.iter_errors(data), key=str):
            msg = ' > '.join([str(e) for e in error.relative_path]) \
                  + ': ' + error.message
            self.logger.error(msg)
            has_error = True

        if has_error:
            raise RuntimeError('The workflow instance has syntax errors.')

    def _semantic_validation(self, data: Dict[str, Any]):
        """Validate the semantics of the JSON workflow execution instance.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        has_error = False

        machine_ids = []
        if 'machines' in data['workflow']:
            for m in data['workflow']['machines']:
                machine_ids.append(m['nodeName'])
        else:
            self.logger.debug('Skipping machines processing.')

        tasks_ids = []
        for j in data['workflow']['jobs']:
            tasks_ids.append(j['name'])
            if 'machine' in j and j['machine'] not in machine_ids:
                self.logger.error('Machine "{}" is not declared in the list of machines.'.format(j['machine']))
                has_error = True

        # since tasks may be declared out of order, their dependencies are only verified here
        for j in data['workflow']['jobs']:
            for p in j['parents']:
                if p not in tasks_ids:
                    self.logger.error(
                        'Parent task "{}}" is not declared in the list of workflow tasks.'.format(p['parentId']))
                    has_error = True

        self.logger.debug('The workflow has {} tasks.'.format(len(tasks_ids)))
        self.logger.debug('The workflow has {} machines.'.format(len(machine_ids)))

        if has_error:
            raise RuntimeError('The workflow instance has semantic errors.')
