#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import jsonschema
import logging
import pathlib
import requests

from logging import Logger
from typing import Any, Dict, Optional


class SchemaValidator:
    """
    Validate JSON files against WfCommons schema (WfFormat). If schema file path
    is not provided, it will look for a local copy of the WfFormat schema, and if
    not available it will fetch the latest schema from the
    `WfFormat schema GitHub <https://github.com/wfcommons/wfformat>`_
    repository.

    :param schema_file_path: JSON schema file path.
    :type schema_file_path: Optional[pathlib.Path]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self, schema_file_path: Optional[pathlib.Path] = None, logger: Optional[Logger] = None) -> None:
        """Create an object of the schema validator class."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.schema = self._load_schema(schema_file_path)

    def validate_instance(self, data: Dict[str, Any]) -> None:
        """
        Perform syntax validation against the schema, and semantic validation.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        self._syntax_validation(data)
        self._semantic_validation(data)

    def _load_schema(self, schema_file_path: Optional[pathlib.Path] = None) -> json:
        """
        Load the schema file. If schema file path is not provided, it will look for
        a local copy of the WfFormat schema, and if not available it will fetch
        the latest schema from the GitHub repository.

        :param schema_file_path: JSON schema file path.
        :type schema_file_path: Optional[pathlib.Path]

        :return: The JSON schema.
        :rtype: json
        """
        if schema_file_path:
            self.logger.info(f'Using schema file: {schema_file_path}')
            return json.loads(open(schema_file_path).read())

        # looking for local copy of schema file
        schema_path = pathlib.Path(f"{pathlib.Path.cwd()}/wfcommons-schema.json")
        if schema_path.exists():
            self.logger.info(f'Using schema file: {schema_path}')
            return json.loads(open(schema_path).read())

        # fetching latest schema file from GitHub repository
        url = 'https://raw.githubusercontent.com/wfcommons/wfformat/master/wfcommons-schema.json'
        response = requests.get(url)
        schema = json.loads(response.content)
        with open(schema_path, 'w') as outfile:
            json.dump(schema, outfile)
        self.logger.info(f"Using latest schema file from GitHub repository (saved local copy into {schema_path}).")
        return schema

    def _syntax_validation(self, data: Dict[str, Any]):
        """
        Validate the JSON workflow execution instance against the schema.

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
        """
        Validate the semantics of the JSON workflow execution instance.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        has_error = False

        machine_ids = []
        if "machines" in data["workflow"]["execution"]:
            for m in data["workflow"]["execution"]["machines"]:
                machine_ids.append(m["nodeName"])
        else:
            self.logger.debug("Skipping machines processing.")

        tasks_ids = []
        for j in data["workflow"]["execution"]["tasks"]:
            tasks_ids.append(j["id"])
            if "machines" in j:
                for m in j["machines"]:
                    if m not in machine_ids:
                        self.logger.error(f"Machine \"{j['machine']}\" is not declared in the list of machines.")
                        has_error = True

        # since tasks may be declared out of order, their dependencies are only verified here
        for j in data["workflow"]["specification"]["tasks"]:
            for p in j["parents"]:
                if p not in tasks_ids:
                    self.logger.error(f"Parent task \"{p}\" is not declared in the list of workflow tasks.")
                    has_error = True

        self.logger.debug(f'The workflow has {len(tasks_ids)} tasks.')
        self.logger.debug(f'The workflow has {len(machine_ids)} machines.')

        if has_error:
            raise RuntimeError('The workflow instance has semantic errors.')
