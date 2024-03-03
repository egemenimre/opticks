# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
import os
from abc import ABC, abstractmethod
from pathlib import Path

from strictyaml import YAML, Map, load

from opticks.utils.yaml_helpers import dict_to_obj_creator


class ImagerComponent(ABC):
    """
    Common base class for the generic imager components.

    Parameters
    ----------
    yaml : YAML
        YAML containing design data
    """

    def __init__(self, yaml: YAML):
        # extract data structures as a class
        class_name = self._params_class_name()
        self.params = dict_to_obj_creator(class_name, yaml.data)(class_name, yaml.data)

    @classmethod
    @abstractmethod
    def schema(cls) -> Map:
        """Schema to be used when converting YAML data to/from a dict."""
        pass

    @classmethod
    @abstractmethod
    def _params_class_name(cls) -> str:
        """Class name after the dict-to-obj conversion."""
        pass

    @classmethod
    def from_yaml_file(cls, file_path: Path):
        """
        Initialise an Imager component from YAML file.

        Parameters
        ----------
        file_path : Path
            Filepath containing design data file (YAML)

        Returns
        -------
        component_data
            Detector component object from the input data

        """

        if not file_path:
            raise ValueError("Path is None. No file to be found.")
        elif not os.path.isfile(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")
        else:
            # Open the file and read its contents
            with open(file_path, "rt") as file:
                # retrieve the component-specific schema dict
                schema = cls.schema()
                data_yaml = load(file.read(), schema, label=file_path)

            return cls(data_yaml)

    @classmethod
    def from_yaml_text(cls, yaml_text: str):
        """
        Initialise an Imager component from YAML text.

        Parameters
        ----------
        yaml_text : str
            Text containing design data (YAML)

        Returns
        -------
        component_data
            Detector component object from the input data

        """

        if not yaml_text:
            raise ValueError("Text content is None.")
        else:
            # Open the text and read its contents

            # retrieve the component-specific schema
            schema = cls.schema()
            data_yaml = load(yaml_text, schema, label=cls.__name__)

            return cls(data_yaml)
