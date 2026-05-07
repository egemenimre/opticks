# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
import os
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict


class ImagerComponent(BaseModel):
    """
    Common base class for the generic imager components.

    Subclasses define their YAML parameters as Pydantic field annotations.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

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
            Component object from the input data
        """

        if not file_path:
            raise ValueError("Path is None. No file to be found.")
        elif not os.path.isfile(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")
        else:
            with open(file_path, "rt") as file:
                return cls.from_yaml_text(file.read())

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
            Component object from the input data
        """

        if not yaml_text:
            raise ValueError("Text content is None.")
        else:
            data = yaml.safe_load(yaml_text)
            return cls(**data)

    def to_yaml_text(self) -> str:
        """
        Serialise the component to a YAML string.

        Allows Unicode characters in the dump.

        Returns
        -------
        str
            YAML text representation of the component
        """
        return yaml.dump(
            self.model_dump(mode="json"),
            default_flow_style=False,
            allow_unicode=True,
        )

    def to_yaml_file(self, file_path: Path) -> None:
        """
        Serialise the component to a YAML file.

        Parameters
        ----------
        file_path : Path
            Filepath to write the YAML data
        """
        with open(file_path, "wt") as file:
            file.write(self.to_yaml_text())
