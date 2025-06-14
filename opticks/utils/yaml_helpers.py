# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for schema and yaml helpers.

"""
from astropy.units import Quantity
from strictyaml import ScalarValidator
from strictyaml.exceptions import YAMLSerializationError

from opticks import u

# add bpp as a unit
u.add_enabled_units([u.def_unit("bpp", u.bit / u.pix)])


def is_quantity(quantity_text: str) -> bool:
    """
    Checks whether the `quantity_text` can be parsed as a valid `Quantity` object.

    Parameters
    ----------
    quantity_text : str
        Text that will be parsed as a `Quantity` object

    Returns
    -------
    is_quantity : bool
        `True` if text can be parsed, `False` otherwise
    """
    try:
        Quantity(quantity_text)
    except ValueError:
        return False
    else:
        return True


class Qty(ScalarValidator):
    """
    Quantity parser for the YAML Schema.
    """

    def validate_scalar(self, chunk) -> Quantity:
        """
        Validator for the `Quantity` object.

        Parameters
        ----------
        chunk : YAMLChunk
            Input YAML Chunk to be parsed

        Returns
        -------
        Quantity
            `Quantity` object parsed from YAML

        """
        val = chunk.contents
        if not is_quantity(chunk.contents):
            chunk.expecting_but_found("when expecting a Quantity")
        else:
            # Only Python 3.6+ supports underscores in numeric literals
            return Quantity(val.replace("_", ""))

    def to_yaml(self, data) -> str:
        """Converts `Quantity` data to YAML."""
        if isinstance(data, Quantity):
            if is_quantity(str(data)):
                return str(f"{data}")
        raise YAMLSerializationError(f"'{data}' is not a Quantity object.")


def _dict_to_obj_init(self, class_name: str, dictionary: dict):
    """
    Converts a (nested) dict to an object.
    """
    for key, value in dictionary.items():
        if isinstance(value, dict):

            # generate the new class from the attributes
            new_class = type(
                key.capitalize(),
                (object,),
                value
                | {
                    "__init__": _dict_to_obj_init,
                    "__str__": lambda x: str(x.__dict__),
                    "__setitem__": lambda x, new_key, new_value: x.__dict__.update(
                        {new_key: new_value}
                    ),
                },
            )

            # instantiate the new class
            value = new_class(key.capitalize(), value)

        setattr(self, key, value)


def dict_to_obj_creator(class_name: str, dictionary: dict):
    """
    Creates a class that keeps the dict as an object.

    Parameters
    ----------
    class_name: str
        Name of the class
    dictionary : dict
        Dictionary to be converted to object
    Returns
    -------
    type
        Class to create the object
    """

    return type(
        class_name,
        (object,),
        {
            "__init__": _dict_to_obj_init,
            "__str__": lambda self: str(self.__dict__),
            # add new item
            "__setitem__": lambda self, new_key, new_value: self.__dict__.update(
                {new_key: new_value}
            ),
        },
    )
