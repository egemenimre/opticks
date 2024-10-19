# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


from pathlib import Path


def process_paths(filepath: Path, *search_dirs: Path):
    """
    Inits a filepath with different alternative locations.

    Searches the file first in the given path, then in the
    current working directory. If it does not exist, tries the
    alternate directories. Returns `None` if the file is not found
    in any of the directories.

    - Nominal path: `current working dir` + `filepath`
    - Alternate paths: `current working dir` + `alternate search dir` + `filepath`

    Parameters
    ----------
    filepath
        File path
    search_dirs
        Alternate search directories

    Returns
    -------
    Path
        Path of the file, `None` if not found
    """
    working_dir = Path.cwd()

    if Path(filepath).exists():
        # check whether the file is at the current working dir
        return filepath.resolve()
    else:
        # build the file path at the current working dir
        file_path = working_dir.joinpath(filepath)
        if file_path.exists():
            return file_path.resolve()
        else:
            # search the remaining directories
            for search_dir in search_dirs:
                file_path = working_dir.joinpath(search_dir).joinpath(filepath)
                if file_path.exists():
                    return file_path.resolve()

    # File is nowhere to be found, return None
    return None
