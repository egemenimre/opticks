# opticks

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/egemenimre/opticks/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/egemenimre/opticks/tree/main)
[![codecov](https://codecov.io/gh/egemenimre/opticks/graph/badge.svg?token=3KRKPXJST3)](https://codecov.io/gh/egemenimre/opticks)
[![Documentation Status](https://readthedocs.org/projects/opticks/badge/?version=latest)](https://opticks.readthedocs.io/en/latest/?badge=latest)

Models and analysis tools for optical system engineering.

Current functionality includes:

- Basic optical system and output electronics sizing analysis (Field of View, Ground Sampling Distance, data rates...)
- Complex aperture definition (Circular, obscurations and more complex shapes)
- PSF and MTF computation with detailed Wavefront Error definition.

See [examples](docs/examples) and [tutorials](docs/tutorials) directories for a few sample notebooks on simple and complex tasks.

The documentation for opticks is here: <https://opticks.readthedocs.io/>

## Requirements

- NumPy and SciPy are used for the underlying mathematical algorithms
- [Astropy](https://www.astropy.org/) provides units and quantity support as well as some other underlying models such as vector mechanics.
- [Prysm](https://github.com/brandondube/prysm/) provides advanced optical design and analysis.
- [Pydantic](https://docs.pydantic.dev/) provides data validation and settings management.
- [PyYAML](https://pyyaml.org/) provides YAML read/write support.
- [Portion](https://github.com/AlexandreDecan/portion) provides interval arithmetic support.
- [Matplotlib](https://matplotlib.org/) provides plotting support.

## Licence

This project is Copyright (c) Egemen Imre and licensed under the terms of the GNU GPL v3+ licence.
