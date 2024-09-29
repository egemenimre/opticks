# opticks

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
- [Pint](https://github.com/hgrecco/pint) provides units and quantity support.
- [Prysm](https://github.com/brandondube/prysm/) provides advanced optical design and analysis.
- [Strictyaml](https://github.com/crdoconnor/strictyaml) provides YAML read/write support.

## Licence

This project is Copyright (c) Egemen Imre and licensed under the terms of the GNU GPL v3+ licence.
