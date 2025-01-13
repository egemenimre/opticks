# opticks: Systems Engineering Tool for Optical Systems

Models and analysis tools for optics system engineering.

Current functionality includes:

- Basic optical system and output electronics sizing analysis (Field of View, Ground Sampling Distance, data rates...)
- Complex aperture definition (Circular, obscurations and more complex shapes)
- PSF and MTF computation with detailed Wavefront Error definition.

Check the [source code on GitHub](https://github.com/egemenimre/opticks).

[Tutorials](tutorials.md) provide information on setting up basic tasks.See also [examples and how-to guides page](how_to_guides.md) for a few sample notebooks to carry out complex tasks.

## Requirements

- NumPy and SciPy are used for the underlying mathematical algorithms
- [Astropy](https://www.astropy.org/) provides units and quantity support as well as some other underlying models such as vector mechanics.
- [Prysm](https://github.com/brandondube/prysm/) provides advanced optical design and analysis.
- [Strictyaml](https://github.com/crdoconnor/strictyaml) provides YAML read/write support.

## Licence

This project is Copyright (c) Egemen Imre and licensed under the terms of the GNU GPL v3+ licence.

```{toctree}
---
maxdepth: 2
caption: Understanding opticks
---
models_index
```

```{toctree}
---
maxdepth: 2
caption: How to start using opticks
---
tutorials
```

```{toctree}
---
maxdepth: 2
caption: How to use opticks for real world tasks
---
how_to_guides
```

```{toctree}
---
maxdepth: 2
caption: References
---
api_index
```
