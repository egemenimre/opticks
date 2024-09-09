# opticks: Sizing Tool for Optical Systems

Notes and models on optics system engineering.

Check the [source code on GitHub](https://github.com/egemenimre/opticks).

See [examples and how-to guides page](how_to_guides.md) for a few sample notebooks.

The documentation for opticks is here: <https://opticks.readthedocs.io/>

## Requirements

- NumPy and SciPy are used for the underlying mathematical algorithms
- [Pint](https://github.com/hgrecco/pint) provides units and quantity support.
- [Prysm](https://github.com/brandondube/prysm/) provides advanced optical design and analysis (optional).
- [Strictyaml](https://github.com/crdoconnor/strictyaml) provides YAML read/write support.

## Licence

This project is Copyright (c) Egemen Imre and licensed under the terms of the GNU GPL v3+ licence.

```{toctree}
---
maxdepth: 2
caption: Understanding opticks
---
models/models_index
```

```{toctree}
---
maxdepth: 2
caption: Tutorials
---
tutorials/yaml_reader.ipynb
```

```{toctree}
---
maxdepth: 2
caption: Examples
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
