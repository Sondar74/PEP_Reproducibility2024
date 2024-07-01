# Pseudo-Energy-Preserving Explicit Runge-Kutta Methods

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://zenodo.org/doi/TODO)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{barrios2024pseudo,
  title={Pseudo-Energy-Preserving Explicit {R}unge-{K}utta Methods},
  author={Barrios, Gabriel and Ketcheson, David I and Ranocha, Hendrik},
  year={2024},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{barrios2024pseudoRepro,
  title={Reproducibility repository for
         "{P}seudo-Energy-Preserving Explicit {R}unge-{K}utta Methods"},
  author={Barrios, Gabriel and Ketcheson, David I and Ranocha, Hendrik},
  year={2024},
  howpublished={\url{https://github.com/Sondar74/PEP_Reproducibility2024}},
  doi={TODO}
}
```

## Abstract

Using a recent characterization of energy-preserving B-series, we derive the
explicit conditions on the coefficients of a Runge-Kutta method that ensure
energy preservation (for Hamiltonian systems) up to a given order in the step
size, which we refer to as the pseudo-energy-preserving (PEP) order.  We study
explicit Runge-Kutta methods with PEP order higher than their classical order.
We provide examples of such methods up to PEP order six, and test them on
Hamiltonian ODE and PDE systems. We find that these methods behave similarly
to exactly energy-conservative methods over moderate time intervals and
exhibit significantly smaller errors, relative to other Runge-Kutta methods
of the same order, for moderately long-time simulations.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/) and Python.
The numerical experiments presented in this article were performed using
Julia vTODO and Python vTODO.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface.

TODO: Describe how to run Julia code:
Then, you need to start Julia in the `code_julia` directory of this
repository and follow the instructions described in the `README.md` file
therein.

To generate all figures of the paper, you can execute the Python code in
the Jupyter notebook `PEP_paper_Reproducibility.ipynb`.



## Authors

- Gabriel Barrios (TODO: affiliation)
- [David I. Ketcheson](https://www.davidketcheson.info) (King Abdullah University of Science and Technology)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
