# Pseudo-Energy-Preserving Explicit Runge-Kutta Methods

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12699752.svg)](https://zenodo.org/doi/10.5281/zenodo.12699752)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@article{barrios2025pseudo,
  title={Pseudo-Energy-Preserving Explicit {R}unge-{K}utta Methods},
  author={Barrios de Leon, Gabriel and Ketcheson, David I and Ranocha, Hendrik},
  journal={ESAIM: Mathematical Modelling and Numerical Analysis (ESAIM: M2AN)},
  volume={59},
  number={2},
  pages={729--748},
  year={2025},
  month={03},
  eprint={2407.15365},
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
  author={Barrios de Leon, Gabriel and Ketcheson, David I and Ranocha, Hendrik},
  year={2024},
  howpublished={\url{https://github.com/Sondar74/PEP_Reproducibility2024}},
  doi={10.5281/zenodo.12699752}
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
to install Python and [Julia](https://julialang.org/).
The numerical experiments presented in this article were performed using
Python v3.10.9  and Julia v1.10.4.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface.

To generate all figures of the paper, you can execute the Python code in
the Jupyter notebook `PEP_paper_Reproducibility.ipynb`.

The Jupyter notebook loads some output data generated by the Julia code
in the `code_julia` directory. To generate the data, you need to start Julia
in the `code_julia` directory of this repository and follow the instructions
described in the `README.md` file therein.



## Authors

- Gabriel Barrios de Léon (Universidad de San Carlos de Guatemala)
- [David I. Ketcheson](https://www.davidketcheson.info) (King Abdullah University of Science and Technology)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
