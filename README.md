# gwtails

This is a Python package to analyze post-merger late-time tails in eccentric binary black hole (BBH) mergers. 

This code is based on the paper "Phenomenology and origin of late-time tails in eccentric binary black hole mergers" by Islam, Faggioli, Khanna, Field, van de Meent and Buonanno. Here, we inspect waveforms produced by merging eccentric binary black holes (BBH) in . These waveforms are all generated using black hole perturbation theory. The trajectories are described by the dynamics of the point-particle orbiting the Kerr black hole using effective-one-body formalism: http://arxiv.org/abs/gr-qc/9811091, http://arxiv.org/abs/gr-qc/0001013, http://arxiv.org/abs/2405.19006. These trajectories are then fed into the time-domain Teukolsky solver developed by Gaurav Khanna and collaborators: http://arxiv.org/abs/gr-qc/0703028, http://arxiv.org/abs/0803.0317, http://arxiv.org/abs/0803.0317, http://arxiv.org/abs/1003.0485, http://arxiv.org/abs/1108.1816, http://arxiv.org/abs/arXiv:2010.04760.

## Getting the package
The latest development version will always be available from the project git repository:
```bash
git clone https://github.com/tousifislam/gwtails
```
## Issue tracker
Known bugs are recorded in the project bug tracker:
https://github.com/tousifislam/gwtails/issues

## License
This code is distributed under the MIT License. Details can be found in the LICENSE file.

## Maintainer
Tousif Islam

## Citation guideline
If you make use of the gwtails framework, please cite the following papers:

```
@article{Islam:2024rhm,
    author = "Islam, Tousif and Faggioli, Guglielmo and Khanna, Gaurav and Field, Scott and van de Meent,  Maarten and Buonanno, Alessandra",
    title = "{Phenomenology and origin of late-time tails in eccentric binary black hole mergers}",
    eprint = "XXXX.XXXXX",
    archivePrefix = "arXiv",
    primaryClass = "gw-qc",
    month = "6",
    year = "2024"
}
```
