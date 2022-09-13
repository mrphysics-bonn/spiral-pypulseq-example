# Spiral Pypulseq Example with ISMRMRD protocol creation

This repository contains a Jupyter Notebook with an example code for a spiral MRI sequence in the Pulseq open source file format [1]. The Pypulseq [2] package for Python is used for the designing pulses and creating the sequence file. The Jupyter Notebook contains a reduced example for better visualization.

The protocol file is written in the ISMRMRD format [3] and can be used for reconstruction of acquired data with the [Python Reco Server](https://github.com/mrphysics-bonn/python-ismrmrd-reco)

## Usage/Installation/Requirements

The Jupyter Notebook can be used inside binder without further installation requirements. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mrphysics-bonn/spiral-pypulseq-example/HEAD?labpath=spiral_example.ipynb)

The Notebook can also be used with a local Python installation containing Jupyter Lab and numpy with the following additional installation requirements:

* Pypulseq: (https://github.com/imr-framework/pypulseq)
* spiraltraj: see README in the spiraltraj folder
* ISMRMRD Python API:
```console
pip install ismrmrd
```
or
```console
git clone https://github.com/ismrmrd/ismrmrd-python
cd ismrmrd-python
git checkout v1.9.3
pip install .
```

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] Ravi, Keerthi, Sairam Geethanath, and John Vaughan. "PyPulseq: A Python Package for MRI Pulse Sequence Design." Journal of Open Source Software 4.42 (2019): 1725., https://github.com/imr-framework/pypulseq

[3] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io
