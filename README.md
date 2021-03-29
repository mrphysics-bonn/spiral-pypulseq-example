# Spiral Pypulseq Example with BART Reconstruction

This repository contains a Jupyter Notebook with an example code for a spiral MRI sequence in the Pulseq open source file format [1]. The Pypulseq [1] package for Python is used for the designing pulses and creating the sequence file.

Additionally, a protocol file is written, which can be used for reconstruction of acquired data.

## Usage/Installation/Requirements

The Jupyter Notebook can be used inside binder without further installation requirements. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mrphysics-bonn/spiral-pypulseq-example/master?filepath=spiral_example.ipynb)

The Notebook can also be used with a local Python installation containing Jupyter Lab and numpy with the following additional installation requirements:

* Pypulseq (https://github.com/imr-framework/pypulseq)
* spiraltraj - see README in the spiraltraj folder
* ISMRMRD Python API:
```console
git clone https://github.com/ismrmrd/ismrmrd-python
cd ismrmrd-python
git checkout 80fecd0
python setup.py install
```

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] Ravi, Keerthi, Sairam Geethanath, and John Vaughan. "PyPulseq: A Python Package for MRI Pulse Sequence Design." Journal of Open Source Software 4.42 (2019): 1725., https://github.com/imr-framework/pypulseq

[3] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io
