# Spiral Pypulseq Example with BART Reconstruction

This repository contains a Jupyter Notebook with an example code for a spiral MRI sequence in the Pulseq open source file format [1]. The Pypulseq [2] package for Python is used for the designing pulses and creating the sequence file.

Additionally, an example reconstruction of spiral MRI data is shown. This data was acquired with an enhanced version of the presented Pulseq sequence. The reconstruction is done with the Berkeley Advanced Reconstruction Toolbox (BART) [3]. Two different reconstructions are perfomed using a k-space trajectory predicted with the Gradient Impulse Response Function (GIRF) and using the nominal k-space trajectory.

## Usage/Installation/Requirements

The Jupyter Notebook can be used inside binder without further installation requirements. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mrphysics-bonn/spiral-pypulseq-example/esmrmb2020?filepath=spiral_example.ipynb)

Reconstruction with the BART toolbox can be quite slow inside binder, if a high number of iterations is chosen. The Notebook can also be used with a local Python installation containing Jupyter Lab and numpy with the following additional installation requirements:

* Pypulseq (https://github.com/imr-framework/pypulseq) 
* BART Toolbox (https://mrirecon.github.io/bart)       
* spiraltraj - see README in the spiraltraj folder

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] Ravi, Keerthi, Sairam Geethanath, and John Vaughan. "PyPulseq: A Python Package for MRI Pulse Sequence Design." Journal of Open Source Software 4.42 (2019): 1725., https://github.com/imr-framework/pypulseq

[3] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart/
