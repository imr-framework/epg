# EPG
Python implementation of the Extended Phase Graph (EPG) algorithm for multi-echo MRI.

A MATLAB version can be find within the same project as https://github.com/imr-framework/epg-matlab

For the theory, see https://www.ncbi.nlm.nih.gov/pubmed/24737382 (Weigel, 2015) or read [this powerpoint presentation](https://github.com/imr-framework/epg-matlab/blob/master/epg/EPG-notes-presentation.pptx).

## Introduction
The extended phase graph algorithm is an alternative to time-domain Bloch simulations. Operating in the Fourier domain, it enables straightforward echo detection and is valuable for multi-echo sequences such as Turbo Spin Echo. The current scripts are capable of simulating the response to regularly spaced RF pulses with arbitrary phase and flip angle with linear gradients and T1, T2 relaxation effects. By “regularly spaced”, we mean that the intervals between pulses are integer multiples of the same Δt. Arbitrary spacing may be approximated with a small Δt and a large number of intervals, but the simulation will be slow.

The code has been translated from the MATLAB version, which was developed by Gehua Tong and Sairam Geethanath based on Weigel's paper (2015).

## Quick Start
After cloning the repository and installing the dependencies, run the following code to simulate and visualize a simple TSE sequence.

```Python

import epg
import psd_library as pl

# Make a predefined sequence
S = pl.make_tse_sequence(alpha=120, N=5, esp=100, use_y90=True, t1t2=(1000,100))
# Initialize & run EPG
epg1 = epg.EPG(S)
epg1.simulate()
# Displays EPG diagram
epg1.display()
# Displays echo strengths 
epg1.plot_echoes()

```




## Using the EPG functions 
Core functions for the EPG algorithm:
* `rf_rotation(phi, alpha)`
* `relaxation(tau, t1, t2, omega)`
* `grad_shift(dk, omega)`

Classes and class methods for running EPG:
* `Sequence(rf, grad, events, time, t1t2=(0,0), name="")`
* `EPG(seq)`
* `EPG.simulate()`
* `EPG.reset()`

Class methods for display & output:
* `Sequence.plot_seq()`
* `EPG.find_echoes()`
* `EPG.plot_echoes()`
* `EPG.display()`
* `EPG.get_sim_history()`

The `psd_library.py` script includes functions for making pre-defined sequences (with custom parameters) for simulation.

## Operators

EPG represents the magnetization in terms of configuration states. At any time, the configuration matrix Ω looks like this:

<p align="center"> <a>
    <img title="eqn1" src="https://github.com/imr-framework/epg-matlab/blob/master/EPG_guide_equations/e1.PNG" width="225">
  </a></p>
	
There are three basic operators in EPG simulation: RF pulse, gradient shift, and T1 & T2 relaxation, represented by the following functions:

* `rf_rotation(phi, alpha)` returns a 3 x 3 rotation matrix
* `relaxation(tau, t1, t2, omega)` returns the new omega after relaxation over interval tau
* `grad_shift(dk, omega)` takes the configuration states and shifts it by an integer dk


### RF Rotation

The rotation operator exchanges the F+, F-, and Z populations within each k. 
The operator is equivalent to a matrix multiplication, Ω=T(Φ,α)Ω, where:

<p align="center"> <a>
    <img title="eqn2" src="https://github.com/imr-framework/epg-matlab/blob/master/EPG_guide_equations/e2.PNG" width="400">
  </a></p>

### Gradient shift
A positive unit gradient (dk = 1) moves all the F populations to a higher k but keeps the Z populations in place. For example:

<p align="center"> <a>
    <img title="eqn3" src="https://github.com/imr-framework/epg-matlab/blob/master/EPG_guide_equations/e3.PNG" width="500">
  </a></p>

With a negative gradient, the populations move in the opposite direction, and as gradients are applied, the matrix grows more columns.

### T1 & T2 Relaxation
T1 and T2 relaxation are represented together in the following operator. Over a time interval τ, the relaxation factors are E1=exp(-τ/T1) and E2=exp(-τ/T2).

For k = 0:

<p align="center"> <a>
    <img title="eqn4" src="https://github.com/imr-framework/epg-matlab/blob/master/EPG_guide_equations/e4.PNG" width="350">
  </a></p>	
  
And for k ≠0,

<p align="center"> <a>
    <img title="eqn5" src="https://github.com/imr-framework/epg-matlab/blob/master/EPG_guide_equations/e5.PNG" width="270">
  </a></p>	

## References

*Weigel, Matthias. "Extended phase graphs: dephasing, RF pulses, and echoes ‐ pure and simple." Journal of Magnetic Resonance Imaging 41.2 (2015): 266-295.*

*Hennig J, Scheffler K. Hyperechoes. Magnetic Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine. 2001 Jul;46(1):6-12.*

