# coreShellPlasmons

****************************************************************************************************************************************
****************************************************************************************************************************************
The code that is used to model the LSP resonances of quasistatic concentric-shell Drude-metal nanoparticles.
****************************************************************************************************************************************
****************************************************************************************************************************************


A short description of each file:


****************************************************************************************************************************************
Directory -- modelCode
****************************************************************************************************************************************

fitFunc\_s.wls:

A Wolfram Language Script file that provides functions that solve Poisson's equation _symbolically_ for a set of concentric spherical shells and use the solutions to calculate various observales. In particular, these functions 1) use the solutions to produce the dielectric-function-dependent polarizability of the particle, 2) plug Drude-model dielectric functions into the polarizability to reveal that the polarizability can be rewritten as the sum of polarizabilities of a set of oscillators, 3) characterize the relevant parameters of these oscillators, and 4) solve for the corresponding absorption cross-sections of each oscillator. The "\_s" suffix denotes that this code is designed to be uploaded to and used on a server or other remote computer. The solution of Poisson's equation becomes prohibitively slow after N = 7 (alpha7lm takes tens of hours to compute on a supercomputer) and so the use of fitFuncConst\_s.wls is recommended for calculations of >7-shell particles.



fitFuncConst\_s.wls:

Similar to fitFunc\_s.wls, but solves the Poisson equation _numerically_. Much faster than fitFunc\_s.wls, but requires that Poisson's equation be solved anew each time a change in dielectric material is required. Integration of this code with a cure-fitting script has not been completed succesfully.



****************************************************************************************************************************************
Directory -- variables
****************************************************************************************************************************************

alphaNlm.wl:

A Wolfram Language object that holds a symbolic polarization variable (alpha) to be uploaded by a Wolfram Language script and used to calculate an absorption cross-section. The numbers N, l, and m stand for the number of shells in the model, the l index of the solution, and the m index of the solution, respectively. When l = 1 and m = 0, the variable describes the polarizability of a set of N z-oriented dipole modes.



****************************************************************************************************************************************
Directory -- notes
****************************************************************************************************************************************

multiShellNotes.tex:

A LaTeX file that contains my notes on how to navigate from first-principles electrostatics to a description of the oscillator parameters of each mode of a multi-shell Drude-model nanosphere.


