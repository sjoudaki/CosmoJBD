# CosmoJBD

We release modified [EFTCAMB](https://github.com/EFTCAMB) and [CosmoMC](https://github.com/cmbant/CosmoMC) modules and associated data files to perform a comprehensive cosmological analysis of Jordan-Brans-Dicke (JBD) gravity. The nonlinear corrections to the Weyl and matter power spectra are included in a modified [HMCODE](https://github.com/alexander-mead/HMcode) that we have calibrated to a hybrid suite of COLA and RAMSES _N_-body simulations. For further information, please [see the readme file](https://github.com/sjoudaki/CosmoJBD/blob/main/readme_cosmojbd).

The underlying theory, cosmological observables and datasets (in particular Planck, KiDS, 2dFLenS, BOSS, Pantheon), treatment of systematic uncertainties, analysis techniques, and fitting pipeline are presented in [Joudaki, Ferreira, Lima, and Winther (2022)](https://arxiv.org/abs/2010.15278). The central chains of the paper are available to [download from here](https://u.pcloud.link/publink/show?code=XZgKscXZgOPsJcYsiY8kDVFbNKIGiLSmoGTk).

A patch with the modifications to RAMSES can be found [here](https://github.com/HAWinther/RamsesPatchApproxMGSolver) and the COLA code used can be found [here](https://github.com/HAWinther/MG-PICOLA-PUBLIC).

# Installation and use

* The code requires OpenMP and Intel Fortran (ifort14 or higher). In order to use the Planck 2018 likelihood, a [separate installation](https://cosmologist.info/cosmomc/readme_planck.html) is needed. 

* Compile the code by executing "make clean", followed by "make" in the source directory. 

* Then run the code with multiple chains and threads by executing "sbatch testjbd.sh" in the main directory.

Please feel free to [contact us](mailto:shahab.joudaki@physics.ox.ac.uk) if you have any questions.
