WASA Preselection
=================
Sources of my software for preselection of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
All files are distributed under GPL license



Required software
=================
    ROOT
framework for calulations

    Pluto
library for Monte Carlo simulation of reaction products. Uses ROOT.

    WMC
software for Monte Carlo simulation of registration of previously simulated particles by WASA detector. Uses ROOT and Pluto.

    RootSorter
framework for data preselection (both WMC-generated and obtained from the measurements). Uses ROOT.

    gnuplot
software for plotting. Is used by some applications performing final analysis.


Needed environment variables
============================
    ROOTSYS
path where ROOT is installed

    WASA_ROOT
path where WMC is installed

    ROOTSORTERSYS
path where RootSorter is installed

    PLUTO_OUTPUT
path where pluto files are stored

    RUNS_DATA
path where data from experiment are located

    RUNS_TMP
path for temporary data

    WMC_DATA
path where data from WMC are located. These data are analysed further by preselection algorithm.

    WASA_OUTPUT_DATA
path to store raw analysis results

    POST_ANALYSIS_DATA
path to store the final results

Directories
===========
    Preselection
Makefile project of software that provides data preselection. Is compiled and runs on wasa00 server. 

    DataFiles
additional information required for data preselection.

    scripts
directory with scripts and additional data for configuring used WASA software and running WASA Monte Carlo and data preselection.



Compiling and running preselection (on wasa00 server)
============================================

    git clone https://github.com/alexkernphysiker/WASA-preselection.git
    cd WASA-preselection
    git submodule update --init --recursive
    
Configure WMC and rootsorter for current experimental setup
    
    cd wasa_scripts/config
    ./config.sh 
    cd ../..

Compile
    
    cd Preselection
    make
    cd ..
    
Go to scripts directory

    cd wasa_scripts
    
Here you can run scripts:

For WMC. (Pluto running application from submodule must have created the files)

    ./QueueWMC <reaction>

For running analysis of the results of WMC work. If reconstruction fits have not been performed

    ./QueueReconstructionSimulation.sh <reaction>

This will create in DataFiles/ directory text files that can be used by fit algorithm
    
    
For running analysis of the results of WMC work. Applying working preselection algorithm.

    ./QueueAnalysis-mc.sh <reaction>

For running analysis of 50 new existing runs of data

    ./QueueAnalysis-data.sh
    

List of implemented reactions
=============================

    He3eta
    He3pi0
    He3pi0pi0
    He3pi0pi0pi0
