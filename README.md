WASA raw data Analysis
======================
Sources of my software for the analysis of raw data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
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


Compiling
=========

    git clone https://github.com/alexkernphysiker/WASA-RawAnalysis.git
    cd WASA-RawAnalysis
    git submodule update --init --recursive
    cd ..
    mkdir WASA-build
    cd WASA-build
    cmake ../WASA-RawAnalysis
    make

Running from build directory
============================

Configure WMC and rootsorter for current experimental setup

    cd config
    ./config.sh 
    cd ..

Run pluto
    

Run WMC

Analyse Monte Carlo

Analyse raw data


List of implemented reactions
=============================

    He3eta
    He3pi0
    He3pi0pi0
    He3pi0pi0pi0
