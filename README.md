WASA raw data Analysis
======================
Sources of my software providing raw data analysis for the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
This is only a part of analysis that somehow requires RootSorter installed on the computer.
All files are distributed under GPL license.


Required software
=================

Framework for calulations

    ROOT

Library for Monte Carlo simulation of reaction products. Uses ROOT.

    Pluto

Software for Monte Carlo simulation of registration of previously simulated particles by WASA detector. Uses ROOT and Pluto.

    WMC

Framework for data preselection (both WMC-generated and obtained from the measurements). Uses ROOT.

    RootSorter

Pluto-less simulation algorithms show debug information in the plots that are created by this software

    gnuplot


Unfortunately, RootSorter cannot be compiler by the newest versions of gcc.
That's why I compile the whole code with gcc that supports only C++11 standard.

Needed environment variables
============================

path where ROOT is installed

    ROOTSYS

path where WMC is installed

    WASA_ROOT

path where RootSorter is installed

    ROOTSORTERSYS

path where resulting files after running simulations are stored

    PLUTO_OUTPUT

path where data from experiment are located

    RUNS_DATA

path where the resulting data after running WMC are stored.

    WMC_DATA

path to store raw analysis results

    WASA_OUTPUT_DATA


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

Running Monte Carlo simulations
===============================

First run simulation for all reactions that are implemented. All these algorithms are implemented in programs named:

    ./sim-*

Then run WMC for all reactions

    ./RunWMC.sh <reaction>

Here is the list of implemented reactions for Monte Carlo

    He3eta
    He3eta_gg (with angulat distribution taken into account)
    He3eta_6g (with angular distribution taken into account)
    He3pi0
    He3pi0pi0
    He3pi0pi0pi0
    pd
    pd_ (with angular distribution taken into account)
    ppn_qf_ (n-spectator)
    bound1-2g (simulation of direct decay of eta coupled by 3He into 2 gammas)
    bound1-6g (simulation of direct decay of eta coupled by 3He into 6 gammas)
    bound2-2g (different p_{Fermi} distribution)
    bound2-6g
    bound3-2g (one more different p_{Fermi} distribution)
    bound3-6g 

After that run analysis for these reaction

    ./RunAnalysis-mc.sh <reaction>

The analysis application will run apropriate analysis algorithms for each reaction

Analyse raw data
================

For the analysis of raw data you should run

    ./RunAnalysis-data.sh <analysis type>

List of analyses for data

    F - forward 3He tracks
    E - pd and ppn reactions
    C - forward 3He + central gammas

These analyses can be provided independently from each other

The result of analysis
======================

I copy all *.root files created in the directory pointed by WASA_OUTPUT_DATA variable to my laptop.
Then I analyse them using the software that is here: https://github.com/alexkernphysiker/WASA-analysis
