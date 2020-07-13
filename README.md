#############################################################

# Ando Lab edits to the pySCA6.0 code #

#############################################################
_README file for pySCA_

_03.2015_ 

_Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds_

_This program is free software distributed under the BSD 3-clause_

_license, please see the file LICENSE for details._
#############################################################

The current version of the Statistical Coupling Analysis (SCA)
analysis is implemented in python. This directory contains the
necessary code for running the SCA calculations, as well examples/tutorials for
the dihydrofolate reductase (DHFR) enzyme family, the S1A serine
proteases, the small G-protein family and the Beta-lactamase enzyme
family. The tutorials are distributed as iPython notebooks; for
details please see: http://ipython.org/notebook.html

#############################################################

For installation instructions, and an introduction to using the
toolbox, please see:

https://github.com/reynoldsk/pySCA

or open the html files included with the pySCA distribution in a browser:

html_docs/index.html

#############################################################

Contents

   Inputs/		  :  Directory containing input files (including those
   		   	     needed for the tutorials)
  
   Outputs/   	   	  :  Directory for output files (empty at install)
  
   html_docs/             :  Directory containing html documentation

   annotate_MSA.py	  :  Python script that annotates alignments
   			     with phylogenetic/taxonomic information   

   scaProcessMSA.py	  :  Python script that conducts some initial
   			     processing of the sequence alignment
  
   scaCore.py             :  Python script that runs the core SCA calculations
  
   scaSectorID.py         :  Python script that defines sectors given the
   			     results of the calculations in scaCore
  
   scaTools.py		  :  The SCA toolbox - contains all functions
   			     needed for the SCA calculations
  
   SCA_DHFR.ipynb	  :  Python notebook example for DHFR
   
   SCA_G.ipynb		  :  Python notebook example for the small G proteins
   
   SCA_betalactamase.ipynb:  Python notebook example for the beta-lactamases

   SCA_S1A.ipynb	  :  Python notebook example for the S1A
                             serine protease
#############################################################
