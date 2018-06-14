
###################
##  README file  ##
###################

iOmicsPASS version 1.0 
Creator: Hiromi WL Koh

To compile, using Terminal/Shell and in this directory, type "make". The compilation requires gcc which can be acquired by installing Xcode for Mac OS X users and Cygwin for Windows users.

Please read manual (iOmicsPASSmanual.pdf) for more information about the tool and how to use it.

Folders:
1) bin 
- holds the executable program after compilation

2) example 
- holds all the example datasets for running iOmicsPASS (see Manual) 

3) Library 
- Contains the dependencies (boost library), required in compiling iOmicsPASS.


4) NetworkFiles 
- holds 3 types of network files for Homo sapiens: (a) Transcription-factor regulatory network, (b) Protein-protein interaction network and (c) Pathway Module file (compiled from Census Pathway Database and Gene Ontology)

5) Rcodes
- contains the R-script for producing a plot of the misclassification error in the cross-validations.

6) src
- holds all the c++ scripts required to compile iOmicsPASS.

7) Windows Binary
-  The windows executable program for iOmicsPASS and an input parameter file for the executable.

