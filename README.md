# MrParse
An aid to decision making in Molecular Replacement

MrParse is a [CCP4](http://www.ccp4.ac.uk) program takes a protein amino acid sequence file and searches for homologs using [PHMMER](http://hmmer.org/) (default) or [HHSEARCH](https://github.com/soedinglab/hh-suite). If supplied with a reflection data file (currently in [MTZ](http://www.ccp4.ac.uk/html/mtzformat.html) format), it can then use [PHASER](https://www.phaser.cimr.cam.ac.uk/index.php/Phaser_Crystallographic_Software) to calculate the Expected Log Likelihood Gains (eLLG) values for the homologs. It also searched the [EBI AlphaFold database](https://alphafold.ebi.ac.uk/) _ab initio_ models. 

It also attempts to classify the sequence according to its secondary structure, and whether any regions are expected to be Coiled-Coil or Transmembrane.

Results are currently displayed in a simple HTML webpage that is rendered using [VUE](https://vuejs.org). The sequence graphics are created using the [PFAM graphics library](https://pfam.xfam.org/generate_graphic), a copy of which is distributed with this code.

### Search Model Finder
The search model finder currently uses [PHMMER](http://hmmer.org/) (distributed with [CCP4](http://www.ccp4.ac.uk)) to search for homologs. The facility to use [HHSEARCH](https://github.com/soedinglab/hh-suite) is almost complete, but is waiting on the HHSEARCH parsing functionality implemented in the GitHub [pull request](https://github.com/biopython/biopython/pull/1965) to be incorporated into the [BioPython](https://biopython.org) release.

### EBI Alphafold database search
The search model finder currently uses [PHMMER](http://hmmer.org/) (distributed with [CCP4](http://www.ccp4.ac.uk)) to search the EBI Alphafold database. This will be replaced with the [3Dbeacons](https://github.com/3D-Beacons) API when it becomes available. 

### Classifiers
* secondary structure classification is currently carried by submitting jobs to the [JPRED](http://www.compbio.dundee.ac.uk/jpred/) server 
* Coiled-Coil classification is carried out with [Deepcoil](https://github.com/labstructbioinf/DeepCoil). This needs to be installed locally.
* Transmembrane classification is carried out with [TMHMM](https://github.com/dansondergaard/tmhmm.py). This needs to be installed locally. 
