# MrParse
Automated Molecular Replacement decision protocol

MrParse is a [CCP4](http://www.ccp4.ac.uk) program takes a protein amino acid sequence file and searches for homologs using [PHMMER](http://hmmer.org/) (default) or [HHSEARCH](https://github.com/soedinglab/hh-suite). If supplied with a reflection data file (currently in [MTZ](http://www.ccp4.ac.uk/html/mtzformat.html) format), it can then use [PHASER](https://www.phaser.cimr.cam.ac.uk/index.php/Phaser_Crystallographic_Software) to calculate the Expected Log Likelihood Gains (eLLG) values for the homologs.

It also attempts to classify the sequence according to its secondary structure, and whether any regions are expected to be Coiled-Coil or Transmembrane.

Results are currently displayed in a simple HTML webpage that is rendered using [VUE](https://vuejs.org). The sequence graphics are created using the [PFAM graphics library](https://pfam.xfam.org/generate_graphic), a copy of which is distributed with this code.

## Notes
### General
The current Biopython shipped with CCP4 is out of date and contains an error that prevents PHMMER log files being parsed correctly. To update the version of Biopython within CCP4, run the command:
```bash
ccp4-python -m pip install --upgrade biopython
```

The latest version of ample is also required. This can be checked out and linked into the current CCP4 installation with the command:
```bash
git clone https://github.com/rigdenlab/ample.git
ample.git/bin/ample_into_ccp4.sh
```

### Search Model Finder
The search model finder currently uses [PHMMER](http://hmmer.org/) (distributed with [CCP4](http://www.ccp4.ac.uk)) to search for homologs. The facility to use [HHSEARCH](https://github.com/soedinglab/hh-suite) is almost complete, but is waiting on the HHSEARCH parsing functionality implemented in the GitHub [pull request](https://github.com/biopython/biopython/pull/1965) to be incorporated into the [BioPython](https://biopython.org) release.

### Classifiers
* secondary structure classification is currently carried by submitting jobs to the [JPRED](http://www.compbio.dundee.ac.uk/jpred/) server 
* Coiled-Coil classification is carried out with [Deepcoil](https://github.com/labstructbioinf/DeepCoil). This needs to be installed locally.
* Transmembrane classification is carried out with [Topcons2](https://github.com/ElofssonLab/TOPCONS2). This functionality is currently broken as the WSDL description file for their online server http://topcons.net/ appears to be missing.

