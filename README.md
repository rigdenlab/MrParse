# MrParse
An aid to decision making in Molecular Replacement

MrParse is a [CCP4](http://www.ccp4.ac.uk) program takes a protein amino acid sequence file and searches for homologs using [PHMMER](http://hmmer.org/) (default) or [HHSEARCH](https://github.com/soedinglab/hh-suite). If supplied with a reflection data file (currently in [MTZ](http://www.ccp4.ac.uk/html/mtzformat.html) format), it can then use [PHASER](https://www.phaser.cimr.cam.ac.uk/index.php/Phaser_Crystallographic_Software) to calculate the Expected Log Likelihood Gains (eLLG) values for the homologs. It also searches the [EBI AlphaFold database](https://alphafold.ebi.ac.uk/) for related models. 

It also attempts to classify the sequence according to its secondary structure, and whether any regions are expected to be Coiled-Coil or Transmembrane.

Results are currently displayed in a simple HTML webpage that is rendered using [VUE](https://vuejs.org). The sequence graphics are created using the [PFAM graphics library](https://pfam.xfam.org/generate_graphic), a copy of which is distributed with this code.

### Simple command line
```bash
mrparse --seqin <PATH TO SEQUENCE FILE>
```

To provide a reflection file and classify the sequence we can provide the following optional flags:
```bash
--hklin <PATH TO MTZ FILE>
--do_classify
```

### Search Model Finder
The search model finder by default uses [PHMMER](http://hmmer.org/) (distributed with [CCP4](http://www.ccp4.ac.uk)) to search for homologs. If installed you can also use [HHSEARCH](https://github.com/soedinglab/hh-suite) using the following flags.
```bash
--search_engine hhsearch
--hhsearch_exe <PATH TO HHSEARCH EXECUTABLE>
--hhsearch_db <PATH TO HHSEARCH PDB70 DATABASE>
```

NOTE: The current Biopython shipped with CCP4 is out of date and contains an error that prevents HHSearch log files being parsed correctly. To update the version of Biopython within CCP4, run the command:
```bash
ccp4-python -m pip install --upgrade biopython==1.76
```

### EBI Alphafold database search
The search model finder currently uses [PHMMER](http://hmmer.org/) (distributed with [CCP4](http://www.ccp4.ac.uk)) to search the EBI Alphafold database. This will be replaced with the [3Dbeacons](https://github.com/3D-Beacons) API when it becomes available. 

### Classifiers
* secondary structure classification is currently carried by submitting jobs to the [JPRED](http://www.compbio.dundee.ac.uk/jpred/) server 
* Coiled-Coil classification is carried out with [Deepcoil](https://github.com/labstructbioinf/DeepCoil). This needs to be installed locally.
* Transmembrane classification is carried out with [TMHMM](https://github.com/dansondergaard/tmhmm.py). This needs to be installed locally. 

Deepcoil and TMHMM executables can be specified with:
```bash
--TMHMM_exe <PATH TO TMHMM EXECUTABLE>
--deepcoil_exe <PATH TO DEEPCOIL EXECUTABLE>
```

