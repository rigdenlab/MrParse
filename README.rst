***********************************************************
MrParse: an aid to decision making in Molecular Replacement
***********************************************************

.. image:: https://readthedocs.org/projects/mrparse/badge/?version=documentation
   :target: https://mrparse.readthedocs.io/en/documentation/?badge=documentation
   :alt: Documentation Status

About
+++++

MrParse is a `CCP4 <http://www.ccp4.ac.uk>`_ program takes a protein amino acid sequence file and searches for homologs using `PHMMER <http://hmmer.org/>`_ (default) or `HHSEARCH <https://github.com/soedinglab/hh-suite>`_. If supplied with a reflection data file (currently in `MTZ <http://www.ccp4.ac.uk/html/mtzformat.html>`_ format), it can then use `PHASER <https://www.phaser.cimr.cam.ac.uk/index.php/Phaser_Crystallographic_Software>`_ to calculate the Expected Log Likelihood Gains (eLLG) values for the homologs. It also searches the `EBI AlphaFold database <https://alphafold.ebi.ac.uk/>`_ for related models.

It also attempts to classify the sequence according to its secondary structure, and whether any regions are expected to be Coiled-Coil or Transmembrane.

Results are currently displayed in a simple HTML webpage that is rendered using `VUE <https://vuejs.org>`_. The sequence graphics are created using the `PFAM graphics library <https://pfam.xfam.org/generate_graphic>`_, a copy of which is distributed with this code.

Installation
++++++++++++

MrParse is distributed with CCP4, although optional software can be installed to get the most out of MrParse. Full details are provided `here <https://mrparse.readthedocs.io/en/documentation/install.html>`_


Simple command line
+++++++++++++++++++

.. code-block:: bash

   mrparse --seqin <PATH TO SEQUENCE FILE>

To provide a reflection file and classify the sequence we can provide the following optional flags:

.. code-block:: bash

   --hklin <PATH TO MTZ FILE>
   --do_classify

Search Model Finder
+++++++++++++++++++
The search model finder by default uses `PHMMER <http://hmmer.org/>`_ (distributed with `CCP4 <http://www.ccp4.ac.uk>`_) to search for homologs. If installed you can also use `HHSEARCH <https://github.com/soedinglab/hh-suite>`_.
Examples of how to use MrParse are provided `here <https://mrparse.readthedocs.io/en/documentation/examples.html>`_


EBI Alphafold database search
+++++++++++++++++++++++++++++
The search model finder currently uses `PHMMER <http://hmmer.org/>`_ (distributed with `CCP4 <http://www.ccp4.ac.uk>`_) to search the EBI Alphafold database. This will be replaced with the `3Dbeacons <https://github.com/3D-Beacons>`_ API when it becomes available.

Classifiers
+++++++++++
* secondary structure classification is currently carried by submitting jobs to the `JPRED <http://www.compbio.dundee.ac.uk/jpred/>`_ server.
* If installed, coiled-Coil classification is carried out with `Deepcoil <https://github.com/labstructbioinf/DeepCoil>`_.
* If installed, transmembrane classification is carried out with `TMHMM <https://github.com/dansondergaard/tmhmm.py>`_.

