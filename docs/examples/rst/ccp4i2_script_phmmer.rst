.. _ccp4i2_script_phmmer:

Performing an MrParse search in CCP4i2 using Phmmer
---------------------------------------------------

The simplest way to run MrParse is with Phmmer. Phmmer is distributed with the CCP4 software suite and therefore is the default search engine used by MrParse.

Input
+++++

MrParse can be found under the 'Bioinformatics including model preparation for Molecular Replacement' menu in the CCP4i2 GUI:

.. figure:: ../images/ccp4i2_mrparse.png
   :width: 50%
   :align: center

Opening MrParse will bring you to the following menu:

.. figure:: ../images/ccp4i2_mrparse_run.png
   :width: 50%
   :align: center

To run MrParse, all you require is a sequence file. If a reflection file is provided, MrParse will also calculate the eLLG for each of the hits identified by Phmmer.
eLLG has been shown to be a better identifier of whether a search model will work in MR than sequence identity and so we recommend providing a reflection file if one is available to you.
Additionally, you can elect to perform sequence classification. If selected, MrParse will attempt to classify the sequence according to it's secondary structure and whether any regions are predicted to be coiled-coil or transmembrane.

.. note::

   To get the most out of MrParse's classification report you may wish to install additional software, see :ref:`installation`.

Output
++++++

When MrParse finished running, an HTML page will pop up showing the results of the search:

.. figure:: ../images/mrparse_results.png
   :width: 50%
   :align: center

The sections of the MrParse report page are highlighted in different colours:

* In red is information on the input reflection file, including resolution, space group and crystal pathology.
* In teal is information about the PDB entries identified by Phmmer and visualisations of the matches.
* In purple is the protein classification report. This includes a secondary structure prediction, a coiled-coil prediction, and a transmembrane prediction.
* In blue is information about the AlphaFold models identified by Phmmer and visualisations of the matches coloured by pLDDT on an orange to blue scale, where orange indicates very low confidence in the model and blue indicates very high confidence in the model.

Individual models can be downloaded by clicking on the name of the model in the report page. They can also be found by navigating to the homologs/models directory within the MrParse run directory.

.. note::

   Future work is planned to better integrate the output of MrParse with CCP4i2.