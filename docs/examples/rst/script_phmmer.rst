.. _script_phmmer:

Performing an MrParse search using Phmmer
-----------------------------------------

The simplest way to run MrParse is with Phmmer. Phmmer is distributed with the CCP4 software suite and therefore is the default search engine used by MrParse.

Input
+++++

The minimum requirement for running MrParse is a sequence file. To run MrParse from the command line with just a sequence file you can use the following command:

.. code-block:: bash

   mrparse --seqin <PATH TO SEQUENCE FILE>

If a reflection file is available for your data, MrParse can calculate the eLLG for each of the hits identified by Phmmer.
eLLG has been shown to be a better indicator of whether a search model will work in MR than sequence identity and so, if available, we recommend providing one to MrParse by adding the following flag:

.. code-block:: bash

   --hklin <PATH TO REFLECTION FILE>


Additionally, you can also elect to perform sequence classification. MrParse will attempt to classify the sequence according to it's secondary structure and whether any regions are predicted to be coiled-coil or transmembrane.
Sequence classification can be enabled by providing the following flag:

.. code-block:: bash

   --do_classify

.. note::

   To get the most out of MrParse's classification report you may wish to install additional software, see :ref:`installation`.

.. note::

  A full list of command line flags for MrParse can be found here: :ref:`mrparse_options`.


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