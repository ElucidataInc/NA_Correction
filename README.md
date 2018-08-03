## **About the package**

Corna package performs Natural Abundance Correction on the output intensity of metabolites, obtained after performing Mass Spectrometry(MS).

**Natural Abundance Correction**- Natural abundance (NA) refers to the abundance of isotopes of a chemical element as naturally found on the planet. While performing analysis,the observed intensity contains contribution from isotopic natural abundance that needs to be corrected. This process is referred as NA Correction.

**Pool Total**- Sum of the intensities of all different number of labeled atoms of the isotope element is called pool total.

**Fractional enrichment**- Normalization of intensities of a metabolite between the range of 0 to 1.

## **Installation of Corna**

**Get the Source Code**

Run this command in your terminal to clone the repository:

*git clone git://github.com/ElucidataInc/NA_Correction.git*

**Environment Setup**

1) Create virtual environment:

*virtualenv env-name*

2) Activate virtual environment:

*source env-name/bin/activate*

3) Install requirements:

*pip install -r requirements.txt*

## **Getting Started**

- LCMS raw intensity file should have the following columns-
Name, Label, Formula, Sample_Name

Below is a sample LCMS *raw_intensity* file-
Inside Demo/data_lcms/testfiles/nacorr_test_1.xlsx

- LCMS/MS raw intensity file should have the following columns-
Original Filename, Component Name, Mass Info, Area

Below is a sample LCMS/MS *raw_intensity* file-
Inside Demo/data_msms/tca_data_mq/raw_intensity_file.csv

LCMS/MS metadata file should have the following columns-
Component Name, Unlabeled Fragment, Formula, Isotopic tracer, Parent Formula

Below is a sample LCMS/MS *metadata* file-
Inside Demo/data_msms/tca_data_mq/metadata.xlsx


**Running the Package-**

1) Use the Demo folder to run the demo.py script. The demo folder contains the
   data foler with data files

2) You can also write your own script by referring to the demo.py script and add
   your data files in the working directory

3) Or Try this Ipython Script-
   Inside Demo_ipython folder, find the ipython notbook for LCMS and LCMS/MS workflows.
   Steps to run :-
    1.) Copy the folder to local system.
    2.) *pip install -r requirements.txt* to install dependencies.
    3.) *jupyter notebook <notebook-name>* 

## **How to Contribute**

1) Check for open issues or open a fresh issue to start a discussion around a feature idea or a bug. 

2) Fork the repository on GitHub to start making your changes to the master branch (or branch off of it).

3) Write a test which shows that the bug was fixed or that the feature works as expected.

4) Send a pull request and bug the maintainer until it gets merged and published.

















