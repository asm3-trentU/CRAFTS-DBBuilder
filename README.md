CRAFTS Lab AI-MS Database Builder (DB)

Program for generating R data table (.RDS) databases for use 
with the NIST/NIJ DART-MS DIT (v3.22).

Developed by: Arun Moorthy - arunmoorthy@trentu.ca

Adapted from code originally developed by Arun and Edward Sisco (edward.sisco@nist.gov) as presented as supplemental material in [Sisco et. al (2021). JASMS 32(3), 685-689](https://pubs.acs.org/doi/10.1021/jasms.0c00416). See original source code at https://github.com/asm3-nist/DART-MS-DBB.

==============================================================================

Notes:

There are two primary scripts in this repository: (1) Basic_DB-BuilderScript.R and (2) Full_DB-BuilderScript.R. As their names suggest, the first script allows a user to build a basic database and the second allows a user to build a more complete database (but requires more information). 

Steps to building a library:
1. Download and extract repository into a folder on your computer. We'll refer to this folder as the CRAFTS-DBBuilder directory.

2. Add a folder of mass spectra (measured or simulated) you'd like to compile into a library to the CRAFTS-DBBuilder directory. The spectra must be split into three subfolders representing measurements of the sample at three different "fragmentation levels". The spectra of the same compound must use the same filename/code in all three folders. At present, the program accepts folders named "1,2,3" or "LOW,MID,HIGH" or "30,60,90". An example folder of spectra is provided in the repository as a model: "Example-Library_DARTM-MassSpectra". 

3. Create a spreadsheet with metadata information for the compounds being included in the library. For a basic library, please include the code/filename without extension, name of the compound, and a target mass-to-charge ratio that can be used by the DIT during library searching. For building a full library, additional information must be included. Example spreadsheets are provided in the repository as a model: "BASIC_DART-ExampleLibrary_Metadata.xlsx" and "FULL_DART-ExampleLibrary_Metadata.xlsx".

4. If creating a full library, install a system appropriate version of [Open Babel](https://openbabel.org/docs/index.html), and then run "Full_DB-BuilderScript.R" in R. Follow the command line prompts and make sure you select the appropriate "master file" and folder for your library. If creating a basic library, simply run "Basic_DB-BuilderScript.R" in R. 

For more information about this project or if you'd like help modifying this code, please contact Arun. To view a video describing the set up and use of the database builder scripts, please see the [CRAFTS Lab youtube channel](https://youtube.com/@craftslab-trentu?si=MC6UjEBhxWv0n4yo). 
