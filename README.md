# streamlit-pacbiohifi
- a streamlit component for reading and plotting the PacBioHiFi sequencing reads.
- new components added such as
   -- filter read length and write filtered fasta reads.
   -- read visualization pre and post filtering.
   -- pattern matching.
   -- specific region extractor.

To run the application, use 
```
conda activate streamlit
streamlit run streamlitpacbiohifi.py 
```
- provide the path and select the options which you need. If you have bigger files then add .toml file in the directory with the size specified. 
- These are indiviual components and a multi page application in development mode, which will give you all access to PacBioHifi from sequence to graphs.
- The app will open like this. You can paste the file path and fetch the names, sequences and plot the lengthe distribution. it supports both the fastq and the fasta PacBioHifi and if the write option is selected, it will write the fasta files also. 

Gaurav Sablok \
University of Potsdam \
Potsdam,Germany


