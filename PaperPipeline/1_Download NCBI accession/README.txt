Folder Name
1. Download NCBI accession


Description

The folder '1. Download NCBI accession' is mainly for downloading virus data according to their NCBI accessions, including 
CDS sequences, genomic sequences, and taxonomy information. Those data are later processed into Ready-to-go data (RTG_Data)
before using them for training models.

Ready-to-go data (RTG_Data) includes:
	AAU.csv          - Amino acids usages (or percentages)
	ATGC.csv         - ATGC percentages
	CDS_LENGTH.csv   - CDS Length data
	CU.csv           - Codon usages (or percentages)
	HOST.csv         - Host labels (This data is extracted for y data in Training)
	HS_CORR_AA.csv   - Amino-acid-wise correlation coefficient to Homo Sapiens RSCU
	HS_CORR.csv      - Correlation coefficient to Homo Sapiens RSCU
	INFO.csv         - Informations of viruses (NCBI accessions)
	PARTITE.csv      - Partite of viruses
	RSCU.csv         - Relative Synonymous Codon Usage
	START_STOP.csv   - Start codons and stop codons related data
	TAXONOMY_RAW.csv - Taxonomy data before encoding
	TAXONOMY.csv     - Taxonomy data after encoding (1, 0, -1)


Origin Data Download:

	The accessions IDs of virus genomes RefSeq were downloaded from NCBI Virus Genomes Resource https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi, 
	and the file is renamed as 'Viral_genome_browser.csv'. 

	The reference human codon usage was downloaded from https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs. 


Script files:

1 Virus_CU_AAU_RSCU.ipynb

	This code file is for downloading data from NCBI website and process the downloaded sequence data into RSCU-related and other
	Conten-CDS related data before using for training models.  

	Contents:

		1. Downloaded csv file from NCBI viral genome browser to ALL_VIRUS_DATABASE.pickle
		2. ALL_VIRUS_DATABASE.pickle to ALL_VIRUS_COMBINED_CDS.pickle
		3. ALL_VIRUS_COMBINED_CDS.pickle to ALL_VIRUS_CODON_USAGE.pickle
		4. ALL_VIRUS_COMBINED_CDS.pickle to ALL_VIRUS_AA_USAGE.pickle
		5. ALL_VIRUS_CODON_USAGE.pickle to ALL_VIRUS_RSCU.pickle
		6. ALL_VIRUS_COMBINED_CDS.pickle to ALL_VIRUS_ATGC.csv
		7. ALL_VIRUS_CODON_USAGE.csv / ALL_VIRUS_RSCU.csv / ALL_VIRUS_AA_USAGE.csv to ALL_VIRUS_Human_Corr.csv
		8. ALL_VIRUS_CODON_USAGE.csv / ALL_VIRUS_RSCU.csv / ALL_VIRUS_AA_USAGE.csv to ALL_VIRUS_Human_Corr_AA.csv


2 Other_details_organise.ipynb

	This code file is for processing the downloaded sequence data (from NCBI) into other data before using for training models.  

	Contents:
	
		1. ALL_VIRUS_DATABASE.pickle to ALL_VIRUS_DATABASE.csv and ALL_VIRUS_info.csv
		2. ALL_VIRUS_DATABASE.pickle to ALL_VIRUS_HOST.csv


3 Single_CDS_CU_AAU_RSCU.ipynb

	This code file is for processing the downloaded sequence data (from NCBI) into single-CDS related databefore using for training models.  

	Contents:
		
		1. ALL_VIRUS_DATABASE.pickle to ALL_VIRUS_ALL_CDS_DETAIL.pickle
		2. ALL_VIRUS_ALL_CDS_DETAIL.pickle to ALL_VIRUS_ALL_CDS.pickle
		3. ALL_VIRUS_ALL_CDS.pickle to ALL_VIRUS_ALL_CDS_SIMPLE.pickle
		4. ALL_VIRUS_ALL_CDS_SIMPLE.pickle to ALL_VIRUS_CDS_Length.csv and ALL_VIRUS_CDS_error.csv
		5. ALL_VIRUS_ALL_CDS_SIMPLE.pickle to ALL_VIRUS_start_stop.csv


4 Organise_data_before_use.ipynb

	This code file is for organising above data into Ready-to-go data (RTG_Data) before using them for training models


5 Download_Taxonomy_data.ipynb

	This code file is for downloading Taxonomy related data from NCBI and organising them into Ready-to-go data (RTG_Data).

	Contents:

		1. Search and Download NCBI Taxonomy IDs of viruses
		2. Manually search NCBI Taxonomy browser to correct the wrong ids
		3. Download Taxonomy Data from NCBI according to NCBI Taxonomy IDs
		4. Encode Taxonomy data into Train-Ready data

		
Authors

	Shuquan Su

		Email:
			Shuquan.Su@student.uts.edu.au

		Address:
			Data Science Institute (DSI), 
			School of Computer Science, 
			Faculty of Engineering and Information Technology (FEIT), 
			University of Technology Sydney (UTS)
			Sydney, New South Wales, Australia

