Name

4. COVID-19 Origin Tracing


Description

The folder '4. COVID-19 Origin Tracing' is mainly for downloading virus data of USA SARS-CoV-2 according to their NCBI accessions, including 
CDS sequences, genomic sequences. Those data are later processed into Ready-to-go data (RTG_Data), and the RTG_Data is later used in 
evalution of HIP. 


Script files:

Mutation_Path_Simulation_Directional.py

    This code file is for simlutation evolution path from a set of start virus and end viruses using HIP as gradient. All possible codon mutations 
    including subsitution, deletion and addtion are applied to the codon composition of viruses to screen out the one generating highest/lowest HIP
    in each iteration. If multiple identical HIP are found, the one has highest RSCU correlation coefficient to end virus RSCU will be chosen. The 
    simulation will terminated when no improvement in previous 200 steps.  


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