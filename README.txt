################################################################################################
  ___                 _                           _  _         _                          _   
 / __|  _  _   _ _   | |_   ___   _ _    _  _    | \| |  ___  | |_  __ __ __  ___   _ _  | |__
 \__ \ | || | | ' \  |  _| / -_) | ' \  | || |   | .` | / -_) |  _| \ V  V / / _ \ | '_| | / /
 |___/  \_, | |_||_|  \__| \___| |_||_|  \_, |   |_|\_| \___|  \__|  \_/\_/  \___/ |_|   |_\_\
        |__/                             |__/                                                 
################################################################################################

The SYNTENY NETWORK tool is a program to score the synteny of candidate 
sequences in relation to a list of known homologous sequences. Therefore, 
evolutionary related candidate sequences will get a high synteny score, whilst
sequences with no homologous gene neighborhood will get a low score.

The workflow was originally developed to work with the output of a GLASSgo sRNA
search as list of trustable hologous sequences. Further the synteny of unknown 
sequences e.g. from a BLAST search can be scored.  To score sequences the header 
in the input FASTA file needs to start with the NCBI accession number followed by 
a : and the starting position of the hit in the genome. 

SYNTENY NETWORK uses the PageRank algorithm combined with a modified version of the
Bellman-Ford algorithm to calculate synteny values. These Values are added to the end 
of sequence headers in the output FASTA files



Sequence Input:

-i --network_file
		fasta file containing sequences used for network construction
		sequences are written into outpath + _network.fasta with the
		synteny value of each sequence added to the header.
		Header of each input sequence needs to have following form:
			"NCBI ACC number":"#-starting nucleotide"-"#-end nucleotide"
		
-t --test_file
		Input is optional. 
		fasta file containing candidate sequences that are checked for 
		network match sequences are written into outpath + _questionable.fasta
		with the synteny value of each sequence added to the header. Header of 
		each input sequence needs to have following form:
			"NCBI ACC number":"#-starting nucleotide"-"#-end nucleotide"
		
-o --outfiles
		path and name where all outfiles are getting stored. Folder is NOT 
		created if it is not present
		
-c --cluster_script
		path to the R synteny extraction script
		
-d --w_dir
		path where temporary files are stored. default is the current directory
		Folder is not created if it is not present
	
-n --network
		produces a network outfile when set to cys or svg.
		svg will produce a svg image of the network used for synteny value 
		calculation
		cys produces a comma seperated text file that can be viewed in cytoscape
		the cys file has the following structure:
			"cluster,connected cluster,PageRank, connection weight"
		both settings will also produce a text file containing each sequence and
		their surrounding clusters
		All these files are written into the outpath

		
-w --synteny_window
		synteny window used for extraction of neighboring proteins. As the
		network will only use 4 proteins in each direction, the synteny window
		should be set between 3000 and 10000 bp. Higher values will increase
		runtime dramatically. default 5000
		- will be removed in further version and replaced by the number of -
			- 		up and downstream proteins that is extracted	-
-- node_normalization
		uses a teleport from the sRNA node that depends on a normalized number 
		of occurrences of the clusters instead of a random teleport if set to True. 
		Default is False
			
Example:

	python SyntenyNetwork2.py -i ./testfiles/sRNA.fasta -t ./testfiles/candidates.fasta  -o synteny
	
	will run the script and add synteny value (SV) for sequences in 
	candidates.fasta and in sRNA.fasta
	The output will be written into the current directory and named:
		synteny_network.fasta 		(sRNA.fasta with SV in header)
		synteny_questionable.fasta 	(candidates.fasta with SV in header) 
		
	python SyntenyNetwork2.py -i sRNA.fasta -t candidates.fasta  -o synteny -n svg
	
	will additionally produce a svg image of the network in the current folder
	named:
		synteny_network.svg
		
