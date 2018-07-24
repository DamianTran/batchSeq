# batchSEQ
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    batchSEQ - a simple batch program for automated downloading and processing of raw sequencing reads
    into sorted BAM files for downstream applications

    Requires the NCBI SRA toolkit, HISAT2, and Samtools

    (c) 2018 The Hope Cancer Lab, McMaster University

   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   
   REQUIREMENTS:
   
   	- Strawberry Perl >= 5.28.0
	- Mingw64 and development toolchains (on Windows)
	- HISAT2
	- Samtools and Htslib
	- NCBI SRA Toolkit (for your OS)
   
   BEFORE RUNNING:

        - In bash, set the PATH environment variable as follows (no spaces after the one following "export"):
                export PATH=$PATH:/c/path/to/samtools/:/c/path/to/SRA/bin/
		ie.    export PATH=$PATH:/c/samtools-develop/:/c/SRA/bin
		-This has to be done once on each new bash terminal startup
        - Make sure a HISAT2 genome index is downloaded from the HISAT site
            - Recommended grch38_snp_tran to capture all possible transcripts
        - If this is not in the C:/hisat2/index/ directory, you must specify it in the command line
            by adding -gen /c/path/to/hisat2/index/grch38_snp_tran/genome
	- Download a RunInfo table from the BioProject of interest from:
		ncbi.nlm.nih.gov/Traces/study/
		-Leave the filename as is (SraRunTable.txt) and move it to the directory
		 where you want to download and process the reads

   PARAMS:

        -dir [directory]: The program will read from the SraRunTable file in the directory
                        specified by -dir and download, process, cleanup, merge, and sort
                        BAM files for every detected sample.
        -acc [file path]: If the SraRunTable is not in the same directory as the one you
                        want this BAM files to be saved and processed in, you will need to specify it in
                        this way.  Otherwise don't use this argument.
        -gen [directory ending in filename base]: The indexes downloaded from the John Hopkins servers
                        contain a common filename such as genome.ht2.1, genome.ht2.2, etc.  Specify the
                        location that this base filename is located so that batchSeq can perform the transcript
                        alignments.
                                                                                                        
    EXECUTION:
	2 ways:

	1: Run by copying into the directory where you'd like to
           create the database, and entering in the bash terminal:
		export PATH=$PATH:/c/path/to/samtools/:/c/path/to/SRA/bin/
		cd /c/the/directory/
		./batchSeq -gen [genome index]

	2: Run by entering in the bash terminal:
		export PATH=$PATH:/c/path/to/samtools/:/c/path/to/SRA/bin/
		cd /c/the/directory/with/batchSeq
		./batchSeq -dir [directory] -gen [genome index]
