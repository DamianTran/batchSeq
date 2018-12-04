/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



    batchSEQ - a simple batch program for automated downloading and processing of raw sequencing reads

    into sorted BAM files for downstream applications



    Requires the NCBI SRA toolkit, HISAT2, and Samtools



    (c) 2018 The Hope Cancer Lab, McMaster University


    DISCLAIMER: This software is provided AS IS with no otherwise implied or guaranteed warranties. The

    end user assumes all potential risks involved with the use of this software, and has read and

    acknowledged the terms included in the original repository's 'licence.txt.'


    Use of this code is free and permissible for personal, academic, and commercial purposes including

    redistribution.


   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



   BEFORE RUNNING:



        - In bash, set the PATH environment variable as follows:

                export PATH=$PATH:"/c/path/to/samtools/":"/c/path/to/SRA/bin/"

        - Make sure a HISAT2 genome index is downloaded from the HISAT site

            - Recommended grch38_snp_tran to capture all possible transcripts

        - If this is not in the C:/hisat2/index/ directory, you must specify it in the command line

            by adding -gen /c/path/to/hisat2/index/grch38_snp_tran/genome



   PARAMS:



        -dir [directory]: The program will read from the SraRunTable file in the directory

                        specified by -dir and download, process, cleanup, merge, and sort

                        BAM files for every detected sample.

        -acc [file path]: If the SraRunTable is not in the same directory as the one you

                        want this BAM files to be saved and processed in, you will need to specify it in

                        this way

        -gen [directory ending in filename base]: The indexes downloaded from the John Hopkins servers

                        contain a common filename such as genome.ht2.1, genome.ht2.2, etc.  Specify the

                        location that this base filename is located so that batchSeq can perform the transcript

                        alignments

                                                                                                        */



#include <stdio.h>

#include <string>

#include <string.h>

#include <stdlib.h>

#include <unistd.h>

#include <iostream>

#include <fstream>

#include <vector>

#include <sstream>

#include <chrono>

#include <thread>



#if defined WIN32 || defined __WIN32

#include <windows.h>

#elif defined __APPLE__

#include <mach-o/dyld.h>

#endif



using namespace std;



int main(int argc, char** argv){



    vector<thread> allThreads; // Thread registry



    unsigned int mP = 1000; // Get the executable directory

    char contextPath[mP];



    #if defined WIN32 || defined __WIN32

    cout << "OS: WINDOWS\n";

    GetCurrentDirectory(mP, contextPath);

    #elif defined __APPLE__

    _NSGetExecutablePath(contextPath, &mP);
	
    cout << "OS: APPLE\n";

    #endif



    string accFileName, // The filename for the accession IDs as an input to this program

                genomeIndex, // The hisat2 genome index to align the reads to

                directory, // The directory for output

                homeDir = contextPath;



    for(size_t i = 0; i < homeDir.size(); ++i){ // Forward slash is universally understood as a path branch delimiter

        if(homeDir[i] == '\\') homeDir[i] = '/';

    }



    #ifdef __APPLE__

    while(homeDir.back() != '/') homeDir.pop_back(); // OSX path will contain the executable name - remove this

    #endif



    directory = homeDir; // By default, the directory is the path of the executable

    genomeIndex = "C:/hisat2/index/grch38/genome"; // Requires hisat2 to be installed in the base drive directory



    unsigned int numThreads = 4; // Default thread count of 4



    for(size_t i = 0; i < argc; ++i){ // Process the arguments

        if(string(argv[i]) == "-acc"){ // We require an argument with the "-acc" prefix for the accession IDs

            if(i+1 < argc){

                if(strchr(argv[i+1], '-') == NULL){

                    accFileName = argv[i+1];

                    ++i;

                }

            }

        }

        else if(string(argv[i]) == "-gen"){ // Obtain the genome index if given

            if(i+1 < argc){

                if(strchr(argv[i+1], '-') == NULL){

                    genomeIndex = argv[i+1];

                    ++i;

                }

            }

        }

        else if(string(argv[i]) == "-dir"){ // Obtain the desired output directory

            if(i+1 < argc){

                if(strchr(argv[i+1], '-') == NULL){

                    directory = argv[i+1];

                    ++i;

                }

            }

        }

        else if(string(argv[i]) == "-threads"){ // Obtain the desired output directory

            if(i+1 < argc){

                if(strchr(argv[i+1], '-') == NULL){

                    try{

                        numThreads = (uint32_t)strtod(argv[i+1], NULL);

                    }catch(...){ }

                    ++i;

                }

            }

        }

    }



    if(accFileName.empty()){

        accFileName = directory + "/SraRunTable.txt"; /* By default, the accession filename

                                            is the default downloaded from the NCBI trace database for a bioproject */

    }



    for(size_t i = 0; i < directory.size(); ++i){

        if(directory[i] == '\\') directory[i] = '/';

    }



    while(directory.back() == '/') directory.pop_back(); // Format to fit algorithms



    system(std::string("cd ").append(homeDir).c_str());



    if(access(accFileName.c_str(), X_OK)){

        cout << "Could not find accession file\n";

        return -1;

    }

    else{

        cout << "Loading accession IDs from " << accFileName << "...\n";

    }

    // Read the accession file

    FILE* accFILE = fopen(accFileName.c_str(), "rb");

    size_t fileSIZE;

    fseek(accFILE, 0, SEEK_END);

    fileSIZE = ftell(accFILE);

    fseek(accFILE, 0, SEEK_SET);



    char fileData[fileSIZE];

    fread(fileData, fileSIZE, sizeof(char), accFILE);

    vector<vector<string>> dataMatrix(1, vector<string>());

    size_t lastIndex = 0;



    for(size_t i = 0; i < fileSIZE; ++i){

        if((fileData[i] == '\t') || (fileData[i] == ',')){

            dataMatrix.back().emplace_back();

            dataMatrix.back().back().assign(fileData + lastIndex, fileData + i);

            lastIndex = i+1;

        }

        else if(fileData[i] == '\n'){

            dataMatrix.back().emplace_back();

            dataMatrix.back().back().assign(fileData + lastIndex, fileData + i);

            dataMatrix.emplace_back();

            lastIndex = i+1;

        }

    }



    // Interpret the table



    vector<string> accIDs,

                  sampleIDs;



    cout << "Detected table with " << dataMatrix.size()-1 << " entries" << endl;



    for(size_t i = 0; i < dataMatrix[0].size(); ++i){ // Check the header for the appropriate columns

        if(dataMatrix[0][i] == "Run"){

            for(size_t j = 1; j < dataMatrix.size(); ++j){

                if(dataMatrix[j].size() > i){

                    accIDs.push_back(dataMatrix[j][i]);

                }

            }

        }

        else if(dataMatrix[0][i] == "Sample_Name"){

            for(size_t j = 1; j < dataMatrix.size(); ++j){

                if(dataMatrix[j].size() > i){

                    sampleIDs.push_back(dataMatrix[j][i]);

                }

            }

        }

    }



    if(accIDs.size() < 1){

	cout << "Could not find accession IDs in " << accFileName << endl;

        return -1;

    }



	cout << "Found accession IDs:\n";

        unsigned int maxPrint = 30, printIndex = 0;

        for(auto& ID : accIDs){

            cout << '\t' << ID;

            if(printIndex < sampleIDs.size()) cout << " from sample " << sampleIDs[printIndex] << '\n';

            ++printIndex;

            if(printIndex == maxPrint) break;

        }

        if(accIDs.size() > maxPrint){

            cout << "\t... " << accIDs.size() - maxPrint << " others\n";

        }

        cout << endl;



        stringstream strbuf;

        string fastq_1_FN,

                    fastq_2_FN,

                    SRA_FN,

                    sam_FN,

                    bam_FN,

                    sample_bam_FN,

                    acc_ID,

                    sample_ID;



        unsigned int accIndex = 0;



        for(; accIndex < accIDs.size(); ++accIndex){ // Download and process the reads one accession at a time



            sample_ID = sampleIDs[accIndex];

            sample_bam_FN = directory + "/" + sample_ID + ".bam";



	    while(!access(sample_bam_FN.c_str(), X_OK)){

                cout << "Detected assembled *.BAM file for " << sample_ID << ": skipping" << endl;

                while(sample_ID == sampleIDs[accIndex]) ++accIndex;

		sample_ID = sampleIDs[accIndex];

                sample_bam_FN = directory + "/" + sample_ID + ".bam";

            }



	    acc_ID = accIDs[accIndex];

	    fastq_1_FN = directory + "/" + acc_ID + "_1.fastq";

            fastq_2_FN = directory + "/" + acc_ID + "_2.fastq";

	    SRA_FN = "~/ncbi/public/sra/" + acc_ID + ".sra";

            sam_FN = directory + "/" + acc_ID + ".sam";

            bam_FN = directory + "/" + acc_ID + ".bam";



	    if(!access(bam_FN.c_str() ,X_OK)){

		if((accIndex == accIDs.size() - 1) || // Assume *.bam files not at the end of the list are complete

			!access(std::string(directory).append("/").append(accIDs[accIndex+1]).append(".bam").c_str(), X_OK)){

                    cout << "Detected fragment *.BAM file for accession " << acc_ID << endl;

		    goto checkMerge;

                }

                else{

                    cout << "Found incomplete fragment " << acc_ID << " for sample " << sample_ID << ": resuming" << endl;

                }

	    }



            cout << "Obtaining read data for " << acc_ID << "... " << endl;

            strbuf << "bash wonderdump.sh -I --split-files --outdir " << directory << " " << acc_ID;

            system(strbuf.str().c_str());

            strbuf.str("");

            cout << endl;



            if(!access(fastq_1_FN.c_str(), X_OK) &&

               !access(fastq_2_FN.c_str(), X_OK)){



                cout << "Aligning paired-end read data into " << acc_ID << ".sam ..." << endl;



                /*

                    HISAT2 perl script with params:

                        -q [fastq input]

                        -p 4 [4 threads]

                        -x [genome index]

                        -1 [forward read file]

                        -2 [reverse read file]

                        -S [output file]

                */



                strbuf << "perl C:/hisat2/hisat2 -q -p " << numThreads << " -x "

                        << genomeIndex << " -1 " << fastq_1_FN << " -2 "

                        << fastq_2_FN << " -S " << sam_FN;



                system(strbuf.str().c_str());

                strbuf.str("");

                cout << endl;

            }



            if(!access(fastq_1_FN.c_str(), X_OK)){ // Clean up FASTA files

                remove(fastq_1_FN.c_str());

            }

            if(!access(fastq_2_FN.c_str(), X_OK)){

                remove(fastq_2_FN.c_str());

            }

	    if(!access(SRA_FN.c_str(), X_OK)){     // Clean up SRA file

		remove(SRA_FN.c_str());

	    }



            cout << "Converting to *.BAM format ... " << endl;

            strbuf << "samtools view -bS --threads " << numThreads << ' ' << sam_FN << " > " << bam_FN;

            system(strbuf.str().c_str());

            strbuf.str("");

            cout << endl;



            if(!access(sam_FN.c_str(), X_OK)){

                remove(sam_FN.c_str());

            }



            if(access(bam_FN.c_str(), X_OK)){

       	         cout << "Failed to produce *.BAM format for accession " << acc_ID << endl;

        	 cout << "Exiting...";

       	         break;

       	    }

            else{

        	 cout << "Successfully produced *.BAM format for accession " << acc_ID << endl;

	    }



	checkMerge:;



        if((accIndex == sampleIDs.size() - 1) ||

                    (sample_ID != sampleIDs[accIndex+1])){ // Merge and sort if all BAM files for this sample have been accounted for



	    allThreads.emplace_back([&, directory, accIDs, accIndex, sampleIDs, sample_ID,

					acc_ID, sample_bam_FN, numThreads](){ // Add a new thread to merge completed fragments

										// While main thread continues on



	    std::stringstream strbuf;



	    std::string sample_sorted_FN = directory + "/" + sample_ID + "_sorted.bam";



            vector<string> parts; // Collect part accession IDs for this sample

            int i = accIndex;

            cout << "Merging:\n";

            while((i >= 0) && (sampleIDs[i] == sample_ID)){

                parts.insert(parts.begin(), directory + "/" + accIDs[i] + ".bam");

                if(access(parts.front().c_str(), X_OK)){

                    cout << "Missing part: " << parts.front() << "\n";

                    cout << "Unable to complete sample " << sample_ID << endl;

                    return -1;

                }

                --i;

            }



            strbuf << "samtools merge --threads " << numThreads << " " << sample_bam_FN;

	    for(auto& part : parts){

		cout << '\t' << part << endl;

		strbuf << " " << part;

	    }

            cout << "Into: " << sample_bam_FN << endl;



            system(strbuf.str().c_str()); // Perform merge using samtools

            strbuf.str("");



            if(!access(sample_bam_FN.c_str(), X_OK)){

                cout << "Successfully created sample *.BAM file:\n" << sample_bam_FN << endl;

                for(auto& part : parts){ // Clean up parts

                    if(!access(part.c_str(), X_OK)){

                        remove(part.c_str());

                    }

                }



                strbuf << "samtools sort " << sample_bam_FN << " -o " << sample_sorted_FN;

                cout << "Sorting: " << sample_bam_FN << endl;

                system(strbuf.str().c_str());

                strbuf.str("");

                cout << endl;



                if(!access(sample_sorted_FN.c_str(), X_OK)){

                    remove(sample_bam_FN.c_str());

                    rename(sample_sorted_FN.c_str(), sample_bam_FN.c_str());

                }

            }

            else{

                cout << "Unable to create sample *.BAM file:\n" << sample_bam_FN << endl;

                return -1;

            }



	    });

        }

    }



    for(auto& thread : allThreads){ // Wait for all merge threads to finish

	thread.join();

    }



    cout << "Processing complete" << endl;



    return 0;

}
