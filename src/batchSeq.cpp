/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    batchSEQ - a simple batch program for automated downloading and processing of raw sequencing reads
    into sorted BAM files for downstream applications

    Requires the NCBI SRA toolkit, HISAT2, and Samtools

    (c) 2018 The Hope Cancer Lab, McMaster University

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

    const size_t mP = 1000; // Get the executable directory
    char contextPath[mP];

    #if defined WIN32 || defined __WIN32
    cout << "OS: WINDOWS\n";
    GetCurrentDirectory(mP, contextPath);
    #elif defined __APPLE__
    _NSGetExecutablePath(path, &mP);
    #endif

    string accFileName, // The filename for the accession IDs as an input to this program
                genomeIndex, // The hisat2 genome index to align the reads to
                directory; // The directory for output

    unsigned int numThreads = 4;

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

    if(directory.empty()){
        directory = contextPath; // Use the directory of the executable if not defined by user
        #ifdef __APPLE__
        while(directory.back != '/') directory.pop_back(); // OSX path will contain the executable name - remove this
        #endif
    }

    for(auto& c : directory){ // Forward slash is universally understood as a path branch delimiter
        if(c == '\\') c = '/';
    }

    if(accFileName.empty())
        accFileName = directory + "/SraRunTable.txt"; /* By default, the accession filename
                                            is the default downloaded from the NCBI trace database for a bioproject */
    if(genomeIndex.empty())
        genomeIndex = "C:/hisat2/index/grch38/genome"; // Requires hisat2 to be installed in the base drive directory

    if(access(accFileName.c_str(), X_OK)){
        cout << "Could not find accession file\n";
        return -1;
    }
    else{
        cout << "Loading accession IDs from " << accFileName << "...\n";
    }

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

    if(accIDs.size() > 0){
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
                    sam_FN,
                    bam_FN,
                    sample_bam_FN,
                    sample_sorted_FN,
                    acc_ID;

        unsigned int accIndex = 0;
        for(; accIndex < accIDs.size(); ++accIndex){ // Pre-check for processed files if resuming
            bam_FN = directory + "/" + accIDs[accIndex] + ".bam";
            sample_bam_FN = directory + "/" + sampleIDs[accIndex] + ".bam";

            if(!access(sample_bam_FN.c_str(), X_OK)){
                if(!access(bam_FN.c_str(), X_OK)){ // Detect if a bam file is present for this sample indicating partial processing
                    remove(bam_FN.c_str()); // Stop at this accession index and reprocess the BAM file
                    break;
                }
            }
        }

        for(; accIndex < accIDs.size(); ++accIndex){ // Download and process the reads one accession at a time

            acc_ID = accIDs[accIndex];

            fastq_1_FN = directory + "/" + acc_ID + "_1.fastq";
            fastq_2_FN = directory + "/" + acc_ID + "_2.fastq";
            sam_FN = directory + "/" + acc_ID + ".sam";
            bam_FN = directory + "/" + acc_ID + ".bam";
            sample_bam_FN = directory + "/" + sampleIDs[accIndex] + ".bam";

            cout << "Obtaining read data for " << acc_ID << "... " << endl;
            strbuf << "fastq-dump -I --split-files " << acc_ID << " --outdir " << directory;
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

                strbuf << "perl C:/hisat2/hisat2 -q -p " << to_string(numThreads) << " -x "
                        << genomeIndex << " -1 " << fastq_1_FN << " -2 "
                        << fastq_2_FN << " -S " << sam_FN;

                system(strbuf.str().c_str());
                strbuf.str("");
                cout << endl;
            }

            if(!access(fastq_1_FN.c_str(), X_OK)){
                remove(fastq_1_FN.c_str());
            }
            if(!access(fastq_2_FN.c_str(), X_OK)){
                remove(fastq_2_FN.c_str());
            }

            cout << "Converting to *.BAM format ... " << endl;
            strbuf << "samtools view -bS " << sam_FN << " > " << bam_FN;
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
            else
                cout << "Successfully produced *.BAM format for accession " << acc_ID << endl;

        }

        if((accIndex == sampleIDs.size() - 1) ||
                    (sampleIDs[accIndex] != sampleIDs[accIndex+1])){ // Merge and sort if all BAM files for this sample have been accounted for

            vector<string> parts; // Collect part accession IDs for this sample
            int i = accIndex;
            cout << "Merging:\n";
            while((i >= 0) && (sampleIDs[i] == sampleIDs[accIndex])){
                parts.push_back(directory + "/" + accIDs[i] + ".bam");
                if(access(parts.back().c_str(), X_OK)){
                    cout << "Missing part: " << parts.back() << "\n";
                    cout << "Unable to complete sample " << sampleIDs[accIndex] << endl;
                    return -1;
                }
                cout << '\t' << parts.back() << '\n';
                --i;
            }
            cout << "Into: " << sample_bam_FN << endl;

            strbuf << "samtools merge " << sample_bam_FN;
            for(auto& p : parts){
                strbuf << " " << p;
            }
            system(strbuf.str().c_str()); // Perform merge using samtools
            strbuf.str("");

            if(!access(sample_bam_FN.c_str(), X_OK)){
                cout << "Successfully created sample *.BAM file:\n" << sample_bam_FN << endl;
                for(auto& p : parts){ // Clean up parts
                    if(!access(p.c_str(), X_OK)){
                        remove(p.c_str());
                    }
                }

                sample_sorted_FN = directory + "/" + sampleIDs[accIndex] + "_sorted.bam";

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
        }

    }
    else{
        cout << "Could not find accession IDs in " << accFileName << endl;
        return -1;
    }

    return 0;
}
