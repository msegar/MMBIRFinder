/*
 * birFinder.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: msegar
 */

#include "defs.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <sys/wait.h>
#include <stdlib.h>

using namespace std;

struct t_read{
	int iLength;
	int iStart;
	int iChromosome;
	string sSequence;
};

// birFinder
int startCandidateReads(){
	string sFinalAlignedFile = confDB.getKey("finalAlignedFile").stringVal;

	if (confDB.getKey("performClustering").boolVal == false)
		return 0;

	cout << "\nDB before: " << vCandidateReads.size() << endl;
	fLogFileOut << "\nDB before: " << vCandidateReads.size() << endl;

	// Let's take the final unaligned reads and find the locations of the matched pair reads
	createSplitReadDatabase(sFinalAlignedFile + ".sam");

	cout << "\nDB after: " << vCandidateReads.size() << endl;
	fLogFileOut << "\nDB after: " << vCandidateReads.size() << endl;

	//createParentReads();

	fLogFileOut << "\nRead database time = " << (int)time(NULL)-time0 << endl;

	return 0;
}

int createSplitReadDatabase(string sAlignedFilename){
	cout << "\nCreating candidate read database starting..." << endl;
	fLogFileOut << "\nCreating candidate read database starting..." << endl;
	ifstream input;
	vector<string> curr;
	char row_delim = '\n';
	char field_delim = '\t';
	string sReadName;
	t_consolidated frag;
	int iChr = 0;
	int chromosome = confDB.getKey("chromosome").intVal;
	int iExludedReads = 0;
	int iBadReads = 0;
	int index = 0;

	fLogFileOut << "Aligned file: " << sProjectDirectory + sAlignedFilename << endl;

	// Read in the aligned SAM file and create a database of anchored half reads and their locations
	// We are treating these as candidate reads
	input.open((sProjectDirectory + sAlignedFilename).c_str());
	for (string row; getline(input, row, row_delim);){
		++index;
		if (index % 10000000 == 0)
			cout << "half-read: " << index << endl;
		// reset vector
		curr.clear();

		istringstream ss(row);
		for(string word; getline(ss, word, field_delim);)
			curr.push_back(word);

		if (row[0]=='@' && row[1]=='S' && row[2]=='Q'){
			// Get the header information
			vReferenceGenome[iChr].samHeader = curr[1].substr(3);
			++iChr;
		} else {
			// We store the read name as the key and the position, length, and sequence as the elements in a struct
			sReadName = curr[0];
			frag.sReadName = curr[0];
			frag.iParentStart = atoi(curr[3].c_str());
			frag.sParentRead = curr[9];
			frag.iParentEnd = frag.iParentStart + frag.sParentRead.length() - 1;
			frag.iFlag = atoi(curr[1].c_str());

			if (sReadName[sReadName.length()-1] == '1'){ // the first half is anchored
				frag.bAnchorLeft = true;
			} else if (sReadName[sReadName.length()-1] == '2'){ // the second half is anchored, first half is unaligned
				frag.bAnchorLeft = false;
			} else {
				fLogFileOut << "\n ** ERROR ** There was an error when searching for the other half of an anchored read" << endl;
				continue;
			}

			for (unsigned int i = 0; i < vReferenceGenome.size(); ++i){
				//cout << curr[2] << " ?= " << vReferenceGenome[i].fastaHeader << endl;
				if (curr[2].find(vReferenceGenome[i].fastaHeader) != string::npos)
					frag.iChromosome = i;
				else if (vReferenceGenome[i].fastaHeader.find(curr[2]) != string::npos) // TODO find better way of searching for fasta header
					frag.iChromosome = i;
			}

			if (frag.iParentStart < 0 || frag.iFlag == 4){
				frag.bBadRead = true;
				++iBadReads;
			} else
				frag.bBadRead = false;

			// limit to single chromosome
			if (chromosome != (frag.iChromosome+1) && chromosome > 0){
				++iExludedReads;
				continue;
			}

			vCandidateReads.push_back(frag);
		}
	}
	input.close();

	fLogFileOut << "Excluded non-chromosome specific reads: " << iExludedReads << endl;
	fLogFileOut << "Bad reads: " << iBadReads << endl;

	cout << "Number candidate reads: " << vCandidateReads.size() << endl;
	fLogFileOut << "Number candidate reads: " << vCandidateReads.size() << endl;

	return 0;
}
// samtools view -S unaligned_1.sam | awk '{OFS="\t"; print ">"$1"-1\n"substr($10,1,length($10)/2)"\n>"$1"-2\n"substr($10,length($10)/2+1,length($10))}' - > filename.fasta

