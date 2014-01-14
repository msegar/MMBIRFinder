/*
 * birConsolidate.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: msegar
 */


#include "defs.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/wait.h>
#include <stdlib.h>
#include <algorithm>
#include <limits.h>
#include <ctime>
#include <cmath>

using namespace std;

// birConsolidate
int startConsolidate(){
	cout << "\nstartConsolidate start..." << endl;
	fLogFileOut << "\nstartConsolidate start..." << endl;
	ofstream output;

	// only perform clustering if user option is specified
	if (confDB.getKey("performClustering").boolVal == true){
		cout << "start candidate read sorting..." << endl;
		fLogFileOut << "start candidate read sorting..." << endl;
		// sort the parent starting locations for each chromosome from smallest to largest
		sort(vCandidateReads.begin(), vCandidateReads.end(), compareStart);

		fLogFileOut << "Before consolidation: " << vCandidateReads.size() << endl;
		cout << "Printing unconsolidated_bir_locations.txt..." << endl;
		output.open((sProjectDirectory + "unconsolidated_bir_locations.txt").c_str());
		for (unsigned int i = 0; i < vCandidateReads.size(); ++i){
			if (vCandidateReads[i].bBadRead)
				continue;
			output << vCandidateReads[i].iChromosome << ", " << vCandidateReads[i].iParentStart << "," << vCandidateReads[i].iParentEnd << endl;
		}
		output.close();

		// Finally, let's consolidated the possible BIR locations to speed up the alignment steps. This is based on overlapping regions
		consolidateLocations();

		fLogFileOut << "\nAfter consolidation: " << vConsolidated.size() << endl;
	}

	// Let's destruct all the elements in the database to save space
	vCandidateReads.clear();
	vConsolidated.clear();

	/*
	 * *************************************************************************************************************************
	 * ******************			THIS IS WHERE THE PROGRAM SPLITS IN TWO. PREVIOUS TO THIS POINT,          ******************
	 * ******************           THE PROGRAM WORKS AS SPECIFIED. THE PROGRAM NOW REQUIRES ANOTHER          ******************
	 * ******************           SCRIPT THAT UTILIZES MYSQL TO SEARCH FOR THE HALF-READS QUICKLY.          ******************
	 * ******************          THE PROGRAM NOW WILL IMPORT THE RESULTS FROM THE SCRIPT AND CONTINUE.      ******************
	 * *************************************************************************************************************************
	 */

	// Find the other size of the half reads
	createParentReads();

	// Make all the reads in a cluster the same length
	getConsensus();

	// Create a consensus read that takes into consideration all the base calls at a specific location
	consolidateBaseCalls();

	// Let's destruct all the elements in the database to save space
	vConsolidated.clear();

	fLogFileOut << "\nConsolidation time = " << (int)time(NULL)-time0 << endl;

	output.open((sProjectDirectory+ "consolidated_bir_locations.txt").c_str());
	for (unsigned int i = 0; i < vCandidateRegions.size(); ++i)
		output << vCandidateRegions[i].iChromosome << ", " << vCandidateRegions[i].iParentStart << "," << vCandidateRegions[i].iParentStart + vCandidateRegions[i].sParentRead.length() << endl;
	output.close();

	return 0;
}


/*
 * ************************************************************
 * consolidateLocations
 *
 * This function takes all the candidate regions and groups them together by the starting location of the unaligned read. If
 * there is any overlap between the reads, then they will be grouped together.
 *
 * A config.txt parameter (minConsolidate) determines how many reads need to be used to consolidate.
 */
void consolidateLocations(){
	cout << "\nConsolidating read locations start..." << endl;
	fLogFileOut << "\nConsolidating read locations start..." << endl;

	// TODO Find better way to consolidate, maybe based on stdev and discard ones that aren't (ie through out outliers)
	int iLastEndLoc = 0;
	unsigned int iMinConsolidate = confDB.getKey("minConsolidate").intVal; // the minimum number of reads to consolidate, otherwise don't add to vector
	vector<t_consolidated> curr;
	int iNumReads = 0;
	int iSkippedCluster = 0;
	bool first = true;
	int iChr = (confDB.getKey("chromosome").intVal == 0 ? 0 : confDB.getKey("chromosome").intVal - 1); // get the user specified chromoome
	bool bOnlyOneChr = (confDB.getKey("chromosome").intVal != 0 ? true :false);
	int sizeCan = vCandidateReads.size();
	int fivePercent = sizeCan / 20;
	string sHalfReadClusters = confDB.getKey("clusterFile").stringVal;

	cout << "Chromsome: " << iChr << endl;
	for (int i = 0; i < sizeCan; ++i){
		if (i % fivePercent == 0)
			cout << "\tcandidateRead: " << i << " of " << sizeCan << " on chromosome " << iChr << " (" << ((i * 100) / sizeCan) << "%)" << endl;
		if (vCandidateReads[i].iChromosome != iChr && vCandidateReads[i].iChromosome != iChr+1)
			continue;
		else if (vCandidateReads[i].iChromosome != iChr && vCandidateReads[i].iChromosome == iChr+1){
			if (bOnlyOneChr) // if the user only wants to cluster one chromosome, we break out here
				break;
			++iChr;
			cout << "\nChromsome: " << iChr << endl;
			if (curr.size() >= iMinConsolidate && curr.size() <= 200){ // TODO make 200 a config variable
				vConsolidated.push_back(curr);
				iNumReads += curr.size();
			} else {
				++iSkippedCluster;
			}
			curr.clear();

			// now instead of copying the code below, we set first to true (since it's the start of a new chromosome) and repeat the loop
			first = true;
			--i;
			continue;
		}

		if (first){
			curr.push_back(vCandidateReads[i]);
			first = false;
			iLastEndLoc = vCandidateReads[i].iParentStart + vCandidateReads[i].sParentRead.length() - 1;
			continue;
		}

		if (vCandidateReads[i].iParentStart > iLastEndLoc){
			if (curr.size() >= iMinConsolidate && curr.size() <= 200){ // TODO make 200 a config variable
				vConsolidated.push_back(curr);
				iNumReads += curr.size();
			} else {
				++iSkippedCluster;
			}
			curr.clear();
			curr.push_back(vCandidateReads[i]);
		} else {
			curr.push_back(vCandidateReads[i]);
		}
		iLastEndLoc = vCandidateReads[i].iParentStart + vCandidateReads[i].sParentRead.length() - 1;

	}

	fLogFileOut << "Number clusters skipped: " << iSkippedCluster << endl;
	cout << "\nNumber reads after consolidation: " << iNumReads << endl;
	fLogFileOut << "Number reads after consolidation: " << iNumReads << endl;
	cout << "Number of clusters: " << vConsolidated.size() << endl;
	fLogFileOut << "Number of clusters: " << vConsolidated.size() << endl;

	// Now the program will print out the consolidate read locations
	// A MySQL and Perl script is required to parse the data to further the program
	// The program will now exit
	fLogFileOut << "\nPrinting half_read_clusters.txt" << endl;
	ofstream output;
	output.open((sProjectDirectory + sHalfReadClusters.c_str()).c_str());
	int size2 = 0;
	for (int i = 0; i < vConsolidated.size(); ++i){
		size2 = vConsolidated[i].size();
		for (int j = 0; j < size2; ++j){
			output << vConsolidated[i].at(j).sReadName << "\t" << vConsolidated[i].at(j).sParentRead << "\t" << vConsolidated[i].at(j).iParentStart << "\t" << vConsolidated[i].at(j).iChromosome << "\t"<< i << endl;
		}
	}
	output.close();

	// print out the half-read clustered locations
	/*fLogFileOut << "Printing half_read_clustered_locations.txt" << endl;
	output.open("half_read_clustered_locations.txt");
	output << "i\tchr\tsize\tmin\tmax\tlength" << endl;
	for (int i = 0; i < iNumReads; ++i){
		int min = 0;
		int max = 0;
		size2 = vConsolidated[i].size();
		for (int j = 0; j < size2; ++j){
			if (j == 0)
				min = vConsolidated[i].at(j).iParentStart;
			else if (vConsolidated[i].at(j).iParentStart < min)
				min = vConsolidated[i].at(j).iParentStart;
			if (max < vConsolidated[i].at(j).iParentEnd)
				max = vConsolidated[i].at(j).iParentEnd;
		}
		output << i << "\t" << vConsolidated[i].at(0).iChromosome << "\t" << size2 << "\t" << min << "\t" << max << "\t" << (max - min + 1) << endl;
	}
	output.close();*/

	fLogFileOut << "\nCluster time = " << (int)time(NULL)-time0 << endl;

	// clean up memory before exiting
	cout << "\nCleaning vCandidateReads..." << endl;
	vCandidateReads.clear();
	cout << "Cleaning vConsolidated..." << endl;
	vConsolidated.clear();

	exit(0);
}

int createParentReads(){
	cout << "\nImporting clusters start... " << endl;
	fLogFileOut << "\nImporting clusters start... " << endl;

	ifstream input;
	input.open((sProjectDirectory + confDB.getKey("mysqlFile").stringVal).c_str());
	// check if file is open
	if (input.is_open() == false){
		cout << "\nThe mysql_results file could not be found" << endl;
		exit(1);
	}

	t_consolidated cons;
	vector<string> curr;
	vector<t_consolidated> vCluster;
	int iCluster = 0;
	int iCurrCluster = 0;
	int iClustersImported = 0;
	int iReadsImported = 0;
	int iChromosome = confDB.getKey("chromosome").intVal;
	bool first = true;
	unsigned int iMinConsolidate = confDB.getKey("minConsolidate").intVal; // the minimum number of reads to consolidate, otherwise don't add to vector
	int iSkipped = 0;

	for (string row; getline(input, row, '\n');){
		istringstream ss(row);
		curr.clear();
		for(string word; getline(ss, word, '\t');)
			curr.push_back(word);

		cons.sParentRead = curr[0];
		cons.sReadName = curr[1].substr(0, curr[1].length()-2);
		cons.iParentStart = atoi(curr[3].c_str());
		cons.iChromosome = atoi(curr[4].c_str());
		iCurrCluster = atoi(curr[5].c_str());

		if (iChromosome != (cons.iChromosome+1) && iChromosome != 0)
			continue;

		// catch reads that happen to have a starting location less than 0
		if (cons.iParentStart <= 0){
			++iSkipped;
			continue;
		}

		if (first){
			iCluster = iCurrCluster;
			first = false;
		}

		if (iCurrCluster == iCluster)
			vCluster.push_back(cons);
		else {
			//if (vCluster.size() >= iMinConsolidate){
				vConsolidated.push_back(vCluster);
				iReadsImported += vCluster.size();
				++iClustersImported;
			//} else
			//	++iSkipped;
			vCluster.clear();
			vCluster.push_back(cons);
			iCluster = iCurrCluster;
		}
	}
	if (vCluster.size() > 0){
		vConsolidated.push_back(vCluster);
		iReadsImported += vCluster.size();
		vCluster.clear();
		++iClustersImported;
	}

	input.close();

	cout << "Reads skipped (start < 0): " << iSkipped << endl;
	fLogFileOut << "\nReads skipped (start < 0): " << iSkipped << endl;
	cout << "Clusters imported: " << iClustersImported << endl;
	fLogFileOut << "Clusters imported: " << iClustersImported << endl;
	cout << "Reads imported: " << iReadsImported << endl;
	fLogFileOut << "Reads imported: " << iReadsImported << endl;
	fLogFileOut << "Cluster import time = " << (int)time(NULL)-time0 << endl;

	ofstream output;
	output.open((sProjectDirectory+ "clustered_locations.txt").c_str());

	output << "i\tj\tstart\tchr\tread" << endl;
	int size = vConsolidated.size();
	int size2 = 0;
	for (int i = 0; i < size; ++i) {
		size2 = vConsolidated[i].size();
		for (int j = 0; j < size2; ++j)
			output << i << ", " << j << ", " << vConsolidated[i].at(j).iParentStart << ", " << vConsolidated[i].at(j).iChromosome << ", " << vConsolidated[i].at(j).sParentRead << endl;
	}

	output.close();

	return 0;
}

/*
 * ************************************************************
 * getConsensus
 *
 * This function takes the grouped (consolidated) reads and makes them all the same length by adding '-' to the beginning and end.
 *
 * For example, if the reads are:
 *
 *   AACGTCGA								  	--AACGTCGA---
 * CGAACGTC			will be transformed into    CGAACGTC-----
 *     CGTCGACGT								----CGTCGACGT
 *
 * This helps for the the next step which is calculate the base calls at each location.
 */
void getConsensus(){
	cout << "\nMaking reads same length..." << endl;
	fLogFileOut << "\nMaking reads same length..." << endl;

	int iLowestStart;
	int iLongest;
	int iCurrStart;
	int diff;
	int size = 0;
	string currString;
	int sizeCons = vConsolidated.size();
	int fivePercent = sizeCons / 20;

	// make all the reads in each consolidation have the same starting point
	for (unsigned int i = 0; i < sizeCons; ++i){
		if (i % fivePercent == 0)
			cout << "consolidated: " << (i+1) << " of " << sizeCons << " (" << ((i * 100) / sizeCons) << "%)" << endl;
		size = vConsolidated[i].size();
		iLowestStart = vConsolidated[i].at(0).iParentStart;
		for (int j = 1; j < size; ++j){
			iCurrStart = vConsolidated[i].at(j).iParentStart;
			if (iLowestStart > iCurrStart)
				iLowestStart = iCurrStart;
		}

		// Now we have the lowest starting point, let's add '-' to the beginning so they all match up
		for (int j = 0; j < size; ++j){
			iCurrStart =  vConsolidated[i].at(j).iParentStart;
			currString = vConsolidated[i].at(j).sParentRead;
			diff = iCurrStart - iLowestStart;
			for (int k = 0; k < diff; ++k){
				currString = "-" + currString;
				++iCurrStart;
			}
			vConsolidated[i].at(j).sParentRead = currString;
			currString = "";
		}

		// Let's repeat only adding '-' to the end
		iLongest = vConsolidated[i].at(0).sParentRead.length();
		for (int j = 1; j < size; ++j){
			iCurrStart = vConsolidated[i].at(j).sParentRead.length();
			if (iLongest < iCurrStart)
				iLongest = iCurrStart;
		}

		for (int j = 0; j < size; ++j){
			iCurrStart = vConsolidated[i].at(j).sParentRead.length();
			currString = vConsolidated[i].at(j).sParentRead;
			for (int k = iCurrStart; k < iLongest; ++k){
				currString = currString + "-";
				++iCurrStart;
			}
			vConsolidated[i].at(j).sParentRead = currString;
			currString = "";
		}
	}
	cout << "exiting..." << endl;
	ofstream output;
	output.open((sProjectDirectory + "clustered_locations2.txt").c_str());

	output << "i\tj\tstart\tchr\tread" << endl;
	int size2 = vConsolidated.size();
	int size3 = 0;
	for (int i = 0; i < size2; ++i) {
		size3 = vConsolidated[i].size();
		for (int j = 0; j < size3; ++j)
			output << i << ", " << j << ", " << vConsolidated[i].at(j).iParentStart << ", " << vConsolidated[i].at(j).iChromosome << ", " << vConsolidated[i].at(j).sParentRead << endl;
		output << endl;
	}

	output.close();
}


/*
 * ************************************************************
 * consolidateBaseCalls
 *
 * This is basically a modified de novo assembly method that takes the base call frequency at each location to create a new read.
 */
void consolidateBaseCalls(){
	cout << "\nCalculating consensus base calls..." << endl;
	fLogFileOut << "\nCalculating consensus base calls..." << endl;

	int iConsolidatedSize = vConsolidated.size();
	int iCurrSize = 0;
	int arr[85];
	int numTot = 0;
	int iReadLength = 0;
	string blank = "";
	string sConsensus = ""; // the consensus string
	int iMinStartPos = 0;
	string sCurrChar;
	int iHighestProb;
	int fivePercent = iConsolidatedSize / 20;
	t_consolidated consensusRead;

	/*
	 * Memory is starting to become an issue. To compensate for this, we will start deleting the vConsolidated cluster after the vCandidateRegions
	 * struct is added. However, erasing in a vector is inefficient since all the memory is copied. There is one way to compensate for this, start from
	 * the back. Since it doesn't matter which cluster we start on, and reallocation won't be affected erasing from the back, let's reverse traverse the
	 * vector.
	 *
	 * However, we need to start rearranging previous methods to start from the back as well. Start with the sorting algorithm and make it so the largest
	 * chromosome is at the front of the vector.
	 */
	// TODO erase starting from the back.

	for (int i = 0; i < iConsolidatedSize; ++i){
		if (i % fivePercent == 0)
			cout << "consolidated: " << i << " of " << iConsolidatedSize << " (" << ((i * 100) / iConsolidatedSize) << "%)" << endl;

		// get the address to make lookups faster
		vector<t_consolidated>& curr = vConsolidated[i];
		iReadLength = curr[0].sParentRead.length();
		sConsensus = "";

		iCurrSize = curr.size();

		iMinStartPos = curr[0].iParentStart;
		for (int j = 1; j < iCurrSize; ++j){
			if (curr[j].iParentStart < iMinStartPos)
				iMinStartPos = curr[j].iParentStart;
		}

		for (int j = 0; j < iReadLength; ++j){
			// initialize array and variables to 0s
			for (int k = 0; k < 85; ++k)
				arr[k] = 0;
			numTot = 0;

			for (int k = 0; k < iCurrSize; ++k){
				++arr[(int) curr[k].sParentRead.at(j)];
				if (curr[k].sParentRead.at(j) != '-')
					++numTot;
			}
			//curr.clear();
			iHighestProb = 0;

			// -s
			if (numTot == 0)
				continue;
			// As
			if ((arr[65] * 10) / numTot > iHighestProb){
				iHighestProb = (arr[65] * 10) / numTot;
				sCurrChar = "A";
			}
			// Cs
			if ((arr[67] * 10) / numTot > iHighestProb){
				iHighestProb = (arr[67] * 10) / numTot;
				sCurrChar = "C";
			}
			// Gs
			if ((arr[71] * 10) / numTot > iHighestProb){
				iHighestProb = (arr[71] * 10) / numTot;
				sCurrChar = "G";
			}
			// Ts
			if ((arr[84] * 10) / numTot > iHighestProb){
				iHighestProb = (arr[84] * 10) / numTot;
				sCurrChar = "T";
			}

			sConsensus += sCurrChar;
			sCurrChar = "-";
		}

		curr.clear();
		if (sConsensus.length() > 900) // TODO THIS NEEDS TO BE FIXED, BUT FOR NOW WE'LL SEE
			continue;
		consensusRead.sParentRead = sConsensus;
		consensusRead.iParentStart = iMinStartPos;
		consensusRead.iChromosome = curr[0].iChromosome;
		consensusRead.bBirCandidateFound = false;
		vCandidateRegions.push_back(consensusRead);
	}

	// sort the consensus reads
	sort(vCandidateRegions.begin(), vCandidateRegions.end(), compareStart);
}


bool compareStart(const t_consolidated &a, const t_consolidated &b){
	if (a.iChromosome < b.iChromosome) return true;
	if (a.iChromosome == b.iChromosome) return a.iParentStart < b.iParentStart;
	return false;
}
