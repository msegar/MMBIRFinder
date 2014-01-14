/*
 * templateFinder.cpp
 *
 *  Created on: Feb 4, 2013
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

// templateFinder
int startTemplateFinder(){
	cout << "\nStarting template finder ..." << endl;
	fLogFileOut << "\nStarting template finder ..." << endl;

	ofstream output;
	output.open((sProjectDirectory+ "final_bir_locs.txt").c_str());
	string sBir;
	string sReference;
	string sBirReversed;
	int iBirStart = 0;
	int iTemplateSearchDistance = confDB.getKey("searchLength").intVal;
	int iMinBirLength = confDB.getKey("minBirLength").intVal;
	//int iMinBirLength = 13;
	int iChr = 0;
	t_alignment_struct tAligned;
	int iSize = vCandidateRegions.size();
	int tenPercent = iSize / 10;

	// statistics variables
	int skipped_short = 0;
	int skipped_stdev = 0;

	// for the candidate regions with the bBirCandidateFound flag set to TRUE...
	for (unsigned int i = 0; i < iSize; ++i){
		if (i % tenPercent == 0)
			cout << "template: " << (i+1) << " of " << iSize << " (" << ((i * 100) / iSize) << "%)" << endl;
		if (!vCandidateRegions[i].bBirCandidateFound)
					continue;

		iChr = vCandidateRegions[i].iChromosome;
		iBirStart = vCandidateRegions[i].iBirStart;
		sBir = vCandidateRegions[i].sBir;
		sBirReversed = getReverseComplement(sBir);
		sReference = vReferenceGenome[iChr].sequence.substr(iBirStart - iTemplateSearchDistance - 1, iTemplateSearchDistance);
		transform(sReference.begin(), sReference.end(),sReference.begin(), ::toupper);

		/*cout << "iChr: " << iChr << endl;
		cout << "iBirStart: " << iBirStart << endl;
		cout << "sBir: " << sBir << endl;
		cout << "sBirReversed: " << sBirReversed << endl;
		cout << "sReference: " << sReference << endl;
		*/
		tAligned = getLocalAlignment(sReference, sBirReversed, 2.0, 5.0);
		//printAlignment(cout, sBirReversed, sReference, iBirStart, tAligned);

		vCandidateRegions[i].sTemplate = tAligned.sAlignedRegionJ + sBirReversed[sBirReversed.length()-1];

		if ((int) vCandidateRegions[i].sTemplate.length() < iMinBirLength){
			++skipped_short;
			continue;
		}
		if ((sBir.length() * 0.8) > (int) vCandidateRegions[i].sTemplate.length()){ // TODO make 0.8 a config variable (80%)
			++skipped_stdev;
			continue;
		}


		vCandidateRegions[i].iTemplateStart = iBirStart - iTemplateSearchDistance - 1 + tAligned.iStartPosJ;
		vCandidateRegions[i].iTemplateEnd = iBirStart - iTemplateSearchDistance - 1 + tAligned.iEndPosJ;

		output << "iChr: " << iChr << endl;
		output << "iBirStart: " << iBirStart << endl;
		output << "INS: " << iBirStart - vCandidateRegions[i].iTemplateEnd + 1 << endl;
		output << "sBir: " << sBir << endl;
		output << "sBirReversed: " << sBirReversed << endl;
		output << "ref: " << sReference << endl;
		output << "bir: " << vCandidateRegions[i].sParentRead << endl;
		output << "TEM: " << vCandidateRegions[i].sTemplate << endl;

		printAlignment(output, sBirReversed, sReference, iBirStart, tAligned);
		output << endl;

		vFinalBirLocs.push_back(vCandidateRegions[i]);
	}
	output.close();

	fLogFileOut << "\nFinal stats: " << endl;
	fLogFileOut << "  skipped (short): " << skipped_short << endl;
	fLogFileOut << "  skipped (stdev): " << skipped_stdev << endl;
	fLogFileOut << "  found: " << vFinalBirLocs.size() << endl;

	// print out the final BIR incidents found and the templates
	//printFinal();

	// print final time
	fLogFileOut << "\nFinal time: " << (int)time(NULL)-time0 << endl;

	// clear memory
	vCandidateRegions.clear();
	vFinalBirLocs.clear();

	return 0;
}


string getReverseComplement(string &str){
	string rev;

	for (unsigned int i = 0; i < str.length(); ++i){
		if (str[i] == 'A')
			rev = 'T' + rev;
		else if (str[i] == 'C')
			rev = 'G' + rev;
		else if (str[i] == 'G')
			rev = 'C' + rev;
		else if (str[i] == 'T')
			rev = 'A' + rev;
		else
			rev = 'N' + rev;
	}

	return rev;
}

