/*
 * birAligner.cpp
 *
 *  Created on: Dec 19, 2012
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

#define STATE_1 0
#define RIGHT_ANCHORED 1
#define LEFT_ANCHORED 2
#define UNKNOWN 3

// birAligner
int startBirFinder(){
	cout << "\nStarting birFinder ..." << endl;
	fLogFileOut << "\nStarting birFinder ..." << endl;

	string sReference;
	string sBirSearchRegion;
	t_alignment_struct tAligned;
	int iFudgeFactor = confDB.getKey("tolerance").intVal;
	int iReadStart = 0; // the starting position of the read
	int iLastPos = 0; // a multi-use variable to keep the location of the aligned reads for RIGHT and LEFT anchored
	int iStartPercentage = 0;
	int iEndPercentage = 0;
	int iMinAlignedLength = confDB.getKey("minAlignedLength").intVal; // 12
	int iChr = 0;
	int size = vCandidateRegions.size();
	int tenPercent = size / 10;

	// statistics variables
	int skipped_short = 0;
	int skipped_pos = 0;
	int skipped_state1 = 0;
	int pos1 = 0;
	int left = 0;
	int right = 0;
	int unknown = 0;
	int iState = 0; // the state

	ofstream output;
	output.open((sProjectDirectory + "consensus_reads.txt").c_str());
	output << "i, chr, start, size, sequence" << endl;
	for (int i = 0; i < size; ++i)
		output << i << ", " << vCandidateRegions[i].iChromosome << ", " << vCandidateRegions[i].iParentStart << ", " << vCandidateRegions[i].sParentRead.length() << ", " << vCandidateRegions[i].sParentRead << endl;
	output.close();

	// Finally, let's repeat what we did at the start to find the BIR loc only this time reinforce it by searching for the template
	for (int i = 0; i < size; ++i){
		if (i % tenPercent == 0)
			cout << "candidate: " << (i+1) << " of " << size << " (" << ((i * 100) / size) << "%)" << endl;
		iReadStart = vCandidateRegions[i].iParentStart;
		iChr = vCandidateRegions[i].iChromosome;
		//cout << "\ni: " << i << "/" << size << ", start: " << iReadStart << ", chr: " << iChr << endl;

		// Lets' take the parent read with a known position and find the corresponding reference sequence
		sBirSearchRegion = vCandidateRegions[i].sParentRead; // the actual parent read sequence
		sReference = vReferenceGenome[iChr].sequence.substr(iReadStart-1, sBirSearchRegion.length());

		// transform sRefernce to upper-case
		transform(sReference.begin(), sReference.end(),sReference.begin(), ::toupper);

		// Perform an alignment
		tAligned = getLocalAlignment(sReference, sBirSearchRegion, 2.0, 5.0);

		//printAlignment(cout, sBirSearchRegion, sReference, iReadStart, tAligned);

		// get the aligned regions percentages
		iStartPercentage = tAligned.iStartPosJ * 100 / sBirSearchRegion.length();
		iEndPercentage = tAligned.iEndPosJ * 100 / sBirSearchRegion.length();

		// if the aligned region starts at the beginning percentage of the region and is a long distance away from the other end, it is RIGHT_ANCHORED
		if (iEndPercentage >= (100 - iFudgeFactor) && iStartPercentage > iFudgeFactor)
			iState = RIGHT_ANCHORED;
		// if the aligned region starts at the end of the region and is a long distance away from the other end, it is LEFT_ANCHORED
		else if (iEndPercentage < (100 - iFudgeFactor) && iStartPercentage <= iFudgeFactor)
			iState = LEFT_ANCHORED;
		// otherwise it is an ideal STATE_1 candidate and requires a FSM to get the coordinates
		else if (iStartPercentage <= iFudgeFactor && (int) iEndPercentage >= (100 - iFudgeFactor))
			iState = STATE_1;
		else
			iState = UNKNOWN;

		switch (iState) {
		case RIGHT_ANCHORED:
			//cout << "\nRight anchored..." << endl;
			iLastPos = tAligned.iStartPosJ;
			sBirSearchRegion = sBirSearchRegion.substr(0, tAligned.iStartPosJ);
			sReference = sReference.substr(0, tAligned.iStartPosI);
			tAligned = tAligned = getLocalAlignment(sReference, sBirSearchRegion, 2.0, 5.0);

			// printAlignment(sBirSearchRegion, sPossibleReference, iReadStart, tAligned);

			if (tAligned.iEndPosJ - tAligned.iStartPosJ < iMinAlignedLength){
				//cout << "\nskipped...." << endl;
				++skipped_short;
				vCandidateRegions[i].bBirCandidateFound = false;
			}
			else if (sBirSearchRegion.length() - tAligned.iEndPosJ > iMinAlignedLength && tAligned.iStartPosJ < iFudgeFactor){
				vCandidateRegions[i].sBir = sBirSearchRegion.substr(tAligned.iEndPosJ + 1, (iLastPos - tAligned.iEndPosJ + 1));
				vCandidateRegions[i].iBirStart = iReadStart + tAligned.iEndPosJ + 1;
				vCandidateRegions[i].iBirEnd = iReadStart + iLastPos - 1;
				vCandidateRegions[i].iBirLength = vCandidateRegions[i].iBirEnd - vCandidateRegions[i].iBirStart + 1;

				/*cout << "\nBIR: " << vCandidateRegions[i].sBir << endl;
				cout << tAligned.iEndPosJ + 1 << " - " << iLastPos - 1 << endl;
				cout << vCandidateRegions[i].iBirStart << ", " << vCandidateRegions[i].iBirEnd << ", " << vCandidateRegions[i].iBirLength << endl;
				cout << "\nPossible2..." << endl;*/

				++right;
				vCandidateRegions[i].bBirCandidateFound = true;
			}
			else{
				++skipped_pos;
				vCandidateRegions[i].bBirCandidateFound = false;
			}
			break;

		case LEFT_ANCHORED:
			//cout << "\nLeft anchored..." << endl;
			iLastPos = tAligned.iEndPosJ;
			sBirSearchRegion = sBirSearchRegion.substr(tAligned.iEndPosJ + 1);
			sReference = sReference.substr(tAligned.iEndPosI + 1);
			tAligned = tAligned = getLocalAlignment(sReference, sBirSearchRegion, 2.0, 5.0);

			// printAlignment(sBirSearchRegion, sPossibleReference, iReadStart, tAligned);

			if (tAligned.iEndPosJ - tAligned.iStartPosJ < iMinAlignedLength){
				//cout << "\nskipped...." << endl;
				++skipped_short;
				vCandidateRegions[i].bBirCandidateFound = false;
			}
			else if (tAligned.iStartPosJ > iMinAlignedLength && (int) sBirSearchRegion.length() - tAligned.iEndPosJ < iFudgeFactor){
				vCandidateRegions[i].sBir = sBirSearchRegion.substr(0, tAligned.iStartPosJ);
				vCandidateRegions[i].iBirStart = iReadStart + iLastPos + 1;
				vCandidateRegions[i].iBirEnd = iReadStart + iLastPos + tAligned.iStartPosJ;
				vCandidateRegions[i].iBirLength = vCandidateRegions[i].iBirEnd - vCandidateRegions[i].iBirStart + 1;

				/*cout << "\nBIR: " << vCandidateRegions[i].sBir << endl;
				cout << iLastPos + 1 << " - " << (iLastPos + tAligned.iStartPosJ + 1) << endl;
				cout << vCandidateRegions[i].iBirStart << ", " << vCandidateRegions[i].iBirEnd << ", " << vCandidateRegions[i].iBirLength << endl;
				cout << "\nPossible2..." << endl;*/

				++left;
				vCandidateRegions[i].bBirCandidateFound = true;
			}
			else{
				++skipped_pos;
				vCandidateRegions[i].bBirCandidateFound = false;
			}
			break;

		case STATE_1:
			//cout << "\nPossible1..." << endl;
			getBirLoc(tAligned, iReadStart, i);
			if (vCandidateRegions[i].bBirCandidateFound)
				++pos1;
			else
				++skipped_state1;
			break;

		default: // THERE IS AN ERROR WITH THE ALIGNMENT
			++unknown;
			//if (bPrintError)
			//printError(sBirSearchRegion, sReference, tAligned);
			break;
		}
		//cout << "\n-------------------------------------------------------------------------------------------------------\n" << endl;
	}
	fLogFileOut << "\nBIR Alignment Stats: " << endl;
	fLogFileOut << "  candidate regions: " << vCandidateRegions.size() << endl;
	fLogFileOut << "  skipped (short): " << skipped_short << endl;
	fLogFileOut << "  skipped (pos): " << skipped_pos << endl;
	fLogFileOut << "  skipped (state1): " << skipped_state1 << endl;
	fLogFileOut << "  state_1: " << pos1 << endl;
	fLogFileOut << "  left: " << left << endl;
	fLogFileOut << "  right: " << right << endl;
	fLogFileOut << "  unknown: " << unknown << endl;

	fLogFileOut << "\nBIR alignment time = " << (int)time(NULL)-time0 << endl;

	return 0;
}



#define STATE_LEFT 0
#define STATE_BIR 1
#define STATE_RIGHT 2

/*
 * ************************************************************
 * getBirLoc
 *
 * Uses a finite state machine to get the start AND end of the BIR location.
 *
 * NOTE: This assumes that the template and read start at the same position
 */
void getBirLoc(t_alignment_struct tAligned, int st, int curr){
	t_consolidated locs;

	string reference = tAligned.sAlignedRegionI;
	string birRegion = tAligned.sAlignedRegionJ;
	int iMissCount = confDB.getKey("missCount").intVal; // 3
	int iHitCount = confDB.getKey("hitCount").intVal; // 4
	int iRefPos = -1;
	int iLastPos = 0;
	int length = birRegion.length();
	int iMisses = 0;
	int iHits = 0;
	int start = 0;
	int end = 0;
	bool isGood = false;

	int nState = STATE_LEFT;

	for (int i = 0; i < length; ++i){
		++iRefPos;

		switch (nState) {
		case STATE_LEFT:
			if (reference[iRefPos] == birRegion[i] || reference[iRefPos == 'N' || birRegion[i] == 'N']){
				++iHits;
				iMisses = 0;
			}
			else{
				if (iMisses == 0)
					iLastPos = i;
				++iMisses;
			}
			if (iMisses == iMissCount){
				nState = STATE_BIR;
				start = iLastPos;
			}
			break;
		case STATE_BIR:
			if (reference[iRefPos] != birRegion[i]){
				++iMisses;
				iHits = 0;
			}
			else{
				if (iHits == 0)
					iLastPos = i;
				++iHits;
			}
			if (iHits == iHitCount){
				nState = STATE_RIGHT;
				end = iLastPos - 1;
				isGood = true;
			}
			break;
		case STATE_RIGHT:
			// do nothing
			break;
		}
	}

	if (isGood){
		vCandidateRegions[curr].bBirCandidateFound = true;
		vCandidateRegions[curr].sBir = birRegion.substr(start, end - start + 1);
		vCandidateRegions[curr].iBirStart = st + start;
		vCandidateRegions[curr].iBirEnd = st + end;
		vCandidateRegions[curr].iBirLength = end - start + 1;
	} else {
		vCandidateRegions[curr].bBirCandidateFound = false;
	}
}
