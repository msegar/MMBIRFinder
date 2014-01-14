/*
 * alignment.cpp
 *
 *  Created on: Dec 19, 2012
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

//double mu = confDB.getKey("mu").doubleVal;
//double delta = confDB.getKey("delta").doubleVal;
double mu = 0.0;
double delta = 0.0;
int arrIPath[1000][1000]; // the backtrack path for i
int arrJPath[1000][1000]; // the backtrack path for j
double matrix[1000][1000];


/*
 * ************************************************************
 * getLocalAlignment
 *
 * This takes in two strings. The first is the reference string, the second is the search string.
 * The function searches the reference for the search string and returns the best LOCAL alignment
 */
t_alignment_struct getLocalAlignment(string seq1, string seq2, double m, double d){
	//cout << seq1.length() << ", " << seq2.length() << endl;
	mu = m;
	delta = d;
	t_alignment_struct tReturn;
	tReturn.iEndPosI = tReturn.iEndPosJ = tReturn.iStartPosI = tReturn.iStartPosJ = 0;
	tReturn.sAlignedRegionI = tReturn.sAlignedRegionJ = "";
	string sBIR = "";
	string sTemplate = "";
	double temp[4];
	int iCase = 0;
	int iLengthS1 = seq1.length();
	int iLengthS2 = seq2.length();
	/*static int arrIPath[seq1.length() + 1][seq2.length() + 1]; // the backtrack path for i
	static int arrJPath[seq1.length() + 1][seq2.length() + 1]; // the backtrack path for j
	static double matrix[seq1.length() + 1][seq2.length() + 1];
	*/

	// initialize matrix to 0s
	for (int i = 0; i <= iLengthS1; ++i){
		for (int j = 0; j <= iLengthS2; ++j){
			matrix[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; ++i)
		temp[i] = 0;

	// Now we're ready for the algorithm!

	for (int i = 1; i < iLengthS1; ++i){
		for (int j = 1; j < iLengthS2; ++j){
			// store the 4 possible values for a local alignment
			temp[0] = matrix[i-1][j-1] + getSimilarityScore(seq1[i-1], seq2[j-1]); // diagonal
			temp[1] = matrix[i-1][j] - delta; // directly above
			temp[2] = matrix[i][j-1] - delta; // to the left
			temp[3] = 0; // 0
			iCase = getMaxArrayValue(temp, 4);
			matrix[i][j] = temp[iCase];
			switch(iCase){
			case 0: // a match or mismatch (diagonal)
				arrIPath[i][j] = i-1;
				arrJPath[i][j] = j-1;
				break;
			case 1: // a deletion in sequence 1
				arrIPath[i][j] = i-1;
				arrJPath[i][j] = j;
				break;
			case 2: // a deletion in sequence 2
				arrIPath[i][j] = i;
				arrJPath[i][j] = j-1;
				break;
			case 3: // (i,j) is the beginning of a sequence (or subsequence)
				arrIPath[i][j] = i;
				arrJPath[i][j] = j;
				break;
			}
		}
	}

	// Ok, now we have our matrix and our backtracing path
	// Let's scan through our matrix and find the maximum score
	double dMax = 0.0;
	int iMaxCoordI = 0; // i coordinate of maximum value
	int iMaxCoordJ = 0; // j coordinate of maximum value
	for (int i = 1; i < iLengthS1; ++i){
		for (int j = 1; j < iLengthS2; ++j){
			if (matrix[i][j] > dMax){
				dMax = matrix[i][j];
				iMaxCoordI = i;
				iMaxCoordJ = j;
			}
		}
	}

	// Now we have our maximum value and the respective coordinates
	// Let's backtrack from that maximum value until we get to 0 (as per a local alignment)
	int iCurrentI = iMaxCoordI;
	int iCurrentJ = iMaxCoordJ;
	int iNextI = arrIPath[iCurrentI][iCurrentJ]; // the next i value to take
	int iNextJ = arrJPath[iCurrentI][iCurrentJ];

	while (((iCurrentI != iNextI) || (iCurrentJ != iNextJ)) && (iNextI != 0) && (iNextJ != 0)){
		if (iNextI == iCurrentI) // A deletion in seq A (the template)
			sTemplate = '-' + sTemplate;
		else // A match or mismatch in seq A
			sTemplate = seq1[iCurrentI - 1] + sTemplate;

		if (iNextJ == iCurrentJ) // A deletion in seq B (the BRI)
			sBIR = '-' + sBIR;
		else // A match or mismatch in seq B
			sBIR = seq2[iCurrentJ - 1] + sBIR;

		iCurrentI = iNextI;
		iCurrentJ = iNextJ;
		iNextI = arrIPath[iCurrentI][iCurrentJ];
		iNextJ = arrJPath[iCurrentI][iCurrentJ];
	}

	sTemplate = seq1[iCurrentI - 1] + sTemplate;
	sBIR = seq2[iCurrentJ - 1] + sBIR;

	// We have our consensus motif!!

	tReturn.sAlignedRegionJ = sBIR;
	tReturn.sAlignedRegionI = sTemplate;
	tReturn.iStartPosJ = iCurrentJ - 1;
	tReturn.iStartPosI = iCurrentI - 1;
	tReturn.iEndPosI = iMaxCoordI - 1;
	tReturn.iEndPosJ = iMaxCoordJ - 1;

	return tReturn;
}

/*
 * ************************************************************
 * getGlobalAlignment
 *
 * This takes in two strings. The first is the reference string, the second is the search string.
 * The function searches the reference for the search string and returns the best GLOBAL alignment
 */
t_alignment_struct getGlobalAlignment(string &seq1, string &seq2, double m, double d){
	mu = m;
	delta = d;
	string sBIR;
	string sTemplate;
	double matrix[seq1.length()+1][seq2.length()+1];
	double temp[3]; // ******* CHANGE TO 4 FOR LOCAL ALIGNEMNT
	int iCase;
	int iLengthS1 = seq1.length();
	int iLengthS2 = seq2.length();

	t_alignment_struct tReturn;


	// initialize matrix to 0s
	for (int i = 0; i <= iLengthS1; ++i){
		for (int j = 0; j < iLengthS2; ++j){
			matrix[i][j] = 0;
		}
	}

	for (int i = 0; i <= iLengthS1; ++i)
		matrix[i][0] = delta*i;
	for (int j = 0; j < iLengthS2; ++j)
		matrix[0][j] = delta*j;

	// Now we're ready for the algorithm!

	for (int i = 1; i <= iLengthS1; ++i){
		for (int j = 1; j <= iLengthS2; ++j){
			// store the 3 possible values for a global alignment
			temp[0] = matrix[i-1][j-1] + getSimilarityScore(seq1[i-1], seq2[j-1]); // diagonal
			temp[1] = matrix[i-1][j] - delta; // directly above
			temp[2] = matrix[i][j-1] - delta; // to the left
			iCase = getMaxArrayValue(temp, 3);
			matrix[i][j] = temp[iCase];
		}
	}

	// Ok, now we have our matrix and our backtracing path
	// Let's scan through our matrix and find the maximum score

	//cout << "dMax: " << dMax << endl;

	// Now we have our maximum value and the respective coordinates
	// Let's backtrack from that maximum value until we get to 0 (as per a local alignment)
	int iCurrI = iLengthS1;
	int iCurrJ = iLengthS2;

	while (iCurrI > 0 || iCurrJ > 0){
		if (iCurrI > 0 && iCurrJ > 0 && matrix[iCurrI][iCurrJ] == (matrix[iCurrI-1][iCurrJ-1] + getSimilarityScore(seq1[iCurrI-1], seq2[iCurrJ-1]))){
			sTemplate = seq1[iCurrI-1] + sTemplate;
			sBIR = seq2[iCurrJ-1] + sBIR;
			--iCurrI;
			--iCurrJ;
		}
		else if (iCurrI > 0 && matrix[iCurrI][iCurrJ] == (matrix[iCurrI-1][iCurrJ] + delta)){
			sTemplate = seq1[iCurrI-1] + sTemplate;
			sBIR = "-" + sBIR;
			--iCurrI;
		}
		else {
			sTemplate = "-" + sTemplate;
			sBIR = seq2[iCurrJ-1] + sBIR;
			--iCurrJ;
		}
	}

	//sTemplate = seq1[iCurrentI - 1] + sTemplate;
	//sBRI = seq2[iCurrentJ - 1] + sBRI;

	// We have our consensus motif!!

	tReturn.sAlignedRegionJ = sBIR;
	tReturn.sAlignedRegionI = sTemplate;
	tReturn.iEndPosJ = 0;
	tReturn.iEndPosI = 0;
	tReturn.iStartPosI = 0;
	tReturn.iStartPosJ = 0;
	return tReturn;
}

/*
 * ************************************************************
 * getSimilarityScore
 *
 * This is a simple function to take in 2 characters and determine if they are equal.
 * If they are then a match of 1.0 is returned
 * Otherwise, a negative mu value is returned to signal a mismatch
 * ************************************************************
 */
double getSimilarityScore(char a, char b){
	double result;

	result = (a == b ? 1.0 : -mu);

	return result;

}

/*
 * ************************************************************
 * getMaxArrayValue
 *
 * Of the four possible locations for a local alignment, this function returns the
 * maximum case (i.e. 1-4) of the four values
 * ************************************************************
 */
int getMaxArrayValue(double array[], int length){
	double max = array[0];
	int iCase = 0; // which of the 4 cases contains the highest

	for (int i = 1; i < length; ++i){
		if (array[i] > max){
			max = array[i];
			iCase = i;
		}
	}
	return iCase;
}
