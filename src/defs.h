/*
 * defs.h
 *
 *  Created on: Dec 12, 2012
 *      Author: msegar
 */

#ifndef DEFS_H_
#define DEFS_H_

#include "ConfigDB.h"
#include "Fragment.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>

using namespace std;

struct t_chromosome {
	string fastaHeader;
	string samHeader;
	string sequence;
};

struct t_consolidated{
	string sReadName;
	string sParentRead;
	int iParentStart;
	int iParentEnd;
	int iBirStart;
	int iBirEnd;
	string sBir;
	int iBirLength;
	int iTemplateStart;
	int iTemplateEnd;
	string sTemplate;
	int iTemplateLength;
	bool bBirCandidateFound;
	int iChromosome;
	bool bAnchorLeft;
	int iFlag;
	bool bBadRead;
};

struct t_alignment_struct{
	int iEndPosI;
	int iStartPosI;
	int iEndPosJ;
	int iStartPosJ;
	string sAlignedRegionI;
	string sAlignedRegionJ;
};

// globals
#ifdef MAIN_CPP_
ConfigDB confDB;
vector<t_consolidated> vCandidateReads; // map the read name to the children fragments
vector<vector <t_consolidated> > vConsolidated; // vector of CONSOLIDATED BIR locations
vector<t_consolidated> vCandidateRegions;
vector<t_consolidated> vFinalBirLocs;
vector<t_chromosome> vReferenceGenome; // the reference genome
int time0;
string sJobId;
string sBaseFileName;
string sProjectDirectory;
ofstream fLogFileOut;

#else
extern ConfigDB confDB;
extern vector<t_consolidated> vCandidateReads; // map the read name to the children fragments
extern vector<vector <t_consolidated> > vConsolidated; // vector of CONSOLIDATED BIR locations
extern vector<t_consolidated> vCandidateRegions;
extern vector<t_consolidated> vFinalBirLocs;
extern vector<t_chromosome> vReferenceGenome; // the reference genome
extern int time0;
extern string sJobId;
extern string sBaseFileName;
extern string sProjectDirectory;
extern ofstream fLogFileOut;

#endif /* MAIN_CPP_ */


// main
void readInReferenceGenome();
int startExecutables();
int processCommandLine(int argc, char* argv[]);
void printUsageAndExit(char *sName);
string setupBaseFileName();
void prepareLogFile(ofstream& fout, string sBaseFileName, string sLogDirName);

// programExecution
int startExecutables();
int executeBwaIndex(string s);
int executeBwaAligner(string s1, string s2, string s3);
int getReads(string s1, string s2, string s3, string s4, bool b1);
int filterOut(string s1, string s2, string s3, string s4);
int convertSAMtoFASTA(string s);

// birFinder
int startCandidateReads();
int createSplitReadDatabase(string s1);
int createParentReads();

// birAligner
int startBirFinder();
int getLastPosition(string &sReference, string &sAligned);
void getBirLoc(t_alignment_struct tAligned, int st, int i);

// birConsolidate
int startConsolidate();
void consolidateLocations();
void getConsensus();
void consolidateBaseCalls();
bool compareStart(const t_consolidated &a, const t_consolidated &b);

// templateFinder
int startTemplateFinder();
string getReverseComplement(string &str);

// alignment
t_alignment_struct getLocalAlignment(string seq1, string seq2, double d1, double d2);
t_alignment_struct getGlobalAlignment(string &seq1, string &seq2, double d1, double d2);
double getSimilarityScore(char a, char b);
int getMaxArrayValue(double array[], int length);

// printOptions
void printConfig(ostream& out);
void printMatches(ostream& out, t_alignment_struct &a);
void printError(ostream& out, string &bir, string &ref, t_alignment_struct tAligned);
void printAlignment(ostream& out, string &bir, string &ref, int readStart, t_alignment_struct tAligned);
void printFinal();
void printReferenceGenomeInfo(ostream& out);
void printCandidateReadsInfo(ostream& out);

#endif /* DEFS_H_ */
