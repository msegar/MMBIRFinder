/*
 * Fragment.h
 *
 *  Created on: Jan 7, 2013
 *      Author: msegar
 */

#ifndef FRAGMENT_H_
#define FRAGMENT_H_

#include <string>

using namespace std;

class Fragment {
private:
	string sAnchoredRead;
	int iAnchoredStart;
	int iFlag;
	string sUnanchoredRead;
	int iUnanchoredStart;
	string sParentRead;
	int iParentStart;
	bool bAnchorLeft;
public:
	Fragment();
	Fragment(string anchor, int anStart, int flag, bool left, string unanchor, int unStart);
	Fragment(string anchor, int anStart, int flag, bool left, string unanchor, int unStart, string parent, int pStart);

	virtual ~Fragment();

	// get functions
	bool isAnchorLeft()		{ return bAnchorLeft; }
	int getAnchoredStart() 	{ return iAnchoredStart; }
	int getFlag()			{ return iFlag; }
	int getUnanchoredStart() { return iUnanchoredStart; }
	bool getAnchorLeft()		{ return bAnchorLeft; }
	int getParentStart() 	{ return iParentStart; }
	string* getAnchoredRead() { return &sAnchoredRead; }
	string* getParentRead() { return &sParentRead; }
	string* getUnanchoredRead() { return &sUnanchoredRead; }

	// set functions
	void setParentRead(string parentRead) { sParentRead = parentRead; }
	void setParentStart(int parentStart) { iParentStart = parentStart; }
	//void setUnanchoredRead(string unanchoredRead) { sUnanchoredRead = unanchoredRead; }
};

#endif /* FRAGMENT_H_ */
