/*
 * Fragment.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: msegar
 */

#include "Fragment.h"

using namespace std;

Fragment::Fragment() {
	sAnchoredRead = "";
	iAnchoredStart = -1;
	iFlag = -1;
	sUnanchoredRead = "";
	iUnanchoredStart = -1;
	sParentRead = "";
	iParentStart = -1;
	bAnchorLeft = true;
}

Fragment::Fragment(string anchor, int anStart, int flag, bool left, string unanchor, int unStart){
	sAnchoredRead = anchor;
	iAnchoredStart = anStart;
	iFlag = flag;
	sUnanchoredRead = unanchor;
	iUnanchoredStart = unStart;
	sParentRead = "";
	iParentStart = -1;
	bAnchorLeft = left;
}

Fragment::Fragment(string anchor, int anStart, int flag, bool left, string unanchor, int unStart, string parent, int pStart){
	sAnchoredRead = anchor;
	iAnchoredStart = anStart;
	iFlag = flag;
	sUnanchoredRead = unanchor;
	iUnanchoredStart = unStart;
	sParentRead = parent;
	iParentStart = pStart;
	bAnchorLeft = left;
}


Fragment::~Fragment() {
	// Auto-generated destructor stub
}
