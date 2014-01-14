/*
 * ConfigDB.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: msegar
 */


#include "ConfigDB.h"

#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

ConfigDB::ConfigDB() {
	keyVal.clear();
	configFileName.clear();
}

ConfigDB::ConfigDB(const char *fileName) {
	initializeDB(fileName);
}

void ConfigDB::initializeDB(const char *fileName) {
	string line; // line reading in
	string key;
	string val;
	string::size_type pos;
	string::size_type posEq;

	keyVal.clear();

	ifstream input;
	input.open(fileName);

	if (input.fail()){
		cout << "Failing here (configDB.cpp)" << endl;
		cout << "Unable to open file (configDB.cpp)" << endl;
		exit(1);
		cout << "***** Something went horribly wrong ***********" << endl;
	}

	while (input.good()){
		getline(input, line);

		// null line or commented line
		posEq = line.find('=');
		if (string::npos == posEq || line[0] == '#')
			continue;

		// trim whitespace before '='
		while ((pos = line.find_first_of(" \t\r\n")) != string::npos && pos < posEq){
			line.erase(pos, 1);
			posEq--;
		}

		// trim whitespace after characters
		while ((pos = line.find_first_of(" \t\r\n", line.length() - 1)) != string::npos){
			line.erase(pos, 1);
		}

		posEq++;
		pos = line.find_first_not_of(" \t\r\n", posEq);

		if (pos == string::npos){
			cout << "Improper format: " << line << endl;
			exit(1);
		}

		if (pos > posEq)
			line.erase(posEq, pos- posEq);

		pos = line.find('=');
		key = line.substr(0, pos);

		//map<string, string>::iterator itr = keyVal.find(key);

		// don't add duplicates
		if (keyVal.find(key) != keyVal.end()){
			continue;
			cout << "uh oh" << endl;
		}

		pos++; // skip '=' sign

		// don't add key with no value
		if (line.length() == pos){
			cout << "Improper format: " << line << endl;
			exit(1);
		}

		keyVal[key] = line.substr(pos);

	}

	input.close();
}

bool ConfigDB::checkKey(const string key){
	if (keyVal.find(key) != keyVal.end())
		return true;
	else
		return false;
}

ConfigDB::Gen_Val ConfigDB::getKey(const string key){
	map<string, string>::iterator iter = keyVal.find(key);

	if (iter == keyVal.end()){
		cout << "Failed key lookup: " + key << endl;
		exit(1);
	}

	ConfigDB::Gen_Val val;

	val.stringVal = iter->second;
	val.intVal = atoi(iter->second.c_str());
	val.ulongVal = val.intVal;
	//val.boolVal = val.intVal;
	val.doubleVal = atof(iter->second.c_str());
	val.charVal = (iter->second)[0];

	if (val.stringVal.compare("true") == 0)
		val.boolVal = true;
	if (val.stringVal.compare("false") == 0)
		val.boolVal = false;

	return val;
}

void ConfigDB::setKey(string key, string val){
	string::size_type pos;

	// trim whitespaces for key and val
	while ((pos = key.find_first_of(" \t\r\n")) != string::npos)
		key.erase(pos, 1);

	while ((pos = val.find_first_of(" \t\r\n")) != string::npos)
		val.erase(pos, 1);

	map<string, string>:: iterator iter = keyVal.find(key);

	// if present, erase
	if (iter != keyVal.end()){
		keyVal.erase(iter);
	}

	keyVal[key] = val;
}

void ConfigDB::setConfigFile(const std::string &name){
	configFileName = name;
}

void ConfigDB::printConfigDB(){
	//for_each(keyVal.begin(), keyVal.end(), printFxn);

	map<string, string>::const_iterator itr;
	for(itr = keyVal.begin(); itr != keyVal.end(); ++itr){
		cout << itr->first << "\t" << itr->second << "\n";
	}
}

ConfigDB::~ConfigDB() {
	// Auto-generated destructor stub
}


