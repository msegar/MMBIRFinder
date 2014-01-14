/*
 * ConfigDB.h
 *
 *  Created on: Dec 12, 2012
 *      Author: msegar
 */

#ifndef CONFIGDB_H_
#define CONFIGDB_H_

#include <string>
#include <map>

class ConfigDB {
public:
	ConfigDB(); // default constructor
	ConfigDB(const char *fileName);
	virtual ~ConfigDB();

	struct Gen_Val {
		bool boolVal;
		double doubleVal;
		unsigned long ulongVal;
		int intVal;
		std::string stringVal;
		char charVal;
	};

	std::string configFileName;

	Gen_Val getKey(const std::string key);
	bool checkKey (const std::string key);
	void setKey (std::string key, std::string val);
	void initializeDB(const char *fileName);
	void setConfigFile(const std::string &name);
	void printConfigDB();

private:
	std::map<std::string, std::string> keyVal;
};


#endif /* CONFIGDB_H_ */
