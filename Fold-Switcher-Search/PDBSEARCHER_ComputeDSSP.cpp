#ifndef ALIGNEDGENOMEANALYZER_CPP
#define ALIGNEDGENOMEANALYZER_cpp

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{
	string inFileName;
	string line;
	
	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
	}

	cout << "Computing DSSP" << endl;
	
	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;

	while(!(infile.eof()))
	{
		getline(infile,line);
		if (line == "")
			break;
			
		string currFileName = line;
		
		string command = "gunzip " + currFileName;
		system((char*)command.c_str());
		
		int len = currFileName.size();
		string root = currFileName.substr(0,len-7);
		
		command = "dssp -i " + root + ".ent -o " + root + ".ss";
		system((char*)command.c_str());
		
		command = "gzip " + root + ".ent";
		system((char*)command.c_str());
	}
	
	infile.close();
}

#endif
