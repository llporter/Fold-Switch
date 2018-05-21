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
	string inFileName,curFileName;	

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
	}

	cout << "Parsing PDB protein sequences" << endl;
	
	string line,command;

	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;

	while (getline(infile,line))
	{	if (line[0] == '>') // we have a new sequence
		{	string PDBName = line.substr(1,4) + "_" + line.substr(6,1);
			curFileName = PDBName + ".pro";
			command = "echo '>" + PDBName + "' > " + curFileName;
			system((char*)command.c_str());
		}
		else
		{	command = "echo " + line + " >> " + curFileName;
			system((char*)command.c_str());
		}
	}

	infile.close();
}

#endif
