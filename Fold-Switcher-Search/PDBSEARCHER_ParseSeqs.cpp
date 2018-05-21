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
const bool FALSE = 0;
const bool TRUE = 1;

int main(int argc, char *argv[])
{
	string inFileName,outFileName,ssFileName;
	string line;
	
	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-out")
			outFileName = argv[++i];
		else if (curWord == "-ss")
			ssFileName  = argv[++i];
	}

	cout << "Making .seqs file from PDB secondary structure download" << endl;
	
	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;
	
	ofstream outfile ((char*)outFileName.c_str());
	ofstream ssfile  ((char*)ssFileName.c_str());

	getline(infile,line);
	while (!(infile.eof()))
	{
		if (line != "" && line[0] == '>' && line.substr(8,3) == "seq")
		do
		{	outfile << line << endl;
			getline(infile,line);
		} while (line != "" && line[0] != '>');
		
		do
		{	ssfile << line << endl;
			getline(infile,line);
		} while (line != "" && line[0] != '>');

	}
	infile.close();
	outfile.close();
	ssfile.close();	
}

#endif
