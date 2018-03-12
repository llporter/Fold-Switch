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

string itos(int n)
{
	int sign;
	string s;

	sign = n;
	n = abs(n);
	do {                          		// generate digits in reverse order
		s += (char(n % 10) + '0');	// get next digit
	} while ((n/=10) > 0);        		// delete it

	if (sign < 0)
		s += '-';

	reverse(s.begin(), s.end());	// This is what the code should look like
                                	// if the string class is compatible with
                                	// the standard C++ string class
	return s;
}

int main(int argc, char *argv[])
{
	string inFileName,outFileName;
	int size;

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-out")
			outFileName = argv[++i];
		else if (curWord == "-size")
			size = atoi(argv[++i]);
	}

	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;
	ofstream outfile ((char*)outFileName.c_str());

	string line;

	int currDir = 1;
	int currIdx = -1;
	
	while (getline(infile,line))
	{	if (line == "")
			break;
		outfile << "PDBSEARCHER_FindFoldSwitchersSingleSequence -ss ss.ss -annotations PDB_annotations.txt -pdb Dir" << itos(currDir) << "/" << line << endl;
		currIdx++;
		if (currIdx == size)
		{
			currIdx = 0;
			currDir++;
		}
	}

	infile.close();
	outfile.close();
}

#endif
