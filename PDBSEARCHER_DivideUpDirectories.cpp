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
	string inFileName,line;
	int size;	

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-size")
			size = atoi(argv[++i]);
	}

	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;

	ofstream pdbfile ("PDBlist.txt");

	int currDir = 1;
	int currIdx = -1;
	string command = "mkdir Dir1";
	system((char*)command.c_str());

	string PDBName = "";
	string outFileName = "";

	while (!(infile.eof()))
	{	getline(infile,line);
		if (line == "")
			break;
		if (line[0] == '>') // we have a new sequence
		{	PDBName = line.substr(1,4) + "_" + line.substr(6,1);
			outFileName = PDBName + ".pro";
			pdbfile << PDBName << endl;
			currIdx++;
			if (currIdx == size)
			{
				currIdx = 0;
				currDir++;
				command = "mkdir Dir" + itos(currDir);
				system((char*)command.c_str());
			}
			command = "echo '>" + PDBName + "' >> Dir" + itos(currDir) + "/" + outFileName;
			system((char*)command.c_str());
		}
		else
		{	command = "echo " + line + " >> Dir" + itos(currDir) + "/" + outFileName;
			system((char*)command.c_str());
		}
	}
	infile.close();
	pdbfile.close();
}

#endif
