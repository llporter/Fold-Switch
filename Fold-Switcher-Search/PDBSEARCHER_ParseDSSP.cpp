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

string toupper (string s)
{
	string ret = "";
	int len = s.size();
	
	for (int i = 0; i < len; i++)
		ret += toupper(s[i]);
		
	return ret;
}

char isSpaceOrBracketOrTabOrControlCharacter(char c)
{
	if ((c == ' ') || (c == '{') || (c == '}') || (c == '\t') || ((int)c == 13))
		return TRUE;
	else
		return FALSE;
}

vector<string> split (const string& s)
{
	vector<string> ret;
	typedef string::size_type string_size;
	string_size i = 0;

	// invariant: we have processed characters in the range [original value of i,i).
	while (i != s.size())
	{
		//ignore leading blanks
		// invariant: characters in the range [i^orig,i^cur) are all spaces or brackets, etc
		while (i != s.size() && isSpaceOrBracketOrTabOrControlCharacter(s[i]))
			++i;

		// find the end of the next word
		string_size j = i;
		// invariant: none of the characters in the range [j^orig,j^cur) is a space or bracket, etc
		while (j != s.size() && !isSpaceOrBracketOrTabOrControlCharacter(s[j]))
			++j;

		// if we found some non-whitespace and non-bracket characters
		if (i != j)
		{
			// copy from s starting at i and taking j-i characters
			ret.push_back(s.substr(i,j-i));
			i = j;
		}
	}
	return ret;
}


int main(int argc, char *argv[])
{
	string inFileName,seqFileName,ssFileName,bothFileName;
	string line,line2;
	vector<string> words;
	
	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-seq")
			seqFileName = argv[++i];
		else if (curWord == "-ss")
			ssFileName  = argv[++i];
		else if (curWord == "-both")
			bothFileName = argv[++i];
	}

	cout << "Parsing DSSP .ss output files" << endl;
	
	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;

	ofstream seqfile   ((char*)seqFileName.c_str());
	ofstream ssfile    ((char*)ssFileName.c_str());
	ofstream bothfile ((char*)bothFileName.c_str());

	while(!(infile.eof()))
	{
		getline(infile,line);
		
		if (line == "")
			break;
		
		string currFileName = line;
		int len = currFileName.size();
		string root = currFileName.substr(0,len-7);
		string PDB  = toupper(root.substr(3,4));
		
		string currSSFileName = root + ".ss";
		ifstream currSsfile ((char*)currSSFileName.c_str());
		if (currSsfile == NULL)
		{	cout << currSSFileName << " not found. Likely that this is not protein." << endl;
			ssfile.close();
			continue;
		}
		
		string prevChain  = "";
		
		while (!(currSsfile.eof()))
		{
			do
			{	getline(currSsfile,line2);
			}	while (line2.size() < 12 || line2.substr(0,12) != "  #  RESIDUE");
			
			string currSeq,currSS,currChain,currAA,currSSRes;
			int currNum,currRes;
			
			do
			{	getline(currSsfile,line2);
				words = split(line2);
				
				if (words.size() < 5)
					break;
					
				cout << "PDB = " << PDB << ", line = " << line << endl;
				
				if (words[1] == "!*")	// end of a subunit
				{
					seqfile << ">" << PDB << ":" << currChain << ":sequence" << endl;
					seqfile << currSeq << endl;
					bothfile << ">" << PDB << ":" << currChain << ":sequence" << endl;
					bothfile << currSeq << endl;

					ssfile << ">" << PDB << ":" << currChain << ":secstr" << endl;
					ssfile << currSS << endl;
					bothfile << ">" << PDB << ":" << currChain << ":secstr" << endl;
					bothfile << currSS << endl;
					continue;		
				}
					
				currNum = atoi(words[0].c_str());
				currRes = atoi(words[1].c_str());
				currChain = words[2];
				currAA = words[3];
				currSSRes = words[4];
				
				if (currChain != prevChain)
				{
					currSeq = "";
					currSS  = "";
				}
				else
				{
					currSeq += currAA;
					currSS  += currSSRes;
				}
			} while (line2 != "");
		}
		currSsfile.close();
	}
	
	infile.close();
	seqfile.close();
	ssfile.close();
	bothfile.close();
}

#endif
