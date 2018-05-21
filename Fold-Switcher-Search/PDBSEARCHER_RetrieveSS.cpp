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

map<int,string> SS;
map<int,string> titles;
map<string,string> annotations;

string toupper (string s)
{
	string ret = "";
	int len = s.size();
	
	for (int i = 0; i < len; i++)
		ret += toupper(s[i]);
		
	return ret;
}

char isSpaceOrBracketOrTabOrControlCharacterOrComma(char c)
{
	if ((c == ' ') || (c == '{') || (c == '}') || (c == '\t') || ((int)c == 13) || (c == ','))
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
		while (i != s.size() && isSpaceOrBracketOrTabOrControlCharacterOrComma(s[i]))
			++i;

		// find the end of the next word
		string_size j = i;
		// invariant: none of the characters in the range [j^orig,j^cur) is a space or bracket, etc
		while (j != s.size() && !isSpaceOrBracketOrTabOrControlCharacterOrComma(s[j]))
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

void retrieveAnnotations (istream& in)
{
	string line;
	
	while (!(in.eof()))
	{
		getline(in,line);
		if (line != "" && line[0] == '>')
		{
			string PDB = toupper(line.substr(1,4));
			string chain = line.substr(6,1);
			string currDomain = PDB + ":" + chain;
			
			annotations[currDomain] = line;
		}
	}
	
	return;
}

void retrieveSS(istream& in)
{
	string line;
	int idx = 0;
	
	getline(in,line);

	while (!(in.eof()))
	{
		if (line != "" && line[0] == '>' && line.substr(8,3) == "sec")
		{	string currRead = "";
			string currName = line.substr(1,6);
		
			do
			{	getline(in,line);
				if (line != "" && line[0] != '>')
					currRead += line;
			} while (line[0] != '>' && !(in.eof()));
			titles[idx] = currName;
			SS[idx] = currRead;
			idx++;
		}
	}
	
	return;
}

string fixSS (string& s1, string& s2)
{
	string ret = "";
	
	int len2 = s2.size();
	
	int idx = 0;
	
	for (int i = 0; i < len2; i++)
	{
		if (s2[i] == '-')
			ret += " ";
		else
		{	ret += s1[idx];
			idx++;
		}
	}
	
	return ret;
}

int main(int argc, char *argv[])
{
	string inFileName,outFileName,ssFileName,annotationsFileName;
	string line;
	vector<string> words;
	
	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-out")
			outFileName = argv[++i];
		else if (curWord == "-ss")
			ssFileName  = argv[++i];
		else if (curWord == "-annotations")
			annotationsFileName = argv[++i];
	}

	cout << "Aligning secondary structure info from homologous proteins" << endl;
	
	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;

	ifstream ssfile ((char*)ssFileName.c_str());
	if (ssfile == NULL)
		cout << "ERROR: couldn't find file " << ssFileName << endl;
	ifstream annotationsfile ((char*)annotationsFileName.c_str());
	if (annotationsfile == NULL)
		cout << "ERROR: couldn't find file " << annotationsFileName << endl;
	
	ofstream outfile ((char*)outFileName.c_str());

	retrieveSS(ssfile);
	cout << "Done with retrieveSS" << endl;
	
	retrieveAnnotations(annotationsfile);
	cout << "Done with retrieveAnnotations" << endl;

	while(!(infile.eof()))
	{
		getline(infile,line);
		if (line.substr(0,6) == "Query_")
		{
			do
			{
				getline(infile,line);
				if (line == "")
					break;
				
				words = split(line);
				
				int currIdx = atoi(words[0].c_str());
				int currStart = atoi(words[1].c_str()) - 1;
				int currEnd = atoi(words[3].c_str());
				string currAlign = words[2];
				
				int currLen = currEnd - currStart + 1;
				
				string tempSS = SS[currIdx].substr(currStart,currLen);
				string fixedSS = fixSS(tempSS,currAlign);
				string currName = titles[currIdx];
				
				if (annotations[currName] == "")
					annotations[currName] = ">" + currName;
								
				outfile << annotations[currName] << endl << currAlign << endl << fixedSS << endl;
				
			} while (line != "" && !(infile.eof()));
		}
	}
	
	infile.close();
	outfile.close();
	ssfile.close();
	annotationsfile.close();
}

#endif
