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
const float minHammDiff = 0.75;
const int   minGraftSize = 12;

map<int,string> SS;
map<int,string> titles;
map<string,string> annotations;

int countSpacesFront (string s)
{
	int ret = 0;
	int len = s.size();
	for (int i = 0; i < len; i++)
	{	if (s[i] == ' ')
			ret++;
		else
			break;
	}
	return ret;
}

int countSpacesRear (string s)
{
	int ret = 0;
	int len = s.size();
	for (int i = len - 1; i >= 0; i--)
	{	if (s[i] == ' ')
			ret++;
		else
			break;
	}
	return ret;
}

string spaces (int len)
{
	string retString="";
	for (int i = 0; i < len; i++)
		retString += " ";
	return retString;
}

string toupper (string s)
{
	string ret = "";
	int len = s.size();
	
	for (int i = 0; i < len; i++)
		ret += toupper(s[i]);
		
	return ret;
}

bool perfectMatch (string s)
{
	int len = s.size();
	
	for (int i = 0; i < len; i++)
		if (s[i] != '.' && s[i] != '-')
			return FALSE;
			
	return TRUE;
}

float differenceScore (char c1, char c2)
{
	if ((c1 == 'E' && c2 == 'B') || (c1 == 'B' && c2 == 'E'))
		return 0.5;
	else if ((c1 == 'S' || c1 == 'T' || c1 == 'H' || c1 == 'G') && (c2 == 'S' || c2 == 'T' || c2 == 'H' || c2 == 'G'))
		return 0.5;
	else
		return 1.0;
}

float hammingDistance (string s1, string s2)
{
	int len1 = s1.size();
	int len2 = s2.size();
	
	if (len1 < len2)
	{	
		int diff = len2 - len1;
		s1 += spaces(diff);
	}
	else if (len2 < len1)
	{	
		int diff = len1 - len2;
		s2 += spaces(diff);
	}

	int ret = 0;
	
	for (int i = 0; i < len1; i++)
		if (s1[i] != s2[i] && s1[i] != ' ' && s2[i] != ' ')
			ret += differenceScore(s1[i],s2[i]);
			
	return ret;
}

int nonBlank (string s1, string s2)
{
	int len1 = s1.size();
	int len2 = s2.size();
	
	if (len1 < len2)
	{	
		int diff = len2 - len1;
		s1 += spaces(diff);
	}
	else if (len2 < len1)
	{	
		int diff = len1 - len2;
		s2 += spaces(diff);
	}

	int ret = 0;
	
	for (int i = 0; i < len1; i++)
		if (s1[i] != ' ' && s2[i] != ' ')
			ret++;
			
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
	//cout << "entering fixSS with " << s1 << " and " << s2 << endl;
	
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
	
	//cout << "returning " << ret << endl;
	
	return ret;
}

int main(int argc, char *argv[])
{
	string inFileName,outFileName,ssFileName,annotationsFileName;
	string line,line2,line3,line4;
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

	cout << "Finding fold-switching proteins" << endl;
	
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

	int PDBIdx = 1;

	getline(infile,line);
	while(!(infile.eof()))
	{
		if (line != "" && line[0] == '>')
		{
			string currRead = "";
			string PDB = line.substr(1,4);
			string chain = line.substr(6,1);
			string currDomain = PDB + "_" + chain;

			cout << "working on " << currDomain << endl;

			do
			{	getline(infile,line);
				if (line != "" && line[0] != '>')
					currRead += line;
			} while (line[0] != '>' && !(infile.eof()));
			
			string command = "echo '>" + currDomain + "' > " + currDomain + ".pro";
			system((char*)command.c_str());
			command = "echo " + currRead + " >> " + currDomain + ".pro";
			system((char*)command.c_str());
			
			command = "/groups/looger/home/loogerl/src/ncbi-blast-2.4.0+/bin/blastp -evalue 1e-10 -num\_descriptions 10000000 -num\_alignments 10000000 -outfmt 3 -db ss.seqs -query " + currDomain + ".pro > " + currDomain + ".log";
			system((char*)command.c_str());
			
			string inFileName2 = currDomain + ".log";
			ifstream infile2 ((char*)inFileName2.c_str());
			if (infile2 == NULL)
				cout << "ERROR: couldn't find file " << inFileName2 << endl;
				
			string outFileName2 = currDomain + ".log2";
			ofstream outfile2 ((char*)outFileName2.c_str());			
			
			bool hitFound = FALSE;
			while(!(infile2.eof()))
			{
				getline(infile2,line2);
				if (line2.substr(0,6) == "Query_")
				{
					words = split(line2);
					string templateName = words[0];
					int templateStart = atoi(words[1].c_str());
					string templateRead = words[2];
					int templateSize = templateRead.size();
					int templateEnd = atoi(words[3].c_str());
					int readTab = line2.find(templateRead);
					//outfile2 << endl << endl << "Template:  " << templateRead << endl;
					do
					{
						getline(infile2,line2);
						
						if (line2 == "")
							break;
				
						words = split(line2);
						if (words.size() < 4)
							continue;
				
						int currIdx = atoi(words[0].c_str());
						int currStart = atoi(words[1].c_str()) - 1;
						int currEnd = atoi(words[3].c_str());
						string currAlign = line2.substr(readTab,templateSize);
						string tempAlign = words[2];
						int hitTabFront = countSpacesFront(currAlign);
						int hitTabRear  = countSpacesRear(currAlign);
						hitFound = TRUE;
					
						int currLen = currEnd - currStart + 1;
				
						string currName = titles[currIdx];
						string tempSS="",fixedSS="";
						if (SS[currIdx].size() >= currEnd)
						{
							tempSS = SS[currIdx].substr(currStart,currLen);
							fixedSS = spaces(hitTabFront) + fixSS(tempSS,tempAlign) + spaces(hitTabRear);
						}
						else
						{
							cout << "ERROR: for " << currDomain << ", idx " << currIdx  << ", name " << currName << ", SS is not long enough." << endl << "line = " << line2 << endl;
						}
							
						if (annotations[currName] == "")
							annotations[currName] = ">" + currName;
						
						outfile2 << annotations[currName] << endl << currAlign << endl << fixedSS << endl;
						
					} while (line2 != "" && !(infile2.eof()));
					outfile2 << endl;
				}
			}
			infile2.close();
			outfile2.close();
			
			if (hitFound == FALSE)
			{	cout << "No hits found for " << currDomain << endl;
				continue;
			}
			
			string inFileName3 = currDomain + ".log2";
			ifstream infile3 ((char*)inFileName3.c_str());
			if (infile3 == NULL)
				cout << "ERROR: couldn't find file " << inFileName3 << endl;
				
			string outFileName3 = currDomain + ".log3";
			ofstream outfile3 ((char*)outFileName3.c_str());			
			
			while (!(infile3.eof()))
			{	getline(infile3,line3);
				if (line3 == "")
				{	outfile3 << endl;
					getline(infile3,line3);
				}
				
				if (line3[0] == '>')
				{
					string currName = line3;
					getline(infile3,line3);
					string currAlign = line3;
					getline(infile3,line3);
					string currSS = line3;
					
					if (perfectMatch(currAlign))
					{	outfile3 << currName << endl << currAlign << endl << currSS << endl;
					}
				}
			}
			infile3.close();
			outfile3.close();

			string inFileName4 = currDomain + ".log3";
			ifstream infile4 ((char*)inFileName4.c_str());
			if (infile4 == NULL)
				cout << "ERROR: couldn't find file " << inFileName4 << endl;
				
			string outFileName4 = currDomain + ".log4";
			ofstream outfile4 ((char*)outFileName4.c_str());			
			
			vector<string> names;
			vector<string> aligns;
			vector<string> ss;
			
			while (!(infile4.eof()))
			{	getline(infile4,line4);
				if (line4 == "")
				{	int numSeqs = ss.size();
					int len = ss[0].size();
					for (int i = 0; i < numSeqs - 1; i++)
					{
						string SSi = ss[i];
						
						for (int j = i+1; j < numSeqs; j++)
						{
							string SSj = ss[j];
							float hammDist = hammingDistance(SSi,SSj);
							int nonBlankResidues = nonBlank(SSi,SSj);
							float hammDiff = (float)(hammDist/(float)(nonBlankResidues));
							if (hammDist > 3 && hammDiff > minHammDiff && nonBlankResidues > minGraftSize)
							{	outfile4 << names[i] << " & " << names[j] << " might fold switch:" << endl;
								outfile4 << SSi << endl << SSj << endl << endl;
								outfile4 << "hammDist = " << hammDist << ", hammDiff = " << hammDiff << endl;
							}
						}
					}
				
					names.clear();
					aligns.clear();
					ss.clear();
					getline(infile4,line4);
				}

				if (line4[0] == '>')
				{
					names.push_back(line4);
					getline(infile4,line4);
					aligns.push_back(line4);
					getline(infile4,line4);
					ss.push_back(line4);
				}
			}
			
			cout << "Done with PDB " << PDBIdx << " " << currDomain << endl;
			PDBIdx++;
		}
	}
	
	infile.close();
	outfile.close();
	ssfile.close();
	annotationsfile.close();
}

#endif
