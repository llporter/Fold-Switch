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
#include <string>
#include <boost/assign/list_of.hpp>

using namespace std;
const bool FALSE = 0;
const bool TRUE = 1;
const float minHammDiff = 0.5;
const int   minGraftSize = 12;
float minHammDist = minGraftSize * minHammDiff;
int numMutationsAllowed;

map<char,float> scorer = boost::assign::map_list_of('*',1.0)('%',0.5)('^',0.25)('-',0.0)('.',0.0);
map<int,string> SS;
map<int,string> titles;
map<string,string> annotations;
map<string,string> aligns;
map<string,string> SSs;
map<string,string> BLASTs;
vector<string> names;
vector<vector<int> > lines;

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

template <class T>
bool isMemberOfVector (vector<T>& vec, T& elem)
{
	int size = vec.size();
	for (int i = 0; i < size; i++)
		if (elem == vec[i])
			return TRUE;
	
	return FALSE;
}

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

string dashes (int len)
{
	string retString="";
	for (int i = 0; i < len; i++)
		retString += "-";
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

bool perfectMatch (string& s)
{
	int len = s.size();
	
	for (int i = 0; i < len; i++)
		if (s[i] != '.' && s[i] != '-' && s[i] != ' ')
			return FALSE;
			
	return TRUE;
}

int mutations (string& s)
{
	int len = s.size();
	int ret = 0;
	for (int i = 0; i < len; i++)
		if (s[i] != '.' && s[i] != '-' && s[i] != ' ')
			ret++;
			
	return ret;
}

float differenceScore (char c1, char c2)
{
	if (c1 == c2)
		return 0.0;
	else if (c1 == '-' || c2 == '-')
		return 0.0;
	else if ((c1 == 'E' && c2 == 'B') || (c1 == 'B' && c2 == 'E'))
		return 0.25;
	else if ((c1 == 'S' || c1 == 'T' || c1 == 'H' || c1 == 'G' || c1 == 'I') && (c2 == 'S' || c2 == 'T' || c2 == 'H' || c2 == 'G' || c2 == 'I'))
		return 0.25;
	else if (c1 == ' ' || c2 == ' ')
		return 0.5;
	else
		return 1.0;
}

pair<float,string> hammingDistance (string s1, string s2)
{
	int len1 = s1.size();
	int len2 = s2.size();
	string diffString = "";
	
	if (len1 < len2)
	{	
		int diff = len2 - len1;
		s1 += dashes(diff);
	}
	else if (len2 < len1)
	{	
		int diff = len1 - len2;
		s2 += dashes(diff);
	}

	float ret = 0.0;
	
	for (int i = 0; i < len1; i++)
	{	float currScore = differenceScore(s1[i],s2[i]);
		ret += currScore;
		if (currScore == 1.0)
			diffString += "*";
		else if (currScore == 0.5)
			diffString += "%";
		else if (currScore == 0.25)
		        diffString += "^";
		else if (s1[i] == '-' || s2[i] == '-')
			diffString += "-";
		else if (s1[i] == s2[i])
			diffString += ".";
	}
			
	return pair<float,string> (ret,diffString);
}

int nonBlank (string s1, string s2)
{
	int len1 = s1.size();
	int len2 = s2.size();

	//cout << "DEBUG: s1 = " << s1 << endl << "s2 = " << s2 << endl << "len1 = " << len1 << ", len2 = " << len2 << endl;
	
	if (len1 < len2)
	{	
		int diff = len2 - len1;
		s1 += dashes(diff);
	}
	else if (len2 < len1)
	{	
		int diff = len1 - len2;
		s2 += dashes(diff);
	}

	int ret = 0;
	
	for (int i = 0; i < len1; i++)
		if (s1[i] != '-' && s2[i] != '-')
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
			ret += "-";
		else
		{	ret += s1[idx];
			idx++;
		}
	}
	
	//cout << "returning " << ret << endl;
	
	return ret;
}

string replacePeriods (string& s1, string& s2)
{
	string ret = "";
	int len = s1.size();
	
	for (int i = 0; i < len; i++)
	{	if (s1[i] != '.')
			ret += s1[i];
		else
			ret += s2[i];
	}
	
	return ret;
}

float diffScore (string& s)
{
	int len = s.size();
	float ret = 0.0;

	for (int i = 0; i < len; i++)
	{	char curr = s[i];
		if (curr == '*')
			ret += 1.0;
		else if (curr == '%')
			ret += 0.5;
	}
	
	return ret;
}

pair<float,int> findMaximalLocalDifference (string& s)
{
	int len = s.size();
	float retFloat = 0.0;
	int retInt = 0;
	
	for (int i = 0; i < len - minGraftSize; i++)
	{
		string curr = s.substr(i,minGraftSize);
		float currDiff = diffScore(curr);
		if (currDiff > retFloat)
		{	retFloat = currDiff;
			retInt = i;
		}
	}
	
	retFloat = retFloat/minGraftSize;
	
	pair<float,int> ret = pair<float,int> (retFloat,retInt);
	
	return ret;
}

void juicy_differences(string&s,vector<string>&v1, vector<int>&v2, vector<int>&v3)
{

  string delimiter = ".";

  size_t pos = s.find(delimiter);
  size_t next_pos = -1;

  //cout << s << endl;

  if (pos == string::npos){

    v1.push_back(s.substr(0,pos));
    v3.push_back(-999);
    v2.push_back(0);
    return;
  }
  else if (s==string(s.size(),'.')){
    v1.push_back("");
    v3.push_back(-998);
    v2.push_back(0);
  }
    

  v1.push_back(s.substr(0,pos));
  if (pos != 0){
    v3.push_back(0);
    //v2.push_back(0);
  }
  //v3.push_back(pos+1);

  int i = 0;
  int delimiter_length = 0;
  int posOld = 0;

  while ((next_pos = s.find(delimiter,pos+1)) != string::npos) {

    delimiter_length++;

    //If there are successive delimiters, pos=next_pos (=pos+1) and continue

    if (next_pos-pos == 1){
      pos = next_pos;
      continue;
      }

    v1.push_back(s.substr(pos+1,next_pos-pos-1));
    v2.push_back(delimiter_length);
    v3.push_back(pos+1);

    posOld = pos;
    pos = next_pos;

    delimiter_length = 0;
  
  }

  if (s.size() > pos+1){
    v1.push_back(s.substr(pos+1,s.size()-1));
    v2.push_back(pos-posOld-v1[v1.size()-2].size());
    v3.push_back(pos+1);
  }
}

template <class T>
void print_vector(vector<T> &v1, string vop = "Vector output: "){

  for  (int i=0; i<v1.size();i++)
    cout << vop << v1[i] << endl;

}

float score_bit(string &s){
  float curScore = 0.0;
  int true_size = 0;

  for(int j=0; j<s.size(); j++){
    curScore += scorer[s[j]];
    if (s[j] != '-'){
      true_size++;
    }
  }
  if (true_size == 0) true_size = 1;

  return curScore;
}

int resize(string &s){

  int ssize = 0;

  for (int i=0;i<s.size();i++){

    if (s[i] != '-') ssize++;

  }

  return ssize;
}

int resize2(string &s){

  int ssize = 0;

  for (int i=0;i<s.size();i++){

    if (s[i] != '-' && s[i] != '.') ssize++;

  }

  return ssize;
}

void score_juicy_bits(vector<string> &v1, vector<float> &v2){

  for (int i=0; i<v1.size(); i++){
    v2.push_back(score_bit(v1[i]));
  }
}


float rescore(float cur_score,int gap_size,float next_score,string &max_string, string &next_string){

  int max_string_size = resize(max_string);
  int next_string_size = resize(next_string);

  int full_length = max_string_size+next_string_size+gap_size;
  float score= (cur_score+next_score)/full_length;

  //cout << max_string_size <<" " <<next_string_size<<" "<<gap_size<<" " << score<<endl;

  return cur_score+next_score-score*gap_size;

}

void attempt_to_extend(string max_string, int max_idx, vector<string> &juicy_strings, \
		       vector<float>&juicy_scores, vector<int> &gap_sizes, vector<int> &good_indices){

  float cur_score = juicy_scores[max_idx];
  int o1 = 1;
  int o2 = 1;
  int o = 0;
  float new_score1 = -1.0;
  float new_score2 = -1.0;
  string string1 = "";
  string string2 = "";
  int ssize1 = 0;
  int ssize2 = 0;

  good_indices.push_back(max_idx);

  if (gap_sizes.size() < juicy_scores.size()) o = 1;

  while (max_idx+o1 < juicy_scores.size() || max_idx-o2 >= 0){

    // print_vector(juicy_strings,"Juicy strings: ");
    //print_vector(gap_sizes, "Gap sizes: ");

    if (max_idx+o1 < juicy_scores.size()){

      new_score1 = rescore(cur_score,gap_sizes[max_idx+o1-o],juicy_scores[max_idx+o1], \
			   max_string,juicy_strings[max_idx+o1]);
      string1 = juicy_strings[max_idx+o1];

      ssize1 = resize(string1)+resize(max_string)+gap_sizes[max_idx+o1-o];
    }
    else{
      new_score1 = -50.0;
      string1 = "";
    }

    if (max_idx-o2 >= 0){

      new_score2 = rescore(cur_score,gap_sizes[max_idx-o2],juicy_scores[max_idx-o2], \
			   max_string,juicy_strings[max_idx-o2]);
      string2 = juicy_strings[max_idx-o2];

      ssize2 = resize(string2)+resize(max_string)+gap_sizes[max_idx-o2];
    }
    else{
      new_score2 = -50.0;
      string2 = "";
    }

    //cout <<"new_score1 : "<<new_score1<<" string1: "<< string1<<endl;
    //cout <<"new_score2 : "<<new_score2<<" string2: "<<string2<<endl;
      
    if (new_score2 < minHammDiff*ssize2 && new_score1 < minHammDiff*ssize1){
      break;
    }

    if (new_score1 >=new_score2){
      good_indices.push_back(max_idx+o1);
      cur_score = new_score1;
      max_string = string1;
      o1++;
    }
    else{
      good_indices.push_back(max_idx-o2);
      cur_score = new_score2;
      max_string = string2;
      o2++;
    }  
  }
} 

void find_optimal_indices(vector<string> &juicy_strings, vector<float>&juicy_scores, \
			 vector<int> &gap_sizes, vector<int> &good_indices){

  //print_vector(juicy_scores,"Juicy scores: ");
  int max_idx = distance(juicy_scores.begin(),max_element(juicy_scores.begin(),juicy_scores.end()));
  //cout << "Max idx: "<<max_idx<<".  Max string: "<<juicy_strings[max_idx]<<endl;
  string max_string = juicy_strings[max_idx];
  attempt_to_extend(max_string,max_idx,juicy_strings,juicy_scores,gap_sizes,good_indices);
  //cout << "Max idx " << max_idx << endl; 
}

void optimal_string(string &optimal_string, vector<string>&juicy_bits, \
		    vector<int> &juicy_spaces,vector<int> &optimal_indices){

  for (int i=0; i<optimal_indices.size(); i++){

    if (juicy_spaces.size() == juicy_bits.size()){
      if (i==0) optimal_string += juicy_bits[optimal_indices[i]];
      else optimal_string += string(juicy_spaces[optimal_indices[i]],'.')+juicy_bits[optimal_indices[i]];
    }

    else{
      if (i < optimal_indices.size()-1){
	optimal_string += juicy_bits[optimal_indices[i]]+string(juicy_spaces[optimal_indices[i]],'.');
      }
      else optimal_string += juicy_bits[optimal_indices[i]];
    }
  }
}

int optimal_idx(vector<int> &optimal_indices, vector<string>&juicy_bits,
		vector<int> &juicy_spaces){

  int idx = 0;

  if (optimal_indices[0] == 0 && juicy_spaces.size() < juicy_bits.size()) return 0;
  else{
    if (optimal_indices[0] == 0) return juicy_spaces[0];

    for (int i=0; i<optimal_indices[0];i++){
      idx += juicy_spaces[i]+juicy_bits[i].size();
      //cout <<"Juicy spaces: "<<juicy_spaces[i]<<" Juicy bits size: " << juicy_bits[i].size()<<endl;
    }

    //cout <<"Optimal idx: "<<idx<<endl;
    //if (juicy_spaces.size() < juicy_bits.size()) idx += juicy_spaces[optimal_indices[0]];
    if (juicy_spaces.size() == juicy_bits.size()) idx += juicy_spaces[optimal_indices[0]];
    //print_vector(optimal_indices,"Optimal Indices: ");
    //cout <<"Optimal idx: "<<idx<<endl;

    return idx;
  }
}

int main(int argc, char *argv[])
{
	string currDomain,outFileName,ssFileName,annotationsFileName;
	string line,line2,line3,line4;
	vector<string> words;
	
	//Process user inputs.

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-pdb")
			currDomain = argv[++i];
		else if (curWord == "-ss")
			ssFileName  = argv[++i];
		else if (curWord == "-annotations")
			annotationsFileName = argv[++i];
		else if (curWord == "-muts")
			numMutationsAllowed = atoi(argv[++i]);
	}

	cout << "Finding fold-switching proteins" << endl;
	
	//open SS file.  Print Error if NULL

	ifstream ssfile ((char*)ssFileName.c_str());
	if (ssfile == NULL)
		cout << "ERROR: couldn't find file " << ssFileName << endl;

	//Open Annotations file.  Print Error if NULL.

	ifstream annotationsfile ((char*)annotationsFileName.c_str());
	if (annotationsfile == NULL)
		cout << "ERROR: couldn't find file " << annotationsFileName << endl;
	
	//Get all SS annotations from .ss file and get PDB names and IDs also

	retrieveSS(ssfile);
	ssfile.close();
	cout << "Done with retrieveSS" << endl;

	//Map annotations to their respective PDBID and Domain
	
	retrieveAnnotations(annotationsfile);
	annotationsfile.close();
	cout << "Done with retrieveAnnotations" << endl;

	//Blast AA sequence of current PDB chain against the rest of the pdb PDB.  Save output as PDBID_Chain.log

	string command = "/groups/looger/home/loogerl/src/ncbi-blast-2.4.0+/bin/blastp -evalue 1e-10 -num\_descriptions 10000000 -num\_alignments 10000000 -outfmt 3 -db ss.seqs -query " + currDomain + ".pro > " + currDomain + ".log";
	
	//Open the .log file
	system((char*)command.c_str());
			
	string inFileName = currDomain + ".log";
	ifstream infile ((char*)inFileName.c_str());
	if (infile == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;
	

	//numLines 
	int numLines=0;

	//Determine number of lines in sequence alignment
	//Convert all BLAST ID numbers for hits to PDBids+chains
			
	while(!(infile.eof()))
	{
		getline(infile,line2);
		if (line2.substr(0,6) == "Query_")
		{
			numLines++;
			do
			{
				getline(infile,line2);
					
				if (line2 == "")
					break;
				
				words = split(line2);
				if (words.size() < 4)
					continue;
				
				int currIdx = atoi(words[0].c_str());				
				string currName = titles[currIdx];
				if (annotations[currName] == "")
					annotations[currName] = ">" + currName;
						
				if (!(isMemberOfVector(names,currName)))
				{	names.push_back(currName);
					annotations[currName] += " (" + itos(currIdx) + ")";
				}
							
						
			} while (line2 != "" && !(infile.eof()));
		}
	}
	infile.close();

	int numNames = names.size();
	vector<vector<string> > linesUsed(numLines);

	ifstream infile2 ((char*)inFileName.c_str());
	if (infile2 == NULL)
		cout << "ERROR: couldn't find file " << inFileName << endl;
	string outFileName2 = currDomain + ".log2";
	ofstream outfile2 ((char*)outFileName2.c_str());			
			
	bool hitFound = FALSE;
	string Query = "";
	int currLine = -1;
	int templateSize;
	
	while(!(infile2.eof()))
	{
		getline(infile2,line2);
		if (line2.substr(0,6) == "Query_")
		{
			currLine++;
			words = split(line2);
			string templateName = words[0];
			// int templateStart = atoi(words[1].c_str());
			string templateRead = words[2];
			templateSize = templateRead.size();
			// int templateEnd = atoi(words[3].c_str());
			int readTab = line2.find(templateRead);
			Query += templateRead;
			//outfile2 << "Template: " << templateRead << endl;
			
			int prevIdx = -1;
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
				string BLAST = line2;
					
				int currLen = currEnd - currStart + 1;
			
				string currName = titles[currIdx];
				linesUsed[currLine].push_back(currName);
				string tempSS="",fixedSS="";
				if (SS[currIdx].size() >= currEnd)
				{
					tempSS = SS[currIdx].substr(currStart,currLen);
					fixedSS = dashes(hitTabFront) + fixSS(tempSS,tempAlign) + dashes(hitTabRear);
				}
				else
					cout << "ERROR: for " << currDomain << ", idx " << currIdx  << ", name " << currName << ", SS is not long enough." << endl << "line = " << line2 << endl;
						
				if (currIdx == prevIdx) // this is a duplicated BLAST line
				{	//cout << "Duplicate: " << currDomain << ", " << currIdx << ", " << currName << endl << "from " << currStart << " to " << currEnd << " = " << currAlign << endl << fixedSS << endl << BLAST << endl;
				}
				else		
				{	aligns[currName] += currAlign;
					SSs[currName] += fixedSS;
					BLASTs[currName] += BLAST + " ";
				}
				
				prevIdx = currIdx;
						
				//outfile2 << annotations[currName] << endl << currAlign << endl << fixedSS << endl << BLAST << endl;
						
			} while (line2 != "" && !(infile2.eof()));
			//outfile2 << endl;
			for (int i = 0; i < numNames; i++)
			{	string currName = names[i];
				if (!(isMemberOfVector(linesUsed[currLine],currName)))
				{	SSs[currName] += dashes(templateSize);
					aligns[currName]  += dashes(templateSize);
					//cout << "adding " << templateSize << " dashes to " << currName << " on line " << currLine << ", curralignsize = " << aligns[currName].size() << endl;
				}
			}
		}
	}
	infile2.close();
			
	outfile2 << "Template:" << endl << Query << endl;
	for (int i = 0; i < numNames; i++)
	{	string currName = names[i];
		outfile2 << annotations[currName] << endl << aligns[currName] << endl << SSs[currName] << endl << endl;
	}
			
	outfile2.close();
			
	if (hitFound == FALSE)
	{	cout << "No hits found for " << currDomain << endl;
		return 0;
	}
	
	string outFileName3 = currDomain + ".log3";
	ofstream outfile3 ((char*)outFileName3.c_str());			

	string currPDB = currDomain;
	if (currPDB.find("/") != string::npos)
		currPDB = currPDB.substr(currPDB.find("/")+1);
	currPDB = toupper(currPDB.substr(0,4)) + ":" + currPDB.substr(5,1);
	
	int i;
	for (i = 0; i < numNames; i++)
	{
		string currNameI = names[i];
		if (currNameI == currPDB)
			break;
	}
			
	string currNameI = names[i];
	string currAnnotationsI = annotations[currNameI];
	string currAlignsI = aligns[currNameI];
	string SSI = SSs[currNameI];
	string BLASTI = BLASTs[currNameI];
	
	//outfile3 << "Template:" << endl;
	//outfile3 << currNameI << endl << replacePeriods(currAlignsI,Query) << endl << SSI << endl << BLASTI << endl << endl;
	
	bool firstHit = TRUE;
	for (int j = 0; j < numNames; j++)
	{
		if (j == i)
			continue;
		string currNameJ = names[j];
		string currAnnotationsJ = annotations[currNameJ];
		string currAlignsJ = aligns[currNameJ];
		string SSJ = SSs[currNameJ];
		string BLASTJ = BLASTs[currNameJ];
		
		int numMutations = mutations(currAlignsJ);
		
		if (numMutations <= numMutationsAllowed)
		{	//outfile3 << currNameJ << endl << replacePeriods(currAlignsJ,Query) << endl << SSJ << endl << BLASTJ << endl;
			pair<float,string> hammDistPair = hammingDistance(SSI,SSJ);
			float hammDist = hammDistPair.first;
			string diffString = hammDistPair.second;
			vector<string> juicy_bits;
			vector<int> juicy_spaces;
			vector<int> juicy_indices;
			vector<float> juicy_scores;
			vector<int> optimal_indices;
			string optimal_diffstring;
			int optimal_diffstring_idx = 0;
			float bit_score = -1.0;

			juicy_differences(diffString,juicy_bits,juicy_spaces,\
					  juicy_indices);

			//print_vector(juicy_bits,"Juicy strings: ");
			//print_vector(juicy_spaces,"Juicy spaces: ");
			//print_vector(juicy_indices,"Juicy indices: ");

			//If there are no "."s in diffString
			if (juicy_indices[0] == -999){
			  if (score_bit(juicy_bits[0]) >= minHammDiff*resize(juicy_bits[0])){
			    optimal_indices.push_back(0);
			    juicy_scores.push_back(score_bit(juicy_bits[0]));
			    optimal_diffstring = juicy_bits[0];
			  }
			  else continue;
			}
			//If diffString is a string of "."s, no structural change; continue.
			else if (juicy_indices[0] == -998){
			  continue;
			}

			else{
			  score_juicy_bits(juicy_bits, juicy_scores);

			  find_optimal_indices(juicy_bits,juicy_scores,juicy_spaces,optimal_indices);

			  sort(optimal_indices.begin(),optimal_indices.end());

			  optimal_string(optimal_diffstring,juicy_bits,juicy_spaces,optimal_indices);

			  if (resize2(optimal_diffstring) < minGraftSize) continue;

			  optimal_diffstring_idx = optimal_idx(optimal_indices, juicy_bits, juicy_spaces); 

			  //cout << "Optimal string: "<<optimal_diffstring<<endl;

			  //print_vector(optimal_indices,"optimal indices: ");
			}

			
    

			int nonBlankResidues = nonBlank(SSI,SSJ);
			float hammDiff = (float)(hammDist/(float)(nonBlankResidues));
			//outfile3 << "hammDist = " << hammDist << ", hammDiff = " << hammDiff << ", nonBlank = " << nonBlankResidues << endl << diffString << endl << endl;
			pair<float,int> maxDiffPair = findMaximalLocalDifference(diffString);
			float maxDiff = maxDiffPair.first;
			int maxDiffIdx = maxDiffPair.second;
			
			if (maxDiff >= minHammDiff)
			{	if (firstHit)
				{
					outfile3 << "Template:" << endl << currNameI << endl << currAnnotationsI << endl << replacePeriods(currAlignsI,Query) << endl << SSI << endl << BLASTI << endl << endl;
					firstHit = FALSE;
				}
				outfile3 << currNameJ << endl << currAnnotationsJ << endl << replacePeriods(currAlignsJ,Query) << endl << SSJ << endl << BLASTJ << endl;
				//outfile3 << "hammDist = " << hammDist << ", hammDiff = " << hammDiff << ", nonBlank = " << nonBlankResidues << ", numMuts = " << numMutations << ", maxDiff = " << maxDiff << ", maxString = " << diffString.substr(maxDiffIdx,minGraftSize) << ", maxSeq = " << Query.substr(maxDiffIdx,minGraftSize) << endl << diffString << endl << endl;
				//cout <<"ODI "<<optimal_diffstring_idx<<" ODS "<< optimal_diffstring.size()<<" QS " \
				    //<<Query.size() <<endl;
				outfile3 << "hammDist = " << hammDist << ", hammDiff = " << hammDiff << ", nonBlank = " << nonBlankResidues << ", numMuts = " << numMutations << ", maxDiff = " << maxDiff << ", maxString = " << optimal_diffstring << ", maxSeq = " << Query.substr(optimal_diffstring_idx,optimal_diffstring.size()) << endl << diffString << endl << endl;
			}
		}
	}
}

#endif
