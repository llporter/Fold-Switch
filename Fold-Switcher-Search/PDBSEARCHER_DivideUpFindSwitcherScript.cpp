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
	string inFileName;
	int size;	

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-in")
			inFileName = argv[++i];
		else if (curWord == "-size")
			size = atoi(argv[++i]);
	}

	string command = "wc -l " + inFileName + " > temp";
	system((char*)command.c_str());
	
	ifstream tempfile ("temp");
	int N_Lines;
	tempfile >> N_Lines;
	tempfile.close();
	command = "rm -f temp";
	system((char*)command.c_str());
	
	int N_Scripts = N_Lines/size;
	int Remainder = N_Lines - N_Scripts * size;
	
	if (N_Scripts == 0)
	{
		int start = 1;
		int stop = N_Lines;
		string scriptName = "script1.com";
		command = "sed -n " + itos(start) + "," + itos(stop) + "p " + inFileName + " > " + scriptName;
		system((char*)command.c_str());
	}
	else
	{
		for (int scriptIdx = 1; scriptIdx <= N_Scripts; scriptIdx++)
		{
			int start = (scriptIdx-1)*size + 1;
			int stop  = scriptIdx*size;
			string scriptName = "script" + itos(scriptIdx) + ".com";
			
			command = "sed -n " + itos(start) + "," + itos(stop) + "p " + inFileName + " > " + scriptName;
			system((char*)command.c_str());
		}
		int start = N_Scripts * size + 1;
		int stop = N_Scripts * size + Remainder;
		string scriptName = "script" + itos(N_Scripts+1) + ".com";
		command = "sed -n " + itos(start) + "," + itos(stop) + "p " + inFileName + " > " + scriptName;
		system((char*)command.c_str());
	}
}

#endif
