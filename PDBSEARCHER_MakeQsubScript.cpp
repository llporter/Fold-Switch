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
	int size;	

	for (int i = 1; i < argc; i++)
	{	string curWord = argv[i];
		if (curWord == "-size")
			size = atoi(argv[++i]);
	}

	string command = "pwd > temp";
	system((char*)command.c_str());
	
	ifstream tempfile ("temp");
	string dir;
	tempfile >> dir;
	tempfile.close();
	command = "rm -f temp";
	system((char*)command.c_str());
	
	ofstream outfile ("qsubscript.com");
	
	for (int i = 1; i <= size; i++)
		outfile << "qsub -cwd -R y -j y -b y -V " + dir + "/script" + itos(i) + ".com" << endl;
	
	outfile.close();
}

#endif
