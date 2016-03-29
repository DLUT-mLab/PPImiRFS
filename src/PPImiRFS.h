#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
// #include <direct.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "Compute.h"
#include "DataGlobal.h"
#include "WPin.h"

using namespace std;

const char* const prog_name = "PPImiRFS";
const char* const prog_version = "PPImiRFS version 1.0 (2014/8/27)";
const char* const log_file_name = "./result/log.txt";


//¥Ú”°∞Ô÷˙–≈œ¢
inline void usage(int exit_value = 0){
	ifstream ihelp("readme.txt");
	while (!ihelp.eof())
	{
		string str;
		getline(ihelp, str);
		cout << str << endl;
	}
	exit(exit_value);
}


string filewppin;
string fileMirFs;
string fileMirPredict;
ofstream olog(log_file_name);
