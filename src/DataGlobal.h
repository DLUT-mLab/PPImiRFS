#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "WPin.h"

using namespace std;

extern ofstream olog;

class DataGlobal
{
public:
	DataGlobal();
	~DataGlobal();
	void InputWppin(const string filePin, WPin &wpin);          //����wppin����
	void GetMirTar(const string, vector<string> &);             //��ȡmiRNA�İл�����Ϣ
};

