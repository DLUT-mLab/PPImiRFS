#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "DataGlobal.h"
#include "Tree.h"
#include "WPin.h"

using namespace std;

class Compute
{
public:
	Compute();
	~Compute();
	void computeMrnaSetFs(WPin &, float &, vector<string> &, vector<string> &);                        //��������miRNA�л��򼯼�Ĺ���������
	void computeMrnaFs(const string, const string, WPin &, float &);                //�������������mRNA��Ĺ���������
	void computeTarNset(vector<string> &, vector<string> &, WPin &);                 //��������miRNA�л���Ľ���
private:
	vector<string> tar;                     //������miRNA�л��򼯵Ĳ���
	vector<string> tarn;                    //������miRNA�л��򼯵Ľ���
};

