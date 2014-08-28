#pragma once
#include <string>
#include <queue>
#include <vector>
#include "EdgeNode.h"
#include "VertexNode.h"
#include "Tree.h"

using namespace std;

class WPin
{
public:
	~WPin();
	WPin(int);
	void ComputeFs(const int, const int, double &fs);
	bool Insertedge(const string &, const string &, float wgh);
	bool InsertVertex(const string &);
	int LocateVertex(const string &);
	bool SetVertex(int v, const string &);
private:
	VertexNode *vertices;                           // �����
	int numVertices;                                            // �������
	int numedges;                                                // ����
	int maxVertices;                                            // ���ɴ�ŵĶ������
};