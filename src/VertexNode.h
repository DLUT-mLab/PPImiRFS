#pragma once
#include <string>
#include "EdgeNode.h"

using namespace std;

class VertexNode
{
public:
	VertexNode();
	~VertexNode();
	void ClearedgeList();           // ɾ�����������ڽӱ�
	bool Appendedge(int, float = 0);  // �����������ڽӱ��м���һ����
	bool Removeedge(int);           // �����������ڽӱ���ɾ��һ����
	string data;                    // ���������
	EdgeNode *edgeList;             // �����ڽӱ���׽ڵ�ָ��
};

