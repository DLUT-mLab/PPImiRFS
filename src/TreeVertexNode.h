#pragma once
#include "TreeEdgeNode.h"

class TreeVertexNode
{
public:
	TreeVertexNode();
	~TreeVertexNode();
	bool AppendTreeEdge(int, float = 0.0);  // �����������ڽӱ��м���һ����
	void ClearedgeList();           //���һ��������ڽӱ��ͷ��ڽӱ���ڵ�
	int level;                    // ���������еĲ���
	int search = 0;
	double dAccumFs = 0.0;         //����ö���Ķ���·�����ۻ�Ȩ�����ֵ
	TreeEdgeNode *edgeList = nullptr;             // �����ڽӱ���׽ڵ�ָ��
};

