#pragma once
#include <stack>
#include "TreeEdgeNode.h"
#include "TreeVertexNode.h"

using namespace std;

class Tree
{
public:
	Tree(const int);
	~Tree();
	bool SetTreeVertex(const int, const int, const double fs);                     // ������ţ����ö���Ĳ���, ����ۻ�Ȩ��
	bool InsertTreeEdge(const int, const int, float = 0.0);       // ����һ����
	int GetFirstAdjTreeVex(int);                                    // ������ţ�ȡ�ö���ĵ�һ���ڽӶ�������
	int GetNextAdjTreeVex(int v, int w);                            // ������ţ�ȡ��v(�����w)����һ���ڽӶ������
	float GetTreeEdge(int, int);                            // ���ݶ�����ţ�ȡ��������֮��ıߵ�Ȩ��
	int GetLevel(const int);                             //���ݶ�����Ż�ȡ�ö���Ĳ���
	double GetAccumFs(const int);                         //���ݶ�����Ż�ȡ�ö�������ۻ�Ȩ��ֵ
	void ComputeMrFs(const int, const int, float &);                //�����������mRNA��Ĺ���������
	bool SetTreeSearchFlag(const int);
	int GetSearchFlag(const int);

private:
	TreeVertexNode *vertices;                                   // �����
	int maxVertices;                                            // ���ɴ�ŵĶ������ 
};

