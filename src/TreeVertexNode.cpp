#include "TreeVertexNode.h"


TreeVertexNode::TreeVertexNode()
{
	level = -1;
}


TreeVertexNode::~TreeVertexNode()
{
	ClearedgeList();
}


/* �����������ڽӱ��м���һ����
����vΪ��������ڽӶ�����ţ�wΪ����ߵ�Ȩֵ*/
bool TreeVertexNode::AppendTreeEdge(int v, float w)
{
	TreeEdgeNode *p = edgeList;
	TreeEdgeNode *q = nullptr;
	// �ҵ����ӱ���ĩ�ڵ㣬ĩ�ڵ��ָ�븳ֵ��q�����������һ���ڵ��adjVex��ֵ��v��ͬ���򷵻�false
	while (p != nullptr)
	{
		if (p->adjTreeVertex == v)
			return false;
		q = p;
		p = p->next;
	}
	// ���ڽӱ��������һ����
	p = new TreeEdgeNode(v, w);
	if (q == nullptr)
		edgeList = p;
	else
		q->next = p;
	return true;
}


/*���һ��������ڽӱ��ͷ��ڽӱ���ڵ�*/
void TreeVertexNode::ClearedgeList()
{
	TreeEdgeNode *p, *q;
	p = edgeList;
	while (p != nullptr)
	{
		q = p->next;
		delete p;
		p = q;
	}
	edgeList = nullptr;
}