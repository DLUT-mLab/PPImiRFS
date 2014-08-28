#include "VertexNode.h"


VertexNode::VertexNode()
{
}


VertexNode::~VertexNode()
{
//	ClearedgeList();
}


/*���һ��������ڽӱ��ͷ��ڽӱ���ڵ�*/
void VertexNode::ClearedgeList()
{
	EdgeNode *p, *q;
	p = edgeList;
	while (p != nullptr)
	{
		q = p->next;
		delete p;
		p = q;
	}
	edgeList = nullptr;
}


/* �����������ڽӱ��м���һ����
����vΪ��������ڽӶ�����ţ�wghΪ����ߵ�Ȩֵ*/
bool VertexNode::Appendedge(int v, float wgh)
{
	EdgeNode *p = edgeList;
	EdgeNode *q = nullptr;
	// �ҵ����ӱ���ĩ�ڵ㣬ĩ�ڵ��ָ�븳ֵ��q�����������һ���ڵ��adjVex��ֵ��v��ͬ���򷵻�false
	while (p != nullptr)
	{
		if (p->adjVertex == v)
			return false;
		q = p;
		p = p->next;
	}
	// ���ڽӱ��������һ����
	p = new EdgeNode(v, wgh);
	if (q == nullptr)
		edgeList = p;
	else
		q->next = p;
	return true;
}


/* �ڶ�����ڽӱ���ɾ��һ����
����vΪҪɾ���ߵ��ڽӶ������*/
bool VertexNode::Removeedge(int v)
{
	EdgeNode *p = edgeList;
	EdgeNode *q = nullptr;
	// �����ڽӱ�������������ߣ���ɾ��
	while (p != nullptr)
	{
		if (p->adjVertex == v)
		{
			if (p == edgeList)
				edgeList = p->next;
			else
				q->next = p->next;
			delete p;
			return true;
		}
		q = p;
		p = p->next;
	}
	return false;
}
