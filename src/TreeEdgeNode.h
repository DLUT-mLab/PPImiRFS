#pragma once
class TreeEdgeNode
{
public:
	~TreeEdgeNode();
	TreeEdgeNode(int = -1, float = 0.0);       // ���캯��
	int adjTreeVertex;                    // �ñ��е��ڽӶ����ڶ�����е����
	float weight;                       // �ñߵ�Ȩֵ
	TreeEdgeNode *next;
};

