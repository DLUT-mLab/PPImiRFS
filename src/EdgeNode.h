#pragma once
class EdgeNode
{
public:
	~EdgeNode();
	EdgeNode(int = -1, float = 0.0);       // ���캯��
	int adjVertex;                    // �ñ��е��ڽӶ����ڶ�����е����
	float weight;                       // �ñߵ�Ȩֵ
	EdgeNode *next;
};
