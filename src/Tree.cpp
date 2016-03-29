#include "Tree.h"


Tree::Tree(const int v) :maxVertices(v)
{
	if (v > 0)
		vertices = new TreeVertexNode[v];
}


Tree::~Tree()
{
	delete[] vertices;
}


/*���ݶ������ţ����ö���Ĳ���
����vΪ������ţ�����level�������ö���Ĳ���*/
bool Tree::SetTreeVertex(const int v, const int level, const double fs)
{
	if (v < 0 || v >= maxVertices)
		return false;
	vertices[v].level = level;
	vertices[v].dAccumFs = fs;
	return true;
}


/*�����в���һ����
����v1Ϊ���׽ڵ㡢v2Ϊ���ӽڵ㣻����wΪ���׺ͺ��ӽڵ���Ȩֵ��Ĭ��ֵΪ0*/
bool Tree::InsertTreeEdge(const int v1, const int v2, float w)
{
	if (v1 == -1 || v2 == -1)
		return false;
	// Ϊ��һ��������ڽӱ�������һ����
	vertices[v1].AppendTreeEdge(v2, w);
	return true;
}


/* ȡ�����е�v������ĵ�һ�����Ӷ������*/
int Tree::GetFirstAdjTreeVex(int v)
{
	if (v < 0 || v >= maxVertices)
		return -1;
	if (vertices[v].edgeList == nullptr)
		return -1;
	else
		return vertices[v].edgeList->adjTreeVertex;
}


/* ������ţ�ȡ�ø��ڵ�v(����ں��ӽڵ�w)����һ���ڽӶ������*/
int Tree::GetNextAdjTreeVex(int v, int w)
{
	if (v < 0 || v >= maxVertices)
		return -1;
	TreeEdgeNode *p = vertices[v].edgeList;
	while (p != nullptr)
	{
		if (p->adjTreeVertex == w)
			break;
		p = p->next;
	}
	if (p == nullptr || p->next == nullptr)
		return -1;
	return p->next->adjTreeVertex;
}


/* ���ݶ�����ţ�ȡ��������֮��ıߵ�Ȩ��*/
float Tree::GetTreeEdge(int v, int w)
{
	float weight = 0.0;
	if (v < 0 || v >= maxVertices)
		return 0.0;
	if (w < 0 || w >= maxVertices)
		return 0.0;
	TreeEdgeNode *edge = vertices[v].edgeList;
	while (edge != nullptr)
	{
		if (edge->adjTreeVertex == w)
		{
			weight = edge->weight;
			break;
		}
		edge = edge->next;
	}
	return weight;
}


/*��ȡ�������*/
int Tree::GetLevel(const int v)
{
	if (v < 0 || v >= maxVertices)
		return false;
	int level = vertices[v].level;
	return level;
}


/*�����������mRNA��Ĺ���������*/

/*��Ȩֵ�۳˷���*/
//void Tree::ComputeMrFs(const int v, const int w, float &fs)
//{
//	float maxfs = 0.0;                                                        //����ͬ��·�����������������ֵ
//	float temp_fs = 1.0;                                                      //�ݴ浱ǰ·����Ȩֵ�˻�
//	stack<int> stack_path;                                                    //���������������ջ
//	stack_path.push(v);
//	int temp_parent = stack_path.top();                                       //�ݴ����Ƚڵ�
//	int temp_child = GetFirstAdjTreeVex(temp_parent);                         //�ݴ溢�ӽڵ�
//
//	/*�������Ӹ��ڵ㵽��������Ҷ�ӽڵ��·������������ջ��*/
//	while (temp_child != -1 && temp_child != w)
//	{
//		//��Ȩֵ�۳˷���
//		temp_fs *= GetTreeEdge(temp_parent, temp_child);
//		stack_path.push(temp_child);
//		temp_parent = stack_path.top();
//		temp_child = GetFirstAdjTreeVex(temp_parent);
//	}
//
//	/*����������Ҷ�ӽڵ���Ŀ��ڵ㣬���ʼ��maxfsֵ*/
//	if (temp_child == w)
//	{
//		//��Ȩֵ�۳˷���
//		maxfs = temp_fs * GetTreeEdge(temp_parent,temp_child);
//	}
//
//	int end = 0;                                                              //������������־λ
//
//	/*�����������ͳ�Ƴ����а���Ŀ��ڵ�·����Ȩ�س˻�������ֵ*/
//	while (!stack_path.empty())
//	{
//		int temp_brother = stack_path.top();                                  //��ǰջβ�ڵ�����������������꣬��ջβ�ڵ��ջ����Ϊ�ֵܽڵ�
//		stack_path.pop();
//		if (stack_path.empty())
//			break;
//		temp_parent = stack_path.top();                                       //��ǰջβ�ڵ�Ϊ�����������������Ƚڵ�
//		temp_child = GetNextAdjTreeVex(temp_parent, temp_brother);            //��ǰ���Ƚڵ�Ľ��������ֵܽڵ����һ�����ӽڵ�
//
//		//��Ȩֵ�۳˷���
//		temp_fs /= GetTreeEdge(temp_parent, temp_brother);                    //��ǰջ�б����·����Ȩ�س˻�
//
//		/*�����һ���������˵�ǰ���Ƚڵ�����������������ϻ��ݣ�ֱ���ҵ�һ��û����������Ƚڵ�*/
//		while (temp_child == -1)
//		{
//			temp_brother = stack_path.top();
//			stack_path.pop();
//
//			/*���������������������Ƴ�ѭ��*/
//			if (stack_path.empty())
//			{
//				end = 1;
//				break;
//			}
//			temp_parent = stack_path.top();
//			temp_child = GetNextAdjTreeVex(temp_parent, temp_brother);
//
//			//��Ȩֵ�۳˷���
//			temp_fs /= GetTreeEdge(temp_parent, temp_brother);
//		}
//
//		if (end == 1)
//			break;
//
//		/*��ǰ���Ƚڵ������������û�б����꣬�������������*/
//		while (temp_child != -1 && temp_child != w)
//		{
//			//��Ȩֵ�۳˷���
//			temp_fs *= GetTreeEdge(temp_parent, temp_child);
//
//			stack_path.push(temp_child);
//			temp_parent = stack_path.top();
//			temp_child = GetFirstAdjTreeVex(temp_parent);
//		}
//
//		/*�ҵ�Ŀ��ڵ㣬��������·����Ŀ��ڵ���ۻ�Ȩ��ֵ�Ƿ����maxfs�������ڣ������maxfsΪ��ǰֵ�����򲻸���*/
//		if (temp_child == w)
//		{
//			//��Ȩֵ�۳˷���
//			maxfs = temp_fs * GetTreeEdge(temp_parent, temp_child) > maxfs ? temp_fs*GetTreeEdge(temp_parent, temp_child) : maxfs;
//		}
//	}
//
//	fs = maxfs;
//}


/*��Ȩֵ�ۼ�ƽ������*/
void Tree::ComputeMrFs(const int v, const int w, float &fs)
{
	float maxfs = 0.0;                                                        //����ͬ��·�����������������ֵ
	float temp_fs = 0.0;                                                      //�ݴ浱ǰ·����Ȩֵ�ۼ�ƽ��ֵ
	stack<int> stack_path;                                                    //���������������ջ
	stack_path.push(v);
	int temp_parent = stack_path.top();                                       //�ݴ����Ƚڵ�
	int temp_child = GetFirstAdjTreeVex(temp_parent);                         //�ݴ溢�ӽڵ�

	//��Ȩֵ�ۼ�ƽ������
	int level = 0;

	/*�������Ӹ��ڵ㵽��������Ҷ�ӽڵ��·������������ջ��*/
	while (temp_child != -1 && temp_child != w)
	{
		//��Ȩֵ�ۼ�ƽ������
		level = GetLevel(temp_child);
		//temp_fs += GetTreeEdge(temp_parent, temp_child) * pow(0.5, level - 1);
		temp_fs += GetTreeEdge(temp_parent, temp_child);

		stack_path.push(temp_child);
		temp_parent = stack_path.top();
		temp_child = GetFirstAdjTreeVex(temp_parent);
	}

	/*����������Ҷ�ӽڵ���Ŀ��ڵ㣬���ʼ��maxfsֵ*/
	if (temp_child == w)
	{
		//��Ȩֵ�ۼ�ƽ������
		level = GetLevel(temp_child);
		//maxfs = temp_fs + GetTreeEdge(temp_parent, temp_child) * pow(0.5, level - 1);
		maxfs = temp_fs + GetTreeEdge(temp_parent, temp_child) > maxfs ? temp_fs + GetTreeEdge(temp_parent, temp_child) : maxfs;
	}

	int end = 0;                                                              //������������־λ

	/*�����������ͳ�Ƴ����а���Ŀ��ڵ�·����Ȩ�س˻�������ֵ*/
	while (!stack_path.empty())
	{
		int temp_brother = stack_path.top();                                  //��ǰջβ�ڵ�����������������꣬��ջβ�ڵ��ջ����Ϊ�ֵܽڵ�
		stack_path.pop();
		if (stack_path.empty())
			break;
		temp_parent = stack_path.top();                                       //��ǰջβ�ڵ�Ϊ�����������������Ƚڵ�
		temp_child = GetNextAdjTreeVex(temp_parent, temp_brother);            //��ǰ���Ƚڵ�Ľ��������ֵܽڵ����һ�����ӽڵ�

		//��Ȩֵ�ۼ�ƽ������
		level = GetLevel(temp_child);
		//temp_fs -= GetTreeEdge(temp_parent, temp_brother)*pow(0.5, level - 1);
		temp_fs -= GetTreeEdge(temp_parent, temp_brother);

		/*�����һ���������˵�ǰ���Ƚڵ�����������������ϻ��ݣ�ֱ���ҵ�һ��û����������Ƚڵ�*/
		while (temp_child == -1)
		{
			temp_brother = stack_path.top();
			stack_path.pop();

			/*���������������������Ƴ�ѭ��*/
			if (stack_path.empty())
			{
				end = 1;
				break;
			}
			temp_parent = stack_path.top();
			temp_child = GetNextAdjTreeVex(temp_parent, temp_brother);

			//��Ȩֵ�ۼ�ƽ������
			level = GetLevel(temp_child);
			//temp_fs -= GetTreeEdge(temp_parent, temp_brother)*pow(0.5, level - 1);
			temp_fs -= GetTreeEdge(temp_parent, temp_brother);
		}

		if (end == 1)
			break;

		/*��ǰ���Ƚڵ������������û�б����꣬�������������*/
		while (temp_child != -1 && temp_child != w)
		{
			//��Ȩֵ�ۼ�ƽ������
			level = GetLevel(temp_child);
			//temp_fs += GetTreeEdge(temp_parent, temp_brother)*pow(0.5, level - 1);
			temp_fs += GetTreeEdge(temp_parent, temp_child);

			stack_path.push(temp_child);
			temp_parent = stack_path.top();
			temp_child = GetFirstAdjTreeVex(temp_parent);
		}

		/*�ҵ�Ŀ��ڵ㣬��������·����Ŀ��ڵ���ۻ�Ȩ��ֵ�Ƿ����maxfs�������ڣ������maxfsΪ��ǰֵ�����򲻸���*/
		if (temp_child == w)
		{

			//��Ȩֵ�ۼ�ƽ������
			level = GetLevel(temp_child);
			//maxfs = temp_fs + GetTreeEdge(temp_parent, temp_child) * pow(0.5, level - 1) > maxfs ? temp_fs*GetTreeEdge(temp_parent, temp_child) : maxfs;
			maxfs = temp_fs + GetTreeEdge(temp_parent, temp_child) > maxfs ? temp_fs + GetTreeEdge(temp_parent, temp_child) : maxfs;
		}
	}

	fs = maxfs / GetLevel(w);
}


double Tree::GetAccumFs(const int ver)
{
	if (ver < 0 || ver >= maxVertices)
		return false;
	double fs = vertices[ver].dAccumFs;
	return fs;
}


bool Tree::SetTreeSearchFlag(const int ver)
{
	if (ver < 0 || ver >= maxVertices)
		return false;
	vertices[ver].search = 1;
	return true;
}


int Tree::GetSearchFlag(const int ver)
{
	if (ver < 0 || ver >= maxVertices)
		return false;
	int flag = vertices[ver].search;
	return flag;
}
