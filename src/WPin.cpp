#include "WPin.h"


/* ���캯��
����sΪ��ͼ���ɴ�ŵĶ������*/
WPin::WPin(int s) : numVertices(0), numedges(0), maxVertices(s)
{
	if (s > 0)
		vertices = new VertexNode[s];
}


/* �������� */
WPin::~WPin()
{
	for (int i = 0; i < numVertices; i++)
		vertices[i].ClearedgeList();
	if (maxVertices != 0)
		delete[] vertices;
}


 /*���ݶ�������ݣ��ҵ������ڶ�����е����
 ����vexΪ���������
 ����ֵΪ�������š��緵��-1��δ�ҵ���صĶ���*/
int WPin::LocateVertex(const string & vex)
{
	for (int i = 0; i < numVertices; i++)
	if (vertices[i].data == vex)
		return i;
	return -1;
}


 /*���ݶ������ţ����ö��������
 ����vΪ������ţ�����vex�������õĶ�������*/
bool WPin::SetVertex(int v, const string &vex)
{
	if (v < 0 || v >= numVertices)
		return false;
	vertices[v].data = vex;
	return true;
}


 /*��ͼ�в���һ������
 ����vexΪҪ����Ķ�������*/
bool WPin::InsertVertex(const string &vex)
{
	// ������ڵĶ�������Ѵ����ֵ���򷵻�
	if (numVertices == maxVertices)
		return false;
	// �ж���ͬ�Ķ����Ƿ���ڣ�����ڣ��򷵻�
	for (int i = 0; i < numVertices; i++)
	if (vertices[i].data == vex)
		return false;
	// ����һ��������
	vertices[numVertices].data = vex;
	vertices[numVertices++].edgeList = nullptr;
	return true;
}


 /*��ͼ�в���һ����
 ����vex1��vex2ΪҪ����ߵ����������ֵ������wghΪ������ߵ�Ȩֵ��Ĭ��ֵΪ0*/
bool WPin::Insertedge(const string &vex1, const string &vex2, float wgh)
{
	// �ҵ����������ڶ�����е���ţ��ֱ�ֵ��v1��v2
	// ����������ֻҪ��һ����ͼ�Ķ������δ�ҵ����򷵻�
	int v1 = LocateVertex(vex1);
	int v2 = LocateVertex(vex2);
	if (v1 == -1 || v2 == -1)
		return false;
	// Ϊ��һ��������ڽӱ�������һ����
	bool a = vertices[v1].Appendedge(v2, wgh);
	// ����ͼ�����������һ������ڽӱ�������һ����
	bool b = vertices[v2].Appendedge(v1, wgh);
	if (a && b)
	{
		numedges++;
		return true;
	}
	else
		return false;
}


//���øĽ��Ĺ������������������������л����Ĺ���������
void WPin::ComputeFs(const int root, const int goal, double &fs)
{

	Tree tree(numVertices);
	queue<int> path_queue;

	int flag = 0;
	int find_level = -1;
	double fs_temp = 1.0;

	if (goal != -1 && root != -1)
	{
		path_queue.push(root);

		tree.SetTreeVertex(root, 0, 1.0);

		while (!path_queue.empty())
		{
			int parent = path_queue.front();
			path_queue.pop();

			int parent_level = tree.GetLevel(parent);
			int search_flag = tree.GetSearchFlag(parent);
			
			if (1 == search_flag || parent_level == find_level)
				continue;
			tree.SetTreeSearchFlag(parent);

			fs_temp = tree.GetAccumFs(parent);

			EdgeNode *p = vertices[parent].edgeList;
			while (p != nullptr)
			{
				int child_level = tree.GetLevel(p->adjVertex);

				if (child_level == -1 || parent_level + 1 == child_level)
				{
					if (p->adjVertex == goal)
					{
						double AccumFs = fs_temp * (p->weight);

						if (flag == 0)
						{
							flag = 1;
							find_level = parent_level + 1;
						}
						if (-1 != child_level)
						{
							double temp = tree.GetAccumFs(p->adjVertex);
							AccumFs = temp >= AccumFs ? temp : AccumFs;
						}

						tree.SetTreeVertex(p->adjVertex, parent_level + 1, AccumFs);
						break;
					}
					else
				    {
						double AccumFs = fs_temp * (p->weight);
						if (-1 != child_level)
						{
							double temp = tree.GetAccumFs(p->adjVertex);
							AccumFs = temp >= AccumFs ? temp : AccumFs;
						}

						path_queue.push(p->adjVertex);
						tree.SetTreeVertex(p->adjVertex, parent_level + 1, AccumFs);
					}
				}
				p = p->next;
			}
		}
	}
	if (1 == flag)
		fs = tree.GetAccumFs(goal);
	else
		fs = 0.0;
}