#include "PPImiRFS.h"

int main(int argc, char* argv[])
{
	/*���������в��������ݲ���ִ�к�������*/

	clock_t start;
	clock_t finish;

	start = clock();

	mkdir("result", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if (argc == 1)
	{
		fprintf(stderr, "Usage:%s [-i value] [-o value] [-p value] [-h]  arglist ... \nparameters are necessary!!!!", argv[0]);
		exit(1);
	}

	for (int i = 1; i < argc; ++i)
	{
		char *para = argv[i];

		if (para[0] == '-')
		{
			switch (para[1])
			{
			case 'i':           //��Ԥ��miRNA�ļ�
				fileMirPredict = argv[++i];
				break;
			case 'p':           //WPPIN�ļ�
				filewppin = argv[++i];
				break;
			case 'o':           //����ļ�
				fileMirFs = argv[++i];
				break;
			case 'v':           //�汾��Ϣ
				cout << prog_version << endl;
				return 0;
			case 'h':           //������Ϣ
				usage();
			default:            //����������ڴ��󣬴�ӡ������Ϣ
				usage();
				break;
			}
		}
		else
		{
			fprintf(stderr, "Usage:%s [-i value] [-o value] [-p value] [-h]  arglist ... ", argv[0]);
			exit(1);
		}
	}
	if (fileMirPredict != "" && filewppin != "" && fileMirFs != "")
	{
		clock_t start1 = clock();

		cout << "inputing datas from input files..." << endl;
		olog << "inputing datas from input files..." << endl;

        string perl_output = "./result/protein_count.txt";
        string perl_script = "./count_protein.pl";

		//ͳ�Ƶ���������
		char cmd[100];
		snprintf(cmd, 100, "%s %s %s %s", perl_script.c_str(), filewppin.c_str(), perl_output.c_str(), log_file_name);
		system(cmd);
		memset(cmd, 0, sizeof(cmd));

		//���뵰������Ŀ
		ifstream icprotein(perl_output);
		int count_protein = 0;
		icprotein >> count_protein;
		icprotein.close();

		//�����Ԥ��miRNA����
		ifstream mirpredict(fileMirPredict.c_str());

		//���ý���ļ���
		ofstream fMirFs(fileMirFs.c_str());

		DataGlobal dg;
		WPin wpin(count_protein);

		//����Ȩ�ص�������������
		snprintf(cmd, 100, "%s", filewppin.c_str());
		dg.InputWppin(cmd, wpin);
		memset(cmd, 0, sizeof(cmd));

		clock_t finish1 = clock();
		cout << "Running time is: " << static_cast<double>(finish1 - start1) / CLOCKS_PER_SEC << "s" << endl;
		olog << "Running time is: " << static_cast<double>(finish1 - start1) / CLOCKS_PER_SEC << "s" << endl;

		cout << endl;
		olog << endl;
		cout << "Computing..." << endl;
		olog << "Computing..." << endl;
		cout << endl;
		olog << endl;

		int predict_num = 0;

		while (!mirpredict.eof())
		{
			clock_t start3 = clock();
			cout << "computing the " << ++predict_num << "pair miRNAs' function similirity..." << endl;
			olog << "computing the " << predict_num << "pair miRNAs' function similirity..." << endl;

			string x;                  //miRNA-a����
			string y;                  //miRNA-b����

			/*��ȡ��Ԥ�������miRNA����*/

			cout << "extract miRNA name..." << endl;
			olog << "extract miRNA name..." << endl;

			mirpredict >> x >> y;
			if ("" == x || "" == y)
			{
				continue;
			}


			if (x == y)
			{
				fMirFs << x << "	" << y << "	" << 1.0 << endl;
				continue;
			}

			//��ȡmiRNA�л�����Ϣ

			cout << "extract miRNA-Target informations..." << endl;
			olog << "extract miRNA-Target informations..." << endl;

			vector<string> tara;            //miRNA-a��Ӧ�İл�����Ϣ������mRNA���֣���Ӧ�İ����У��ɽӽ���ֵ
			vector<string> tarb;            //miRNA-b��Ӧ�İл�����Ϣ������ͬ��
			dg.GetMirTar(x, tara);
			dg.GetMirTar(y, tarb);

			if (tara.size() == 0 || tarb.size() == 0)
			{
				fMirFs << x << "	" << y << "	" << 0.0 << endl;
				continue;
			}

			Compute com;

			com.computeTarNset(tara, tarb, wpin);

			//��������miRNA�л��򼯼�Ĺ���������

			cout << "computing the function similirity of the target sets of miRNAs..." << endl;
			olog << "computing the function similirity of the target sets of miRNAs..." << endl;

			float mrnafs = 0.0;
			com.computeMrnaSetFs(wpin, mrnafs, tara, tarb);

			fMirFs << x << "	" << y << "	" << mrnafs << endl;

			clock_t finish3 = clock();
			cout << "Running time is: " << static_cast<double>(finish3 - start3) / CLOCKS_PER_SEC << "s" << endl;
			cout << endl;
			olog << "Running time is: " << static_cast<double>(finish3 - start3) / CLOCKS_PER_SEC << "s" << endl;
			olog << endl;
		}

		fMirFs.close();
		mirpredict.close();

	}
	else
	{
		fprintf(stderr, "Usage:%s [-i value] [-r value] [-p value] [-h]  arglist ... \nparameters are enough!!!!", argv[0]);
		exit(1);
	}

	finish = clock();
	cout << "Whole Running time is: " << static_cast<double>(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	olog << "Whole Running time is: " << static_cast<double>(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
