#define _CRT_SECURE_NO_WARNINGS

#include<locale.h>
#include<fstream>
#include<sstream>
#include<time.h>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#undef max
#undef min
#include "matrix_new.cpp"
#include "filtering.h"
#include "logger.h"
#include "texture_filtering.h"
#undef max
#undef min
#include<Windows.h>
using namespace std;

#define Random(a, b) (double)rand() / (RAND_MAX + 1)*((b)-(a)) + (a)
#define CALC_Times_Matrix 5
#define CALC_Times_Filter 1
int threads = 1;
long Matrix<int>::Destroed = 0;
long Matrix<int>::Created = 0;
bool Matrix<int>::DEBUG = false;

long Matrix<double>::Destroed = 0;
long Matrix<double>::Created = 0;
bool Matrix<double>::DEBUG = false;

double calcTime(double*& times, int count)
{
	double sred = 0;
	for (int i = 0; i < count; i++)
	{
		sred += times[i];
	}
	sred /= count;
	double otcl = 0;
	for (int i = 0; i < count; i++)
	{
		otcl += abs(sred - times[i]);
	}
	otcl /= count;
	double l = sred - otcl * 2;
	double r = sred + otcl * 2;
	double time = 0;
	int tcount = 0;
	for (int i = 0; i < count; i++)
		if (times[i] >= l && times[i] <= r)
		{
			time += times[i];
			tcount++;
		}
	return time / tcount;
}

double Speed(double& tsync, double& tparal)
{
	return tsync / tparal;
}

double Efficiency(double& speed, int &streams)
{
	return speed / streams;
}

//template<class T>
//void Test()
//{
//	Matrix<T> a = Matrix<T>(4, 4, [](long n, long m)->double {return 2; });
//	Matrix<T> b = Matrix<T>(a);
//	cout << a << endl;
//	cout << b << endl;
//	Method m = Method::Cilk_for_index;
//	cout << a.Sum(m) << endl;
//	cout << a.Add(b, m) << endl;
//	cout << a.Mul(m) << endl;
//	cout << a.Max(b, m) << endl;
//	//cout<< a.Add(b, Method::Cilk_spawn);
//}

namespace one
{
	template<class T>
	void Zadanie()
	{
		string opnames[4] = { "Summation ", "Addition", "Multiplication", "Maximum" };
		string menames[7] = { "Sequentional", "Omp For", "Omp Sections", "Cilk For", "Cilk Spawn", "Cilk Index", "Cilk For+Index" };
		double posled[4][4];//4 метода 4 набора данных
		string temp;
		stringstream str;
		double time, spd, eff;
		double *times = new double[CALC_Times_Matrix];
		int k = 0;
		bool hasData = false;
		Matrix<T> a,b;
		char buf[20];
		int c = 0;
		logg << "Method;AsyncMode;Data;Thread;Time;Speed\n";
		for (int op = 0; op < 4; op++)//операци€
		{
			for (int m = 0; m < 7; m++)//метод распараллеливани€
			{
				logg << opnames[op] << ";" << menames[m] << ";\n";
				for (int d = 1000, k=0; d <= 8000; d += 1750, k++)
				{
					hasData = false;
					for (threads = 2; threads < 5; threads++)
					{
						omp_set_num_threads(threads);
						__cilkrts_end_cilk();
						sprintf_s(buf, "%d", threads);
						__cilkrts_set_param("nworkers", buf);
						cilk_for(int i = 0; i < 20; i++)
							c += 0;
						//cout << endl << __cilkrts_get_nworkers() << endl;
						if (m == 0)
							logg << ";;" << d << ";" << 1 << ";";
						else
							logg << ";;" << d << ";" << threads << ";";
						for (int try_ = 0; try_ < CALC_Times_Matrix; try_++)
						{
							if (!hasData)
							{
								a = Matrix<T>::Ones(d, d);
								b = Matrix<T>::E(d, d);
								hasData = true;
							}
							times[try_] = -omp_get_wtime();
							switch (op)
							{
							case 0:
							{
								a.Sum((Method)m, threads);
								break;
							}
							case 1:
							{
								a.Add(b, (Method)m, threads);
								break;
							}
							case 2:
							{
								a.Mul((Method)m, threads);
								break;
							}
							case 3:
							{
								a.Max(b, (Method)m, threads);
								break;
							}
							default:
								break;
							}
							times[try_] += omp_get_wtime();
						}
						time = calcTime(times, CALC_Times_Matrix);
						logg << doubleToString(time) << ";";
						if (m == 0)
						{
							posled[op][k] = time;
							logg << "\n";
							break;
						}
						else
						{
							spd = Speed(posled[op][k], time);
							logg << doubleToString(spd) << ";\n";
						}
					}
				}
				logg.Flush();
			}
		}
	}
}

namespace two
{
	Median* medianM;
	Gauss* gaussM;
	void Initialize()
	{
		if (medianM == nullptr)
		{
			medianM = new Median[3];
			medianM[0] = medianFiltering;
			medianM[1] = medianFilteringAsyncSort;
			medianM[2] = medianFilteringCilkFor;

			gaussM = new Gauss[3];
			gaussM[0] = LineFilteringGauss;
			gaussM[1] = LineFilteringGaussParal;
			gaussM[2] = LineFilteringGaussCilk1Dx2;
		}
	}

	void Zadanie()
	{
		//путь до папки с картинками
		string dir = "C:\\Users\\Public\\Documents\\TestImage\\";
		logg = Log("Table2.csv", true, true);
		//названи€ картинок
		string fnames[3] = { "1", "2", "3"};
		string method_names[2] = { "Median ", "Gauss "};
		string parallel_method[3] = { "Sequentional", "Omp Async Sort", "Cilk For + Index"};
		double posled[2][3];//2 метода 3 набора данных
		string temp;
		stringstream str;
		double time, spd, eff;
		double *times = new double[CALC_Times_Filter];
		int k = 0, c=0;
		char buf[5];
		BITMAPFILEHEADER head;
		BITMAPINFOHEADER info;
		RGBQUAD **RGB, **result;
		Initialize();
		logg << "Method;Image;Window;Thread;Time;Speed;\n";
		for (int m = 0; m < 2; m++)//метод фильтрации
			for (int p = 0; p < 3; p++)// параллельный метод
			{
				logg << method_names[m] << parallel_method[p] << ";\n";
					for (int img = 0; img < 3; img++)// номер картинки
					{
						str = stringstream();
						str << dir << fnames[img] << ".bmp";
						BMPRead(RGB, head, info, str.str().c_str());
						result = new RGBQUAD*[info.biHeight];
						for (int i = 0; i < info.biHeight; i++)
						{
							result[i] = new RGBQUAD[info.biWidth];
						}

						for (int window = 3; window <= 11; window += 4)
						{

							for (threads = 2; threads < 5; threads++)
							{
								if (m == 0)
									logg << ";" << fnames[img] << ";" << window << ";" << 1 << ";";
								else
									logg << ";" << fnames[img] << ";" << window << ";" << threads << ";";
								omp_set_num_threads(threads);
								__cilkrts_end_cilk();
								sprintf_s(buf, "%d", threads);
								__cilkrts_set_param("nworkers", buf);
								cilk_for(int i = 0; i < 20; i++)
									c += 0;
								for (int try_ = 0; try_ < CALC_Times_Filter; try_++)
								{
									times[try_] = -omp_get_wtime();
									if (m == 0)
										medianM[p](RGB, info.biHeight, info.biWidth, window, window, result, ShellSort);
									else
										gaussM[p](RGB, info.biHeight, info.biWidth, window, window, result);
									times[try_] += omp_get_wtime();
								}
								time = calcTime(times, CALC_Times_Filter);
								logg << time << ";";
								if (threads == 2)
								{
									str = stringstream();
									str << dir << "lab5\\" << fnames[img] << "(" << method_names[m] << parallel_method[p] << ")[" << window << "].bmp";
									BMPWrite(result, head, info, str.str().c_str());
								}

								if (m == 0)
								{
									posled[m][img] = time;
									logg << "\n";
									break;
								}
								else
								{
									spd = Speed(posled[m][img], time);
									logg << spd << ";\n";
								}

							}

						}

						for (int i = 0; i < info.biHeight; i++)
						{
							delete[] result[i];
							delete[] RGB[i];
						}
						delete[] result;
						delete[] RGB;
					}
			}
	}

	void Test()
	{
		string dir = "C:\\Users\\Public\\Documents\\TestImage\\";
		Initialize();
		stringstream str;
		BITMAPFILEHEADER head;
		BITMAPINFOHEADER info;
		RGBQUAD **RGB, **result;
		str = stringstream();
		str << dir << "1" << ".bmp";
		BMPRead(RGB, head, info, str.str().c_str());
		result = new RGBQUAD*[info.biHeight];
		for (int i = 0; i < info.biHeight; i++)
			result[i] = new RGBQUAD[info.biWidth];
		gaussM[2](RGB, info.biHeight, info.biWidth, 8, 8, result);
		str = stringstream();
		str << dir << "lab5\\" << "1"<< ".bmp";
		BMPWrite(result, head, info, str.str().c_str());
	}
}

namespace three
{
	txFilter* methods;

	void Initialize()
	{
		if (methods == nullptr)
		{
			methods = new txFilter[5];
			methods[0] = textureFilter;
			methods[1] = textureFilterOmpOutside;
			methods[2] = textureFilterCilkOutside;
			methods[3] = textureFilterOmpOutsideVec;
			methods[4] = textureFilterCilkOutsideVec;
		}
	}

	void Zadanie()
	{
		string dir = "C:\\Users\\Public\\Documents\\TestImage\\";
		logg = Log("Table3.csv", true, true);
		string fnames[2] = { "3", "4" };
		string method_names[5] = { "Texture Sequence", "Textute Omp For", "Texture Cilk For" ,"Texture Omp For + Index", "Texture Cilk For + Index"};
		double posled[3][3];//1 метод 3 набора данных 3 режима
		string temp;
		double time, spd;
		int k = 0;
		int c = 0;
		BITMAPFILEHEADER head;
		BITMAPINFOHEADER info;
		RGBQUAD **RGB;
		char buf[5];
		stringstream str;
		float **M, **U, **R, **E; // матрицы статистик
		Initialize();
		logg << "Method;Image;Window;Thread;Time;Speed;\n";
		for (int m = 0; m < 5; m++)//метод фильтрации
		{
			logg << method_names[m] << ";\n";
				for (int img = 0; img < 2; img++)// номер картинки
				{
					str = stringstream();
					str << dir << fnames[img] << ".bmp";
					BMPRead(RGB, head, info, str.str().c_str());
					M = new float*[info.biHeight];
					U = new float*[info.biHeight];
					R = new float*[info.biHeight];
					E = new float*[info.biHeight];
					for (int i = 0; i<info.biHeight; i++)
					{
						M[i] = new float[info.biWidth]();
						U[i] = new float[info.biWidth]();
						R[i] = new float[info.biWidth]();
						E[i] = new float[info.biWidth]();
					}
					for (int window = 2, k=0; window <= 10; window += 4, k++)
					{
						for (threads = 2; threads < 5; threads++)
						{
							if (m == 0)
								logg << ";" << fnames[img] << ";" << window << ";" << 1 << ";";
							else
								logg << ";" << fnames[img] << ";" << window << ";" << threads << ";";
							omp_set_num_threads(threads);
							__cilkrts_end_cilk();
							sprintf_s(buf, "%d", threads);
							__cilkrts_set_param("nworkers", buf);
							cilk_for(int i = 0; i < 20; i++)
								c += 0;

							time = -omp_get_wtime();
							methods[m](RGB, info.biHeight, info.biWidth, window, window, M, U, R, E);
							time += omp_get_wtime();

							logg << time << ";";
							if (m == 0)
							{
								posled[img][k] = time;
								logg << "\n";
								break;
							}
							else
							{
								spd = Speed(posled[img][k], time);
								logg << spd << ";\n";
							}

						}
						str = stringstream();
						//папка с файлами \ текстуры \ название файла \ название метода \ размер окна \ характеристика
						str << dir << "texture\\" << fnames[img] <<"-M-"<< "(" << method_names[m] << ")["<<window<<"].bmp";
						formImage(head, info, M, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
						str = stringstream();
						str << dir << "texture\\" << fnames[img] << "-U-" << "(" << method_names[m] << ")[" << window << "].bmp";
						formImage(head, info, U, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
						str = stringstream();
						str << dir << "texture\\" << fnames[img] << "-R-" << "(" << method_names[m] << ")[" << window << "].bmp";
						formImage(head, info, R, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
						str = stringstream();
						str << dir << "texture\\" << fnames[img] << "-E-" << "(" << method_names[m] << ")[" << window << "].bmp";
						formImage(head, info, E, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);

					}
					for (int i = 0; i < info.biHeight; i++)
						delete[] RGB[i];
					delete[] RGB;
				}
		}
	}

	void Test()
	{
		BITMAPFILEHEADER head;
		BITMAPINFOHEADER info;
		RGBQUAD **RGB, **result;
		string dir = "C:\\Users\\Public\\Documents\\TestImage\\";
		stringstream str;
		str = stringstream();
		str << dir << "1" << ".bmp";
		BMPRead(RGB, head, info, str.str().c_str());
		float **M, **U, **R, **E; // матрицы статистик
		M = new float*[info.biHeight];
		U = new float*[info.biHeight];
		R = new float*[info.biHeight];
		E = new float*[info.biHeight];
		for (int i=0;i<info.biHeight;i++)
		{
				M[i] = new float[info.biWidth]();
				U[i] = new float[info.biWidth]();
				R[i] = new float[info.biWidth]();
				E[i] = new float[info.biWidth]();
		}
		textureFilterCilkOutside(RGB, info.biHeight, info.biWidth, 6, 6, M, U, R, E);

		str = stringstream();
		str << dir <<"texture\\" <<"M.bmp";
		formImage(head, info, M, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
		str = stringstream();
		str << dir << "texture\\" << "U.bmp";
		formImage(head, info, U, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
		str = stringstream();
		str << dir << "texture\\" << "R.bmp";
		formImage(head, info, R, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
		str = stringstream();
		str << dir << "texture\\" << "E.bmp";
		formImage(head, info, E, info.biHeight, info.biWidth, str.str().c_str(), 0.2, 0.8);
	}
}


void menu()
{
	int choise = -1;
	while (choise != 0)
	{
		cout << "[0]¬ыход \n[1]«адание 1\n[2]«адание 2\n[3]«адание 3\n";
		cin >> choise;
		switch (choise)
		{
		case 0:
		{
			return;
		}
		case 1:
		{
			logg = Log("Table1.csv", true, true);
			one::Zadanie<int>();
			one::Zadanie<double>();
			logg.Close();
			system("pause");
			break;
		}
		case 2:
		{
			two::Zadanie();
			system("pause");
			break;
		}
		case 3:
		{
			three::Zadanie();
			system("pause");
			break;
		}
		default:
			break;
		}
		system("cls");
	}
}


int main()
{
	setlocale(0, "");
	//omp_set_num_threads(4);
	//__cilkrts_set_param("nworkers", "4");
	menu();
	//three::Test();
	//two::Test();
	//Test<int>();
	//two::Test();
	//Test<int>();
	system("pause");
	return 0;
}