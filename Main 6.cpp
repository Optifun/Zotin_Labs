#define _CRT_SECURE_NO_WARNINGS
#include "../../Zotin_Labs/BMPFileRW.h"
#include "../../Zotin_Labs/special_points.h"
#include <locale.h>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;
string doubleToString(double number)
{
	char* buf = new char[256];
	int size = sprintf(buf, "%f", number);
	string num = string(buf, size);
	delete[] buf;
	return num;
}

float Round(float number, int digits) {
	long long temp = floor(number*powf(10, digits));
	return temp / pow(10, digits);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Hello World!\n";
	RGBQUAD **RGB, **RGBresult;
	BITMAPFILEHEADER head;
	BITMAPINFOHEADER info;
	string str3 = "444.bmp";
	string str2 = "44.bmp";
	string str = "4.bmp";
	BMPRead(RGB, head, info, str.c_str());
	RGBresult = new RGBQUAD*[info.biHeight];
	for (int i = 0; i < info.biHeight; i++)
		RGBresult[i] = new RGBQUAD[info.biWidth];
	BMPRead(RGB, head, info, str.c_str());
	SpecialPointsPosled(RGB, info.biHeight, info.biWidth, 10, RGBresult);	
	BMPWrite(RGBresult, head, info, str2.c_str());
	//SpecialPointsParalOmp(RGB, info.biHeight, info.biWidth, 10, RGBresult);
	//BMPWrite(RGBresult, head, info, str3.c_str());

	//Файл для результата
	std::ofstream out;
	out.open("ResultData.csv", ios::out);

	//Количество видов функций:
	int FunctArrSize = 1;
	//Суммарное количество вложенных функций
	int FunctArrArrSize[] { 5 };
	//Задание массива указателей (первая последовательная, остальные - параллельные реализации)
	void*** FunctArr = new void**[FunctArrSize];
	FunctArr[0] = new void*[FunctArrArrSize[0]]{ SpecialPointsPosled , SpecialPointsParalOmp, SpecialPointsParalCilkFor, SpecialPointsParalOmpCilk, SpecialPointsParalCilkForIndex};
	//Задание массива имён функций
	string ** FunctNameArr = new string*[FunctArrSize];
	FunctNameArr[0] = new string[FunctArrArrSize[0]]{ "Особые точки последов." , "Особ. точки парал. OMP" , "Особ. точки парал. CILKFOR" , "Особ. точки парал. OMPCILK" , "Особ. точки парал. CILKINDEX"};
	//Количество наборов данных
	int NDFunctionSetSize = 4;

	//Тут необходимо задать параметры функции
	typedef double*(*typeOfFunct)(RGBQUAD** &RGB, int height, int width, double threshold, RGBQUAD** &RGBresult);
	typeOfFunct Xxx;
	//Временные переменные
	double* Times;

	double tmpTime, posledTime, paralTime, time_Start, summ1;
	for (int k = 0; k < FunctArrSize; k++)
	{
		out << "Название;V;S;Time1;Time2;Time3;Time4;Time5\n";
		cout << "Название\tV\tS\tTime1\tTime2\tTime3\tTime4\tTime5\n";

		for (int j = 0; j < NDFunctionSetSize; j++)
		{
			//Действия установки НД
			str = to_string(j + 1) + ".bmp";
			BMPRead(RGB, head, info, str.c_str());
			//
			out << FunctNameArr[k][0].c_str() << " НД" << j + 1 << ";";
			cout << FunctNameArr[k][0].c_str() << " НД" << j + 1 << "\t";
			//
			Xxx = (typeOfFunct)FunctArr[k][0];
			double summ = 0;
			for (int sr = 0; sr < 1; sr++)
			{
				time_Start = omp_get_wtime();
				//Тут необходимо задать параметры функции
				Times = Xxx(RGB, info.biHeight, info.biWidth, 10, RGBresult);
				summ += omp_get_wtime() - time_Start;
			}
			posledTime = summ / 1;
			out << " " << Round(posledTime, 3) << ";" << Round(Times[0], 5) << ";" << Times[1] << ";" << Round(Times[0], 5) << ";" << Round(Times[0], 5) << ";" << Round(Times[0], 5) << ";";
			cout << Round(posledTime, 3) << "\t" << Round(Times[0], 5) << "\t" << Times[1] << "\t" << Round(Times[0], 5) << "\t" << Round(Times[0], 5) << "\t" << Round(Times[0], 5) << "\t";
			out << "\n";
			cout << "\n";
			for (int i = 1; i < FunctArrArrSize[k]; i++) {
				
				for (int g = 2; g < 5; g++)
				{
					out << FunctNameArr[k][i].c_str() << " Thr: " << g << ";";
					cout << FunctNameArr[k][i].c_str() << " Thr: " << g << "\t";
					__cilkrts_end_cilk();
					__cilkrts_set_param("nworkers", doubleToString(g).c_str());
					int hjk = 0;
					cilk_for(int i = 0; i < 20; i++)
						hjk++;
					omp_set_num_threads(g);
					Xxx = (typeOfFunct)FunctArr[k][i];
					summ = 0;
					for (int sr = 0; sr < 1; sr++)
					{
						time_Start = omp_get_wtime();
						//Тут необходимо задать параметры функции
						Times = Xxx(RGB, info.biHeight, info.biWidth, 10, RGBresult);
						summ += omp_get_wtime() - time_Start;
					}
					paralTime = summ / 1;
					out << " " << Round(paralTime, 3) << ";" << " " << Round(posledTime / paralTime, 3) << ";" << Round(Times[0], 5) << ";" << Round(Times[1], 5) << ";" << Round(Times[2], 5) << ";" << Round(Times[3], 5) << ";" << Round(Times[4], 5) << ";";
					cout << Round(paralTime, 3) << "\t" << Round(posledTime / paralTime, 3) << "\t" << Round(Times[0], 5) << "\t" << Round(Times[1], 5) << "\t" << Round(Times[2], 5) << "\t" << Round(Times[3], 5) << "\t" << Round(Times[4], 5) << "\t";
					out << "\n";
					cout << "\n";
				}
				
			}
			out << "\n";
			cout << "\n";
		}
	}
	out.close();
	system("PAUSE");
}