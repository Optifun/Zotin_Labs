#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_opmul.h>
#include <stdio.h>
#include <omp.h>
#include <string>
#include <list>

using namespace std;

typedef struct HarrisPoint {
	double R;
	int Y;
	int X;
} HarrisPoint;

#pragma region Posled
void ShellSortPosled(HarrisPoint* Arr, int Size) {
	for (int step = Size / 2; step > 0; step /= 2)
		for (int i = step; i < Size; i++)
			for (int j = i - step; j >= 0 && Arr[j].R < Arr[j + step].R; j -= step)
			{
				HarrisPoint tmp = Arr[j];
				Arr[j] = Arr[j + step];
				Arr[j + step] = tmp;
			}
}

//Формирование матрицы коэффициентов для фильтрации Гаусса
double** GetGaussMatrixPosled(int RH, int RW, double q) {
	double** Result = new double*[RH * 2 + 1];
	for (int i = 0; i < RH * 2 + 1; i++)
		Result[i] = new double[RW * 2 + 1];
	double SUM = 0;
	for (int Y = -RH; Y <= RH; Y++)
		for (int X = -RW; X <= RW; X++) {
			double CF = (1 / (2 * 3.14159265358979323846 * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));
			Result[Y + RH][X + RW] = CF;
			SUM += CF;
		}
	for (int Y = -RH; Y <= RH; Y++)
		for (int X = -RW; X <= RW; X++)
			Result[Y + RH][X + RW] /= SUM;
	return Result;
}
void LineFilteringGaussPosled(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrixPosled(RH, RW, RW / 3.0); //Сигма тут
	RGBresult = new RGBQUAD*[height];
	for (int Y = 0; Y < height; Y++)
	{
		RGBresult[Y] = new RGBQUAD[width];
		for (int X = 0; X < width; X++)
		{
			double rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -RH; DY <= RH; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -RW; DX <= RW; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					double tmp = CoefMatrix[DY + RH][DX + RW];
					rgbBlue += RGB[KY][KX].rgbBlue * tmp;
					rgbGreen += RGB[KY][KX].rgbGreen * tmp;
					rgbRed += RGB[KY][KX].rgbRed * tmp;
				}
			}
			if (rgbBlue < 0)	rgbBlue = 0;
			if (rgbBlue > 255)	rgbBlue = 255;
			if (rgbGreen < 0)	rgbGreen = 0;
			if (rgbGreen > 255)	rgbGreen = 255;
			if (rgbRed < 0)		rgbRed = 0;
			if (rgbRed > 255)	rgbRed = 255;
			RGBresult[Y][X].rgbBlue = rgbBlue;
			RGBresult[Y][X].rgbGreen = rgbGreen;
			RGBresult[Y][X].rgbRed = rgbRed;
		}
	}
	for (int i = 0; i < RW; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

void TestAfterDiffPosled(double** &Grad, int height, int width) {
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			if (Grad[Y][X] > 255) Grad[Y][X] = 255.0;
			if (Grad[Y][X] < 0) Grad[Y][X] = 0.0;
		}
}
void DiffPosled(double** &I, int height, int width, double** &Iresult, double** &CoefMatrixX) {
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			double Summ = 0;;
			for (int DY = Y - 1; DY <= Y + 1; DY++)
			{
				int KY = DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = X - 1; DX <= X + 1; DX++)
				{
					int KX = DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					Summ += I[KY][KX] * CoefMatrixX[DY - Y + 1][DX - X + 1];
				}
			}
			Iresult[Y][X] = Summ;
		}
}

void SelectPixel(RGBQUAD** &RGB, int height, int width, int Y, int X) {
	for (int DY = Y - 1; DY <= Y + 1; DY++)
	{
		int KY = DY;
		if (KY < 0)
			KY = 0;
		if (KY > height - 1)
			KY = height - 1;
		for (int DX = X - 1; DX <= X + 1; DX++)
		{
			int KX = DX;
			if (KX < 0)
				KX = 0;
			if (KX > width - 1)
				KX = width - 1;
			if (KX == X && KY == Y)
				continue;
			else {
				RGB[KY][KX].rgbBlue = 0;
				RGB[KY][KX].rgbRed = 255;
				RGB[KY][KX].rgbGreen = 0;
			}
		}
	}
}

double* SpecialPointsPosled(RGBQUAD** &RGB, int height, int width, double threshold, RGBQUAD** &RGBresult) {
	//Копирование RGB в RGBResult
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			RGBresult[Y][X].rgbBlue = RGB[Y][X].rgbBlue;
			RGBresult[Y][X].rgbGreen = RGB[Y][X].rgbGreen;
			RGBresult[Y][X].rgbRed = RGB[Y][X].rgbRed;
			RGBresult[Y][X].rgbReserved = RGB[Y][X].rgbReserved;
		}

	double* Times = new double[5];
	double TmpTimes;
////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Фильтрация методом Гаусса изображение
	RGBQUAD** RGBFiltered = new RGBQUAD*[height];
	for (int i = 0; i < height; i++)
		RGBFiltered[i] = new RGBQUAD[width];
	LineFilteringGaussPosled(RGB, height, width, 5, 5, RGBFiltered);

	//Массив яркости	
	double** I = new double*[height];
	for (int Y = 0; Y < height; Y++) {
		I[Y] = new double[width];
		for (int X = 0; X < width; X++)
			I[Y][X] = RGBFiltered[Y][X].rgbRed * 0.299 + RGBFiltered[Y][X].rgbGreen * 0.587 + RGBFiltered[Y][X].rgbBlue * 0.144;
	}
	Times[0] = omp_get_wtime() - TmpTimes;
////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Объявление вспомогательных массивов
	double** DiffX1 = new double*[height];
	for (int Y = 0; Y < height; Y++)	
		DiffX1[Y] = new double[width];

	double** DiffY1 = new double*[height];
	for (int Y = 0; Y < height; Y++)
		DiffY1[Y] = new double[width];

	double** DiffXY = new double*[height];
	for (int Y = 0; Y < height; Y++)
		DiffXY[Y] = new double[width];

	double** DiffX = new double*[height];
	for (int Y = 0; Y < height; Y++)
		DiffX[Y] = new double[width];

	double** DiffY = new double*[height];
	for (int Y = 0; Y < height; Y++)
		DiffY[Y] = new double[width];

	//Формирование матриц коэффициентов для поиска градиентов
	double** CoefX = new double*[3];
	CoefX[0] = new double[3] { 1, 0, -1 };
	CoefX[1] = new double[3] { 1, 0, -1 };
	CoefX[2] = new double[3] { 1, 0, -1 };
	double** CoefY = new double*[3];
	CoefY[0] = new double[3] { 1, 1, 1 };
	CoefY[1] = new double[3] { 0, 0, 0 }; 
	CoefY[2] = new double[3] { -1, -1, -1 };
		
	//Вычисление градиентов яркости (свёртка)
	DiffPosled(I, height, width, DiffX1, CoefX);
	DiffPosled(I, height, width, DiffY1, CoefY);
	
	DiffPosled(DiffX1, height, width, DiffX, CoefX);
	DiffPosled(DiffX1, height, width, DiffXY, CoefY);
	DiffPosled(DiffY1, height, width, DiffY, CoefY);
	
	TestAfterDiffPosled(DiffX, height, width);
	TestAfterDiffPosled(DiffY, height, width);
	TestAfterDiffPosled(DiffXY, height, width);

	Times[1] = omp_get_wtime() - TmpTimes;
////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Формирование массива отклика угла
	double max = 0, min = 0;
	long Counter = 0;
	double** R = new double*[height];
	for (int Y = 0; Y < height; Y++) {
		R[Y] = new double[width];
		for (int X = 0; X < width; X++) {
			double A = DiffX[Y][X];
			double B = DiffY[Y][X];
			double C = DiffXY[Y][X];
			double tmp = (A * B - C * C) - 0.04 * ((A + B) * (A + B));
			if (tmp > 0) {
				if (tmp > max)
					max = tmp;
				R[Y][X] = tmp;
				Counter++;
			}
			else R[Y][X] = 0;
		}
	}

	//Очистка лишнего
	for (int i = 0; i < height; i++)
		delete[] DiffX1[i];
	delete[] DiffX1;
	for (int i = 0; i < height; i++)
		delete[] DiffY1[i];
	delete[] DiffY1;
	for (int i = 0; i < 3; i++)
		delete[] CoefY[i];
	delete[] CoefY;
	for (int i = 0; i < 3; i++)
		delete[] CoefX[i];
	delete[] CoefX;
	for (int i = 0; i < height; i++)
		delete[] I[i];
	delete[] I;
	for (int i = 0; i < height; i++)
		delete[] DiffY[i];
	delete[] DiffY;
	for (int i = 0; i < height; i++)
		delete[] DiffX[i];
	delete[] DiffX;
	for (int i = 0; i < height; i++)
		delete[] DiffXY[i];
	delete[] DiffXY;
	for (int i = 0; i < height; i++)
		delete[] RGBFiltered[i];
	delete[] RGBFiltered;

	Times[2] = omp_get_wtime() - TmpTimes;
////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Составление списка точек Харриса
	max *= ((double)threshold / 100);
	HarrisPoint* HarrisPoints = new HarrisPoint[Counter];
	Counter = 0;
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
			if (R[Y][X] > max) {
				HarrisPoints[Counter].Y = Y;
				HarrisPoints[Counter].X = X;
				HarrisPoints[Counter].R = R[Y][X];
				Counter++;
			}
	
	//Сортировка
	ShellSortPosled(HarrisPoints, Counter);
	
	//Отсечение по порогу точек
	list<HarrisPoint> ToPrint;
	for (int i = 0; i < Counter; i++) {
		bool Dobavit = true;
		for each (HarrisPoint Elem in ToPrint)
		{
			if (pow(HarrisPoints[i].Y - Elem.Y, 2) + pow(HarrisPoints[i].X - Elem.X, 2) < 20 * 20) {
				Dobavit = false;
				break;
			}
		}
		if (Dobavit) ToPrint.push_back(HarrisPoints[i]);
	}
	
	Times[3] = omp_get_wtime() - TmpTimes;
////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Визуализация особоых точек
	for each (HarrisPoint Elem in ToPrint)
		SelectPixel(RGBresult, height, width, Elem.Y, Elem.X);
	Times[4] = omp_get_wtime() - TmpTimes;
	//Очистка лишнего
	for (int i = 0; i < height; i++)
		delete[] R[i];
	delete[] R;
	delete[] HarrisPoints;
////////////////////////////////////////////////////////////
	return Times;
}
#pragma endregion



#pragma region RapalOmp
void ShellSortParalOmp(HarrisPoint* Arr, int Size) {
	for (int step = Size / 2; step > 0; step /= 2)
#pragma omp parallel for shared(Arr) firstprivate(Size, step) schedule(guided)
		for (int i = step; i < Size; i++)
			for (int j = i - step; j >= 0 && Arr[j].R < Arr[j + step].R; j -= step)
			{
				HarrisPoint tmp = Arr[j];
				Arr[j] = Arr[j + step];
				Arr[j + step] = tmp;
			}
}

//Формирование матрицы коэффициентов для фильтрации Гаусса параллельное
double** GetGaussMatrixParalOmp(int RH, int RW, double q) {
	double** Result = new double*[RH * 2 + 1];
	for (int i = 0; i < RH * 2 + 1; i++)
		Result[i] = new double[RW * 2 + 1];
	double SUM = 0;
#pragma omp parallel for reduction(+:SUM) firstprivate(RH, RW, q) shared(Result) schedule(static, RH)
	for (int Y = -RH; Y <= RH; Y++)
		for (int X = -RW; X <= RW; X++) {
			double CF = (1 / (2 * 3.14159265358979323846 * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));
			Result[Y + RH][X + RW] = CF;
			SUM += CF;
		}
#pragma omp parallel for firstprivate(RH, RW, q, SUM) shared(Result) schedule(static, RH)
	for (int Y = -RH; Y <= RH; Y++)
		for (int X = -RW; X <= RW; X++) {
			Result[Y + RH][X + RW] /= SUM;
		}
	return Result;
}
void LineFilteringGaussParalOmp(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrixParalOmp(RH, RW, RW / 3.0); //Сигма тут
	RGBresult = new RGBQUAD*[height];
#pragma omp parallel for firstprivate(RH, RW, height, width) shared(RGB, RGBresult) schedule(static, RH * 2 + 1)
	for (int Y = 0; Y < height; Y++)
	{
		RGBresult[Y] = new RGBQUAD[width];
		for (int X = 0; X < width; X++)
		{
			double rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -RH; DY <= RH; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -RW; DX <= RW; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					double tmp = CoefMatrix[DY + RH][DX + RW];
					rgbBlue += RGB[KY][KX].rgbBlue * tmp;
					rgbGreen += RGB[KY][KX].rgbGreen * tmp;
					rgbRed += RGB[KY][KX].rgbRed * tmp;
				}
			}
			if (rgbBlue < 0)	rgbBlue = 0;
			if (rgbBlue > 255)	rgbBlue = 255;
			if (rgbGreen < 0)	rgbGreen = 0;
			if (rgbGreen > 255)	rgbGreen = 255;
			if (rgbRed < 0)		rgbRed = 0;
			if (rgbRed > 255)	rgbRed = 255;
			RGBresult[Y][X].rgbBlue = rgbBlue;
			RGBresult[Y][X].rgbGreen = rgbGreen;
			RGBresult[Y][X].rgbRed = rgbRed;
		}
	}
	for (int i = 0; i < RW; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

void TestAfterDiffParalOmp(double** &Grad, int height, int width) {
#pragma omp parallel for firstprivate(height, width) shared(Grad) schedule(guided)
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			if (Grad[Y][X] > 255) Grad[Y][X] = 255.0;
			if (Grad[Y][X] < 0) Grad[Y][X] = 0.0;
		}
}
void DiffParalOmp(double** &I, int height, int width, double** &Iresult, double** &CoefMatrixX) {
#pragma omp parallel for firstprivate(height, width) shared(I, Iresult, CoefMatrixX) schedule(static, 3)
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			double Summ = 0;;
			for (int DY = Y - 1; DY <= Y + 1; DY++)
			{
				int KY = DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = X - 1; DX <= X + 1; DX++)
				{
					int KX = DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					Summ += I[KY][KX] * CoefMatrixX[DY - Y + 1][DX - X + 1];
				}
			}
			Iresult[Y][X] = Summ;
		}
}

void SelectPixelParalOmp(RGBQUAD** &RGB, int height, int width, int Y, int X) {
#pragma omp parallel for firstprivate(height, width, Y, X) shared(RGB)
	for (int DY = Y - 1; DY <= Y + 1; DY++)
	{
		int KY = DY;
		if (KY < 0)
			KY = 0;
		if (KY > height - 1)
			KY = height - 1;
		for (int DX = X - 1; DX <= X + 1; DX++)
		{
			int KX = DX;
			if (KX < 0)
				KX = 0;
			if (KX > width - 1)
				KX = width - 1;
			if (KX == X && KY == Y)
				continue;
			else {
				RGB[KY][KX].rgbBlue = 0;
				RGB[KY][KX].rgbRed = 255;
				RGB[KY][KX].rgbGreen = 0;
			}
		}
	}
}

double* SpecialPointsParalOmp(RGBQUAD** &RGB, int height, int width, double threshold, RGBQUAD** &RGBresult) {
#pragma omp parallel for firstprivate(height, width) shared(RGB, RGBresult) schedule(guided)
	//Копирование RGB в RGBResult
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			RGBresult[Y][X].rgbBlue = RGB[Y][X].rgbBlue;
			RGBresult[Y][X].rgbGreen = RGB[Y][X].rgbGreen;
			RGBresult[Y][X].rgbRed = RGB[Y][X].rgbRed;
			RGBresult[Y][X].rgbReserved = RGB[Y][X].rgbReserved;
		}

	double* Times = new double[5];
	double TmpTimes;
	////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Фильтрация методом Гаусса изображение
	RGBQUAD** RGBFiltered = new RGBQUAD*[height];
	for (int i = 0; i < height; i++)
		RGBFiltered[i] = new RGBQUAD[width];
	LineFilteringGaussParalOmp(RGB, height, width, 5, 5, RGBFiltered);

	//Массив яркости	
	double** I = new double*[height];
	for (int Y = 0; Y < height; Y++) {
		I[Y] = new double[width];
		for (int X = 0; X < width; X++)
			I[Y][X] = RGBFiltered[Y][X].rgbRed * 0.299 + RGBFiltered[Y][X].rgbGreen * 0.587 + RGBFiltered[Y][X].rgbBlue * 0.144;
	}
	Times[0] = omp_get_wtime() - TmpTimes;
	////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Объявление вспомогательных массивов
	double** DiffX1;
	double** DiffY1;
	double** DiffXY;
	double** DiffX;
	double** DiffY;
#pragma omp parallel sections
	{
#pragma omp section 
		{
			DiffX1 = new double*[height];
			for (int Y = 0; Y < height; Y++)
				DiffX1[Y] = new double[width];
		}
#pragma omp section
		{
			DiffY1 = new double*[height];
			for (int Y = 0; Y < height; Y++)
				DiffY1[Y] = new double[width];
		}
#pragma omp section 
		{
			DiffXY = new double*[height];
			for (int Y = 0; Y < height; Y++)
				DiffXY[Y] = new double[width];
		}
#pragma omp section
		{
			DiffX = new double*[height];
			for (int Y = 0; Y < height; Y++)
				DiffX[Y] = new double[width];
		}
#pragma omp section
		{
			DiffY = new double*[height];
			for (int Y = 0; Y < height; Y++)
				DiffY[Y] = new double[width];
		}
	}


	//Формирование матриц коэффициентов для поиска градиентов
	double** CoefX = new double*[3];
	CoefX[0] = new double[3]{ 1, 0, -1 };
	CoefX[1] = new double[3]{ 1, 0, -1 };
	CoefX[2] = new double[3]{ 1, 0, -1 };
	double** CoefY = new double*[3];
	CoefY[0] = new double[3]{ 1, 1, 1 };
	CoefY[1] = new double[3]{ 0, 0, 0 };
	CoefY[2] = new double[3]{ -1, -1, -1 };

	//Вычисление градиентов яркости (свёртка)
	DiffParalOmp(I, height, width, DiffX1, CoefX);
	DiffParalOmp(I, height, width, DiffY1, CoefY);

	DiffParalOmp(DiffX1, height, width, DiffX, CoefX);
	DiffParalOmp(DiffX1, height, width, DiffXY, CoefY);
	DiffParalOmp(DiffY1, height, width, DiffY, CoefY);

	TestAfterDiffParalOmp(DiffX, height, width);
	TestAfterDiffParalOmp(DiffY, height, width);
	TestAfterDiffParalOmp(DiffXY, height, width);

	Times[1] = omp_get_wtime() - TmpTimes;
	////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Формирование массива отклика угла
	double max = 0, min = 0;
	long Counter = 0;
	double** R = new double*[height];
#pragma omp parallel for firstprivate(height, width) shared(DiffX, DiffY, DiffXY, R) schedule(guided)
	for (int Y = 0; Y < height; Y++) {
		R[Y] = new double[width];
		for (int X = 0; X < width; X++) {
			double A = DiffX[Y][X];
			double B = DiffY[Y][X];
			double C = DiffXY[Y][X];
			double tmp = (A * B - C * C) - 0.04 * ((A + B) * (A + B));
			if (tmp > 0) {
				if (tmp > max)
					max = tmp;
				R[Y][X] = tmp;
				Counter++;
			}
			else R[Y][X] = 0;
		}
	}

	//Очистка лишнего
#pragma omp parallel sections firstprivate(height)
	{
#pragma omp section 
		{
			for (int i = 0; i < height; i++)
				delete[] DiffX1[i];
			delete[] DiffX1;
			for (int i = 0; i < height; i++)
				delete[] DiffY1[i];
			delete[] DiffY1;
		}
#pragma omp section
		{
			for (int i = 0; i < 3; i++)
				delete[] CoefY[i];
			delete[] CoefY;
			for (int i = 0; i < 3; i++)
				delete[] CoefX[i];
			delete[] CoefX;
		}
#pragma omp section 
		{

			for (int i = 0; i < height; i++)
				delete[] I[i];
			delete[] I;
			for (int i = 0; i < height; i++)
				delete[] DiffY[i];
			delete[] DiffY;
		}
#pragma omp section
		{
			for (int i = 0; i < height; i++)
				delete[] DiffX[i];
			delete[] DiffX;
			for (int i = 0; i < height; i++)
				delete[] DiffXY[i];
			delete[] DiffXY;
			for (int i = 0; i < height; i++)
				delete[] RGBFiltered[i];
			delete[] RGBFiltered;
		}
	}

	Times[2] = omp_get_wtime() - TmpTimes;
	////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Составление списка точек Харриса
	max *= ((double)threshold / 100);
	HarrisPoint* HarrisPoints = new HarrisPoint[Counter];
	Counter = 0;

	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
			if (R[Y][X] > max) {
				HarrisPoints[Counter].Y = Y;
				HarrisPoints[Counter].X = X;
				HarrisPoints[Counter].R = R[Y][X];
				Counter += 1;
			}

	//Сортировка
	ShellSortParalOmp(HarrisPoints, Counter);

	//Отсечение по порогу точек
	list<HarrisPoint> ToPrint;
	for (int i = 0; i < Counter; i++) {
		bool Dobavit = true;
		for each (HarrisPoint Elem in ToPrint)
		{
			if (pow(HarrisPoints[i].Y - Elem.Y, 2) + pow(HarrisPoints[i].X - Elem.X, 2) < 20 * 20) {
				Dobavit = false;
				break;
			}
		}
		if (Dobavit) ToPrint.push_back(HarrisPoints[i]);
	}

	Times[3] = omp_get_wtime() - TmpTimes;
	////////////////////////////////////////////////////////////
	TmpTimes = omp_get_wtime();
	//Визуализация особоых точек
	for each (HarrisPoint Elem in ToPrint)
		SelectPixelParalOmp(RGBresult, height, width, Elem.Y, Elem.X);
	Times[4] = omp_get_wtime() - TmpTimes;
	//Очистка лишнего
	for (int i = 0; i < height; i++)
		delete[] R[i];
	delete[] R;
	delete[] HarrisPoints;
	////////////////////////////////////////////////////////////
	return Times;
}
#pragma endregion