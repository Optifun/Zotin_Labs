#pragma once
#include"filtering.h"
#include<list>
#include<string>
using namespace std;

//Реализации формирования гистограмм
#pragma region formHist

// Формирует гистограмму на основе карты яркости в рамке с радиусами RH, RW на позиции (x,y)
float* formHist(BYTE** BrMap, int height, int width, int x, int y, int RH, int RW)
{
	int index = 0;
	int coordX;
	int coordY;
	//инициализирую нулями
	float* hist = new float[256]();
	//прохожу по рамке
	for (int Y = -RH; Y <= RH; Y++)
	{
		coordY = y + Y;
		for (int X = -RW; X <= RW; X++)
		{
			coordX = x + X;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= width)
				coordX = width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= height)
				coordY = height - 1;
			//инкрементирует элемент гистограммы, соответствующий яркости текущего пикселя
			hist[BrMap[coordY][coordX]] += 1;

		}
	}

	//определю количество пикселей в рамке
	int size = (RH * 2 + 1)*(RW * 2 + 1);
	//нормирую гистограмму
	for (int i = 0; i < 256; i++)
		hist[i] /= size;
	return hist;
}

// Формирует гистограмму на основе карты яркости в рамке с радиусами RH, RW на позиции (x,y)
//Параллельный вариант с использованием циклов из ОМП
float* formHistOMP(BYTE** BrMap, int height, int width, int x, int y, int RH, int RW)
{
	int index = 0;
	int coordX;
	int coordY;
	//инициализирую нулями
	float* hist = new float[256]();
	//прохожу по рамке
	#pragma omp parallel for shared(BrMap, hist) firstprivate(x, y, width, height, RH, RW) schedule(dynamic, 3)
	for (int Y = -RH; Y <= RH; Y++)
	{
		coordY = y + Y;
		for (int X = -RW; X <= RW; X++)
		{
			coordX = x + X;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= width)
				coordX = width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= height)
				coordY = height - 1;
			//инкрементирует элемент гистограммы, соответствующий яркости текущего пикселя
			#pragma omp critical
			{
			hist[BrMap[coordY][coordX]] += 1;
			}

		}
	}

	//определю количество пикселей в рамке
	int size = (RH * 2 + 1)*(RW * 2 + 1);
	//нормирую гистограмму
	for (int i = 0; i < 256; i++)
		hist[i] /= size;
	return hist;
}

#pragma endregion

//Реализации получения метрик
#pragma region getMetrics

//Создаёт метрики по текущей гистограмме
void getMetrics(float &m2, float &u, float &r, float &e, float* &hist)
{
	float m = 0;
	for (int i = 0; i < 256; i++)
		m += hist[i] * i;

	for (int i = 0; i < 256; i++)
	{
		m2 += pow((i - m), 2)*hist[i];
		e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
		u += pow(hist[i], 2);
	}
	r = 1 - (1 / (1 + m2));
	e *= -1;
}

//Создаёт метрики по текущей гистограмме
//Параллельный вариант с использованием циклов и секций из ОМП
void getMetricsOmp(float &m2, float &u, float &r, float &e, float* &hist)
{
	float m = 0;
	#pragma omp parallel for reduction(+:m) schedule(dynamic, 45)
	for (int i = 0; i < 256; i++)
		m += hist[i] * i;

	#pragma omp parallel sections reduction(+:m2, e, u)
	{
		#pragma omp section
		{
			for (int i = 0; i < 256; i++)
				m2 += pow((i - m), 2)*hist[i];
		}

		#pragma omp section
		{
			for (int i = 0; i < 256; i++)
				e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
		}

		#pragma omp section
		{
			for (int i = 0; i < 256; i++)
				u += pow(hist[i], 2);
		}
	}
	r = 1 - (1 / (1 + m2));
	e *= -1;
}

#pragma endregion

//Реализации текстурных признаков
#pragma region textureFilter

//На вход изображение и размер рамки
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
void textureFilter(Bitmap &image, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	M = new float*[image.height];
	U = new float*[image.height];
	R = new float*[image.height];
	E = new float*[image.height];
	BYTE **Brightness = new BYTE*[image.height];
	float *hist;//256
	for (int y = 0; y < image.height; y++)
	{
		M[y] = new float[image.width]();
		U[y] = new float[image.width]();
		R[y] = new float[image.width]();
		E[y] = new float[image.width]();
		Brightness[y] = new BYTE[image.width];
	}

	for (int y = 0; y < image.height; y++)
		for (int x = 0; x < image.width; x++)
			//перевожу цветную картинку в карут яркости
			Brightness[y][x] = image.map[y][x].rgbRed*0.299 + image.map[y][x].rgbGreen*0.587 + image.map[y][x].rgbBlue*0.114;

	for (int y = 0; y < image.height; y++)
		for (int x = 0; x < image.width; x++)
		{
			//получаю гистограмму для окна
			hist = formHist(Brightness, image.height, image.width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
			delete[] hist;
		}
	//возвращаю метрики наверх
}

//На вход изображение и размер рамки
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Параллельная реализация, использует распараллеливание операций при обработке рамки
void textureFilterOmpInside(Bitmap &image, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	M = new float*[image.height];
	U = new float*[image.height];
	R = new float*[image.height];
	E = new float*[image.height];
	BYTE **Brightness = new BYTE*[image.height];
	float *hist;//256
	#pragma omp parallel for shared(M, U, R, E, Brightness) schedule(dynamic, 50)
	for (int y = 0; y < image.height; y++)
	{
		M[y] = new float[image.width]();
		U[y] = new float[image.width]();
		R[y] = new float[image.width]();
		E[y] = new float[image.width]();
		Brightness[y] = new BYTE[image.width];
	}
	#pragma omp parallel for shared(Brightness, image) schedule(dynamic, 50)
	for (int y = 0; y < image.height; y++)
		for (int x = 0; x < image.width; x++)
			//перевожу цветную картинку в карут яркости
			Brightness[y][x] = image.map[y][x].rgbRed*0.299 + image.map[y][x].rgbGreen*0.587 + image.map[y][x].rgbBlue*0.114;

	for (int y = 0; y < image.height; y++)
		for (int x = 0; x < image.width; x++)
		{
			//получаю гистограмму для окна
			hist = formHistOMP(Brightness, image.height, image.width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetricsOmp(M[y][x], U[y][x], R[y][x], E[y][x], hist);
			delete[] hist;
		}
	//возвращаю метрики наверх
}

#pragma endregion

//Формирует изображение по массиву метрик T[,]
//fname - имя выходного файла
//min<t1<t2<max - пороговые значения
void formImage(BITMAPFILEHEADER head, BITMAPINFOHEADER info, float** T, int Height, int Width, string fname, float t1, float t2, float _min = 0, float _max = 8000)
{
	RGBQUAD** out = new RGBQUAD*[Height];
	for (int y = 0; y < Height; y++)
	{
		out[y] = new RGBQUAD[Width]();
		for (int x = 0; x < Width; x++)
		{
			//Если в первом диапазоне, то крашу в зеленый
			if (T[y][x] > t2 && T[y][x] < _max)
				out[y][x].rgbGreen = 255;

			//Если во втором диапазоне, то крашу в желтый
			if (T[y][x] > t1 &&T[y][x] < t2)
			{
				out[y][x].rgbRed = 255;
				out[y][x].rgbGreen = 255;
			}

			//Если в первом диапазоне, то крашу в красный
			if (T[y][x] > _min && T[y][x] < t1)
				out[y][x].rgbRed = 255;
		}
	}
	//Сохраняю в файл
	BMPWrite(out, head, info, fname.c_str());
}