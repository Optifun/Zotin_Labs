#pragma once
#include"filtering.h"
#include<list>
#include<string>
#include<cilk\reducer_vector.h>
#include<vector>
using namespace std;

typedef void(*txFilter)(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E);

//Реализации формирования гистограмм
#pragma region formHist

// Формирует гистограмму на основе карты яркости в рамке с радиусами RH, RW на позиции (x,y)
vector<float> formHist(BYTE** &BrMap, int height, int width, int x, int y, int RH, int RW)
{
	int index = 0;
	int coordX;
	int coordY;
	//инициализирую нулями
	vector<float> hist = vector<float>(256);
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
vector<float> formHistOMP(BYTE** &BrMap, int height, int width, int x, int y, int RH, int RW)
{
	int index = 0;
	int coordX;
	int coordY;
	//инициализирую нулями
	vector<float> hist = vector<float>(256);
	//прохожу по рамке
	#pragma omp parallel for shared(BrMap, hist) firstprivate(x, y, width, height, RH, RW) schedule(dynamic, 35)
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


// Формирует гистограмму на основе карты яркости в рамке с радиусами RH, RW на позиции (x,y)
vector<float> formHistCilkIndex(BYTE** &BrMap, int height, int width, int x, int y, int RH, int RW)
{
	int index = 0;
	int lenx = (RW * 2 + 1);
	int leny = (RH * 2 + 1);
	int temp;
	//инициализирую нулями
	vector<float> hist = vector<float>(256);
	int *KXs = new int[lenx];

	for (int X = -RW; X <= RW; X++)
	{
		temp = (x + X < 0) ? 0 : x + X;
		KXs[index] = (temp > width - 1) ? width - 1 : temp;
		index++;
	}
	index = 0;

	//прохожу по рамке
	for (int Y = -RH; Y <= RH; Y++)
	{
		temp = (y + Y < 0) ? 0 : y + Y;
		temp = (temp > height - 1) ? height - 1 : temp;
		hist[ BrMap[temp] [KXs[0:lenx]] ] += 1;
	
	}
	//инкрементирует элемент гистограммы, соответствующий яркости текущего пикселя

	//определю количество пикселей в рамке
	int size = lenx*leny;
	//нормирую гистограмму
	hist.data()[0:255] /= size;
	return hist;
}

#pragma endregion

//Реализации получения метрик
#pragma region getMetrics

//Создаёт метрики по текущей гистограмме
void getMetrics(float &m2, float &u, float &r, float &e, vector<float> &hist)
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
void getMetricsOmp(float &m2, float &u, float &r, float &e, vector<float> &hist)
{
	float m = 0;
	#pragma omp parallel for reduction(+:m) schedule(dynamic, 45)
	for (int i = 0; i < 256; i++)
		m += hist[i] * i;
	float step = 256 / 4;
	int mod = 256 % 4;
	int t1 = 0;
	int t2 = step * 1 + mod;
	int t3 = t2 + step;
	int t4 = t3 + step;
	int t5 = t4 + step;
	#pragma omp parallel sections reduction(+:m2, e, u) shared(t1, t2, t3, t4, t5)
	{
		#pragma omp section
		{
			for (int i = t1; i < t2; i++)
			{
				m2 += pow((i - m), 2)*hist[i];
				e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
				u += pow(hist[i], 2);
			}
		}

		#pragma omp section
		{
			for (int i = t2; i < t3; i++)
			{
				m2 += pow((i - m), 2)*hist[i];
				e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
				u += pow(hist[i], 2);
			}
		}

		#pragma omp section
		{
			for (int i = t3; i < t4; i++)
			{
				m2 += pow((i - m), 2)*hist[i];
				e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
				u += pow(hist[i], 2);
			}
		}

		#pragma omp section
		{
			for (int i = t4; i < t5; i++)
			{
				m2 += pow((i - m), 2)*hist[i];
				e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
				u += pow(hist[i], 2);
			}
		}
	}
	r = 1 - (1 / (1 + m2));
	e *= -1;
}

void getMetricsCilk(float &m2, float &u, float &r, float &e, vector<float> &hist)
{
	cilk::reducer<cilk::op_add<float>> m(0);
	cilk::reducer<cilk::op_add<float>> _e(0);
	cilk::reducer<cilk::op_add<float>> _u(0);
	cilk::reducer<cilk::op_add<float>> _m2(0);
	cilk_for(int i = 0; i < 256; i++)
		*m += hist[i] * i;
	float moment = m.get_value();
	cilk_for(int i = 0; i < 256; i++)
	{
		*_m2 += pow((i - moment), 2)*hist[i];
		*_e += (hist[i] != 0) ? hist[i] * log2(hist[i]) : 0;
		*_u += pow(hist[i], 2);
	}
	m2 = _m2.get_value();
	r = 1 - (1 / (1 + m2));
	u = _u.get_value();
	e = _e.get_value() * -1;
}

#pragma endregion

//Реализации текстурных признаков
#pragma region textureFilter

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
void textureFilter(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	vector<float> hist;//256
	BYTE **Brightness = new BYTE*[height]();

	for (int y = 0; y < height; y++)
	{
		Brightness[y] = new BYTE[width]();
		//перевожу цветную картинку в карту яркости
		for (int x = 0; x < width; x++)
			Brightness[y][x] = image[y][x].rgbRed*0.299 + image[y][x].rgbGreen*0.587 + image[y][x].rgbBlue*0.114;
	}

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			//получаю гистограмму для окна
			hist = formHist(Brightness, height, width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
	//возвращаю метрики наверх
}

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Параллельная реализация, использует распараллеливание операций при обработке рамки
void textureFilterOmpInside(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	BYTE **Brightness = new BYTE*[height];
	vector<float> hist;//256

	#pragma omp parallel for shared(Brightness, image) schedule(dynamic, 50)
	for (int y = 0; y < height; y++)
	{
		Brightness[y] = new BYTE[width];

		for (int x = 0; x < width; x++)
			//перевожу цветную картинку в карут яркости
			Brightness[y][x] = image[y][x].rgbRed*0.299 + image[y][x].rgbGreen*0.587 + image[y][x].rgbBlue*0.114;
	}

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			//получаю гистограмму для окна
			hist = formHistOMP(Brightness, height, width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetricsOmp(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
	//возвращаю метрики наверх
}

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Параллельная реализация, использует распараллеливание внешнего цикла на Omp For
void textureFilterOmpOutside(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	BYTE **Brightness = new BYTE*[height];

	#pragma omp parallel for shared(Brightness, image) schedule(dynamic, 50)
	for (int y = 0; y < height; y++)
	{
		Brightness[y] = new BYTE[width];

		for (int x = 0; x < width; x++)
			//перевожу цветную картинку в карут яркости
			Brightness[y][x] = image[y][x].rgbRed*0.299 + image[y][x].rgbGreen*0.587 + image[y][x].rgbBlue*0.114;
	}

	#pragma omp parallel for shared(Brightness, image, M, U, R, E) schedule(dynamic, 50)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			//получаю гистограмму для окна
			vector<float> hist = formHist(Brightness, height, width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
	//возвращаю метрики наверх
}

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Параллельная реализация, использует распараллеливание внешнего цикла на Omp For
void textureFilterOmpOutsideVec(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	BYTE **Brightness = new BYTE*[height];

#pragma omp parallel for shared(Brightness, image) schedule(dynamic, 50)
	for (int y = 0; y < height; y++)
	{
		Brightness[y] = new BYTE[width];

		for (int x = 0; x < width; x++)
			//перевожу цветную картинку в карут яркости
			Brightness[y][x] = image[y][x].rgbRed*0.299 + image[y][x].rgbGreen*0.587 + image[y][x].rgbBlue*0.114;
	}

#pragma omp parallel for shared(Brightness, image, M, U, R, E) schedule(dynamic, 50)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			//получаю гистограмму для окна
			vector<float> hist = formHistCilkIndex(Brightness, height, width, x, y, rh, rw);
			//получаю метрики для текущего положения окна
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
	//возвращаю метрики наверх
}

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Паралелльный вариант с внешним распараллеливанием cilk_for
void textureFilterCilkOutside(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	BYTE **Brightness = new BYTE*[height];
	cilk_for(int y = 0; y < height; y++)
	{
		Brightness[y] = new BYTE[width];
		//перевожу цветную картинку в карут яркости
		for (int x = 0; x < width; x++)
			Brightness[y][x] = image[y][x].rgbRed*0.299 + image[y][x].rgbGreen*0.587 + image[y][x].rgbBlue*0.114;
	}

	cilk_for(int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			vector<float> hist = formHist(Brightness, height, width, x, y, rh, rw);
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
}

//На вход изображение и размер рамки
//Указатели M,U,R,E должны быть инициализированы
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
//Паралелльный вариант с внешним распараллеливанием cilk_for и cilk_index
void textureFilterCilkOutsideVec(RGBQUAD **&image, int height, int width, int rh, int rw, float **&M, float **&U, float **&R, float **&E)
{
	BYTE **Brightness = new BYTE*[height];
	cilk_for(int y = 0; y < height; y++)
	{
	Brightness[y] = new BYTE[width];
	Brightness[y][0:width] = image[y][0:width].rgbRed*0.299 + image[y][0:width].rgbGreen*0.587 + image[y][0:width].rgbBlue*0.114;
	}

	cilk_for(int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			vector<float> hist = formHistCilkIndex(Brightness, height, width, x, y, rh, rw);
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
		}
}

#pragma endregion

//Оценивает карту метрик, возвращает минимальное и максимальное значение признаков
//map - массив признаков
//Height, Width - высота и ширина карты
void getInterval(float** &map, int Height, int Width, float& _min, float& _max)
{
	_min = map[0][0];
	_max = map[0][0];
	for (int y = 0; y < Height; y++)
		for (int x = 0; x < Width; x++)
		{
			_max = max(_max, map[y][x]);
			_min = min(_min, map[y][x]);
		}
}

//Функция формирующая картинку по карте признаков
//head, info - служебные поля BMP файла
//T - карта признаков, Height, Width - высота и ширина карты
//fname - название файла в который будет сохранен результат
//t1, t2 - пороги значений, принимают значение в диапазонах (0, 0.5) (0.5, 1) соответственно
void formImage(BITMAPFILEHEADER head, BITMAPINFOHEADER info, float** &T, int Height, int Width, string fname, float t1, float t2)
{
	float Max=0, Min=0;
	getInterval(T, Height, Width, Min, Max);
	float T1 = (Max - Min)*t1 + Min;
	float T2 = (Max - Min)*t2 + Min;
	//std::cout << "Min = " << Min << endl;
	//std::cout << "Max = " << Max << endl;
	//std::cout << "T1 = " << T1 << endl;
	//std::cout << "T2 = " << T2 << endl;
	RGBQUAD** out = new RGBQUAD*[Height];
	for (int y = 0; y < Height; y++)
	{
		out[y] = new RGBQUAD[Width]();
		for (int x = 0; x < Width; x++)
		{
			//Если в первом диапазоне, то крашу в зеленый
			if (T[y][x] >= T2 && T[y][x] <= Max)
				out[y][x].rgbGreen = 255;

			//Если во втором диапазоне, то крашу в желтый
			if (T[y][x] >= T1 && T[y][x] < T2)
			{
				out[y][x].rgbRed = 255;
				out[y][x].rgbGreen = 255;
			}

			//Если в первом диапазоне, то крашу в красный
			if (T[y][x] >= Min && T[y][x] < T1)
				out[y][x].rgbRed = 255;
		}
	}
	//Сохраняю в файл
	BMPWrite(out, head, info, fname.c_str());
	for (int i = 0; i < Height; i++)
		delete[] out[i];
	delete[] out;
}