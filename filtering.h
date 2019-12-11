#pragma once
#include<omp.h>
#include"BMPFileRW.h"
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>
#include<cilk/reducer_opadd.h>
#include<cilk/reducer_max.h>
#include<cilk/reducer_opmul.h>
#include"sorting.h"

using namespace std;
//Сигнатура метода сортирующего байты
typedef BYTE*(*ByteSortingMethod)(BYTE* arr, long length, IntComparer compare);

//Класс битмап изображения, хранящий поле пикселей ширину и высоту
class Bitmap
{
public:
	Bitmap(RGBQUAD** &image, int _w, int _h)
	{
		width = _w;
		height = _h;
		map = new RGBQUAD*[height];
		for (int i = 0; i < height; i++)
		{
		map[i] = new RGBQUAD[width]();
		memcpy(map[i], image[i], width * sizeof(RGBQUAD));
		}
	}
	Bitmap(Bitmap& t)
	{
		width = t.width;
		height = t.height;
		map = new RGBQUAD*[height];
		for (int i = 0; i < height; i++)
		{
			map[i] = new RGBQUAD[width];
			memcpy(map[i], t.map[i], width * sizeof(RGBQUAD));
		}
	}

	~Bitmap()
	{
		for (int i = 0; i < height; i++)
			delete[] map[i];
		delete[] map;
	}
	int width;
	int height;
	RGBQUAD **map;
};

#pragma region getMedial

//заполнение медиального массива
//image - исходная картинка
//(x,y) - центр рамки
//RH, RW - радиусы рамки по высоте(height) и ширине(width)
RGBQUAD* getMedial(RGBQUAD **&image, int width, int height, int x, int y, int RH, int RW)
{
	int index = 0;
	RGBQUAD* barray = new RGBQUAD[(2 * RW + 1) * (2 * RH + 1)];
	int coordX;
	int coordY;
	for (int dy = -RH; dy <= RH; dy++)
	{
		coordY = y + dy;
		for (int dx = -RW; dx <= RW; dx++) 
		{
			coordX = x + dx;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= width)
				coordX = width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= height)
				coordY = height - 1;
			barray[index] = image[coordY][coordX];
			index++;
		}
	}
	return barray;
}

//заполнение медиального массива с Cilk For
//image - исходная картинка
//(x,y) - центр рамки
//RH, RW - радиусы рамки по высоте(height) и ширине(width)
RGBQUAD* getMedialCilkFor(RGBQUAD **&image, int width, int height, int x, int y, int RH, int RW)
{
	RGBQUAD* barray = new RGBQUAD[(2 * RW + 1) * (2 * RH + 1)];
	cilk_for (int dy = -RH; dy <= RH; dy++)
	{
		int coordY = y + dy;
		for (int dx = -RW; dx <= RW; dx++)
		{
			int coordX = x + dx;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= width)
				coordX = width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= height)
				coordY = height - 1;

			// индекс = строка * кол-во элем в строке + столбец
			int index = (RH + dy) *(2 * RH + 1) + (RW + dx);
			barray[index] = image[coordY][coordX];
		}
	}
	return barray;
}

#pragma endregion

#pragma region sortRGB

//Сортировка массива РГБ
RGBQUAD* sortRGB(RGBQUAD* arr, long length, ByteSortingMethod sort)
{
	BYTE *red = new BYTE[length];
	BYTE *blue = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *red1;
	BYTE *blue1;
	BYTE *green1;
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		blue[i] = arr[i].rgbBlue;
		green[i] = arr[i].rgbGreen;
	}

	red1 = sort(red, length, AscInt);
	blue1 = sort(blue, length, AscInt);
	green1 = sort(green, length, AscInt);
	delete[] red;
	delete[] blue;
	delete[] green;
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = {blue1[i], green1[i],red1[i],  0};
	delete[] red1;
	delete[] green1;
	delete[] blue1;
	return narr;
}

//сортировка массива РГБ с Cilk spawn + vectorization
RGBQUAD* sortRGBVec(RGBQUAD* arr, long length, ByteSortingMethod sort)
{
	BYTE *red = new BYTE[length];
	BYTE *blue = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *red1;
	BYTE *blue1;
	BYTE *green1;
	red[0:length] = arr[0:length].rgbRed;
	blue[0:length] = arr[0:length].rgbBlue;
	green[0:length] = arr[0:length].rgbGreen;

	cilk_spawn[&red, length, &red1, sort]() {
	red1 = sort(red, length, AscInt);
	}();
	cilk_spawn[&blue, length, &blue1, sort]() {
		blue1 = sort(blue, length, AscInt);
	}();
	green1 = sort(green, length, AscInt);
	RGBQUAD* narr = new RGBQUAD[length];
	cilk_sync;
	for (int i = 0; i < length; i++)
		narr[i] = { blue1[i], green1[i], red1[i], 0 };
	delete[] red;
	delete[] blue;
	delete[] green;
	delete[] red1;
	delete[] green1;
	delete[] blue1;
	return narr;
}

//сортировка массива РГБ с Omp sections+for
RGBQUAD* sortRGBAsync(RGBQUAD* arr, long length, ByteSortingMethod sort)
{
	BYTE *red = new BYTE[length];
	BYTE *blue = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *red1;
	BYTE *blue1;
	BYTE *green1;
	#pragma omp parallel for shared(arr, red, blue, green)
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		blue[i] = arr[i].rgbBlue;
		green[i] = arr[i].rgbGreen;
	}
#pragma omp parallel sections shared(red, blue, green, length, red1, blue1, green1)
	{
#pragma omp section
		{
		red1 = sort(red, length, AscInt);
		}
#pragma omp section
		{
		blue1 = sort(blue, length, AscInt);
		}
#pragma omp section
		{
		green1 = sort(green, length, AscInt);
		}
	}
	delete[] red;
	delete[] green;
	delete[] blue;
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = { blue1[i], green1[i],red1[i],  0 };
	delete[] red1;
	delete[] green1;
	delete[] blue1;
	return narr;
}

#pragma endregion

#pragma region medianFiltering

//медианная фильтрация
//RGB - исходное изображение
//RH, RW - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//Возвращает RGBQUAD -  изображение с примененным на нём медианную фильтрацию
void medianFiltering(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, ByteSortingMethod method)
{
	RGBresult = new RGBQUAD*[height];
	RGBQUAD *temp1, *temp2;
	int size = (2 * RH + 1) * (2 * RW + 1);
	for (int y = 0; y < height; y++)
	{
		RGBresult[y] = new RGBQUAD[width];
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, RH, RW); //заполняю медиальный массив
			temp2 = sortRGB(temp1, size, method); // сортирую каждую из компонент
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация c распараллеливанием сортировки по компонентам
//RGB - исходное изображение
//RH, RW - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//Возвращает RGBQUAD -  изображение с примененным на нём медианную фильтрацию
void medianFilteringAsyncSort(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, ByteSortingMethod method)
{
	RGBresult = new RGBQUAD*[height];
	RGBQUAD *temp1, *temp2;
	int size = (2 * RH + 1) * (2 * RW + 1);
	for (int y = 0; y < height; y++)
	{
		RGBresult[y] = new RGBQUAD[width];
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, RH, RW); //заполняю медиальный массив
			temp2 = sortRGBAsync(temp1, size, method);
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация Omp For
//RGB - исходное изображение
//RH, RW - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//Возвращает RGBQUAD -  изображение с примененным на нём медианную фильтрацию
void medianFilteringAsync(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, ByteSortingMethod method)
{
	RGBresult = new RGBQUAD*[height]; // на выходе картинка с примененным фильтром
	RGBQUAD *temp1, *temp2;
	int size = (2 * RH + 1) * (2 * RW + 1);
	#pragma omp parallel for private(temp1, temp2) shared(RGB, RGBresult) schedule(static, RH)
	for (int y = 0; y < height; y++)
	{
		RGBresult[y] = new RGBQUAD[width];
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, RH, RW); //заполняю медиальный массив
			temp2 = sortRGB(temp1, size, method);
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация c Cilk For
//image - исходное изображение
//wHeight, wWidth - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//Возвращает изображение с примененным на нём медианную фильтрацию
void medianFilteringCilkFor(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, ByteSortingMethod method)
{
	RGBresult = new RGBQUAD*[height];
	int size = (2 * RH + 1) * (2 * RW + 1);
	cilk_for (int y = 0; y < height; y++)
	{
		RGBresult[y] = new RGBQUAD[width];
		for (int x = 0; x < width; x++)
		{
			//распараллеливание силк фор здесь даёт замедление
			RGBQUAD *temp1 = getMedial(RGB, width, height, x, y, RH, RW);
			RGBQUAD *temp2 = sortRGBVec(temp1, size, method);
			RGBresult[y][x] = temp2[size / 2];
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация Omp Sections
//image - исходное изображение
//wHeight, wWidth - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//Возвращает изображение с примененным на нём медианную фильтрацию
RGBQUAD** medianFilteringAsyncSec(Bitmap &image, int wHeight, int wWidth, ByteSortingMethod method, int sections)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp1, *temp2;
	float step = image.width / sections;
	int mod = image.width % sections;
	int t1 = 0;
	int t2 = step * 1+mod;
	int t3 = t2 + step;
	int t4 = t3+step;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		#pragma omp parallel sections private(temp1,temp2) shared(image, out)
		{
			#pragma omp section
			{
				for (int x = t1; x < t2; x++)
				{
					//в окне H x W ложу пиксели в массив temp
					temp1 = getMedial(image.map, image.width, image.height, x, y, wHeight, wWidth); //заполняю медиальный массив
					temp2 = sortRGB(temp1, wHeight * wWidth, method);
					out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
					delete[] temp1;
					delete[] temp2;
				}
			}
#pragma omp section
			{
				if (sections>=2)
				for (int x = t2; x < t3; x++)
				{
					//в окне H x W ложу пиксели в массив temp
					temp1 = getMedial(image.map, image.width, image.height, x, y, wHeight, wWidth); //заполняю медиальный массив
					temp2 = sortRGB(temp1, wHeight * wWidth, method);
					out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
					delete[] temp1;
					delete[] temp2;
				}
			}
#pragma omp section
			{
				if (sections>=3)
				for (int x = t3; x < t4; x++)
				{
					//в окне H x W ложу пиксели в массив temp
					temp1 = getMedial(image.map, image.width, image.height, x, y, wHeight, wWidth); //заполняю медиальный массив
					temp2 = sortRGB(temp1, wHeight * wWidth, method);
					out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
					delete[] temp1;
					delete[] temp2;
				}
			}
#pragma omp section
			{
				if (sections>=4)
				for (int x = t4; x < t4+step; x++)
				{
					//в окне H x W ложу пиксели в массив temp
					temp1 = getMedial(image.map, image.width, image.height, x, y, wHeight, wWidth); //заполняю медиальный массив
					temp2 = sortRGB(temp1, wHeight * wWidth, method);
					out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
					delete[] temp1;
					delete[] temp2;
				}
			}
		}

	}
	return out;
}

#pragma endregion

#pragma region LinerFiltering

//Линейный средний фильтр последовательный; RH, RW - размеры рангов скользящего окна
void LineFilteringSred(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
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
					rgbBlue += RGB[KY][KX].rgbBlue;
					rgbGreen += RGB[KY][KX].rgbGreen;
					rgbRed += RGB[KY][KX].rgbRed;
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue / ((RH * 2 + 1) * (RW * 2 + 1));
			RGBresult[Y][X].rgbGreen = rgbGreen / ((RH * 2 + 1) * (RW * 2 + 1));
			RGBresult[Y][X].rgbRed = rgbRed / ((RH * 2 + 1) * (RW * 2 + 1));
		}
}

//Линейный средний фильтр параллельный; RH, RW - размеры рангов скользящего окна
void LineFilteringSredParal(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
#pragma omp parallel for firstprivate(RH, RW, height, width) shared(RGB, RGBresult) schedule(static, RH * 2 + 1)
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
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
					rgbBlue += RGB[KY][KX].rgbBlue;
					rgbGreen += RGB[KY][KX].rgbGreen;
					rgbRed += RGB[KY][KX].rgbRed;
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue / ((RH * 2 + 1) * (RW * 2 + 1));
			RGBresult[Y][X].rgbGreen = rgbGreen / ((RH * 2 + 1) * (RW * 2 + 1));
			RGBresult[Y][X].rgbRed = rgbRed / ((RH * 2 + 1) * (RW * 2 + 1));
		}
}

#pragma endregion

#pragma region GausMatrix

//Формирование матрицы коэффициентов для фильтрации Гаусса
double** GetGaussMatrix(int RH, int RW, double q) {
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

//Формирование матрицы коэффициентов для фильтрации Гаусса параллельное
double** GetGaussMatrixParal(int RH, int RW, double q) {
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

//Формирование матрицы коэффициентов для фильтрации Гаусса
double** GetGaussMatrixCilk(int RH, int RW, double q) 
{
	double** Result = new double*[RH * 2 + 1];
	cilk_for (int i = 0; i < RH * 2 + 1; i++)
		Result[i] = new double[RW * 2 + 1];

	cilk::reducer<cilk::op_add<double>> SUM = 0;
	cilk_for (int Y = -RH; Y <= RH; Y++)
		for (int X = -RW; X <= RW; X++) {
			double CF = (1 / (2 * 3.14159265358979323846 * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));
			Result[Y + RH][X + RW] = CF;
			*SUM += CF;
		}
	double sum = SUM.get_value();
	cilk_for (int Y = -RH; Y <= RH; Y++)
			Result[Y + RH][0:2*RW + 1] /= sum;
	return Result;
}
#pragma endregion

#pragma region GausFiltering

//Линейный фильтр Гаусса последовательный; RH, RW - размеры рангов скользящего окна
//Возвращает RGBresult указывающий на выходную картинку
void LineFilteringGauss(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrix(RH, RW, RW / 3.0); //Сигма тут
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

//Линейный фильтр Гаусса параллельеный; RH, RW - размеры рангов скользящего окна
//Возвращает RGBresult указывающий на выходную картинку
void LineFilteringGaussParal(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrix(RH, RW, RW / 3.0); //Сигма тут
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


//Линейный фильтр Гаусса последовательный; RH, RW - размеры рангов скользящего окна
//Возвращает RGBresult указывающий на выходную картинку
void LineFilteringGaussCilk(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrixCilk(RH, RW, RW / 3.0); //Сигма тут
	RGBresult = new RGBQUAD*[height];
	cilk_for(int Y = 0; Y < height; Y++)
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

	cilk_for (int i = 0; i < RW; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

#pragma endregion

//Пример использования
//int main() {
//	cout << "Hello World!\n";
//	RGBQUAD **RGB, **RGBresult;
//	BITMAPFILEHEADER head;
//	BITMAPINFOHEADER info;
//	string str = "1.bmp";
//	string str2 = "11.bmp";
//	BMPRead(RGB, head, info, "1.bmp");
//	RGBresult = new RGBQUAD*[info.biHeight];
//	for (int i = 0; i < info.biHeight; i++)
//		RGBresult[i] = new RGBQUAD[info.biWidth];
//	LineFilteringGauss(RGB, info.biHeight, info.biWidth, 5, 5, RGBresult);
//	BMPWrite(RGBresult, head, info, str2.c_str());
//}