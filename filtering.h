#pragma once
#include<omp.h>
#include"BMPFileRW.h"
#include"sorting.h"

//—игнатура метода сортирующего байты
typedef BYTE*(*ByteSortingMethod)(BYTE* arr, long length, IntComparer compare);

// ласс битмап изображени€, хран€щий поле пикселей ширину и высоту
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
		map[i] = new RGBQUAD[width];
		for (int j = 0; j < width; j++)
			map[i][j] = image[i][j];
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

//заполнение медиального массива
RGBQUAD* getMedial(Bitmap &image, int x, int y, int rh, int rw)
{
	int index = 0;
	RGBQUAD* barray = new RGBQUAD[rw * rh];
	int coordX;
	int coordY;
	for (int dy = -rh / 2; dy < rh / 2 + rh % 2; dy++)
	{
		coordY = y + dy;
		for (int dx = -rw / 2; dx < rw / 2 + rw % 2; dx++)
		{
			coordX = x + dx;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= image.width)
				coordX = image.width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= image.height)
				coordY = image.height - 1;
			barray[index] = image.map[coordY][coordX];
			index++;

		}
	}
	return barray;
}

//—ортировка массива –√Ѕ
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

//ѕараллельна€ сортировка массива –√Ѕ
RGBQUAD* sortRGBAsync(RGBQUAD* arr, long length, ByteSortingMethod sort)
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
#pragma omp parallel sections shared(red, blue, green, length) lastprivate(red1, blue1, green1)
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

//медиальна€ фильтраци€
RGBQUAD** medialFiltering(Bitmap &image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp1, *temp2;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
			temp2 = sortRGB(temp1, wHeight * wWidth, method); // сортирую каждую из компонент
			out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
	return out;
}

//медиальна€ фильтраци€ c распараллеливанием сортировки по компонентам
RGBQUAD** medialFilteringAsyncSort(Bitmap image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp1, *temp2;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
			temp2 = sortRGBAsync(temp1, wHeight * wWidth, method);
			out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
	return out;
}

//медиальна€ фильтраци€ for
RGBQUAD** medialFilteringAsync(Bitmap &image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height]; // на выходе картинка с примененным фильтром
	RGBQUAD *temp1, *temp2;
	#pragma omp parallel for private(temp1, temp2) shared(image, out) schedule(static, wHeight)
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
			temp2 = sortRGB(temp1, wHeight * wWidth, method);
			out[y][x] = temp2[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
	return out;
}

//медиальна€ фильтраци€ sections
RGBQUAD** medialFilteringAsyncSec(Bitmap &image, int wHeight, int wWidth, ByteSortingMethod method, int sections)
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
		#pragma omp parallel sections private(temp1,temp2, y) shared(image, out)
		{
			#pragma omp section
			{
				for (int x = t1; x < t2; x++)
				{
					//в окне H x W ложу пиксели в массив temp
					temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
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
					temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
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
					temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
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
					temp1 = getMedial(image, x, y, wHeight, wWidth); //заполн€ю медиальный массив
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