#pragma once
#include<omp.h>
#include"BMPFileRW.h"
#include"sorting.h"

typedef BYTE*(*ByteSortingMethod)(BYTE* arr, long length, IntComparer compare);

class Bitmap
{
public:// 73.3 155.2 //75.21 159.21
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
	int width;
	int height;
	RGBQUAD **map;
};

//заполнение медиального массива
RGBQUAD* getMedial(Bitmap image, int x, int y, int rh, int rw)
{
	int index = 0;
	RGBQUAD* barray = new RGBQUAD[2 * rw * 2 * rh];
	int coordX;
	int coordY;
	for (int dy = -rh; dy < rh; dy++)
	{
		coordY = y + dy;
		for (int dx = -rw; dx < rw; dx++)
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

RGBQUAD* sortRGB(RGBQUAD* arr, long length, ByteSortingMethod method)
{
	BYTE *red = new BYTE[length];
	BYTE *blue = new BYTE[length];
	BYTE *green = new BYTE[length];
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		blue[i] = arr[i].rgbBlue;
		green[i] = arr[i].rgbGreen;
	}

	red = method(red, length, Ascending);
	blue = method(blue, length, Ascending);
	green = method(green, length, Ascending);
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = {blue[i], green[i],red[i],  0};
	return narr;
}

//медиальная фильтрация
RGBQUAD** medialFiltering(Bitmap image, int wHeight, int wWidth)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp = getMedial(image, x, y, wHeight/2, wWidth/2); //заполняю медиальный массив
			//temp = BubbleEven(temp, wHeight * wWidth, [](RGBQUAD one, RGBQUAD two)->bool {return one.rgbBlue+one.rgbGreen+one.rgbRed > two.rgbBlue+two.rgbGreen+two.rgbRed; }); //сортирую массив
			temp = sortRGB(temp, wHeight * wWidth, ShellSort);
			out[y][x] = temp[wHeight * wWidth / 2]; // вытаскиваю срединный элемент
		}
	}
	return out;
}