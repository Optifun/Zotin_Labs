#pragma once
#include<omp.h>
#include"BMPFileRW.h"
#include"sorting.h"

//��������� ������ ������������ �����
typedef BYTE*(*ByteSortingMethod)(BYTE* arr, long length, IntComparer compare);

//����� ������ �����������, �������� ���� �������� ������ � ������
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
	int width;
	int height;
	RGBQUAD **map;
};

//���������� ����������� �������
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

//���������� ������� ���
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

	red = method(red, length, AscInt);
	blue = method(blue, length, AscInt);
	green = method(green, length, AscInt);
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = {blue[i], green[i],red[i],  0};
	return narr;
}

//������������ ���������� ������� ���
RGBQUAD* sortRGBAsync(RGBQUAD* arr, long length, ByteSortingMethod method)
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
#pragma omp parallel sections shared(red, blue, green, length)
	{
#pragma omp section
		{
	red = method(red, length, AscInt);
		}
#pragma omp section
		{
	blue = method(blue, length, AscInt);
		}
#pragma omp section
		{
	green = method(green, length, AscInt);
		}
	}
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = { blue[i], green[i],red[i],  0 };
	return narr;
}

//���������� ����������
RGBQUAD** medialFiltering(Bitmap image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//� ���� H x W ���� ������� � ������ temp
			temp = getMedial(image, x, y, wHeight/2, wWidth/2); //�������� ���������� ������
			//temp = BubbleEven(temp, wHeight * wWidth, [](RGBQUAD one, RGBQUAD two)->bool {return one.rgbBlue+one.rgbGreen+one.rgbRed > two.rgbBlue+two.rgbGreen+two.rgbRed; }); //�������� ������
			temp = sortRGB(temp, wHeight * wWidth, method);
			out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
		}
	}
	return out;
}

//���������� ���������� c ������������������ ���������� �� �����������
RGBQUAD** medialFilteringAsyncSort(Bitmap image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height];
	RGBQUAD *temp;
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//� ���� H x W ���� ������� � ������ temp
			temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
																	//temp = BubbleEven(temp, wHeight * wWidth, [](RGBQUAD one, RGBQUAD two)->bool {return one.rgbBlue+one.rgbGreen+one.rgbRed > two.rgbBlue+two.rgbGreen+two.rgbRed; }); //�������� ������
			temp = sortRGBAsync(temp, wHeight * wWidth, method);
			out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
		}
	}
	return out;
}

//���������� ���������� for
RGBQUAD** medialFilteringAsync(Bitmap image, int wHeight, int wWidth, ByteSortingMethod method)
{
	RGBQUAD **out = new RGBQUAD*[image.height]; // �� ������ �������� � ����������� ��������
	RGBQUAD *temp;
	#pragma omp parallel for firstprivate(method, wHeight, wWidth) shared(image, out)
	for (int y = 0; y < image.height; y++)
	{
		out[y] = new RGBQUAD[image.width];
		for (int x = 0; x < image.width; x++)
		{
			//� ���� H x W ���� ������� � ������ temp
			temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
			temp = sortRGB(temp, wHeight * wWidth, method); //�������� ������
			out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
		}
	}
	return out;
}

//���������� ���������� sections
//RGBQUAD** medialFilteringAsyncSort(Bitmap image, int wHeight, int wWidth, ByteSortingMethod method, int sections)
//{
//	RGBQUAD **out = new RGBQUAD*[image.height];
//	RGBQUAD *temp;
//	float step = image.width / sections;
//	int mod = image.width % sections;
//	int t1 = 0;
//	int t2 = step * 1+mod;
//	int t3 = t2 + step;
//	int t4 = t3+step;
//	for (int y = 0; y < image.height; y++)
//	{
//		out[y] = new RGBQUAD[image.width];
//		#pragma omp parallel sections private(temp, y) shared(image, out)
//		{
//			#pragma omp section
//			{
//				for (int x = t1; x < t2; x++)
//				{
//					temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
//					temp = sortRGB(temp, wHeight * wWidth, method);
//					out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
//				}
//			}
//#pragma omp section
//			{
//				if (sections>=2)
//				for (int x = t2; x < t3; x++)
//				{
//					temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
//					temp = sortRGB(temp, wHeight * wWidth, method);
//					out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
//				}
//			}
//#pragma omp section
//			{
//				if (sections>=3)
//				for (int x = t3; x < t4; x++)
//				{
//					temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
//					temp = sortRGB(temp, wHeight * wWidth, method);
//					out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
//				}
//			}
//#pragma omp section
//			{
//				if (sections>=4)
//				for (int x = t4; x < t4+step; x++)
//				{
//					temp = getMedial(image, x, y, wHeight / 2, wWidth / 2); //�������� ���������� ������
//					temp = sortRGB(temp, wHeight * wWidth, method);
//					out[y][x] = temp[wHeight * wWidth / 2]; // ���������� ��������� �������
//				}
//			}
//		}
//
//	}
//	return out;
//}