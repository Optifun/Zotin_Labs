#pragma once
#include"filtering.h"
#include<list>

float* formHist(Bitmap &image, BYTE** BrMap, int x, int y, int RH, int RW)
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

			if (coordX >= image.width)
				coordX = image.width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= image.height)
				coordY = image.height - 1;
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

//На вход изображение и размер рамки
//Выход:
//Моменты второго порядка, однородность, относительная гладкость, ентропия
//Moments, Uniform, Relative, Enthropy
void textureFilter(Bitmap &image, int rh, int rw, float **M, float **U, float **R, float **E)
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
			Brightness[y][x] = image.map[y][x].rgbRed*0.299 + image.map[y][x].rgbGreen*0.587 + image.map[y][x].rgbBlue*0.114;
	for (int y = 0; y < image.height; y++)
		for (int x = 0; x < image.width; x++)
		{
			hist = formHist(image, Brightness, x, y, rh, rw);
			getMetrics(M[y][x], U[y][x], R[y][x], E[y][x], hist);
			delete[] hist;
		}
}

void formImage(float** T, int Height, int Width, string fname, float _min, float t1, float t2, float _max)
{

}