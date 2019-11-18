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

//Ћинейный фильтр последовательный; RH, RW - размеры рангов скольз€щего окна
void LineFilteringSred(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -RH; DY < RH; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -RW; DX < RW; DX++)
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

//‘ормирование матрицы коэффициентов дл€ фильтрации √аусса
double** GetGaussMatrix(int RH, int RW, int q) {
	double** Result = new double*[RH * 2 + 1];
	for (int i = 0; i < RH * 2 + 1; i++)
		Result[i] = new double[RW * 2 + 1];
	double SUM = 0;
	double CF;
	for (int Y = -RH; Y < RH; Y++)
		for (int X = -RW; X < RW; X++) {
			CF = (1 / (2 * 3.14159265358979323846 * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));

			SUM += CF;
			Result[Y + RH][X + RW] = CF;
		}
	for (int Y = -RH; Y < RH; Y++)
		for (int X = -RW; X < RW; X++) {
			Result[Y + RH][X + RW] /= SUM;
		}
	return Result;
}
//Ћинейный фильтр √аусса последовательный; RH, RW - размеры рангов скольз€щего окна
void LineFilteringGauss(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GetGaussMatrix(RH, RW, 10); //—игма тут
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -RH; DY < RH; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -RW; DX < RW; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					rgbBlue += RGB[KY][KX].rgbBlue * CoefMatrix[DY + RH][DX + RW];
					rgbGreen += RGB[KY][KX].rgbGreen * CoefMatrix[DY + RH][DX + RW];
					rgbRed += RGB[KY][KX].rgbRed * CoefMatrix[DY + RH][DX + RW];
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue;
			RGBresult[Y][X].rgbGreen = rgbGreen;
			RGBresult[Y][X].rgbRed = rgbRed;
		}
	for (int i = 0; i < RW; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

//ѕример использовани€
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