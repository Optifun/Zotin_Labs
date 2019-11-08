#pragma once
#include <Windows.h>
#include <iostream>
#include <fstream>

void BMPWrite(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *);
void BMPRead(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *);


void BMPRead(RGBQUAD** &rgb, BITMAPFILEHEADER &header, \
	BITMAPINFOHEADER &bmiHeader, const char* fin)
{
    // Открываем файл
	std::ifstream InFile(fin, std::ios::binary);
	// Считываем заголовок файла
	InFile.read((char*)(&header), sizeof(BITMAPFILEHEADER));
	// Считываем заголовочную часть изображения
	InFile.read((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // Выделяем память под массив RGB хранящий структуры RGBQUAD
    rgb = new RGBQUAD*[bmiHeader.biHeight];
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        rgb[i] = new RGBQUAD[bmiHeader.biWidth];
    }
    // Считываем данные изображения в массив структур RGB 
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
			InFile.read((char*)(&rgb[i][j]), 3); // .rgbBlue .rgbGreen .rgbRed
        }
    }
    // Закрываем файл
	InFile.close();
	
}

void BMPWrite(RGBQUAD** &rgb, BITMAPFILEHEADER &header, \
	BITMAPINFOHEADER &bmiHeader, const char *fout)
{
    // Открываем файл для записи изображения в формат BMP
	std::ofstream OutFile(fout, std::ios::binary);
	//// Записываем заголовок файла
	OutFile.write((char*)(&header), sizeof(BITMAPFILEHEADER));
	//// Записываем заголовочную часть изображения
	OutFile.write((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // Записываем данные изображения из массив структур RGB в файл 
    for (int i = 0; i < bmiHeader.biHeight; i++)
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
			OutFile.write((char*)&(rgb[i][j]), 3);// .rgbBlue .rgbGreen .rgbRed;
        }
    // закрываем файл
	OutFile.close();
}
