#include<stdlib.h>
#include<stdio.h>

struct BITMAPFILEHEADER {
	unsigned int    bfType;
	unsigned long   bfSize;
	unsigned int    bfReserved1;
	unsigned int    bfReserved2;
	unsigned long   bfOffBits; 
};

struct BITMAPINFOHEADER {
	unsigned int    biSize;
	int             biWidth;
	int             biHeight;
	unsigned short  biPlanes;
	unsigned short  biBitCount;
	unsigned int    biCompression;
	unsigned int    biSizeImage;
	int             biXPelsPerMeter;
	int             biYPelsPerMeter;
	unsigned int    biClrUsed;
	unsigned int    biClrImportant;
};

struct RGBQUAD {
	int   rgbBlue;
	int   rgbGreen;
	int   rgbRed;
	int   rgbReserved; 
	//не знаю как иначе сравнивать пиксели
	bool operator>(RGBQUAD& other)
	{
		int sum1 = rgbBlue + rgbGreen + rgbRed;
		int sum2 = other.rgbBlue + other.rgbGreen + other.rgbRed;
		return sum1 > sum2;
	}
};

static unsigned short read_u16(FILE *fp); 
static unsigned int read_u32(FILE *fp); 
static int read_s32(FILE *fp);

static void write_u16(unsigned short input, FILE *fp); 
static void write_u32(unsigned int input, FILE *fp); 
static void write_s32(int input, FILE *fp);

void BMPWrite(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *); 
void BMPRead(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *);

static unsigned short read_u16(FILE *fp) { 
	unsigned char b0, b1;  
	b0 = fgetc(fp);   
	b1 = fgetc(fp);  
	return ((b1 << 8) | b0);
}

static unsigned int read_u32(FILE *fp) {
	unsigned char b0, b1, b2, b3;
	b0 = fgetc(fp);     
	b1 = fgetc(fp);
	b2 = fgetc(fp);     
	b3 = fgetc(fp);
	return ((((((b3 << 8) | b2) << 8) | b1) << 8) | b0);
} 

static int read_s32(FILE *fp) {
	unsigned char b0, b1, b2, b3;     
	b0 = fgetc(fp);    
	b1 = fgetc(fp); 
	b2 = fgetc(fp);
	b3 = fgetc(fp);
	return ((int)(((((b3 << 8) | b2) << 8) | b1) << 8) | b0); 
}

static void write_u16(unsigned short input, FILE *fp)
{
	fputc(input, fp); 
	fputc(input >> 8, fp);
} 

static void write_u32(unsigned int input, FILE *fp) {
	fputc(input, fp); 
	fputc(input >> 8, fp); 
	fputc(input >> 16, fp);
	fputc(input >> 24, fp);
}

static void write_s32(int input, FILE *fp) {
	fputc(input, fp);  
	fputc(input >> 8, fp);   
	fputc(input >> 16, fp);    
	fputc(input >> 24, fp);
}

void BMPRead(RGBQUAD** &rgb, BITMAPFILEHEADER &header, BITMAPINFOHEADER &bmiHeader, const char* fin) 
{     // Открываем файл  
	FILE * pFile = fopen(fin, "rb");  
	// Считываем заголовок файла  
	header.bfType = read_u16(pFile); 
	header.bfSize = read_u32(pFile);   
	header.bfReserved1 = read_u16(pFile);   
	
	header.bfReserved2 = read_u16(pFile);  
	header.bfOffBits = read_u32(pFile);   
	// Считываем заголовочную часть изображения   
	bmiHeader.biSize = read_u32(pFile);
	bmiHeader.biWidth = read_s32(pFile); 
	bmiHeader.biHeight = read_s32(pFile);   
	bmiHeader.biPlanes = read_u16(pFile);   
	bmiHeader.biBitCount = read_u16(pFile);  
	bmiHeader.biCompression = read_u32(pFile); 
	bmiHeader.biSizeImage = read_u32(pFile);    
	bmiHeader.biXPelsPerMeter = read_s32(pFile);  
	bmiHeader.biYPelsPerMeter = read_s32(pFile);   
	bmiHeader.biClrUsed = read_u32(pFile);   
	bmiHeader.biClrImportant = read_u32(pFile);    
	/* Выделяем память под массив RGB хранящий структуры RGBQUAD */ 
	rgb = new RGBQUAD*[bmiHeader.biHeight];     
	for (int i = 0; i < bmiHeader.biHeight; i++)   
	{        
		
		rgb[i] = new RGBQUAD[bmiHeader.biWidth]; 
	}    
	/* Считываем данные изображения в массив структур RGB */    
	for (int i = 0; i < bmiHeader.biHeight; i++)    
	{         
		for (int j = 0; j < bmiHeader.biWidth; j++) 
		{            
			rgb[i][j].rgbBlue = fgetc(pFile);       
			rgb[i][j].rgbGreen = fgetc(pFile);           
			rgb[i][j].rgbRed = fgetc(pFile);        
		}    
	}    
	// Закрываем файл  
	fclose(pFile);  
} 

	void BMPWrite(RGBQUAD** &rgb, BITMAPFILEHEADER &header,  BITMAPINFOHEADER &bmiHeader, const char *fout) 
	{
			/* Открываем файл для записи изображения в формат BMP */    
		FILE * oFile = fopen(fout, "wb");    
		// Записываем заголовок файла    
		write_u16(header.bfType, oFile); 
		write_u32(header.bfSize, oFile);   
		write_u16(header.bfReserved1, oFile); 
		write_u16(header.bfReserved2, oFile); 
		write_u32(header.bfOffBits, oFile);   
		// Записываем заголовочную часть изображения   
		write_u32(bmiHeader.biSize, oFile);  
		write_s32(bmiHeader.biWidth, oFile);   
		write_s32(bmiHeader.biHeight, oFile);   
		write_u16(bmiHeader.biPlanes, oFile);   
		write_u16(bmiHeader.biBitCount, oFile); 
		write_u32(bmiHeader.biCompression, oFile); 
		write_u32(bmiHeader.biSizeImage, oFile);    
		write_s32(bmiHeader.biXPelsPerMeter, oFile);  
		write_s32(bmiHeader.biYPelsPerMeter, oFile);   
		write_u32(bmiHeader.biClrUsed, oFile);   
		write_u32(bmiHeader.biClrImportant, oFile);   
		/* Записываем данные изображения из массива структур RGB в файл */   
		for (int i = 0; i < bmiHeader.biHeight; i++)    
			for (int j = 0; j < bmiHeader.biWidth; j++)   
			{            
				
				fputc(rgb[i][j].rgbBlue, oFile);     
				fputc(rgb[i][j].rgbGreen, oFile);     
				fputc(rgb[i][j].rgbRed, oFile);       
			}        
		// закрываем файл  
		fclose(oFile); 
	}