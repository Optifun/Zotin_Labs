#pragma once
#include<fstream>
#include<iostream>
#include<sstream>
#include<iomanip>

using namespace std;
class Log
{
public:
	Log()
	{
		Fout = false;
		Cout = true;
	}

	Log(string path, bool file, bool console = true)
	{
		this->file = ofstream(path);
		Cout = console;
		Fout = file;
	}
	
	//Запись строки в файл
	Log& operator<<(string str)
	{
		if (Cout)
			cout << str;
		if (Fout)
			file << str;
		return *this;
	}

	//Записть дробного числа в файл
	Log& operator<<(double num)
	{
		if (Cout)
			cout << num;
		if (Fout)
			file << num;
		return *this;
	}

	//Запись целого числа в файл
	Log& operator<<(int num)
	{
		if (Cout)
			cout << num;
		if (Fout)
			file << num;
		return *this;
	}

	//Закрывает файл с сохранением содержимого
	void Close()
	{
		file.close();
	}

	//Выгружает данные буфера в файл
	void Flush()
	{
		file.flush();
	}

	//Log operator=(const Log& log)
	//{
	//	file.close();
	//	file = ofstream(log.filename);
	//	Cout = log.Cout;
	//	Fout = log.Fout;
	//	filename = log.filename;
	//	return *this;
	//}
	//
	//~Log()
	//{
	//	Close();
	//}

private:
	ofstream file;
	string filename;
	bool Cout;
	bool Fout;
};

//глобальная переменная логгера
Log logg;