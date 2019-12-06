#pragma once 
#include <iostream>
#include <math.h>
#include <omp.h>
#include <string>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_opmul.h>
using namespace std;

//Сигнатура [](long, long)->double {return ...;}
typedef double(*fillF)(long, long);
#define MaxOf(a, b) ((a)< (b)) ? (b) : (a)

static string doubleToString(double number)
{
	char buf[256];
	int size = sprintf_s(buf, "%f", number);
	string num = string(buf, size);
	return num;
}

template<class T>
class Filler
{
public:
	Filler(fillF func)
	{
		generator = func;
	}

	T operator()(long x, long y)
	{
		return generator(x, y);
	}
private:
	fillF generator;
};

// Метод параллелизации алгоритма
enum Method {
	//Последовательный вариант
	Sequence,
	//Параллельный цикл на ОМП
	ParallelFor,
	//Параллельные секции на ОМП
	ParallelSections,
	//Параллельный цикл на Силке
	Cilk_for,
	//Спавн параллельных веток на Силке
	Cilk_spawn,
	//Индексная нотация(векторизация) на Силке
	Cilk_index,
	//Параллельный цикл + векторизация на Силке
	Cilk_for_index};

//Класс матрицы
template<class Type>
class Matrix
{
public:

	//Количество созданных матриц за время выполнения программы
	static long Created;
	
	//Количество удаленных матриц за время выполнения программы
	static long Destroed;
	
	//Включение режима отладки
	//Нумерует каждую созданную матрицу и выводит сообщения о создании и удалении каждой матрицы
	static bool DEBUG;

//Конструкторы матриц
#pragma region Constructors

	//Возвращает матрицу размерностью NxM заполненную единицами
	static Matrix<Type> Ones(long _n, long _m)
	{
		return Matrix<Type>(_n, _m, Filler<Type>([](long, long)->double {return 1; }));
	}
	
	//Возвращает матрицу размерностью NxM заполненную нулями
	static Matrix<Type> Zeros(long _n, long _m)
	{
		return Matrix<Type>(_n, _m, Filler<Type>([](long, long)->double {return 0; }));
	}
	
	//Возвращает единичную матрицу размерностью NxM
	static Matrix<Type> E(long _n, long _m)
	{
		auto t = min(_n, _m);
		return Matrix<Type>(t, t, Filler<Type>([](long x, long y)->double {return (x == y) ? 1 : 0; }));
	}

	//Возвращает матрицу размерностью 1х1 с элементом равным нулю
	Matrix()
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = 1;
		m = 1;
		arr = allocateArr(1, 1);
		arr[0][0] = 0;
	}

	//Возвращает матрицу размерностью NxM заполненную нулями
	Matrix(long _n, long _m)
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = _n;
		m = _m;
		arr = allocateArr(n, m);
		fillMatrix(Filler<Type>([](long, long) ->double { return 0; }));
	}

	//Возвращает копию матрицы
	Matrix(Matrix<Type>& M)
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = M.n;
		m = M.m;
		copyArr(M.arr, true);
	}

	//Возвращает матрицу размерностью NxM, элементы которой заполнены используя функцию сигнатуры fillF
	Matrix(long _n, long _m, Filler<Type> filler)
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = _n;
		m = _m;
		arr = allocateArr(n, m);
		fillMatrix(filler);
	}

	//Высвобождает всю выделенную матрицей память
	~Matrix()
	{
		if (DEBUG)
			cout << "D:[" << num <<"]";
		Destroed++;
		for (long i = 0; i < n; i++)
		{
			//for (long j = 0; j < m; j++)
			//	delete[] arr[i][j];
			delete[] arr[i];
		}
		delete[] arr;
		//cout << "Умир" << endl;
	}

	static void operator delete(void* mem, size_t length)
	{
		delete(mem);
	}
	
	static void operator delete[](void* mem, size_t length)
	{
		delete(mem);
	}

#pragma endregion

// Обращения по индексу к матрице
#pragma region Selectors


	//Возращает iю строку матрицы
	Type* getRow(long _i)
	{
		Type * arr1 = new Type[m];
		for (long j = 0; j < m; j++)
			arr1[j] = arr[_i][j];
		return arr1;
	}

	//Возращает jй столбец матрицы
	Type* getCol(long _j)
	{
		Type * arr1 = new Type[n];
		for (long j = 0; j < n; j++)
			arr1[j] = arr[j][_j];
		return arr1;
	}

	//Возвращает элемент матрицы
	Type* getElement(long i, long j)
	{
		if (i >= n || j >= m|| i<0 || j<0)
			throw new exception("Индекс поиска больше размерности матрицы");
		return &arr[i][j];
	}

	//Возвращает iю строку матрицы, представленную в виде столбца
	Matrix<Type> operator[] (long i)
	{
		double **arr1;
		if (m == 1)
		{
			arr1 = allocateArr(1, 1);
			**arr1 = arr[i][0];
			return Matrix(1, 1, arr1);
		}
		if (n == 1)
		{
			arr1 = allocateArr(1, 1);
			**arr1 = arr[0][i];
			return Matrix(1, 1, arr1);
		}
		arr1 = allocateArr(m, 1);
		for (long j = 0; j < m; j++)
		{
			arr1[j][0] = arr[i][j];
		}
		return Matrix(m, 1, arr1);
	}

	Matrix<Type> operator=(Matrix<Type> & M)
	{
		bool tr=false;
		if (arr != NULL &&(n != M.n || m != M.m))
		{
			tr = true;
			if (n>1)
				for (long i = 0; i < n; i++)
					delete[] arr[i];
			delete[] arr;
		}
		m = M.m;
		n = M.n;
		//arr = allocateArr(n, m);
		copyArr(M.arr, tr);
		return *this;

	}

	//Сравнение на равенство матриц
	bool operator == (Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			return false;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				if (arr[i][j] != M.arr[i][j])
					return false;
		return true;
	}

#pragma endregion

// Перемножение матриц
#pragma region Multiplication

	//Быстрое умножение матриц
	Matrix operator* (Matrix &_m)
	{
		if (m != _m.n)
			throw new exception("Были умножены матрицы неверной размерности");
		Matrix a = Matrix(n, _m.m);
		Matrix M = _m.T();
		for (long i = 0; i < a.n; i++)
		{
			for (long j = 0; j < a.m; j++)
			{
				for (long p = 0; p < m; p++)
				{
					auto a1 = arr[i][p];
					auto a2 = M.arr[j][p];
					a.arr[i][j] += a1*a2;
				}
			}
		}
		return a;
	}

#pragma endregion

	// Сумма элементов матрицы
	Type Sum(Method m = Method::Sequence, int _threads = 1);

	// Сложение матриц
	Matrix<Type> Add(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1);
	
	// Разность матриц
	Matrix<Type> Sub(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1);

	// Умножение элементов матриц
	Type Mul(Method m = Method::Sequence, int _threads = 1);

	// Вычисление максимального элемента матриц
	Type Max(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1);


	friend ostream& operator<<(ostream& stream, Matrix<Type>& M)
	{
		stream << "\n";
		if (M.n == 1 && M.m == 1)
		{
			stream << **M.arr;
			return stream;
		}
		if (M.n == 1)
		{
			stream << "[";
			long i = 0;
			for (i = 0; i < M.m - 1; i++)
				stream << M.arr[0][i] << ", ";
			stream << M.arr[0][i]<< " ]";
			return stream;
		}
		else
		{
			stream << "[[";
			long i = 0;
			long j = 0;
			for (i = 0; i < M.n-1; i++)
			{
				for (j = 0; j < M.m - 1; j++)
					stream << M.arr[i][j] << ", ";
				stream << M.arr[i][j] << " ],\n";
			}
			for (j = 0; j < M.m - 1; j++)
				stream << M.arr[i][j] << ", ";
			stream << M.arr[i][j] << " ]]";
			return stream;
		}

	}

	//Количество строк матрицы
	long N()
	{
		return n;
	}

	//Количество столбцов матрицы
	long M()
	{
		return m;
	}

private:
//Сложение элементов матрицы
#pragma region Summation
	// Последовательное сложение элементов матрицы
	Type SumSeq()
	{
		Type S = 0;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

	// Сложение элементов с параллельным циклом
	Type SumOmpFor()
	{
		Type S = 0;
		#pragma omp parallel for reduction(+:S)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

	// Сложение элементов с параллельными секциями
	Type SumOmpSec(int _threads)
	{
		Type S = 0;
		float step = n / _threads;
		int mod = n % _threads;
		int t1 = 0;
		int t2 = step * 1 + mod;
		int t3 = t2 + step;
		int t4 = t3 + step;
#pragma omp parallel sections reduction(+:S) firstprivate(t1, t2, t3, t4, step)
		{
#pragma omp section
			{
				for (long i = t1; i < t2; i++)
					for (long j = 0; j < m; j++)
						S += arr[i][j];
			}
#pragma omp section
			{
				if (_threads>1)
				for (long i = t2; i < t3; i++)
					for (long j = 0; j < m; j++)
						S += arr[i][j];
			}
#pragma omp section
			{
				if (_threads>2)
				for (long i = t3; i < t4; i++)
					for (long j = 0; j < m; j++)
						S += arr[i][j];
			}
#pragma omp section
			{
				if (_threads>3)
				for (long i = t4; i < t4 + step; i++)
					for (long j = 0; j < m; j++)
						S += arr[i][j];
			}
		}
		return S;
	}
	// Сложение элементов с параллельным циклом на Силке
	Type SumCilkFor()
	{
		//редуктор на сложение равный нулю
		cilk::reducer_opadd<Type> S(0);
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S.get_value();
	}

	// Сложение элементов с векторизацией
	Type SumCilkVec()
	{
		Type S;
		// прохожусь по i = 0:n, по j = 0:m, применяю редукцию на сложение
		S = __sec_reduce_add(arr[0:n][0:m]);
		return S;
	}

	// Сложение элементов со спавном параллельных веток
	Type SumCilkSpawn()
	{
		//редуктор на сложение равный нулю
		cilk::reducer_opadd<Type> S(0);
		for (long i = 0; i < n; i++)
		{
			//Спав новой ветки
			cilk_spawn [&S, i, this]() {
				for (long j = 0; j < m; j++)
					S += arr[i][j];
			}();
		}
		cilk_sync; //синхронизирую ветки здесь
		return S.get_value();
	}
	// Сложение элементов с параллельным циклом и векторизацией
	Type SumCilkForVec()
	{
		cilk::reducer_opadd<Type> S(0);
		cilk_for (long i = 0; i < n; i++)
			S += __sec_reduce_add(arr[i][0:m]);
		return S.get_value();
	}
#pragma endregion

// Сложение двух матриц
#pragma region Addition

	//Обычное сложение матриц
	Matrix<Type> AddSeq (Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}
	// Сложение двух матриц с параллельным циклом
	Matrix<Type> AddOmpFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)	
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		#pragma omp parallel for shared(res, arr, _m)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}

	//Сложение двух матриц с параллельными секциями
	Matrix<Type> AddOmpSec(Matrix<Type> &_m, int _threads)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		float step = n / _threads;
		int mod = n % _threads;
		int t1 = 0;
		int t2 = step * 1 + mod;
		int t3 = t2 + step;
		int t4 = t3 + step;
		#pragma omp parallel sections shared(res, arr, _m) firstprivate(t1, t2, t3, t4, step)
		{
			#pragma omp section
			{
				for (long i = t1; i < t2; i++)
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] + _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>1)
				for (long i = t2; i < t3; i++)
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] + _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>2)
				for (long i = t3; i < t4; i++)
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] + _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>3)
				for (long i = t4; i < t4 + step; i++)
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] + _m.arr[i][j];
			}
		}
		return res;
	}

	// Сложение двух матриц с параллельным циклом на Силке
	Matrix<Type> AddCilkFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}

	// Сложение двух матриц со спавном веток
	Matrix<Type> AddCilkSpawn(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				cilk_spawn [&res, &_m, i, this]() {
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] + _m.arr[i][j];
				}();
		cilk_sync;
		return res;
	}

	// Сложение двух матриц с векторизацией
	Matrix<Type> AddCilkVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		res.arr[0:n][0:m] = arr[0:n][0:m] + _m.arr[0:n][0:m];
		return res;
	}

	// Сложение двух матриц с параллельным циклом и векторизацией
	Matrix<Type> AddCilkForVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for(long i = 0; i < n; i++)
				res.arr[i][0:m] = arr[i][0:m] + _m.arr[i][0:m];
		return res;
	}

#pragma endregion

// Перемножение элементов матрицы
#pragma region Multiplication

	// Последовательное перемножение матриц
	Type MulSeq()
	{
		Type Prod = 1;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				Prod *= arr[i][j];
		return Prod;
	}

	// Перемножение матриц с параллельным циклом
	Type MulOmpFor()
	{
		Type Prod = 1;
		#pragma omp parallel for reduction(*:Prod)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				Prod *= arr[i][j];
		return Prod;
	}

	// Перемножение матриц с параллельными секциями
	Type MulOmpSec(int _threads)
	{
		Type Prod = 1;
		float step = n / _threads;
		int mod = n % _threads;
		int t1 = 0;
		int t2 = step * 1 + mod;
		int t3 = t2 + step;
		int t4 = t3 + step;
		#pragma omp parallel sections reduction(*:Prod) shared(arr) firstprivate(t1, t2, t3, t4, step)
		{
			#pragma omp section
			{
				for (long i = t1; i < t2; i++)
					for (long j = 0; j < m; j++)
						Prod *= arr[i][j];
			}

			#pragma omp section
			{
				if (_threads>1)
				for (long i = t2; i < t3; i++)
					for (long j = 0; j < m; j++)
						Prod *= arr[i][j];
			}

			#pragma omp section
			{
				if (_threads>2)
				for (long i = t3; i < t4; i++)
					for (long j = 0; j < m; j++)
						Prod *= arr[i][j];
			}

			#pragma omp section
			{
				if (_threads>3)
				for (long i = t4; i < t4 + step; i++)
					for (long j = 0; j < m; j++)
						Prod *= arr[i][j];
			}

		}
		return Prod;
	}

	// Перемножение матриц с параллельным циклом на Силке
	Type MulCilkFor()
	{
		cilk::reducer<cilk::op_mul<Type>> Prod(1);
		//Type Prod = 1;
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				*Prod *= arr[i][j];
		return Prod.get_value();
	}

	// Перемножение матриц со спавном новых веток
	Type MulCilkSpawn()
	{
		cilk::reducer<cilk::op_mul<Type>> Prod(1);
		for (long i = 0; i < n; i++)
			cilk_spawn[&Prod, i, this]() {
				for (long j = 0; j < m; j++)
					*Prod *= arr[i][j];
			}();
		cilk_sync;
		return Prod.get_value();
	}

	// Перемножение матриц с векторизацией
	Type MulCilkVec()
	{
		Type Prod = 1;
		Prod = __sec_reduce_mul(arr[0:n][0:m]);
		return Prod;
	}

	// Перемножение матриц с параллельным циклом и векторизацией
	Type MulCilkForVec()
	{
		cilk::reducer<cilk::op_mul<Type>> Prod(1);
		cilk_for(long i = 0; i < n; i++)
			*Prod *= __sec_reduce_mul(arr[i][0:m]);
		return Prod.get_value();
	}
#pragma endregion

// Поиск максимального элемента двух матриц
#pragma region Maximum

	// Нахождение максимального элемента двух матриц
	Type MaxSeq(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		Type m1 = arr[0][0], m2 = M.arr[0][0];
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
			{
					m1 = (m1 < arr[i][j]) ? arr[i][j] : m1;;
					m2 = (m2 < M.arr[i][j]) ? arr[i][j] : m2;
			}
		return MaxOf(m1, m2);
	}

	// Нахождение максимального элемента двух матриц с параллельным циклом
	Type MaxOmpFor(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		Type m1 = arr[0][0], m2 = M.arr[0][0];
		#pragma omp parallel for shared(m1, m2, arr, M)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
			{
				m1 = (m1 < arr[i][j]) ? arr[i][j] : m1;;
				m2 = (m2 < M.arr[i][j]) ? M.arr[i][j] : m2;
			}
		return MaxOf(m1, m2);
	}

	// Нахождение максимального элемента двух матриц с параллельными секциями
	Type MaxOmpSec(Matrix<Type> &M, int _threads)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		Type m1 = arr[0][0], m2 = M.arr[0][0];
		float step = n / _threads;
		int mod = n % _threads;
		int t1 = 0;
		int t2 = step * 1 + mod;
		int t3 = t2 + step;
		int t4 = t3 + step;
		#pragma omp parallel sections reduction(max:m1, m2) shared(arr) firstprivate(t1, t2, t3, t4, step)
		{
			#pragma omp section
			{
				for (long i = t1; i < t2; i++)
					for (long j = 0; j < m; j++)
					{
						m1 = MaxOf(m1, arr[i][j]);
						m2 = MaxOf(m2, M.arr[i][j]);
					}

			}
			#pragma omp section
			{
				if (_threads>1)
				for (long i = t2; i < t3; i++)
					for (long j = 0; j < m; j++)
					{
						m1 = MaxOf(m1, arr[i][j]);
						m2 = MaxOf(m2, M.arr[i][j]);
					}

			}

			#pragma omp section
			{
				if (_threads>2)
				for (long i = t3; i < t4; i++)
					for (long j = 0; j < m; j++)
					{
						m1 = MaxOf(m1, arr[i][j]);
						m2 = MaxOf(m2, M.arr[i][j]);
					}

			}

			#pragma omp section
			{
				if (_threads>3)
				for (long i = t4; i < t4 + step; i++)
					for (long j = 0; j < m; j++)
					{
						m1 = MaxOf(m1, arr[i][j]);
						m2 = MaxOf(m2, M.arr[i][j]);
					}

			}
		}
		return MaxOf(m1, m2);
	}

	// Нахождение максимального элемента двух матриц с параллельным циклом на силке
	Type MaxCilkFor(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		cilk::reducer<cilk::op_max<Type>> m1(arr[0][0]);
		cilk::reducer<cilk::op_max<Type>> m2(M.arr[0][0]);
		cilk_for(long i = 0; i < n; i++)
		{
			for (long j = 0; j < m; j++)
			{
				m1->calc_max(arr[i][j]);
				m2->calc_max(M.arr[i][j]);
			}
		}
		return MaxOf(m1.get_value(), m2.get_value());
	}

	// Нахождение максимального элемента двух матриц со спавном веток
	Type MaxCilkSpawn(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		cilk::reducer<cilk::op_max<Type>> m1(arr[0][0]);
		cilk::reducer<cilk::op_max<Type>> m2(M.arr[0][0]);
		for (long i = 0; i < n; i++)
			cilk_spawn[&m1, &m2, &M, i, this]() {
				for (long j = 0; j < m; j++)
				{
					m1->calc_max(arr[i][j]);
					m2->calc_max(M.arr[i][j]);
				}
			}();
		cilk_sync;
		return MaxOf(m1.get_value(), m2.get_value());
	}
	
	// Нахождение максимального элемента двух матриц с векторизацией
	Type MaxCilkVec(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		Type m1 = arr[0][0];
		Type m2 = M.arr[0][0];

		m1 = __sec_reduce_max(arr[0:n][0:m]);
		m2 = __sec_reduce_max(M.arr[0:n][0:m]);

		return MaxOf(m1, m2);
	}
	
	// Нахождение максимального элемента двух матриц с параллельным циклом и векторизацией
	Type MaxCilkForVec(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("Матрицы разной размерности");
		Type m1 = arr[0][0];
		Type m2 = M.arr[0][0];
		cilk_for(long i = 0; i < n; i++)
		{
			m1 = __sec_reduce_max(arr[i][0:m]);
			m2 = __sec_reduce_max(M.arr[i][0:m]);
		}
		return MaxOf(m1, m2);
	}

#pragma endregion

//
// не знаю зачем написал, засоряю код
//
// Разность двух матриц
#pragma region Substraction
	//Обычное вычитание матриц
	Matrix<Type> SubSeq(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

	Matrix<Type> SubOmpFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		#pragma omp parallel for shared(res, arr, _m)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

	Matrix<Type> SubOmpSec(Matrix<Type> &_m, int _threads)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		float step = n / _threads;
		int mod = n % _threads;
		int t1 = 0;
		int t2 = step * 1 + mod;
		int t3 = t2 + step;
		int t4 = t3 + step;
		#pragma omp parallel sections shared(res, arr, _m) firstprivate(t1, t2, t3, t4, step)
		{
			#pragma omp section
			{
				for (long i = t1; i < t2; i++)
					for (long j = 0; j < m; j++)
						res.arr[i][j] = arr[i][j] - _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>1)
					for (long i = t2; i < t3; i++)
						for (long j = 0; j < m; j++)
							res.arr[i][j] = arr[i][j] - _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>2)
					for (long i = t3; i < t4; i++)
						for (long j = 0; j < m; j++)
							res.arr[i][j] = arr[i][j] - _m.arr[i][j];
			}
			#pragma omp section
			{
				if (_threads>3)
					for (long i = t4; i < t4 + step; i++)
						for (long j = 0; j < m; j++)
							res.arr[i][j] = arr[i][j] - _m.arr[i][j];
			}
		}
		return res;
	}

	Matrix<Type> SubCilkFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for(long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

	Matrix<Type> SubCilkSpawn(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				cilk_spawn[&res, &_m, i, this]() {
				for (long j = 0; j < m; j++)
					res.arr[i][j] = arr[i][j] - _m.arr[i][j];
			}();
			cilk_sync;
			return res;
	}

	Matrix<Type> SubCilkVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		res.arr[0:n][0:m] = arr[0:n][0:m] - _m.arr[0:n][0:m];
		return res;
	}

	Matrix<Type> SubCilkForVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("Недопустимо сложение матриц разной размерности");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for(long i = 0; i < n; i++)
			res.arr[i][0:m] = arr[i][0:m] + _m.arr[i][0:m];
		return res;
	}

#pragma endregion

	void fillMatrix(Filler<Type> func)
	{
		Type t = 0;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
			{
				t = func((i + 1), (j + 1));
				arr[i][j] = t;

			}
	}

	Matrix(long _n, long _m, Type**& arr1)
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = _n;
		m = _m;
		copyArr(arr1, true);
	}

	void copyArr(Type**& arr1,  bool alloc)
	{
		if (alloc)
			arr = new Type*[n];
		for (long i = 0; i < n; i++)
		{
			if (alloc)
				arr[i] = new Type[m];
			//memcpy(arr[i], arr1[i], m * sizeof(double*)*sizeof(double));
			for (long j = 0; j < m; j++)
			{
				//arr[i][j] = new double[1];
				arr[i][j] = arr1[i][j];
			}
		}
	}

	Type** allocateArr(long n, long m)
	{
		Type ** arr1 = new Type*[n];
		for (long i = 0; i < n; i++)
		{
			arr1[i] = new Type[m];
			//for (int j = 0; j < m; j++)
			//	arr1[i][j] = new double[1];
		}
		return arr1;
	}

	long n;
	long m;
	long num; // порядковый номер матрицы
	Type **arr;
};


// Сложение элементов матрицы
#pragma region Summation
template<class Type>
Type Matrix<Type>::Sum(Method m = Method::Sequence, int _threads = 1)
{
	omp_set_num_threads(_threads);
	__cilkrts_set_param("nworkers", doubleToString(_threads).c_str());
	__cilkrts_init();
	switch (m)
	{
	case Method::Sequence:
	{
		return SumSeq();
	}
	case Method::ParallelFor:
	{
		return SumOmpFor();
	}
	case Method::ParallelSections:
	{
		return SumOmpSec(_threads);
	}
	case Method::Cilk_for:
	{
		return SumCilkFor();
	}
	case Method::Cilk_index:
	{
		return SumCilkVec();
	}
	case Method::Cilk_spawn:
	{
		return SumCilkSpawn();
	}
	case Method::Cilk_for_index:
	{
		return SumCilkForVec();
	}
	default:
		break;
	}
}

#pragma endregion

// Сложение двух матриц
#pragma region Addition
template<class Type>
Matrix<Type> Matrix<Type>::Add(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1)
{
	omp_set_num_threads(_threads);
	__cilkrts_set_param("nworkers", doubleToString(_threads).c_str());
	__cilkrts_init();
	switch (m)
	{
	case Method::Sequence:
	{
		return AddSeq(M);
	}
	case Method::ParallelFor:
	{
		return AddOmpFor(M);
	}
	case Method::ParallelSections:
	{
		return AddOmpSec(M, _threads);
	}
	case Method::Cilk_for:
	{
		return AddCilkFor(M);
	}
	case Method::Cilk_index:
	{
		return AddCilkVec(M);
	}
	case Method::Cilk_spawn:
	{
		return AddCilkSpawn(M);
	}
	case Method::Cilk_for_index:
	{
		return AddCilkForVec(M);
	}
	default:
		break;
	}
}

#pragma endregion

// Умножение элементов матрицы
#pragma region Multiplication
template<class Type>
Type Matrix<Type>::Mul(Method m = Method::Sequence, int _threads = 1)
{
	omp_set_num_threads(_threads);
	__cilkrts_set_param("nworkers", doubleToString(_threads).c_str());
	__cilkrts_init();
	switch (m)
	{
	case Method::Sequence:
	{
		return MulSeq();
	}
	case Method::ParallelFor:
	{
		return MulOmpFor();
	}
	case Method::ParallelSections:
	{
		return MulOmpSec(_threads);
	}
	case Method::Cilk_for:
	{
		return MulCilkFor();
	}
	case Method::Cilk_index:
	{
		return MulCilkVec();
	}
	case Method::Cilk_spawn:
	{
		return MulCilkSpawn();
	}
	case Method::Cilk_for_index:
	{
		return MulCilkForVec();
	}
	default:
		break;
	}
}
#pragma endregion

//
// не знаю зачем написал, засоряю код
//
// Разность двух матриц
#pragma region Substraction
template<class Type>
Matrix<Type> Matrix<Type>::Sub(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1)
{
	omp_set_num_threads(_threads);
	__cilkrts_set_param("nworkers", doubleToString(_threads).c_str());
	__cilkrts_init();
	switch (m)
	{
	case Method::Sequence:
	{
		return SubSeq(M);
	}
	case Method::ParallelFor:
	{
		return SubOmpFor(M);
	}
	case Method::ParallelSections:
	{
		return SubOmpSec(M, _threads);
	}
	case Method::Cilk_for:
	{
		return SubCilkFor(M);
	}
	case Method::Cilk_index:
	{
		return SubCilkVec(M);
	}
	case Method::Cilk_spawn:
	{
		return SubCilkSpawn(M);
	}
	case Method::Cilk_for_index:
	{
		return SubCilkForVec(M);
	}
	default:
		break;
	}
}

#pragma endregion

// Поиск максимального элемента
#pragma region Maximum
template<class Type>
Type Matrix<Type>::Max(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1)
{
	omp_set_num_threads(_threads);
	__cilkrts_set_param("nworkers", doubleToString(_threads).c_str());
	__cilkrts_init();
	switch (m)
	{
	case Method::Sequence:
	{
		return MaxSeq(M);
	}
	case Method::ParallelFor:
	{
		return MaxOmpFor(M);
	}
	case Method::ParallelSections:
	{
		return MaxOmpSec(M, _threads);
	}
	case Method::Cilk_for:
	{
		return MaxCilkFor(M);
	}
	case Method::Cilk_index:
	{
		return MaxCilkVec(M);
	}
	case Method::Cilk_spawn:
	{
		return MaxCilkSpawn(M);
	}
	case Method::Cilk_for_index:
	{
		return MaxCilkForVec(M);
	}
	default:
		break;
	}
}
#pragma endregion
