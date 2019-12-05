#pragma once 
#include<iostream>
#include<math.h>
#include<Windows.h>
#include <omp.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
using namespace std;

//��������� [](long, long)->double {return ...;}
typedef double(*fillF)(long, long);

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

enum Method {Sequence, ParallelFor, ParallelSections, Cilk_for, Cilk_spawn, Cilk_index, Cilk_for_index};
//����� �������
template<class Type>
class Matrix
{
public:

	//���������� ��������� ������ �� ����� ���������� ���������
	static long Created;
	
	//���������� ��������� ������ �� ����� ���������� ���������
	static long Destroed;
	
	//��������� ������ �������
	//�������� ������ ��������� ������� � ������� ��������� � �������� � �������� ������ �������
	static bool DEBUG;

#pragma region Constructors


	//���������� ������� ������������ 1�1 � ��������� ������ ����
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

	//���������� ������� ������������ NxM ����������� ������
	Matrix(long _n, long _m)
	{
		num = Created;
		if (DEBUG)
			cout << "C:[" << num << "]";
		Created++;
		n = _n;
		m = _m;
		arr = allocateArr(n, m);
		fillMatrix([](long, long) ->double { return 0; });
	}

	//���������� ����� �������
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

	//���������� ������� ������������ NxM, �������� ������� ��������� ��������� ������� ��������� fillF
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

	//������������ ��� ���������� �������� ������
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
		//cout << "����" << endl;
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

#pragma region Selectors


	//��������� i� ������ �������
	Type* getRow(long _i)
	{
		Type * arr1 = new Type[m];
		for (long j = 0; j < m; j++)
			arr1[j] = arr[_i][j];
		return arr1;
	}

	//��������� j� ������� �������
	Type* getCol(long _j)
	{
		Type * arr1 = new Type[n];
		for (long j = 0; j < n; j++)
			arr1[j] = arr[j][_j];
		return arr1;
	}

	//���������� ������� �������
	Type* getElement(long i, long j)
	{
		if (i >= n || j >= m|| i<0 || j<0)
			throw new exception("������ ������ ������ ����������� �������");
		return &arr[i][j];
	}

	//���������� i� ������ �������, �������������� � ���� �������
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

	//��������� �� ��������� ������
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

#pragma region Multiplication

	//������� ��������� ������
	Matrix operator* (Matrix &_m)
	{
		if (m != _m.n)
			throw new exception("���� �������� ������� �������� �����������");
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

	//����������� ������� ����� ������������ ������
	Matrix<Type> mulAsync(Matrix<Type> &_m)
	{
		if (m != _m.n)
			throw new exception("���� �������� ������� �������� �����������");
		Matrix<Type> a = Matrix<Type>(n, _m.m);
		Matrix<Type> M = _m.T();
		#pragma omp parallel for shared(M)
		for (long i = 0; i < a.n; i++)
		{
			for (long j = 0; j < a.m; j++)
			{
				for (long p = 0; p < m; p++)
				{
					a.arr[i][j] +=arr[i][p] * M.arr[j][p];
				}
			}
		}
		return a;
	}

	//������� ��������� ������
	Matrix<Type> dot (Matrix<Type> &_m)
	{
		if (m != _m.n)
			throw new exception("���� �������� ������� �������� �����������");
		Matrix<Type> a = Matrix<Type>(n, _m.m);//���� ����������
		for (long i = 0; i < a.n; i++)
		{
			for (long j = 0; j < a.m; j++)
			{
				for (long p = 0; p < m; p++)
				{
					auto a1 = arr[i][p];
					auto a2 = _m.arr[p][j];
					a.arr[i][j] += a1*a2;
				}
			}
		}
		return a;
	}

	//����������� ������� ����� ������������ ������
	Matrix<Type> dotAsync(Matrix<Type> &_m)
	{
		if (m != _m.n)
			throw new exception("���� �������� ������� �������� �����������");
		Matrix<Type> a = Matrix<Type>(n, _m.m);//���� ����������
		Type** arr1 = a.arr;
		Type ** arr2 = _m.arr;
		long M = a.m;
		long N = a.n;
		Type a1, a2;
		#pragma omp parallel for private(a1, a2)
		for (long i = 0; i < N; i++)
		{
			for (long j = 0; j < M; j++)
			{
				for (long p = 0; p < m; p++)
				{
					a1 = arr[i][p];
					a2 = arr2[p][j];
					arr1[i][j] += a1*a2;
				}
			}
		}
		return a;
	}

#pragma endregion

#pragma region Summation

	Type Sum(Method m, int _threads)
	{
		switch (m)
		{
		case Method::Sequence:
		{
			return SumSeq();
		}
		case Method::ParallelFor:
		{
			return SumOmpFor(_threads);
		}
		case Method::ParallelSections:
		{
			return SumOmpSec(_threads);
		}
		case Method::Cilk_for:
		{
			return SumCilkFor(_threads);
		}
		case Method::Cilk_index:
		{
			return SumCilkVec(_threads);
		}
		case Method::Cilk_spawn:
		{
			return SumCilkSpawn(_threads);
		}
		case Method::Cilk_for_index:
		{
			return SumCilkForVec(_threads);
		}
		default:
			break;
		}
	}

	//������� �������� ������
	Matrix<Type> operator+ (Matrix<Type> &_m)
	{
		if (m != _m.m|| n!=_m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}

#pragma endregion

#pragma region Substraction

	//������� �������� ������
	Matrix<Type> operator- (Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

#pragma endregion

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

	//���������� ����� �������
	long N()
	{
		return n;
	}

	//���������� �������� �������
	long M()
	{
		return m;
	}

private:
#pragma region Summation

	Type SumSeq()
	{
		Type S = 0;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

	Type SumOmpFor(int _threads)
	{
		Type S = 0;
		omp_set_num_threads(_threads);
		#pragma omp parallel for reduction(+:S)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

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

	Type SumCilkFor(int _threads)
	{
		cilk::reducer_opadd<Type> S(0);
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S.get_value();
	}

	Type SumCilkVec(int _threads)
	{
		Type S;
		S = __sec_reduce_add(arr[0:n][0:m]);
		return S;
	}

	Type SumCilkSpawn(int _threads)
	{
		cilk::reducer_opadd<Type> S(0);
		for (long i = 0; i < n; i++)
		{
			cilk_spawn [&S, i, this]() {
				for (long j = 0; j < m; j++)
					S += arr[i][j];
			}();
		}
		cilk_sync;
		return S.get_value();
	}

	Type SumCilkForVec(int _threads)
	{
		cilk::reducer_opadd<Type> S(0);
		for (long i = 0; i < n; i++)
			S += __sec_reduce_add(arr[i][0:m]);
		return S.get_value();
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
	long num; // ���������� ����� �������
	Type **arr;
};