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

//��������� [](long, long)->double {return ...;}
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

// ����� �������������� ���������
enum Method {
	//���������������� �������
	Sequence,
	//������������ ���� �� ���
	ParallelFor,
	//������������ ������ �� ���
	ParallelSections,
	//������������ ���� �� �����
	Cilk_for,
	//����� ������������ ����� �� �����
	Cilk_spawn,
	//��������� �������(������������) �� �����
	Cilk_index,
	//������������ ���� + ������������ �� �����
	Cilk_for_index};

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

//������������ ������
#pragma region Constructors

	//���������� ������� ������������ NxM ����������� ���������
	static Matrix<Type> Ones(long _n, long _m)
	{
		return Matrix<Type>(_n, _m, Filler<Type>([](long, long)->double {return 1; }));
	}
	
	//���������� ������� ������������ NxM ����������� ������
	static Matrix<Type> Zeros(long _n, long _m)
	{
		return Matrix<Type>(_n, _m, Filler<Type>([](long, long)->double {return 0; }));
	}
	
	//���������� ��������� ������� ������������ NxM
	static Matrix<Type> E(long _n, long _m)
	{
		auto t = min(_n, _m);
		return Matrix<Type>(t, t, Filler<Type>([](long x, long y)->double {return (x == y) ? 1 : 0; }));
	}

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
		fillMatrix(Filler<Type>([](long, long) ->double { return 0; }));
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

// ��������� �� ������� � �������
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

// ������������ ������
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

#pragma endregion

	// ����� ��������� �������
	Type Sum(Method m = Method::Sequence, int _threads = 1);

	// �������� ������
	Matrix<Type> Add(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1);
	
	// �������� ������
	Matrix<Type> Sub(Matrix<Type> &M, Method m = Method::Sequence, int _threads = 1);

	// ��������� ��������� ������
	Type Mul(Method m = Method::Sequence, int _threads = 1);

	// ���������� ������������� �������� ������
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
//�������� ��������� �������
#pragma region Summation
	// ���������������� �������� ��������� �������
	Type SumSeq()
	{
		Type S = 0;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

	// �������� ��������� � ������������ ������
	Type SumOmpFor()
	{
		Type S = 0;
		#pragma omp parallel for reduction(+:S)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S;
	}

	// �������� ��������� � ������������� ��������
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
	// �������� ��������� � ������������ ������ �� �����
	Type SumCilkFor()
	{
		//�������� �� �������� ������ ����
		cilk::reducer_opadd<Type> S(0);
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				S += arr[i][j];
		return S.get_value();
	}

	// �������� ��������� � �������������
	Type SumCilkVec()
	{
		Type S;
		// ��������� �� i = 0:n, �� j = 0:m, �������� �������� �� ��������
		S = __sec_reduce_add(arr[0:n][0:m]);
		return S;
	}

	// �������� ��������� �� ������� ������������ �����
	Type SumCilkSpawn()
	{
		//�������� �� �������� ������ ����
		cilk::reducer_opadd<Type> S(0);
		for (long i = 0; i < n; i++)
		{
			//���� ����� �����
			cilk_spawn [&S, i, this]() {
				for (long j = 0; j < m; j++)
					S += arr[i][j];
			}();
		}
		cilk_sync; //������������� ����� �����
		return S.get_value();
	}
	// �������� ��������� � ������������ ������ � �������������
	Type SumCilkForVec()
	{
		cilk::reducer_opadd<Type> S(0);
		cilk_for (long i = 0; i < n; i++)
			S += __sec_reduce_add(arr[i][0:m]);
		return S.get_value();
	}
#pragma endregion

// �������� ���� ������
#pragma region Addition

	//������� �������� ������
	Matrix<Type> AddSeq (Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}
	// �������� ���� ������ � ������������ ������
	Matrix<Type> AddOmpFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)	
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		#pragma omp parallel for shared(res, arr, _m)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}

	//�������� ���� ������ � ������������� ��������
	Matrix<Type> AddOmpSec(Matrix<Type> &_m, int _threads)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
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

	// �������� ���� ������ � ������������ ������ �� �����
	Matrix<Type> AddCilkFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] + _m.arr[i][j];
		return res;
	}

	// �������� ���� ������ �� ������� �����
	Matrix<Type> AddCilkSpawn(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
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

	// �������� ���� ������ � �������������
	Matrix<Type> AddCilkVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		res.arr[0:n][0:m] = arr[0:n][0:m] + _m.arr[0:n][0:m];
		return res;
	}

	// �������� ���� ������ � ������������ ������ � �������������
	Matrix<Type> AddCilkForVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for(long i = 0; i < n; i++)
				res.arr[i][0:m] = arr[i][0:m] + _m.arr[i][0:m];
		return res;
	}

#pragma endregion

// ������������ ��������� �������
#pragma region Multiplication

	// ���������������� ������������ ������
	Type MulSeq()
	{
		Type Prod = 1;
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				Prod *= arr[i][j];
		return Prod;
	}

	// ������������ ������ � ������������ ������
	Type MulOmpFor()
	{
		Type Prod = 1;
		#pragma omp parallel for reduction(*:Prod)
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				Prod *= arr[i][j];
		return Prod;
	}

	// ������������ ������ � ������������� ��������
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

	// ������������ ������ � ������������ ������ �� �����
	Type MulCilkFor()
	{
		cilk::reducer<cilk::op_mul<Type>> Prod(1);
		//Type Prod = 1;
		cilk_for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				*Prod *= arr[i][j];
		return Prod.get_value();
	}

	// ������������ ������ �� ������� ����� �����
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

	// ������������ ������ � �������������
	Type MulCilkVec()
	{
		Type Prod = 1;
		Prod = __sec_reduce_mul(arr[0:n][0:m]);
		return Prod;
	}

	// ������������ ������ � ������������ ������ � �������������
	Type MulCilkForVec()
	{
		cilk::reducer<cilk::op_mul<Type>> Prod(1);
		cilk_for(long i = 0; i < n; i++)
			*Prod *= __sec_reduce_mul(arr[i][0:m]);
		return Prod.get_value();
	}
#pragma endregion

// ����� ������������� �������� ���� ������
#pragma region Maximum

	// ���������� ������������� �������� ���� ������
	Type MaxSeq(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
		Type m1 = arr[0][0], m2 = M.arr[0][0];
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
			{
					m1 = (m1 < arr[i][j]) ? arr[i][j] : m1;;
					m2 = (m2 < M.arr[i][j]) ? arr[i][j] : m2;
			}
		return MaxOf(m1, m2);
	}

	// ���������� ������������� �������� ���� ������ � ������������ ������
	Type MaxOmpFor(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
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

	// ���������� ������������� �������� ���� ������ � ������������� ��������
	Type MaxOmpSec(Matrix<Type> &M, int _threads)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
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

	// ���������� ������������� �������� ���� ������ � ������������ ������ �� �����
	Type MaxCilkFor(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
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

	// ���������� ������������� �������� ���� ������ �� ������� �����
	Type MaxCilkSpawn(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
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
	
	// ���������� ������������� �������� ���� ������ � �������������
	Type MaxCilkVec(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
		Type m1 = arr[0][0];
		Type m2 = M.arr[0][0];

		m1 = __sec_reduce_max(arr[0:n][0:m]);
		m2 = __sec_reduce_max(M.arr[0:n][0:m]);

		return MaxOf(m1, m2);
	}
	
	// ���������� ������������� �������� ���� ������ � ������������ ������ � �������������
	Type MaxCilkForVec(Matrix<Type> &M)
	{
		if (m != M.m || n != M.n)
			throw new exception("������� ������ �����������");
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
// �� ���� ����� �������, ������� ���
//
// �������� ���� ������
#pragma region Substraction
	//������� ��������� ������
	Matrix<Type> SubSeq(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		for (long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

	Matrix<Type> SubOmpFor(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
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
			throw new exception("����������� �������� ������ ������ �����������");
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
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		cilk_for(long i = 0; i < n; i++)
			for (long j = 0; j < m; j++)
				res.arr[i][j] = arr[i][j] - _m.arr[i][j];
		return res;
	}

	Matrix<Type> SubCilkSpawn(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
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
			throw new exception("����������� �������� ������ ������ �����������");
		Matrix<Type> res = Matrix<Type>(n, m);
		res.arr[0:n][0:m] = arr[0:n][0:m] - _m.arr[0:n][0:m];
		return res;
	}

	Matrix<Type> SubCilkForVec(Matrix<Type> &_m)
	{
		if (m != _m.m || n != _m.n)
			throw new exception("����������� �������� ������ ������ �����������");
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
	long num; // ���������� ����� �������
	Type **arr;
};


// �������� ��������� �������
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

// �������� ���� ������
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

// ��������� ��������� �������
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
// �� ���� ����� �������, ������� ���
//
// �������� ���� ������
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

// ����� ������������� ��������
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
