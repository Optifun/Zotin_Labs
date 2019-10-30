﻿#pragma once
#include<omp.h>
#include<queue>
typedef bool(*IntComparer)(int, int);
typedef bool(*DoubleComparer)(double, double);

typedef int*(*SortingInt)(int* arr, long length, IntComparer compare);//Алгоритм сортировки целых чисел
typedef double*(*SortingDouble)(double* arr, long length, DoubleComparer compare);//Алгоритм сортировки чисел с плавающей точкой

const IntComparer AscInt = [](int one, int two)->bool { return one > two; }; //сортирует по возрастанию
const IntComparer DescInt = [](int one, int two)->bool { return one < two; }; //сортирует по убыванию
const DoubleComparer AscDouble = [](double one, double two)->bool { return one > two; }; //сортирует по возрастанию
const DoubleComparer DescDouble = [](double one, double two)->bool { return one < two; }; //сортирует по убыванию

#define Random(a, b) (double)rand() / (RAND_MAX + 1)*((b)-(a)) + (a)

//Не советую трогать, изучать. Лучше не стоит..
namespace magic
{
template<class T>
class Comparator
{
public:
	bool operator()(const T& one, const T& two)
	{
		return ASC == order_ ? (two < one) : (one < two);
	}
	//по возрастанию
	static const Comparator& Asc()
	{
		return asc_;
	}
	//по убыванию
	static const Comparator& Desc()
	{
		return desc_;
	}
private:
	enum Order_ { ASC, DESC } const order_;
	Comparator(Order_ order) : order_(order) {}; // первичная инициализация
	static Comparator asc_;
	static Comparator desc_;
};
template<class T> Comparator<T> Comparator<T>::asc_(Comparator<T>::ASC); //инициализация
template<class T> Comparator<T> Comparator<T>::desc_(Comparator<T>::DESC);
}

//Проверяет последовательность на упорядоченноть
template<class T>
bool isOrdered(T* arr, long length, IntComparer order)
{
	bool ord = true;
	for (int i = 0; i < length - 1; i++)
		if (order(arr[i], arr[i + 1]))
			return false;
	return true;
}

//Заполняет последовательность случайными числами
template<class T>
void Randomize(T *arr, long length)
{
	for (long i = 0; i < length; i++)
		arr[i] = Random(-50, 50);
}

//Находит максимальный элемент в последовательности аргументов
template<class T>
T Max(long count, T args, ...)
{
	T* arg = &args;
	T _max = args;
	for (; count; count--, arg++)
	{
		if (*arg > _max)
			_max = *arg;
	}
	return _max;
}

//Обмен значениями
template<class T>
void Swap(T& one, T& two)
{
	T temp = one;
	one = two;
	two = temp;
}

//Сортировка обменами
template<class T, typename TComparer>
T* Bubble(T* arr, long length, TComparer compare)
{
	bool sorted = false;
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	while (!sorted)
	{
		sorted = true;
		for (long i = 0; i < length-1; i++)
		{
			if (compare(narr[i], narr[i + 1]))
			{
				Swap(narr[i], narr[i + 1]);
				sorted = false;
			}
		}
	}
	return narr;

}

//Сортировка чет-нечет
template<class T, typename TComparer>
T* BubbleEven(T* arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i=0; !sorted; i++)
	{
		sorted = true;
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				sorted = false;
			}
	}
	return narr;
}

//Сортировка чет-нечет параллельная
template<class T, typename TComparer>
T* BubbleEvenAsync(T* arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i = 0; !sorted; i++)
	{
		sorted = true;
		#pragma omp parallel for shared(i, sorted, narr) 
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				sorted = false;
			}
	}
	return narr;
}

//Сортировка Шелла
template<class T, typename TComparer>
T* ShellSort(T *arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	long step, i, j, tmp;

	for (step = length / 2; step > 0; step /= 2)
		// Перечисление элементов, которые сортируются на определённом шаге
		for (i = step; i < length; i++)
			// Перестановка элементов внутри подсписка, пока i-тый не будет отсортирован
			for (j = i - step; j >= 0 && compare(narr[j] , narr[j + step]); j -= step)
				Swap(narr[j], narr[j + step]);
	return narr;
}

//Сортировка Шелла параллельная
template<class T, typename TComparer>
T* ShellSortAsync(T *arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	long step, i, j, tmp;

	for (step = length / 2; step > 0; step /= 2)
		// Перечисление элементов, которые сортируются на определённом шаге
		#pragma omp parallel for shared(narr)
		for (i = step; i < length; i++)
			// Перестановка элементов внутри подсписка, пока i-тый не будет отсортирован
			for (j = i - step; j >= 0 && compare(narr[j], narr[j + step]); j -= step)
				Swap(narr[j], narr[j + step]);
	return narr;
}

//Быстрая сортировка
template<class T>
void quickSortR(T* a, long N) {
	long i = 0, j = N - 1;
	T  p;
	p = a[N >> 1];
	do {
		while (a[i] < p) i++;
		while (a[j] > p) j--;

		if (i <= j) {
			Swap(a[i], a[j]);
			i++; j--;
		}
	} while (i < j);
	//if (j > 0) cout << "( 0 : " << j << " )";
	//if (N > i) cout << "( " << i << " : " << N-1 << " )";
	//cout << endl;
	if (j > 0)
		quickSortR(a, j + 1);
	if (N > i)
		quickSortR(a + i, N - i);
}

//template<class T>
int* quickSortAsync(int* arr, long length) {
	int* narr = new int[length];
	memcpy(narr, arr, sizeof(int)*length);

	// очередь границ интервалов
	std::queue<int> iterators = std::queue<int>();
	long i = 0, j = length - 1;
	int p; // опорный элемент
	long b = 0, e = length-1; //начало и конец массива(включительно)
	iterators.push(b);
	iterators.push(e);
	for (int partLen = length; partLen > 0; partLen /= 2) // дроблю массив на подмассивы
	{
		
		int count = length / partLen; //count = [1, 2, ..., length]
		//#pragma omp parallel for shared(narr, odd) firstprivate(i, j, part, partLen)
		for (int part = 0; part < count; part++)
		{
			// достаю из очереди начало и конец массива
			b = i = iterators.front(); iterators.pop();
			e = j = iterators.front(); iterators.pop();
			if (b >= e)
				continue;
			p = narr[(e + b + 1) / 2]; // выбираю опорный элемент

			//вывожу массив
			for (int _i=b; _i <= e; _i++)
				std::cout << narr[_i] << " ";
			std::cout << std::endl;
			std::cout << "i = " << i << " j = " << j << std::endl;
			std::cout <<"p = " << p << std::endl;

			//стандартный алгоритм
			do {
				while (narr[i] < p) i++;
				while (narr[j] > p) j--;

				if (i <= j) {
					std::cout << narr[i] << "<->" << narr[j] << std::endl;
					Swap(narr[i], narr[j]);
					for (int _i = b; _i <= e; _i++)
						std::cout << narr[_i] << " ";
					std::cout << std::endl;
					i++; j--;
				}
			} while (i < j);


			if (j > b) // если j не достигнула начала
			{
				iterators.push(b); // i = b
				std::cout << "(" << b << " : ";
				auto t = e - j; // избыток, чтобы не уйти за границы массива
				if (t < 0)
				{
					iterators.push(j); // j = j
					std::cout << j + t << ")";
				}
				else
				{
					iterators.push(j);
					std::cout << j << ")";
				}

			}
			else // иначе создаю пустой интервалл
			{
				std::cout << "(" << b << " : ";
				std::cout << b - 1 << ")";
				iterators.push(b);
				iterators.push(b-1);
			}

			if (e >= i) //если i не достигла конца
			{
				//возможно прибавлять b это очень лишнее действие
				auto t = e - (b + i + 1);
				if (t < 0)
				{
					iterators.push(i +1 + t);
					std::cout << "(" << b + i +1 + t << " : ";
				}
				else
				{
					iterators.push(b + i +1);
					std::cout << "(" << b + i +1 << " : ";
				}
				iterators.push(e);
				std::cout << e << ")\n";
			}
			else
			{
				std::cout << "(" << b+i << " : ";
				std::cout << b+i - 1 << ")\n";
				iterators.push(i);
				iterators.push(i-1);
			}
		}
		std::cout << std::endl;
	}
	return narr;
}