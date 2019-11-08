#pragma once
#include<omp.h>
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

//Сортировка пузырьком
template<class T, typename TComparer>
T* Bubble(T* arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i = 0; i<length; i++)
	{
		//sorted = true;
		for (long j = 0; j < length - i - 1; j ++)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				//sorted = false;
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
	for (long i=0; i<length; i++)
	{
		//sorted = true;
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				//sorted = false;
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
	for (long i = 0; i<length; i++)
	{
		//sorted = true;
		#pragma omp parallel for shared(i, sorted, narr,compare) 
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				//sorted = false;
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
	long step, i, j;

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
	long step, i, j;

	for (step = length / 2; step > 0; step /= 2)
		// Перечисление элементов, которые сортируются на определённом шаге
		#pragma omp parallel for shared(narr, compare, step) private(j)
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
	if (j > 0)
		quickSortR(a, j + 1);
	if (N > i)
		quickSortR(a + i, N - i);
}

//Быстрая сортировка распараллеленая секциями
template<class T>
void quickSortAsync2(T* a, long N, int threads) {

	long i = 0, j = N - 1;
	T temp, p;
	p = a[N >> 1];
	if (threads < 1)
		threads = 1;
	do {
	#pragma omp parallel sections shared(a, i, j) num_threads(threads)
		{
		#pragma omp section
			{ while (a[i] < p) i++; }
		#pragma omp section
			{ while (a[j] > p) j--; }

		}
		if (i <= j) {
			Swap(a[i], a[j]);
			i++; j--;
		}
	} while (i < j);

	#pragma omp parallel sections shared(j, i, a) num_threads(threads)
	{
		#pragma omp section
		{
			if (j > 0) quickSortAsync2(a, j + 1, threads-1);
		}
		#pragma omp section
		{
			if (N > i) quickSortAsync2(a + i, N - i, threads-2);
		}
	}
}

////Быстрая сортировка для 1-4 потоков
//template<class T>
//void quickSortOneThrd(T* a, long N) {
//
//	long i = 0, j = N - 1;
//	T temp, p;
//	p = a[N >> 1];
//	do {
//		while (a[i] < p) i++;
//		while (a[j] > p) j--;
//
//		if (i <= j) {
//			Swap(a[i], a[j]);
//			i++; j--;
//		}
//	} while (i < j);
//
//	if (j > 0) quickSortOneThrd(a, j + 1);
//	if (N > i) quickSortOneThrd(a + i, N - i);
//}
//
//template<class T>
//void quickSortTwoThrd(T* a, long N) {
//
//	long i = 0, j = N - 1;
//	T temp, p;
//	p = a[N >> 1];
//	do {
//		while (a[i] < p) i++;
//		while (a[j] > p) j--;
//
//		if (i <= j) {
//			temp = a[i]; a[i] = a[j]; a[j] = temp;
//			i++; j--;
//		}
//	} while (i < j);
//
//#pragma omp parallel sections
//	{
//#pragma omp section
//		{
//			if (j > 0) quickSortOneThrd(a, j + 1);
//		}
//#pragma omp section
//		{
//			if (N > i) quickSortOneThrd(a + i, N - i);
//		}
//	}
//}
//
//template<class T>
//void quickSortFourThrd(T* a, long N) {
//
//	long i = 0, j = N - 1;
//	T temp, p;
//	p = a[N >> 1];
//	do {
//		while (a[i] < p) i++;
//		while (a[j] > p) j--;
//
//		if (i <= j) {
//			temp = a[i]; a[i] = a[j]; a[j] = temp;
//			i++; j--;
//		}
//	} while (i < j);
//
//	T* a1 = a;
//	long N1 = j + 1;
//	long i1 = 0, j1 = N1 - 1;
//	if (j > 0) {
//		T temp1, p1;
//		p1 = a1[N1 >> 1];
//		do {
//			while (a1[i1] < p1) i1++;
//			while (a1[j1] > p1) j1--;
//
//			if (i1 <= j1) {
//				temp1 = a1[i1]; a1[i1] = a1[j1]; a1[j1] = temp1;
//				i1++; j1--;
//			}
//		} while (i1 < j1);
//	}
//
//	T* a2 = a + i;
//	long N2 = N - i;
//	long i2 = 0, j2 = N2 - 1;
//	if (N > i) {
//		T temp2, p2;
//		p2 = a2[N2 >> 1];
//		do {
//			while (a2[i2] < p2) i2++;
//			while (a2[j2] > p2) j2--;
//
//			if (i2 <= j2) {
//				temp2 = a2[i2]; a2[i2] = a2[j2]; a2[j2] = temp2;
//				i2++; j2--;
//			}
//		} while (i2 < j2);
//	}
//
//#pragma omp parallel sections
//	{
//#pragma omp section 
//		{
//			if (j1 > 0 && j > 0) quickSortOneThrd<T>(a1, j1 + 1);
//		}
//#pragma omp section
//		{
//			if (N1 > i1 && j > 0) quickSortOneThrd<T>(a1 + i1, N1 - i1);
//		}
//#pragma omp section 
//		{
//			if (j2 > 0 && N > i) quickSortOneThrd<T>(a2, j2 + 1);
//		}
//#pragma omp section
//		{
//			if (N2 > i && N > i) quickSortOneThrd<T>(a2 + i2, N2 - i2);
//		}
//	}
//}
//
////Быстрая сортировка (main function)
//template<class T>
//void quickSortAsync(T* a, long N) {
//	switch (omp_get_max_threads())
//	{
//	case 1:
//		quickSortOneThrd<T>(a, N);
//		break;
//	case 2:
//		quickSortTwoThrd<T>(a, N);
//		break;
//	case 3:
//		quickSortTwoThrd<T>(a, N);
//		break;
//	default:
//		quickSortFourThrd<T>(a, N);
//		break;
//	}
//}