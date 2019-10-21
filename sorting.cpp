#pragma once
#include<omp.h>
//#include<Windows.h>

//template<typename Type>
typedef bool(*IntComparer)(int, int);
typedef bool(*DoubleComparer)(double, double);

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

template<class T>
void Swap(T& one, T& two)
{
	T temp = one;
	one = two;
	two = temp;
}

template<class T, typename TComparer>
T* Bubble(T* arr, long length, TComparer copare)
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

template<class T, typename TComparer>
T* BubbleEvenAsync(T*& arr, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i = 0; !sorted; i++)
	{
		sorted = true;
		#pragma omp for shared(i, sorted, narr) 
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (compare(narr[j], narr[j + 1]))
			{
				Swap(narr[j], narr[j + 1]);
				sorted = false;
			}
	}
	return narr;
}

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

template<class T>
T* quickSortR(T* &arr, long length) {
	// На входе - массив a[], a[N] - его последний элемент.

	long i = 0, j = length - 1; 		// поставить указатели на исходные места
	T* narr = new T[length];
	memcpy(narr, array, sizeof(T)*length);
	T temp, p;

	p = narr[length >> 1];		// центральный элемент
						// процедура разделения
	do {
		while (narr[i] < p) i++;
		while (narr[j] > p) j--;

		if (i <= j) {
			Swap(narr[i], narr[j]);
			i++; j--;
		}
	} while (i <= j);

	// рекурсивные вызовы, если есть, что сортировать 
	if (j > 0) quickSortR(narr, j);
	if (length > i) quickSortR(narr + i, length - i);
}

template<class T>
T* quickSortAsync(T* arr, long length) {
	// На входе - массив a[], a[N] - его последний элемент.

	long i = 0, j = length - 1; 		// поставить указатели на исходные места
	T* narr = new T[length];
	memcpy(narr, array, sizeof(T)*length);
	T temp, p;
	p = arr[length >> 1];		// центральный элемент
								// процедура разделения
	for (int partLen = length / 2; partLen > 1; partLen /= 2) // дроблю массив на подмассивы
	{//parallel for
		int count = length / partLen;
		for (int part = 0; part < count; part++)
		{
			p = partLen / 2 + part*partLen;
			do {
				while (narr[i] < p) i++;
				while (narr[j] > p) j--;

				if (i <= j) {
					Swap(narr[i], narr[j]);
					i++; j--;
				}
			} while (i <= j);
		}
	}
	// рекурсивные вызовы, если есть, что сортировать 
	if (j > 0) quickSortR(arr, j);
	if (length > i) quickSortR(arr + i, length - i);
}