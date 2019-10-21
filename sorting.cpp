#pragma once
#include<omp.h>
#include<Windows.h>

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
T* ShellSort(T *array, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, array, sizeof(T)*length);
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
T* ShellSortAsync(T *array, long length, TComparer compare)
{
	T* narr = new T[length];
	memcpy(narr, array, sizeof(T)*length);
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