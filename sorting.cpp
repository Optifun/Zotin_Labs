#pragma once
#include<omp.h>
#include<Windows.h>

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

template<class T>
T* Bubble(T* arr, long length)
{
	bool sorted = false;
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	while (!sorted)
	{
		sorted = true;
		for (long i = 0; i < length-1; i++)
		{
			if (narr[i] > narr[i + 1])
			{
				Swap(narr[i], narr[i + 1]);
				sorted = false;
			}
		}
	}
	return narr;

}

//Сортировка чет-нечет
template<class T>
T* BubbleEven(T* arr, long length=-1)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i=0; !sorted; i++)
	{
		sorted = true;
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (narr[j] > narr[j + 1])
			{
				Swap(narr[j], narr[j + 1]);
				sorted = false;
			}
	}
	return narr;
}

template<class T>
T* BubbleEvenAsync(T*& arr, long length)
{
	T* narr = new T[length];
	memcpy(narr, arr, sizeof(T)*length);
	bool sorted = false;
	for (long i = 0; !sorted; i++)
	{
		sorted = true;
		#pragma omp for shared(i, sorted, narr) 
		for (long j = (i + 1) % 2; j < length - 1; j += 2)
			if (narr[j] > narr[j + 1])
			{
				Swap(narr[j], narr[j + 1]);
				sorted = false;
			}
	}
	return narr;
}