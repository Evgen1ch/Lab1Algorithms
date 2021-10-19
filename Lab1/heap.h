#pragma once


inline size_t left(size_t i)
{
	return i * 2 + 1;
}

inline size_t right(size_t i)
{
	return i * 2 + 2;
}

template <typename T>
void heapify(T arr[], const size_t n, const size_t index)
{
	size_t l = left(index);
	size_t r = right(index);

	size_t largest = index;
	if (l < n && arr[l] > arr[index])
		largest = l;
	if (r < n && arr[r] > arr[largest])
		largest = r;
	if (largest != index)
	{
		std::swap(arr[index], arr[largest]);
		heapify(arr, n, largest);
	}
}

template <typename T>
void buildHeap(T arr[], const int n)
{
	for (int i = n / 2 - 1; i >= 0; i--)
		heapify(arr, n, i);
}

template <typename T>
void heapSort(T arr[], const int n)
{
	buildHeap(arr, n);
	for (int i = n - 1; i >= 0; --i)
	{
		swap(arr[0], arr[i]);
		heapify(arr, i, 0);
	}
}
