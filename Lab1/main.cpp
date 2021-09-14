#include <algorithm>
#include <iostream>
#include <fstream>
#include <thread>
#include <random>
#include <memory>
#include <numeric>
#include <ctime>

#include "Timer.h"

using std::cout;
using std::endl;
using std::swap;
using std::ifstream;
using std::ofstream;
using std::string;

template <typename T> void heapify(T arr[], size_t n, size_t index);
template <typename T> void buildHeap(T arr[], int n);
template <typename T> void heapSort(T arr[], int n);
template <typename T> void printArray(T arr[], int n);
int intRand(const int& min, const int& max);
std::shared_ptr<int> getRandomArray(size_t size);
void writeCSV(const std::vector<uint32_t>& counts, const std::vector<std::vector<double>*>& results, const std::vector<double>& averages, const string& filename);
string createDateFilename();

int main()
{
	const int experimentsCount = 10;
	const uint32_t startSize = 1000;
	const uint32_t endSize = 100000;
	const uint32_t step = 1000;

	const uint32_t rows = endSize / step;

	std::vector<uint32_t> counts(rows);
	std::vector<std::vector<double>*> results(rows);
	std::vector<double> averages;
	for (size_t i = 0; i < rows; i++)
	{
		results[i] = new std::vector<double>(experimentsCount);
	}
	for (uint32_t i = 0, j = startSize; j <= endSize; j += step, i++)
	{
		counts[i] = j;
	}

	for (int i = 0; i < experimentsCount; ++i)
	{
		cout << "Experiment #" << i+1 << endl;
		for (size_t j = startSize, k = 0; j <= endSize; j+=step, k++)
		{
			Timer timer;
			auto arr = getRandomArray(j);
			timer.start();
			heapSort(arr.get(), j);
			timer.end();
			//cout << '\t' << j << " elements. " << "Elapsed time: " << timer.getDuration() << "s" << endl;
			(*results[k])[i] = timer.getDuration();
		}
	}

	for (size_t i = 0; i < rows; i++)
	{
		averages[i] = std::accumulate(results[i]->begin(), results[i]->end(), 0.0);
		averages[i] /= experimentsCount;
	}

	writeCSV(counts, results, averages, "../../Lab1PythonPart/" + createDateFilename());
	return 0;
}


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
		swap(arr[index], arr[largest]);
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


template <typename T>
void printArray(T arr[], const int n)
{
	try
	{
		cout << "[";
		for (int i = 0; i < n - 1; ++i)
		{
			cout << arr[i] << ", ";
		}
		cout << arr[n - 1] << "]" << endl;
	}
	catch(const std::exception& ex)
	{
		cout << ex.what() << endl;
	}
}

int intRand(const int& min, const int& max) {
	static thread_local std::mt19937 generator(static_cast<int32_t>(time(nullptr)));
	const std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}

std::shared_ptr<int> getRandomArray(size_t  size)
{
	int* result = new int[size];
	for (size_t i = 0; i < size; ++i)
	{
		result[i] = intRand(INT_MIN, INT_MAX);
	}
	return std::shared_ptr<int>(result);
}

void writeCSV(const std::vector<uint32_t>& counts, const std::vector<std::vector<double>*>& results,
	const std::vector<double>& averages, const string& filename)
{
	try
	{
		ofstream file(filename);
		const size_t cols = results[0]->size();
		for (size_t i = 0; i < counts.size(); ++i)
		{
			file << counts[i] << ",";
			for (size_t j = 0; j < cols; ++j)
			{
				file << (*results[i])[j] << ",";
			}
			file << averages[i] << endl;
		}
		file.close();
	}
	catch(const std::exception& ex)
	{
		cout << ex.what() << endl;
	}
	
}

string createDateFilename()
{
	time_t t = time(nullptr);

	tm localTime{};
	localtime_s(&localTime, &t);

	return std::to_string(localTime.tm_mday) + "_" +
		std::to_string(localTime.tm_mon) + "_" +
		std::to_string(localTime.tm_year - 100) + "_" +
		std::to_string(localTime.tm_hour) + "_" +
		std::to_string(localTime.tm_min) + ".txt";
}
