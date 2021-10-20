#include <algorithm>
#include <iostream>
#include <fstream>
#include <thread>
#include <random>
#include <memory>
#include <numeric>
#include <ctime>

#include <glm/glm.hpp>

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
void writeCSV(const std::vector<uint32_t>& counts, const std::vector<double>& averages, const string& filename);
void writeABC(double a, double b, double c, const std::string& filename);
string createDateFilename();

double n2logn2(double n);
double n2logn(double n);
double nlogn(double n);
double n2(double n);


uint32_t experimentsCount = 10;
uint32_t startSize = 100;
uint32_t endSize = 10000;
uint32_t step = 100;


void proceed_args(int argc, char* argv[], string& filenameRes, string& filenameABC)
{
	char* startSizeArg = argv[1];
	char* endSizeArg = argv[2];
	char* stepArg = argv[3];
	try
	{
		char* endptr;
		startSize = strtoul(startSizeArg, &endptr, 10);
		endSize = strtoul(endSizeArg, &endptr, 10);
		step = strtoul(stepArg, &endptr, 10);
		if (argc > 4)
		{
			char* expCountArg = argv[4];
			experimentsCount = strtoul(expCountArg, &endptr, 10);
		}
		
		if(argc > 5)
			filenameRes = string(argv[5]);
		if (argc > 6)
			filenameABC = string(argv[6]);
			
	}
	catch(std::range_error& ex){
		cout << "Error occurred: " <<ex.what() << endl;
	}
	
}

int main(int argc, char* argv[])
{
	string filenameResults = createDateFilename();
	string filenameABC = "abc-" + createDateFilename();
	if(argc < 4)
	{
		cout << "Need more arguments. Expected 3, given " << argc - 1 << endl;
		return -1;
	}
	proceed_args(argc, argv, filenameResults, filenameABC);

	if(endSize <= startSize)
	{
		cout << "End size cannot be less or equal to start size" << endl;
		return -1;
	}
	if(step > endSize - startSize)
	{
		cout << "Step cannot be greater than difference between start size and end size" << endl;
		return -1;
	}
	

	const uint32_t rows = endSize / step;

	std::vector<uint32_t> counts(rows);
	std::vector<std::vector<double>*> results(rows);
	std::vector<double> averages(rows);
	for (size_t i = 0; i < rows; i++)
	{
		results[i] = new std::vector<double>(experimentsCount);
	}
	for (uint32_t i = 0, j = startSize; j <= endSize; j += step, i++)
	{
		counts[i] = j;
	}

	for (uint32_t i = 0; i < experimentsCount; ++i)
	{
		Timer iterationTimer;
		cout << "Experiment #" << i + 1 << ". ";
		iterationTimer.start();
		for (size_t j = startSize, k = 0; j <= endSize; j += step, k++)
		{
			Timer timer;
			auto arr = getRandomArray(j);
			timer.start();
			heapSort(arr.get(), j);
			timer.end();
			(*results[k])[i] = timer.getDuration();
		}
		iterationTimer.end();
		cout << "Elapsed time: " << iterationTimer.getDuration() << endl;
	}

	for (size_t i = 0; i < rows; i++)
	{
		averages[i] = std::accumulate(results[i]->begin(), results[i]->end(), 0.0);
		averages[i] /= experimentsCount;
	}



	glm::mat<2, 2, double, glm::qualifier::defaultp> mat(0.0);
	glm::mat<2, 2, double, glm::qualifier::defaultp> mat1(0.0);
	glm::mat<2, 2, double, glm::qualifier::defaultp> mat2(0.0);


	for (size_t i = 0; i < results.size(); i++)
	{
		mat[0][0] += n2logn2(counts[i]);
		mat[0][1] += nlogn(counts[i]);
		mat[1][0] += nlogn(counts[i]);
		mat[1][1] += 1;

		const double t1 = averages[i] * nlogn(counts[i]);
		const double t2 = averages[i];

		mat1[0][0] += t1;
		mat1[0][1] += t2;
		mat1[1][0] += nlogn(counts[i]);
		mat1[1][1] += 1;

		mat2[0][0] += n2logn2(counts[i]);
		mat2[0][1] += nlogn(counts[i]);
		mat2[1][0] += t1;
		mat2[1][1] += t2;
	}

	const double det = determinant(mat);
	const double det1 = determinant(mat1);
	const double det2 = determinant(mat2);

	const double a = det1 / det;
	const double b = det2 / det;

	glm::mat<3, 3, double, glm::qualifier::defaultp> mat3_0(0.0);
	glm::mat<3, 3, double, glm::qualifier::defaultp> mat3_1(0.0);
	glm::mat<3, 3, double, glm::qualifier::defaultp> mat3_2(0.0);
	glm::mat<3, 3, double, glm::qualifier::defaultp> mat3_3(0.0);

	for (size_t i = 0; i < results.size(); i++)
	{
		mat3_0[0][0] += n2logn2(counts[i]);
		mat3_0[0][1] += n2logn(counts[i]);
		mat3_0[0][2] += nlogn(counts[i]);
		mat3_0[1][0] += n2logn(counts[i]);
		mat3_0[1][1] += n2(counts[i]);
		mat3_0[1][2] += counts[i];
		mat3_0[2][0] += nlogn(counts[i]);
		mat3_0[2][1] += counts[i];
		mat3_0[2][2] += 1;

		const double t1 = averages[i] * nlogn(counts[i]);
		const double t2 = averages[i] * counts[i];
		const double t3 = averages[i];
		
		mat3_1[0][0] += t1;
		mat3_1[0][1] += t2;
		mat3_1[0][2] += t3;
		mat3_1[1][0] += n2logn(counts[i]);
		mat3_1[1][1] += n2(counts[i]);
		mat3_1[1][2] += counts[i];
		mat3_1[2][0] += nlogn(counts[i]);
		mat3_1[2][1] += counts[i];
		mat3_1[2][2] += 1;

		mat3_2[0][0] += n2logn2(counts[i]);
		mat3_2[0][1] += n2logn(counts[i]);
		mat3_2[0][2] += nlogn(counts[i]);
		mat3_2[1][0] += t1;
		mat3_2[1][1] += t2;
		mat3_2[1][2] += t3;
		mat3_2[2][0] += nlogn(counts[i]);
		mat3_2[2][1] += counts[i];
		mat3_2[2][2] += 1;

		mat3_3[0][0] += n2logn2(counts[i]);
		mat3_3[0][1] += n2logn(counts[i]);
		mat3_3[0][2] += nlogn(counts[i]);
		mat3_3[1][0] += n2logn(counts[i]);
		mat3_3[1][1] += n2(counts[i]);
		mat3_3[1][2] += counts[i];
		mat3_3[2][0] += t1;
		mat3_3[2][1] += t2;
		mat3_3[2][2] += t3;
	}

	const double det3_0 = determinant(mat3_0);
	const double det3_1 = determinant(mat3_1);
	const double det3_2 = determinant(mat3_2);
	const double det3_3 = determinant(mat3_3);

	const double A = det3_1 / det3_0;
	const double B = det3_2 / det3_0;
	const double C = det3_3 / det3_0;

	writeCSV(counts, averages, filenameResults);
	writeABC(a, b, 0, filenameABC);
	writeABC(A, B, C, filenameABC+"_c");
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
	std::random_device dev;
	static thread_local std::mt19937 generator(dev());
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


void writeCSV(const std::vector<uint32_t>& counts, const std::vector<double>& averages, const string& filename)
{
	try
	{
		ofstream file(filename);
		for (size_t i = 0; i < counts.size(); ++i)
		{
			file << counts[i] << ", " << averages[i] << endl;
		}
		file.close();
	}
	catch(const std::exception& ex)
	{
		cout << ex.what() << endl;
	}
	
}

void writeABC(double a, double b, double c, const string& filename)
{
	try
	{
		ofstream file(filename);
		file << a << ", " << b << ", " << c;
		file.close();  
	}
	catch (const std::exception& ex)
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

inline double n2logn2(double n)
{
	return (n * n * log2(n) * log2(n));
}

inline double n2logn(double n)
{
	return n * n * log2(n);
}

inline double nlogn(double n)
{
	return (n * log2(n));
}

inline double n2(double n)
{
	return n * n;
}


