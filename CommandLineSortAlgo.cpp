#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <time.h>
#include <chrono>
#include <math.h>
#include "SortAlgorithms.h"
using namespace std;
using namespace std::chrono;

template <class T>
void HoanVi(T& a, T& b)
{
	T x = a;
	a = b;
	b = x;
}

//-------------------------------------------------

// Hàm phát sinh mảng dữ liệu ngẫu nhiên
void GenerateRandomData(int a[], int n)
{
	srand((unsigned int)time(NULL));

	for (int i = 0; i < n; i++)
	{
		a[i] = rand() % n;
	}
}

// Hàm phát sinh mảng dữ liệu có thứ tự tăng dần
void GenerateSortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}
}

// Hàm phát sinh mảng dữ liệu có thứ tự ngược (giảm dần)
void GenerateReverseData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = n - 1 - i;
	}
}

// Hàm phát sinh mảng dữ liệu gần như có thứ tự
void GenerateNearlySortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}
	srand((unsigned int)time(NULL));
	for (int i = 0; i < 10; i++)
	{
		int r1 = rand() % n;
		int r2 = rand() % n;
		HoanVi(a[r1], a[r2]);
	}
}

// Decide data type 
int DataType(char* order) {
	if (_strcmpi(order, "-rand") == 0)
		return 0;
	if (_strcmpi(order, "-nsorted") == 0)
		return 3;
	if (_strcmpi(order, "-sorted") == 0)
		return 1;
	if (_strcmpi(order, "-rev") == 0)
		return 2;
	else
		return 4; // unknown data type!
}

void GenerateData(int a[], int n, int dataType)
{
	switch (dataType)
	{
	case 0:	// ngẫu nhiên
		GenerateRandomData(a, n);
		break;
	case 1:	// có thứ tự
		GenerateSortedData(a, n);
		break;
	case 2:	// có thứ tự ngược
		GenerateReverseData(a, n);
		break;
	case 3:	// gần như có thứ tự
		GenerateNearlySortedData(a, n);
		break;
	default:
		printf("Error: unknown data type!\n");
		break;
	}
}


//Read file for command 1 and command 2
int* ReadFile(const char* filename, int& n) {
	ifstream fin(filename);
	if (!fin) {
		cout << endl << "Cannot open file " << filename << " !!!" << endl;
		return NULL;
	}
	/*cout << endl << "Opened successfully file " << filename << " !!!" << endl;*/
	fin >> n; //Get the size of array
	int* a = new int[n];
	for (int i = 0; i < n; i++) {
		fin >> a[i];
	}
	fin.close();
	return a;
}

//Write file for command 2 and 5
void WriteFile(const char* filename, int* a, int n) {
	ofstream fout(filename);
	if (!fout) {
		cout << endl << "Cannot open file " << filename << " !!!" << endl;
		return;
	}
	/*cout << endl << "Opened successfully file " << filename << " !!!" << endl;*/
	fout << n << endl;
	for (int i = 0; i < n; ++i) {
		fout << a[i] << " ";
	}
	fout.close();
}
// Parse algorithm and mode
void ParseAlgoAndMode(char* Algo, const char* Mode, int*& a, int n, unsigned long long& count_compare, duration<double, milli>& timeCount) {
	if (_strcmpi(Algo, "selection-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compSelectionSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeSelectionSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "insertion-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compInsertionSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeInsertionSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "bubble-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compBubbleSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeBubbleSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "shaker-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compShakerSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeShakerSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "shell-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compShellSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeShellSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "heap-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compHeapSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeHeapSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "merge-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mod
			compMergeSort(a, 0, n - 1, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeMergeSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "quick-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compQuickSort(a, 0, n - 1, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeQuickSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "counting-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compCountingSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeCountingSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "radix-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compRadixSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeRadixSort(a, n, count_compare, timeCount);
		}
		return;
	}
	else if (_strcmpi(Algo, "flash-sort") == 0) {
		if (_strcmpi(Mode, "-comp") == 0) {
			//Run algorithm with comparison mode
			compFlashSort(a, n, count_compare);
		}
		else {
			//Run algorithm with running time mode
			timeFlashSort(a, n, count_compare, timeCount);
		}
	}
}

void BothMode(char* Algo, int*& a, int n, unsigned long long& count_compare, duration<double, milli>& timeCount) {
	if (_strcmpi(Algo, "selection-sort") == 0) {
		//Run algorithm with both mode
		timeSelectionSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "insertion-sort") == 0) {
		//Run algorithm with both mode
		timeInsertionSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "bubble-sort") == 0) {
		//Run algorithm with both mode
		timeBubbleSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "shaker-sort") == 0) {
		//Run algorithm with both mode
		timeShakerSort(a, n, count_compare, timeCount);;
		return;
	}
	else if (_strcmpi(Algo, "shell-sort") == 0) {
		//Run algorithm with both mode
		timeShellSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "heap-sort") == 0) {
		//Run algorithm with both mode
		timeHeapSort(a, n, count_compare, timeCount);;
		return;
	}
	else if (_strcmpi(Algo, "merge-sort") == 0) {
		//Run algorithm with both mode
		timeMergeSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "quick-sort") == 0) {
		//Run algorithm with both mode
		timeQuickSort(a, n, count_compare, timeCount);;
		return;
	}
	else if (_strcmpi(Algo, "counting-sort") == 0) {
		//Run algorithm with both mode
		timeCountingSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "radix-sort") == 0) {
		//Run algorithm with both mode
		timeRadixSort(a, n, count_compare, timeCount);
		return;
	}
	else if (_strcmpi(Algo, "flash-sort") == 0) {
		//Run algorithm with both mode
		timeFlashSort(a, n, count_compare, timeCount);
		return;
	}
}

void DisplayOfCommandLine3(unsigned long long comparison[], duration<double, milli> time0, duration<double, milli> time1, duration<double, milli> time2, duration<double, milli> time3, char* mode) {
	//Random data
	cout << endl << "Input order: Randomize";
	cout << endl << "------------------------";
	cout << endl << "Running time (if required): ";
	if (_strcmpi(mode, "-time") == 0 || _strcmpi(mode, "-both") == 0)
		cout << time0.count() << " ms";
	cout << endl << "Comparisons (if required): ";
	if (_strcmpi(mode, "-comp") == 0 || _strcmpi(mode, "-both") == 0)
		cout << comparison[0];
	cout << endl;
	//Nearly sorted data
	cout << endl << "Input order: Nearyly Sorted";
	cout << endl << "------------------------";
	cout << endl << "Running time (if required): ";
	if (_strcmpi(mode, "-time") == 0 || _strcmpi(mode, "-both") == 0)
		cout << time1.count() << " ms";
	cout << endl << "Comparisons (if required): ";
	if (_strcmpi(mode, "-comp") == 0 || _strcmpi(mode, "-both") == 0)
		cout << comparison[1];
	cout << endl;
	//Sorted data
	cout << endl << "Input order: Sorted";
	cout << endl << "------------------------";
	cout << endl << "Running time (if required): ";
	if (_strcmpi(mode, "-time") == 0 || _strcmpi(mode, "-both") == 0)
		cout << time2.count() << " ms";
	cout << endl << "Comparisons (if required): ";
	if (_strcmpi(mode, "-comp") == 0 || _strcmpi(mode, "-both") == 0)
		cout << comparison[2];
	cout << endl;
	//Reversed data
	cout << endl << "Input order: Reversed";
	cout << endl << "------------------------";
	cout << endl << "Running time (if required): ";
	if (_strcmpi(mode, "-time") == 0 || _strcmpi(mode, "-both") == 0)
		cout << time3.count() << " ms";
	cout << endl << "Comparisons (if required): ";
	if (_strcmpi(mode, "-comp") == 0 || _strcmpi(mode, "-both") == 0)
		cout << comparison[3];
}


// Release memomry of array
void FreeData(int*& a) {
	delete[]a;
	a = NULL;
}

void AlgorithmMode(int argc, char* argv[]) {
	int mark = 0; // Mark if 
	if (strstr(argv[3], "txt") != NULL) { //Command line 1: Command line 1: [Execution file] -a [Algorithm] [Given input] [Output parameter(s)]
		int n = 0; // size of array a
		int* a = ReadFile(argv[3], n);
		if (a == NULL) {
			cout << endl << "Failed to generate array!!" << endl;
			return;
		}
		unsigned long long comparision = 0;
		duration<double, milli> time;
		ParseAlgoAndMode(argv[2], argv[4], a, n, comparision, time);
		cout << endl << "Algorithm: " << argv[2];
		cout << endl << "Input file: " << argv[4];

		cout << endl << "------------------------";
		cout << endl << "Running time (if required): ";
		if (_strcmpi(argv[4], "-time") == 0 || _strcmpi(argv[4], "-both") == 0)
			cout << time.count() << " ms";
		cout << endl << "Comparisons (if required): ";
		if (_strcmpi(argv[4], "-comp") == 0 || _strcmpi(argv[4], "-both") == 0)
			cout << comparision;

		WriteFile("output.txt", a, n);
		FreeData(a);
	}
	else  { //Command line 2 or 3
		int n = atoi(argv[3]); // Size of an array
		if (argc == 5) { // Command line 3: [Execution file] -a [Algorithm] [Input size] [Output parameter(s)]  (Run a sorting algorithm on ALL data arrangements of a specified size).
			int* randomArray = new int[n];
			int* sortedArray = new int[n];
			int* reverseArray = new int[n];
			int* nearlySortedArray = new int[n];
			GenerateRandomData(randomArray, n);
			GenerateSortedData(sortedArray, n);
			GenerateReverseData(reverseArray, n);
			GenerateNearlySortedData(nearlySortedArray, n);
			unsigned long long comparison[4] = { 0 };
			duration<double, milli> time0;				// 0: random, 1: nearly sorted
			duration<double, milli> time1;				// 2: Sorted, 3: Reversed
			duration<double, milli> time2;
			duration<double, milli> time3;
			cout << endl << "Algorithm: " << argv[2];
			cout << endl << "Input size: " << n;
			if (_strcmpi(argv[4], "-both") == 0) { //Both options
				BothMode(argv[2], randomArray, n, comparison[0], time0);			// for random data
				BothMode(argv[2], nearlySortedArray, n, comparison[1], time1);	// for nearly sorted data
				BothMode(argv[2], sortedArray, n, comparison[2], time2);			// for sorted data
				BothMode(argv[2], reverseArray, n, comparison[3], time3);		// for reversed data
				DisplayOfCommandLine3(comparison, time0, time1, time2, time3, argv[4]);
			}
			else {
				ParseAlgoAndMode(argv[2], argv[4], randomArray, n, comparison[0], time0);		// for random data;
				ParseAlgoAndMode(argv[2], argv[4], nearlySortedArray, n, comparison[1], time1);  // for nearly sorted data
				ParseAlgoAndMode(argv[2], argv[4], sortedArray, n, comparison[2], time2);		// for sorted data
				ParseAlgoAndMode(argv[2], argv[4], reverseArray, n, comparison[3], time3);		// for reversed data
				DisplayOfCommandLine3(comparison, time0, time1, time2, time3, argv[4]);
				// 0: random, 1: nearly sorted
				// 2: Sorted, 3: Reversed
			}
			WriteFile("input_1.txt", randomArray, n);
			WriteFile("input_2.txt", nearlySortedArray, n);
			WriteFile("input_3.txt", sortedArray, n);
			WriteFile("input_4.txt", reverseArray, n);

			FreeData(randomArray);
			FreeData(nearlySortedArray);
			FreeData(sortedArray);
			FreeData(reverseArray);
		} // end if (argc == 5) { // Command line 3: [Execution file] -a [Algorithm] [Input size] [Output parameter(s)]   (Run a sorting algorithm on ALL data arrangements of a specified size).

		else if (argc == 6) { // Command line 2: [Execution file] -a [Algorithm] [Input size] [Input order] [Output parameter(s)]
			int dataType = DataType(argv[4]);
			if (dataType == 4) {
				cout << endl << "Entered wrong command line!!" << endl;
				return;
			}
			int* a = new int[n];
			unsigned long long comparison = 0;
			duration<double, milli> time;
			GenerateData(a, n, dataType);
			WriteFile("input.txt", a, n);
			ParseAlgoAndMode(argv[2], argv[5], a, n, comparison, time);
			cout << endl << "Algorithm: " << argv[2];
			cout << endl << "Input size: " << n;
			cout << endl << "Input order: " << argv[4];
			cout << endl << "------------------------";
			cout << endl << "Running time (if required): ";
			if (_strcmpi(argv[5], "-time") == 0 || _strcmpi(argv[5], "-both") == 0)
				cout << time.count() << " ms";
			cout << endl << "Comparisons (if required): ";
			if (_strcmpi(argv[5], "-comp") == 0 || _strcmpi(argv[5], "-both") == 0)
				cout << comparison;

			WriteFile("output.txt", a, n);
			FreeData(a);
		}
		else {
			cout << endl << "Entered wrong command line!!" << endl;
		}
	}
}

// Duplicate array
int* DuplicateArray(int* src, int n) {
	int* b = new int[n];
	for (int i = 0; i < n; i++) {
		b[i] = src[i];
	}
	return b;
}
void CompareMode(int argc, char* argv[]) {
	if (argc == 5) { // Command line 4: [Execution file] -c [Algorithm 1] [Algorithm 2] [Given input]
		int n = 0; // size of array a
		int* a = ReadFile(argv[4], n);
		if (a == NULL) {
			cout << endl << "Failed to generate array!!" << endl;
			return;
		}
		int* b = DuplicateArray(a, n);
		unsigned long long comparison[2] = { 0 };
		duration<double, milli> time0;
		duration<double, milli> time1;
		BothMode(argv[2], a, n, comparison[0], time0);
		BothMode(argv[3], b, n, comparison[1], time1);
		cout << endl << "Algorithm: " << argv[2] << " | " << argv[3];
		cout << endl << "Input file: " << endl << argv[4] << endl;
		cout << endl << "Input size: " << n;
		cout << endl << "------------------------";
		cout << endl << "Running time: " << time0.count() << " ms | " << time1.count() << " ms";
		cout << endl << "Comparisons: " << comparison[0] << " | " << comparison[1];

		FreeData(a);
		FreeData(b);
	}
	else if (argc == 6) { // Command 5: [Execution file] -c [Algorithm 1] [Algorithm 2] [Input size] [Input order],  Run two sorting algorithms on the data generated automatically
		int n = atoi(argv[4]);
		int dataType = DataType(argv[5]);
		if (dataType == 4) {
			cout << endl << "Entered wrong command line!!" << endl;
			return;
		}
		int* a = new int[n];
		GenerateData(a, n, dataType);
		WriteFile("input.txt", a, n);
		int* b = DuplicateArray(a, n);
		unsigned long long comparison[2] = { 0 };
		duration<double, milli> time0;
		duration<double, milli> time1;
		BothMode(argv[2], a, n, comparison[0], time0);
		BothMode(argv[3], b, n, comparison[1], time1);
		cout << endl << "Algorithm: " << argv[2] << " | " << argv[3];
		cout << endl << "Input size: " << n;
		cout << endl << "Input order: " << argv[5];
		cout << endl << "------------------------";
		cout << endl << "Running time: " << time0.count() << " ms | " << time1.count() << " ms";
		cout << endl << "Comparisons: " << comparison[0] << " | " << comparison[1];
		FreeData(a);

	}
	else {
		cout << endl << "Entered wrong command line!!" << endl;
	}
}
void GenerateTestArray(int**& a, int n) {
	a = new int* [11];
	for (int i = 0; i < 11; i++)
		a[i] = new int[n];
	GenerateRandomData(a[0], n);
	//GenerateReverseData(a[0], n);
	/*GenerateNearlySortedData(a[0], n);*/
	for (int i = 1; i < 11; i++) {
		a[i] = DuplicateArray(a[0], n);
	}
}
void DeleteTestArray(int**& a, int n) {
	for (int i = 0; i < 11; i++) {
		delete[]a[i];
		a[i] = NULL;
	}
	delete a;
	a = NULL;
}
int main(int argc, char* argv[]) {
	// Command line
	if (_strcmpi(argv[1], "-a") == 0) { //This is Algorithm mode
		if (argc <= 4 || argc >= 7) {
			cout << endl << "You entered the wrong command line!!!" << endl;
			return 1;
		}
		cout << endl << "ALGORITHM MODE";
		AlgorithmMode(argc, argv);

	}
	else if (_strcmpi(argv[1], "-c") == 0) { //This is Comparision mode
		if (argc <= 4 || argc >= 7) {
			cout << endl << "You entered the wrong command line!!!" << endl;
			return 1;
		}
		cout << endl << "COMPARE MODE";
		CompareMode(argc, argv);
	}
	return 1;
	
	//FOR TESTING

	/*int n = 50000;
	int* a = new int[n];
	GenerateRandomData(a, n);
	WriteFile("input.txt", a, n);*/

	/*vector<int> dataSize(6);
	dataSize[0] = 500000;
	vector <unsigned long long> comparison(11, 0);
	vector <duration <double, milli>> time(11);
	cout << endl << "=====Data size: " << dataSize[0] << " =====" << endl;
	int** a; GenerateTestArray(a, dataSize[0]);
	int i = 0;
	int j = 0;
	cout << endl << "Generate random successfully !!!" << endl;
	timeSelectionSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "SelSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;


	timeInsertionSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "InsSort: " << comparison[j] << " | " << time[j].count() << endl;;
	j++;


	timeBubbleSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "BubbleSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;


	timeShakerSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "ShakerSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeShellSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "ShellSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeHeapSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "HeapSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeMergeSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "MergeSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeQuickSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "QSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeCountingSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "CountSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeRadixSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "RadixSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;

	timeFlashSort(a[j], dataSize[i], comparison[j], time[j]);
	cout << "FlashSort: " << comparison[j] << " | " << time[j].count() << endl;
	j++;
	cout << "Run " << j - 1 << " successfully!!" << endl;


	DeleteTestArray(a, dataSize[i]);
	comparison.clear();
	time.clear();
	comparison.~vector();
	time.~vector();

	cout << endl << endl;*/
}
