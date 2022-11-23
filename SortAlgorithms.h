#pragma once
#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
using namespace std;
using namespace std::chrono;

int* createTemp(int* a, int n) {
	int* temp = new int[n];
	for (int i = 0; i < n; i++)
	{
		temp[i] = a[i];
	}
	return temp;
}


//For Testing
//Input array
void inputArray(int*& a, int n) {
	cout << "==============INPUT ARRAY==============\n";
	for (int i = 0; i < n; i++) {
		cout << "a[" << i << "]: ";
		cin >> a[i];
	}
}
//Print array
void printArray(int* a, int n) {
	cout << "\n==============OUTPUT ARRAY==============\n";
	for (int i = 0; i < n; i++) {
		cout << a[i] << " ";
	}
}

//Selection Sort, Insertion Sort, Bubble Sort, Shaker Sort, Shell Sort, 
//Heap Sort, Merge Sort, Quick Sort, Counting Sort, Radix Sort, and Flash Sort.


//##################################################
////Sorting functions
//Selection sort
//Find the smallest ele from unsorted array and swap with i
void compSelectionSort(int* a, int n, unsigned long long& comparisonNum) {
	for (int i = 0; ++comparisonNum && i < n; i++) {
		//Loop thru unsorted array
		int minIndex = i; //Set min index
		for (int j = i + 1; ++comparisonNum && j < n; j++) {
			if (++comparisonNum && a[minIndex] > a[j]) {
				minIndex = j;
			}
		}
		//Only swap if min is not at i
		if (++comparisonNum && minIndex != i)
			swap(a[minIndex], a[i]);
	}

}
//Time complex O(n^2)
void timeSelectionSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	compSelectionSort(a, n, comparisonNum);

	auto end = high_resolution_clock::now();
	//translate to ms
	timeCount = duration_cast<microseconds>(end - start);
	/* printArray(a, n);*/
}
//End
//##################################################


//##################################################
//Insertion sort
void compInsertionSort(int* a, int n, unsigned long long& comparisonNum) {
	int i, key, j;
	for (i = 1; ++comparisonNum && i < n; i++)

	{
		key = a[i];
		j = i - 1;
		/*Move element which is greater than key to the value key is later a position compared to its original position */
		while (++comparisonNum && j >= 0 && ++comparisonNum && a[j] > key) {
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = key;
	}
}

void timeInsertionSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	// Run sort
	compInsertionSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();
	timeCount = duration_cast<microseconds>(end - start);
	//printArray(a, n);
}
//End
//##################################################


//##################################################
//Bubble sort
//Compare 2 ele at the time and swap them -> increment to next eles
//Highest ele will go to the right
void compBubbleSort(int* a, int n, unsigned long long& comparisonNum) {
	for (int i = 0; ++comparisonNum && i < n - 1; i++) {
		//Go thru 
		for (int j = 0; ++comparisonNum && j < n - i - 1; j++) {

			if (++comparisonNum && a[j] > a[j + 1]) {
				swap(a[j], a[j + 1]);
			}
		}
	}
}
//Time complex O(n^2)
void timeBubbleSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	//Run sort
	compBubbleSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();

	timeCount = duration_cast<microseconds>(end - start);
	/*printArray(a, n);*/

}
//End
//##################################################



//##################################################
//Shaker sort
void compShakerSort(int* a, int n, unsigned long long& comparisonNum)
{
	int left = 0;
	int right = n - 1;
	int k = 0;
	while (++comparisonNum&& left < right)
	{
		for (int i = left; ++comparisonNum && i < right; i++)
		{
			if (++comparisonNum && a[i] > a[i + 1])
			{
				swap(a[i], a[i + 1]);
				k = i;
			}
		}
		right = k;
		for (int i = right; ++comparisonNum && i > left; i--)
		{
			if (++comparisonNum && a[i] < a[i - 1])
			{
				swap(a[i], a[i - 1]);
				k = i;
			}
		}
		left = k;
	}
}

void timeShakerSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	// Run sort
	compShakerSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();
	timeCount = duration_cast<microseconds>(end - start);
	/*printArray(a, n);*/
}
// End
//##################################################



//##################################################
//Shell sort
void compShellSort(int*& a, int n, unsigned long long& comparisonNum) {
	for (int gap = n / 2; ++comparisonNum && gap > 0; gap = gap / 2) {
		for (int i = gap; ++comparisonNum && i < n; i++) {
			// add a[i] to the elements that have been gap sorted
			// save a[i] in temp and make a hole at position i
			int temp = a[i];

			// shift earlier gap-sorted elements up until the correct
			// location for a[i] is found           
			int j;
			for (j = i; (++comparisonNum && j >= gap) && (++comparisonNum && a[j - gap] > temp); j = j - gap) {
				a[j] = a[j - gap];
			}

			//  put temp (the original a[i]) in its correct location
			a[j] = temp;
		}
	}
}

void timeShellSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	compShellSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();

	timeCount = duration_cast<microseconds>(end - start);

	/* printArray(a, n);*/
}
//End
//##################################################




//##################################################
//Heap sort
// To heapify a subtree rooted with node i which is 
// an index in a[]. n is size of heap 
void heapify(int* a, int n, int i, unsigned long long& comparisonNum) {
	int largest = i; // Initialize largest as root 
	int l = 2 * i + 1; // left = 2*i + 1 
	int r = 2 * i + 2; // right = 2*i + 2 

	// If left child is larger than root 

	if ((++comparisonNum && l < n) && (++comparisonNum && a[l] > a[largest]))
		largest = l;

	// If right child is larger than largest so far 
	if ((++comparisonNum && r < n) && (++comparisonNum && a[r] > a[largest]))
		largest = r;

	// If largest is not root 
	if (++comparisonNum && largest != i)
	{
		swap(a[i], a[largest]);

		// Recursively heapify the affected sub-tree 
		heapify(a, n, largest, comparisonNum);
	}
}
// main function to do heap sort 
void compHeapSort(int* a, int n, unsigned long long& comparisonNum) {
	// Build heap (reaange aay) 
	for (int i = n / 2 - 1; ++comparisonNum && i >= 0; i--)
		heapify(a, n, i, comparisonNum);

	// One by one extract an element from heap 
	for (int i = n - 1; ++comparisonNum && i > 0; i--)
	{
		// Move current root to end 
		swap(a[0], a[i]);

		// call max heapify on the reduced heap 
		heapify(a, i, 0, comparisonNum);
	}
}
// Time complex O(nlogn)
void timeHeapSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	// Run sort
	compHeapSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();
	timeCount = duration_cast<microseconds>(end - start);
	/*printArray(a, n);*/
}
//End
//##################################################



//##################################################
//Merge sort
void mergeArrays(int* a, int l, int m, int r, unsigned long long& comparisonNum) {
	//First half
	int x = m - l + 1;
	//Second half
	int y = r - (m + 1) + 1;
	//create 2 halves array
	int* lArray = new int[x];
	int* rArray = new int[y];

	//Put both sides of a into left and right array
	for (int i = 0; ++comparisonNum && i < x; i++)
		lArray[i] = a[l + i];

	for (int j = 0; ++comparisonNum && j < y; j++)
		rArray[j] = a[m + 1 + j];

	int il = 0;
	int ir = 0;
	//Merged array will start at each l
	int ia = l;
	//Merge 2 arrays by adding the lesser element from each arrays
	while ((++comparisonNum && il < x) && (++comparisonNum && ir < y)) {
		if (++comparisonNum && lArray[il] < rArray[ir]) {
			a[ia] = lArray[il];
			ia++;
			il++;
		}

		else {
			a[ia] = rArray[ir];
			ia++;
			ir++;
		}
	}

	//Put remaining in a
	while (++comparisonNum && il < x) {
		a[ia] = lArray[il];
		ia++;
		il++;
	}
	while (++comparisonNum && ir < y) {
		a[ia] = rArray[ir];
		ia++;
		ir++;
	}

	delete[] lArray;
	delete[] rArray;
}

void compMergeSort(int* a, int l, int r, unsigned long long& comparisonNum) {
	if (++comparisonNum && l < r) {
		//Find middle point
		int m = (l + r) / 2;
		//Call recursive
		compMergeSort(a, l, m, comparisonNum);
		compMergeSort(a, m + 1, r, comparisonNum);
		mergeArrays(a, l, m, r, comparisonNum);
	}

}

void timeMergeSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	compMergeSort(a, 0, n - 1, comparisonNum);
	auto end = high_resolution_clock::now();

	timeCount = duration_cast<microseconds>(end - start);

	//printArray(a, n);
}
//End
//##################################################



//##################################################
//Quick sort
/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
aay, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int partition(int* a, int low, int high, int pi, unsigned long long& comparisonNum)
{
	swap(a[pi], a[high]);;
	int i = low;
	int j = high - 1;
	while (++comparisonNum && i <= j) {
		while (a[i] < a[pi])
			i++;
		while (++comparisonNum && a[j] >= a[pi])
			j--;
		if (i < j)
			swap(a[i], a[j]);
	}
	swap(a[high], a[i]);
	return i;
}
/* The main function that implements QuickSort
a[] --> aay to be sorted,
low --> Starting index,
high --> Ending index */
void compQuickSort(int* a, int low, int high, unsigned long long& comparisonNum)
{
	if (++comparisonNum && low < high)
	{
		/* pi is partitioning index, a[p] is now
		at right place */
		int pi = (low + high) / 2;
		int m = partition(a, low, high, pi, comparisonNum);

		// Separately sort elements before  
		// partition and after partition  
		compQuickSort(a, low, m - 1, comparisonNum);
		compQuickSort(a, m + 1, high, comparisonNum);
	}
}

void timeQuickSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	// Run sort
	compQuickSort(a, 0, n - 1, comparisonNum);
	auto end = high_resolution_clock::now();
	timeCount = duration_cast<microseconds>(end - start);
	//printArray(a, n);
}
//End
//##################################################




//##################################################
// Counting sort
// Function that sort the given a
void compCountingSort(int* a, int n, unsigned long long& comparisonNum)
{
	int* output = new int[n]; // The output will have sorted an array
	int max = a[0];
	int min = a[0];

	int i;
	for (i = 1; ++comparisonNum && i < n; i++)
	{
		if (++comparisonNum && a[i] > max)
			max = a[i]; // Maximum value in array
		else if (++comparisonNum && a[i] < min)
			min = a[i]; // Minimum value in array
	}

	int k = max - min + 1; // Size of count array

	int* count_array = new int[k]; // Create a count_array to store count of each individual a value
	for (i = 0; ++comparisonNum && i < k; i++)
		count_array[i] = 0;

	for (i = 0; ++comparisonNum && i < n; i++)
		count_array[a[i] - min]++; // Store count of each individual a value

	/* Change count_array so that count_array now contains actual
	 position of a values in output array */
	for (i = 1; ++comparisonNum && i < k; i++)
		count_array[i] += count_array[i - 1];

	// Populate output array using count_array and a array
	for (i = 0; ++comparisonNum && i < n; i++)
	{
		output[count_array[a[i] - min] - 1] = a[i];
		count_array[a[i] - min]--;
	}

	for (i = 0; ++comparisonNum && i < n; i++)
		a[i] = output[i]; // Copy the output array to a, so that a now contains sorted values

	delete[] output;
	delete[] count_array;
}

void timeCountingSort(int* a, int n, unsigned long long& comparisonNum, chrono::duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	// Run sort
	compCountingSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();
	timeCount = duration_cast<microseconds>(end - start);
	/*printArray(a, n);*/


}
// End
//##################################################



//##################################################
//Radix sort
int maxDigits(int* a, int n, unsigned long long& comparisonNum) {
	int max = a[0];

	for (int i = 1; ++comparisonNum && i < n; i++) {
		if (++comparisonNum && a[i] > max)
			max = a[i];
	}

	int count = 0;

	while (++comparisonNum && max != 0) {
		count++;
		max = max / 10;
	}

	return count;
}

void compRadixSort(int* a, int n, unsigned long long& comparisonNum) {
	// Count max digits
	int maxDig = maxDigits(a, n, comparisonNum);

	vector < vector<int> > hash(10);
	// To mark if the nth position of hash table contains value or not
	for (int i = 1; ++comparisonNum && i <= maxDig; i++) { // i=1, array is sorted by increasing order of units; i=2,.....of dozens,.... 
		for (int j = 0; ++comparisonNum && j < n; j++) {
			int remain = (a[j] / (int)pow(10, i - 1)) % 10; // get the digits of unit, dozen, hundred,...
			hash[remain].push_back(a[j]); // if this array needs sorting by decreasing, let hash[9-main].
		}
		int count = 0; // Assign hash table to array

		for (int k = 0; ++comparisonNum && k < 10; k++) {
			for (int l = 0; ++comparisonNum && l < hash[k].size(); l++) {
				a[count++] = hash[k][l];
			}
			hash[k].clear();
		}
	}
}

void timeRadixSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	compRadixSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();

	timeCount = duration_cast<microseconds>(end - start);

	/*printArray(a, n);*/

}
//End
//##################################################




//##################################################
//Flash sort
void compFlashSort(int* a, int n, unsigned long long& comparisonNum) {
	int minVal = a[0];
	int max = 0;
	//Calculate number of bucket
	int m = int(0.45 * n);

	int* l = new int[m];

	for (int i = 0; ++comparisonNum && i < m; i++)
		l[i] = 0;
	//Find minValue and max location
	for (int i = 1; ++comparisonNum && i < n; i++) {
		if (++comparisonNum && a[i] < minVal)
			minVal = a[i];
		if (++comparisonNum && a[i] > a[max])
			max = i;
	}
	if (++comparisonNum && a[max] == minVal)
		return;
	//
	double d = (double)(m - 1) / (a[max] - minVal);
	for (int i = 0; ++comparisonNum && i < n; i++) {
		int k = int(d * (a[i] - minVal));
		l[k]++;
	}
	for (int i = 1; ++comparisonNum && i < m; i++)
		l[i] += l[i - 1];
	swap(a[max], a[0]);

	int move = 0;
	int j = 0, k = m - 1, t = 0;
	int flash = 0;
	//Put elements in correct buckets
	while (++comparisonNum && move < n - 1) {

		while (++comparisonNum && j > l[k] - 1) {
			j++;
			k = int(d * (a[j] - minVal));
		}
		flash = a[j];

		if (++comparisonNum && k < 0)
			break;

		while (++comparisonNum && j != l[k]) {
			k = int(d * (flash - minVal));
			int holder = a[t = --l[k]];
			a[t] = flash;
			flash = holder;
			move++;
		}
	}

	compInsertionSort(a, n, comparisonNum);
	delete[] l;
}

void timeFlashSort(int* a, int n, unsigned long long& comparisonNum, duration<double, milli>& timeCount) {
	auto start = high_resolution_clock::now();
	compFlashSort(a, n, comparisonNum);
	auto end = high_resolution_clock::now();

	timeCount = duration_cast<microseconds>(end - start);

	/*printArray(a, n);*/

}
//End
//##################################################



//####################################################################################################

//MAIN FUNCTION
//int main() {
//    //CODE HERE
//
//    return 1;