/*Requirements
    2 Mode comparison (so sánh 2 algo) và algorithm (chạy 1 algo)
    1. Đọc a lấy dữ liệu -> sort -> chọn option (time hay comparison or both)
    2. Chạy algo với size, order mình nhập
    3. Chạy toàn bộ algo với size mình nhập
    4. Compare 2 algo với file a
    5. Compare 2 algo với size, order mình nhập
*/
/*Convention
    Initialize functions
    kiểu dữ liệu truyền vào các algo là int* &a,
        unsigned int &comparisonNum, với biến thời gian &timeCount
    Gom phần helper của sort algo trong dấu //###############
    Comment code
*/
/*Structure code
    Ham main
        Truyền tham chiếu cho biến comparisonNum, timeCount
        Tao vector/mang chứa số comparison và time chạy cho từng thuật toán
        Write file từng dòng là từng tên của algo và số comp và time
*/

/*Tri's Part
Insertion
Heap
Quick
Shaker
Counting
*/

#include <iostream>
#include <time.h>
#include <chrono>
using namespace std;

//Create temp
int* createTemp(int* a, int n) {
    int* temp = new int[n];
    for (int i = 0; i < n; i++)
    {
        temp[i] = a[i];
    }
    return temp;
}

// Hàm phát sinh mảng dữ liệu ngẫu nhiên
void GenerateRandomData(int*& a, int n)
{
    srand((unsigned int)time(NULL));

    for (int i = 0; i < n; i++)
    {
        a[i] = rand() % n;
    }
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


//##################################################
//Insertion sort
void insertionSortFunction(int* a, int n, unsigned long long& comparisonNum) {
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


void insertionSort(int* a, int n, unsigned long long& comparisonNum, std::chrono::duration<double>& timeCount) {
    int* temp = new int[n];
    GenerateRandomData(temp, n);
    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    // Run sort
    insertionSortFunction(temp, n, comparisonNum);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeCount = end - start;
    printArray(temp, n);
    cout << "Completed at "  << fixed << timeCount.count() << "ms.";
    delete[]temp;
}
//End
//##################################################

//##################################################
//Heap sort
// To heapify a subtree rooted with node i which is 
// an index in a[]. n is size of heap 
void heapify(int* a, int n, int i, unsigned long long &comparisonNum)
{
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
void heapSortFunction(int* a, int n, unsigned long long& comparisonNum)
{
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
void heapSort(int* a, int n, unsigned long long& comparisonNum, std::chrono::duration<double>& timeCount) {
    int* temp = new int[n];
    GenerateRandomData(temp, n);
    printArray(temp, n);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    // Run sort
    heapSortFunction(temp, n, comparisonNum);
    auto end = std::chrono::steady_clock::now();
    timeCount = end - start;
    printArray(temp, n);
    cout << "Completed at " << fixed << timeCount.count() << "ms.";
    delete[]temp;
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
int partition(int* a, int low, int high, unsigned long long &comparisonNum)
{
    int pivot = a[high]; // pivot  
    int i = (low - 1); // Index of smaller element  

    for (int j = low; ++comparisonNum && j <= high - 1; j++)
    {
        // If current element is smaller than the pivot  
        if (++comparisonNum && a[j] < pivot)
        {
            i++; // increment index of smaller element  
            swap(a[i], a[j]);
        }
    }
    swap(a[i + 1], a[high]);
    return (i + 1);
}

/* The main function that implements QuickSort
a[] --> aay to be sorted,
low --> Starting index,
high --> Ending index */
void quickSortFunction(int *a, int low, int high, unsigned long long& comparisonNum)
{
    if (low < high)
    {
        /* pi is partitioning index, a[p] is now
        at right place */
        int pi = partition(a, low, high, comparisonNum);

        // Separately sort elements before  
        // partition and after partition  
        quickSortFunction(a, low, pi - 1, comparisonNum);
        quickSortFunction(a, pi + 1, high, comparisonNum);
    }
}

void quickSort(int* a, int n, unsigned long long& comparisonNum, std::chrono::duration<double>& timeCount) {
    int* temp = new int[n];
    GenerateRandomData(temp, n);
    chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    // Run sort
    quickSortFunction(temp, 0, n-1, comparisonNum);
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    timeCount = end - start;
    printArray(temp, n);
    cout << "Completed at " << fixed << timeCount.count() << "ms.";
    delete[]temp;
}
//End
//##################################################

//##################################################
//Shaker sort
void shakerSortFunction(int *a, int n, unsigned long long &comparisonNum)
{
    int left = 0;
    int right = n - 1;
    int k = 0;
    while (++comparisonNum && left < right)
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


void shakerSort(int* a, int n, unsigned long long& comparisonNum, chrono::duration<double>& timeCount) {
    int* temp = new int[n];
    GenerateRandomData(temp, n);
    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    // Run sort
    quickSortFunction(temp, 0, n - 1, comparisonNum);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeCount = end - start;
    printArray(temp, n);
    cout << "Completed at " << fixed << timeCount.count() << "ms.";
    delete[]temp;
}
// End
//##################################################

//##################################################
// Counting sort
// Function that sort the given a
void countingSortFunction(int* a, int n, unsigned long long& comparisonNum)
{
    int *output = new int[n]; // The output will have sorted an array
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

    int *count_array = new int [k]; // Create a count_array to store count of each individual a value
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

void countingSort(int* a, int n, unsigned long long& comparisonNum, chrono::duration<double>& timeCount) {
    int* temp = new int[n];
    GenerateRandomData(temp, n);
    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    // Run sort
    countingSortFunction(temp, n, comparisonNum);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeCount = end - start;
    printArray(temp, n);
    cout << "Completed at " << fixed << timeCount.count() << "ms.";
    delete[]temp;
}
// End
//##################################################







int main() {
    int n = 0;
    cout << "Input: ";
    cin >> n;

    unsigned long long count = 0;
    std::chrono::duration<double> time;
    int* a = new int[n];
    //inputArray(a, n);
    countingSort(a, n, count, time);
    cout << "comparison Num" << count;
    delete[]a;

    return 0;
}
