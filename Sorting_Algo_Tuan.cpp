/*Requirements
    2 Mode comparison (so sánh 2 algo) và algorithm (chạy 1 algo)
    1. Đọc input lấy dữ liệu -> sort -> chọn option (time hay comparison or both)
    2. Chạy algo với size, order mình nhập
    3. Chạy toàn bộ algo với size mình nhập
    4. Compare 2 algo với file input
    5. Compare 2 algo với size, order mình nhập
*/
/*Convention
    Initialize functions
    Viết hàm <NAME>SortFunction để sort và đếm comparison
    Viết hàm <NAME>Sort để chạy cái hàm trên và đếm thời gian
    kiểu dữ liệu truyền vào các algo là int* &a, 
        unsigned int &comparisonNum, với biến thời gian &timeCount
    Truyền mảng a vào thì phải truyền mảng temp để sort (không thay đổi mảng a)
    Gom phần helper của sort algo trong dấu //###############
    Comment code
*/
/*Structure code
    Ham main
        Truyền tham chiếu cho biến comparisonNum, timeCount
        Tao vector/mang chứa số comparison và time chạy cho từng thuật toán
        Write file từng dòng là từng tên của algo và số comp và time
*/

/*Tuan's Part
Selection
Bubble
Merge
Radix
Shell
Flash
*/

#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
using namespace std;
using namespace std::chrono;

////Initialize functions

////Helper functions
void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
}
//Create temp
int* createTemp(int *a, int n) {
    int* temp = new int[n];
    for (int i = 0; i < n; i++)
    {
        temp[i] = a[i];
    }
    return temp;
}
//For Testing
//Input array
void inputArray(int* &a, int n) {
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
////Sorting functions
//Selection sort
//Find the smallest ele from unsorted array and swap with i
void comp_SelectionSort(int *a, int n, unsigned long long &comparisonNum) {
    comparisonNum = 0;
    for (int i = 0; ++comparisonNum && i < n; i++) {
        //Loop thru unsorted array
        int minIndex = i; //Set min index
        for (int j = i + 1;++comparisonNum && j < n; j++) {
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
void time_SelectionSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    comp_SelectionSort(temp, n, comparisonNum);
    
    auto end = high_resolution_clock::now();
    //translate to ms
    timeCount = duration_cast<microseconds>(end - start);
    printArray(temp, n);

}
//End
//##################################################


//##################################################
//Bubble sort
//Compare 2 ele at the time and swap them -> increment to next eles
//Highest ele will go to the right
void comp_BubbleSort(int* a, int n, unsigned long long &comparisonNum) {
    comparisonNum = 0;
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
void time_BubbleSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    //Run sort
    comp_BubbleSort(temp, n, comparisonNum);
    auto end = high_resolution_clock::now();  

    timeCount = duration_cast<microseconds>(end - start);
    printArray(temp, n);

}
//End
//##################################################


//##################################################
//Merge sort
void mergeArrays(int* a, int l, int m, int r, unsigned long long &comparisonNum) {
    //First half
    int x = m - l + 1;
    //Second half
    int y =  r - (m + 1) + 1;
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
            comparisonNum++;
            a[ia] = rArray[ir];
            ia++;
            ir++;
        }
    }

    //Put remaining in a
    while(++comparisonNum && il < x) {
        a[ia] = lArray[il];
        ia++;
        il++;
    }
    while(++comparisonNum && ir < y) {
        a[ia] = rArray[ir];
        ia++;
        ir++;
    }

    delete[] lArray;
    delete[] rArray;
}

void comp_MergeSort(int* a, int l, int r, unsigned long long &comparisonNum) {
    comparisonNum = 0;
    if (++comparisonNum && l < r) {
        //Find middle point
        int m = (l + r)/2;
        //Call recursive
        comp_MergeSort(a, l, m, comparisonNum);
        comp_MergeSort(a, m + 1, r, comparisonNum);
        mergeArrays(a, l, m, r, comparisonNum);
    }
    
}

void time_MergeSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    comp_MergeSort(temp, 0, n - 1, comparisonNum);
    auto end = high_resolution_clock::now();  

    timeCount = duration_cast<microseconds>(end - start);
        
    printArray(temp, n);

}
//End
//##################################################


//##################################################
//Radix sort
int maxDigits(int* a, int n, unsigned long long &comparisonNum) {
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
void comp_RadixSort(int* a, int n, unsigned long long &comparisonNum) {
    comparisonNum = 0;
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

		for (int k = 0;++comparisonNum && k < 10; k++) {
			for (int l = 0;l < hash[k].size();l++){
				a[count++] = hash[k][l];
			}
			hash[k].clear();
		}
	}
}

void time_RadixSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    comp_RadixSort(temp, n ,comparisonNum);
    auto end = high_resolution_clock::now();  

    timeCount = duration_cast<microseconds>(end - start);
        
    printArray(temp, n);
}
//End
//##################################################


//##################################################
//Shell sort
void comp_ShellSort(int* a, int n, unsigned long long &comparisonNum) {
    comparisonNum = 0;
    for (int gap = n / 2; ++comparisonNum && gap > 0; gap = gap / 2) {
        for (int i = gap; ++comparisonNum && i < n; i++) {
            // add a[i] to the elements that have been gap sorted
            // save a[i] in temp and make a hole at position i
            int temp = a[i];
 
            // shift earlier gap-sorted elements up until the correct
            // location for a[i] is found           
            int count = 0;
            for (int j = i; (++comparisonNum && j >= gap) && (++comparisonNum && a[j - gap] > temp); j = j - gap) {
                a[j] = a[j - gap];
                count = j;
            }
                
            //  put temp (the original a[i]) in its correct location
            a[count] = temp;
        }
    }
}

void time_ShellSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    comp_ShellSort(temp, n ,comparisonNum);
    auto end = high_resolution_clock::now();  

    timeCount = duration_cast<microseconds>(end - start);
        
    printArray(temp, n);
}
//End
//##################################################


//##################################################
//Flash sort
//https://steemit.com/popularscience/@krishtopa/the-fastest-sort-method-flashsort
//https://www.youtube.com/watch?v=CAaDJJUszvE 
void comp_FlashSort(int* a, int n, unsigned long long &comparisonNum) {
    comparisonNum = 0;
    int max = 0, min = a[0];
    //Estimate num of buckets
    int m = 0.45 * n;
    int* l = new int(m);
    //find Min max
    for (int i = 1; ++comparisonNum && i < n; i++) {
        if (++comparisonNum && a[i] < min) {
            min = a[i];
        }
        if (++comparisonNum && a[i] > a[max]) {
            max = i;
        }
    }
    if (++comparisonNum && max == min)
        return;
    //Count num of ele in buckets
    for (int i = 0; ++comparisonNum && i < n; i++) {
        int k = floor((m - 1) * (a[i] - min) / (max - min));
        l[k]++;
    }
    //Find last index of each buckets
    for (int i = 1; ++comparisonNum && i < m; i++) {
        l[i] = l[i] + l[i - 1];
    }
    
    //Permutation
    int holder = a[0];
    int move = 0;
    int flash = 0;

    int k = 0, t = 0, j = 0;

    while (++comparisonNum && move < n - 1) {
        flash = a[j];
        while (++comparisonNum && j != l[k]) {
            k = floor((m - 1) * (a[j] - min) / (max - min));
            l[k]--;
            holder = a[t - l[k]];
            a[t] = holder;
            flash = holder;
            move++;
        }
        
    }
    //insertion sort
    for (int j = 1; ++comparisonNum && j < n; j++) {
        holder = a[j];
        int i = j - 1;
        while(++comparisonNum && i >= 0 && a[i] > holder)
            a[i + 1] = a[i--];
        
        a[i + 1] = holder;
    }

}

void time_FlashSort(int* a, int n, unsigned long long &comparisonNum, duration<double>& timeCount) {
    int* temp = createTemp(a, n);
    auto start = high_resolution_clock::now();  
    comp_FlashSort(temp, n ,comparisonNum);
    auto end = high_resolution_clock::now();  

    timeCount = duration_cast<microseconds>(end - start);
        
    printArray(temp, n);
}
//End
//##################################################


int main() {
    int n = 0;
    cout << "Input: ";
    cin >> n;
    
    unsigned long long count = 0;
    duration<double> time;
    int* a = new int[n];
    inputArray(a, n);
    time_RadixSort(a, n, count, time);

    return 1;
}




//Tham khao
/*
measure time:
https://www.geeksforgeeks.org/measure-execution-time-function-cpp/ 
https://github.com/leduythuccs/Sorting-Algorithms/blob/master/sorting-methods/flash_sort.h 
*/