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
#include <time.h>
using namespace std;

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
    cout << "==============INPUT ARRAY==============";
    for (int i = 0; i < n; i++) {
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


////Sorting functions
//Selection sort
//Find the smallest ele from unsorted array and swap with i
void selectionSortFunction(int *a, int n, unsigned int &comparisonNum) {
    for (int i = 0; i < n; i++) {
        //Loop thru unsorted array
        int minIndex = i; //Set min index
        for (int j = i + 1; j < n; j++) {
            if (a[minIndex] > a[j]) {
                minIndex = j;
            } 
        }
        //Only swap if min is not at i
        if (minIndex != i)
            swap(a[minIndex], a[i]);
    }
    
}
//Time complex O(n^2)
void selectionSort(int* a, int n, unsigned int &comparisonNum, clock_t &timeCount) {
    int* temp = createTemp(a, n);
    clock_t start = clock();
    //Run sort
    selectionSortFunction(temp, n, comparisonNum);
    
    clock_t end = clock();
    //translate to ms
    timeCount = (float)(end - start)/CLOCKS_PER_SEC;
}

//Bubble sort
//Compare 2 ele at the time and swap them -> increment to next eles
//Highest ele will go to the right
void bubbleSortFunction(int* a, int n, unsigned int &comparisonNum) {
    for (int i = 0; i < n - 1; i++) {
        //Go thru 
        for (int j = 0; j < n - i - 1; j++) {

            if (a[j] > a[j + 1]) {
                swap(a[j], a[j + 1]);
            }   
        }   
    }
}
//Time complex O(n^2)
void bubbleSort(int* a, int n, unsigned int &comparisonNum, clock_t &timeCount) {
    int* temp = createTemp(a, n);
    clock_t start = clock();
    //Run sort
    bubbleSortFunction(temp, n, comparisonNum);
    clock_t end = clock();

    timeCount = (float)(end - start)/CLOCKS_PER_SEC;
}

void mergeArrays(int* &a, int l, int m, int r, unsigned int &comparisonNum) {
    int nl = 0;
    int nr = 0;
    int na = 0;
    //First half
    int x = m - l + 1;
    //Second half
    int y = r - m;
    //create 2 halves array
    int* lArray = new int[x];
    int* rArray = new int[y];


    for (int i = 0; i < x; i++)
        lArray[i] = a[l + i];
    for (int i = 0; i < y; i++)
        rArray[i] = a[m + l + i];
    //Merge 2 arrays
    while (nl < x && nr < y) {
        if (lArray[nr] <= rArray[nl]) {
            a[na++] = lArray[nl++];
        }

        else if(lArray[nr] > rArray[nl]) {
            a[na++] = rArray[nr++];
        }
    }

    //Put remaining in a
    while(nl < x) {
        a[na++] = lArray[nl++];
    }
    while(nr < x) {
        a[na++] = rArray[nr++];
    }

}

void mergeSortFunction(int* a, int l, int r, unsigned int &comparisonNum) {
    if (l < r) {
        int m = (l + r)/2;

        mergeSortFunction(a, l, m, comparisonNum);
        mergeSortFunction(a, m + 1, r, comparisonNum);
        mergeArrays(a, l, m, r, comparisonNum);
    }
    
}

void mergeSort(int* a, int n, unsigned int &comparisonNum, clock_t &timeCount) {
    int* temp = createTemp(a, n);
    clock_t start = clock();
    mergeSortFunction(temp, 0, n - 1, comparisonNum);
    clock_t end = clock();

    timeCount = (float)(end - start)/CLOCKS_PER_SEC;
}



