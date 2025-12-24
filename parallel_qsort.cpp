#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

void quicksort_serial(std::vector<int>& arr, int left, int right) {
    if (left >= right) return;
    int pivot = arr[(left + right) / 2];
    int i = left, j = right;
    while (i <= j) {
        while (arr[i] < pivot) ++i;
        while (arr[j] > pivot) --j;
        if (i <= j) std::swap(arr[i++], arr[j--]);
    }
    quicksort_serial(arr, left, j);
    quicksort_serial(arr, i, right);
}

// TODO: Implement quicksort_threaded (spawn threads conditionally)
void quicksort_threaded(std::vector<int>&arr, int left, int right, int min_arr_size, int remaining_threads){
	//Base case, array is either empty (left==right) or already sorted
	if (left >= right) 
		return;
    	
	int size = right - left;
	//if the size of the partitioned array is smaller than a certain size, just sort serially
    	if (size < min_arr_size || remaining_threads <= 1) {
        	quicksort_serial(arr, left, right);
        	return;
    	}
	int pivot = arr[(left + right) / 2];
        int i = left, j = right;

    // Partition
	while (i <= j) {
        while (arr[i] < pivot) ++i;
        while (arr[j] > pivot) --j;
        if (i <= j) std::swap(arr[i++], arr[j--]);
    }

	int next_threads=remaining_threads/2;
	
	std::thread t(
        quicksort_threaded,//passing quicksort function
        std::ref(arr),  //passing array by reference
        i,             
        right,          
        next_threads,   //how many threads to create
       	min_arr_size
    );
    quicksort_threaded(arr, left, j, next_threads, min_arr_size);

    t.join();

}
// TODO: Implement quicksort_openmp (#pragma omp task / taskwait)
void quicksort_omp(std::vector<int>& arr, int left, int right, int min_arr_size)
{
    if (left >= right) return;

    int size = right - left;
    if (size < min_arr_size) {
        quicksort_serial(arr, left, right);
        return;
    }

    int pivot = arr[(left + right) / 2];
    int i = left, j = right;

    // Partition
    while (i <= j) {
        while (arr[i] < pivot) ++i;
        while (arr[j] > pivot) --j;
        if (i <= j) std::swap(arr[i++], arr[j--]);
    }

    // Create tasks for parallel work
    #pragma omp task shared(arr) //sorting left array
    quicksort_omp(arr, left, j, min_arr_size);

    #pragma omp task shared(arr) //sorting right array
    quicksort_omp(arr, i, right, min_arr_size);

    #pragma omp taskwait   //waiting for both subtasks to complete

}   

// g++ -O3 -std=c++20 -pthread parallel_qsort.cpp -o parallel_qsort.x
// g++ -O3 -std=c++20 -fopenmp parallel_qsort.cpp -o parallel_qsort.x


int main() {

    std::ofstream csv("par_qsort_results.csv");
    csv << "num_threads,serial_time,threaded_time, omp_time, threaded_speedup, omp_speedup, threaded_eff, omp_eff\n";
    const int N = 1'000'000;
    std::vector<int> arr(N);
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(0, 1'000'000);
for (int num_threads: {1,2,4,8,16}){
    int max_threads=num_threads;
    omp_set_num_threads(num_threads);
    for (int& x : arr) x = dist(gen);
    //Timing the sorting of the array serially
   auto start = std::chrono::high_resolution_clock::now();
    quicksort_serial(arr, 0, N - 1);
    auto end = std::chrono::high_resolution_clock::now();
    double serial_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Serial quicksort time: " << serial_time << " s\n";
    for (int& x : arr) x = dist(gen); //reshuffle the array
	
	int min_arr_size=20'000;
	
	start=std::chrono::high_resolution_clock::now();
	quicksort_threaded(arr, 0, N - 1, min_arr_size, max_threads);

       
     end= std::chrono::high_resolution_clock::now();
    double threaded_time=std::chrono::duration <double>(end-start).count();
    std::cout<<"Threaded quicksort time with "<<num_threads<<" threads:"<<threaded_time<<"s\n";
	

// reshuffle
for (int& x : arr) x = dist(gen);

start = std::chrono::high_resolution_clock::now();

#pragma omp parallel //Creating team of threads that execute the below code
{
    #pragma omp single
    quicksort_omp(arr, 0, N - 1, min_arr_size); //only one thread calls quicksort_omp, the others wait for tasks
}

end = std::chrono::high_resolution_clock::now();
double omp_time=std::chrono::duration<double>(end - start).count();

std::cout << "OpenMP quicksort time with "<<num_threads<<" threads:"
          << omp_time
          << " s\n";
double threaded_speedup=serial_time/threaded_time;
double omp_speedup=serial_time/omp_time;
double threaded_eff=threaded_speedup/num_threads;
double omp_eff=omp_speedup/num_threads;
csv << num_threads << ","
    << serial_time << ","
    << threaded_time << ","
    << omp_time << ","
    << threaded_speedup << ","
    << omp_speedup << ","
    << threaded_eff << ","
    << omp_eff
    << "\n";

}
}	
