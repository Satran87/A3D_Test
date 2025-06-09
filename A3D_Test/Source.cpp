//https://www.geeksforgeeks.org/external-sorting/
//https://en.wikipedia.org/wiki/External_sorting
//http://qaru.site/questions/13923/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
//https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
// C++ program to implement external sorting using 
// merge sort 
#define _CRT_SECURE_NO_WARNINGS

//#define  NO_ASYNC

const unsigned _110Mb = 115343360;

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <ctime>

#ifndef NO_ASYNC
#include <future>
#endif
using namespace std;
unsigned RealMaxLineSize = 11;//+1 íà /r
class my_array
{
public:
	my_array():currentIndex(0)
	{
		const auto baseSize = 2500000;
		arr = new unsigned[baseSize];
		size = baseSize;
	}

	//Óâåëè÷èòü åñëè íàäî
	void check()
	{
		if (size-10< currentIndex)
		{
			const auto arrTmp = new unsigned[size+100000]{};
			std::copy_n(arr, currentIndex+1, arrTmp);
			delete[]arr;
			arr = arrTmp;
			size += 100000;
		}
	}
	void add(unsigned value,long index)
	{
		arr[index] = value;
		currentIndex = index;
	}
	unsigned * get_raw() const
	{
		return arr;
	}
	~my_array()
	{
		delete[]arr;
	}
private:
	long currentIndex;
	long size;
	unsigned *arr;
};

bool fileExist(const char* filename)
{
		if (FILE *file = fopen(filename, "r")) 
		{
			fclose(file);
			return true;
		}
			return false;
}

struct MinHeapNode
{
	// The element to be stored 
	unsigned element;

	// index of the array from which the element is taken 
	long long i;
};

// Prototype of a utility function to swap two min heap nodes 
void swap(MinHeapNode* x, MinHeapNode* y);

// A class for Min Heap 
class MinHeap
{
	MinHeapNode* harr; // pointer to array of elements in heap 
	long long heap_size;	 // size of min heap 

public:
	// Constructor: creates a min heap of given size 
	MinHeap(MinHeapNode a[], long long size);

	// to heapify a subtree with root at given index 
	void MinHeapify(long long);

	// to get index of left child of node at index i 
	long long left(long long i) { return (2 * i + 1); }

	// to get index of right child of node at index i 
	long long right(long long i) { return (2 * i + 2); }

	// to get the root 
	MinHeapNode getMin() { return harr[0]; }

	// to replace root with new node x and heapify() 
	// new root 
	void replaceMin(MinHeapNode x)
	{
		harr[0] = x;
		MinHeapify(0);
	}
};

// Constructor: Builds a heap from a given array a[] 
// of given size 
MinHeap::MinHeap(MinHeapNode a[], long long size)
{
	heap_size = size;
	harr = a; // store address of array 
	long long i = (heap_size - 1) / 2;
	while (i >= 0)
	{
		MinHeapify(i);
		i--;
	}
}

// A recursive method to heapify a subtree with root 
// at given index. This method assumes that the 
// subtrees are already heapified 
void MinHeap::MinHeapify(long long i)
{
	long long l = left(i);
	long long r = right(i);
	long long smallest = i;
	if (l < heap_size && harr[l].element < harr[i].element)
		smallest = l;
	if (r < heap_size && harr[r].element < harr[smallest].element)
		smallest = r;
	if (smallest != i)
	{
		swap(&harr[i], &harr[smallest]);
		MinHeapify(smallest);
	}
}

// A utility function to swap two elements 
void swap(MinHeapNode* x, MinHeapNode* y)
{
	MinHeapNode temp = *x;
	*x = *y;
	*y = temp;
}

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(unsigned arr[], long l, long m, long r)
{
	long i, j, k;
	long n1 = m - l + 1;
	long n2 = r - m;

	/* create temp arrays */
	auto L = new unsigned[n1];
	auto R = new unsigned[n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++)
		L[i] = arr[l + i];
	for (j = 0; j < n2; j++)
		R[j] = arr[m + 1 + j];

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray 
	j = 0; // Initial index of second subarray 
	k = l; // Initial index of merged subarray 
	while (i < n1 && j < n2)
	{
		if (L[i] <= R[j])
			arr[k++] = L[i++];
		else
			arr[k++] = R[j++];
	}

	/* Copy the remaining elements of L[], if there
	are any */
	while (i < n1)
		arr[k++] = L[i++];

	/* Copy the remaining elements of R[], if there
	are any */
	while (j < n2)
		arr[k++] = R[j++];
	delete[]L;
	delete[]R;
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(unsigned arr[], long l, long r)
{
	if (l < r)
	{
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		long m = l + (r - l) / 2;

		// Sort first and second halves 
		mergeSort(arr, l, m);
		mergeSort(arr, m + 1, r);

		merge(arr, l, m, r);
	}
}


// Merges k sorted files. Names of files are assumed 
// to be 1, 2, 3, ... k 
void mergeFiles(char *output_file, const long long k)
{
	auto in=new std::fstream[k];
	for (long long i = 0; i < k; ++i)
		// Open output files in read mode. 
		in[i] = std::fstream(to_string(i), std::ios::binary | istream::in);

	// FINAL OUTPUT FILE 
	
	std::fstream out = std::fstream(output_file, std::ios::binary | istream::out);

	// Create a min heap with k heap nodes. Every heap node 
	// has first element of scratch output file 
	const auto harr = new MinHeapNode[k];
	long long i;
	std::string tmpStr;
	tmpStr.reserve(RealMaxLineSize);
	for (i = 0; i < k; ++i)
	{
		// break if no output file is empty and 
		// index i will be no. of input files 
		if (std::getline(in[i], tmpStr))
			harr[i].element = std::stoul(tmpStr, nullptr, 10);
		else
			break;

		harr[i].i = i; // Index of scratch output file 
	}
	MinHeap hp(harr, i); // Create the heap 

	long long count = 0;

	// Now one by one get the minimum element from min 
	// heap and replace it with next element. 
	// run till all filled input files reach EOF 
	while (count != i)
	{
		// Get the minimum element and store it in output file 
		MinHeapNode root = hp.getMin();
		out << root.element<< "\n\r";

		// Find the next element that will replace current 
		// root of heap. The next element belongs to same 
		// input file as the current min element. 
		if (std::getline(in[root.i], tmpStr))
			root.element = std::stoul(tmpStr, nullptr, 10);
		else
		{
			root.element = 0xffffffff;//UINT_MAX
			++count;
		}

		// Replace root with next element of input file 
		hp.replaceMin(root);
	}

	// close input and output files 
	for (long long i1 = 0; i1 < k; ++i1)
	{
		in[i1].close();
	}
	delete[]in;
	delete[]harr;
	out.close();
}
void CreateSortChunkFile( unsigned* arr, long long next_output_file, long i)
{
	auto out =std::fstream(std::to_string(next_output_file), std::ios::binary | istream::out);
	mergeSort(arr, 0, i - 1);
	for (int j = 0; j < i; ++j)
		out << arr[j] << "\n\r";
#ifndef NO_ASYNC
	delete[] arr;
#endif
	out.close();
}

// Using a merge-sort algorithm, create the initial runs 
// and divide them evenly among the output files 
long long createInitialRuns(char *input_file, int run_size)
{
	// For big input file 
	auto in = std::fstream(input_file, std::ios::binary | istream::in);


	// allocate a dynamic array large enough 
	// to accommodate runs of size run_size 
	auto run_sizeTmp = run_size;
	my_array arrMy;

	bool more_input = true;
	long long next_output_file = 0;

	long i;
	std::string line;
#ifndef NO_ASYNC
	unsigned MaxThread = 4;// std::thread::hardware_concurrency();
	if (MaxThread == 0)
		MaxThread = 1;
	auto myThreads = new future<void>[MaxThread];
#endif
	std::string tmpStr;
	tmpStr.reserve(RealMaxLineSize);
	unsigned short ttt = 0;
	while (more_input)
	{
		i = 0;
		// write run_size elements into arr from input file 
		while (0 < run_sizeTmp)
		{
			if(std::getline(in, tmpStr))
			{
				run_sizeTmp = run_sizeTmp - tmpStr.size();
				arrMy.add(std::stoul(tmpStr, nullptr, 10),i);
				arrMy.check();
			}
			else 
			{
				more_input = false;
				break;
			}
			++i;
		}
		run_sizeTmp = run_size;
		// sort array using merge sort 
#ifdef NO_ASYNC
		CreateSortChunkFile(arrMy.get_raw(), next_output_file, i);
#else
		auto arrTmp = new unsigned[i];
		std::copy_n(arrMy.get_raw(), i, arrTmp);
		unsigned myCounter = 0;
		while (true)
		{
			if(myThreads[myCounter].valid())
			{
				auto status=myThreads[myCounter].wait_for(0ms);
				if(status!=future_status::ready)//Åùå ðàáîòàåò
					++myCounter;
				else
					break;
			}
			else
				break;
			if(myCounter>MaxThread-1)
			{
				std::cout << "here";
				myCounter = 0;
			}
				
		}
		if (ttt < myCounter)
			ttt = myCounter;
		myThreads[myCounter] = std::async(std::launch::async, CreateSortChunkFile, arrTmp, next_output_file, i);
#endif
		++next_output_file;
	}
	std::cout << "MaxThread used " << ttt+1<<std::endl;
#ifndef NO_ASYNC
	for(unsigned myCounter1=0;myCounter1<MaxThread;++myCounter1)
	{
		if (myThreads[myCounter1].valid())
			myThreads[myCounter1].wait();
	}
#endif
	in.close();
	delete[] myThreads;
	return next_output_file;
}

// For sorting data stored on disk 
void externalSort(char* input_file, char *output_file, int run_size)
{
	// read the input file, create the initial runs, 
	// and assign the runs to the scratch output files 
	const auto totalFiles=createInitialRuns(input_file, run_size);

	// Merge the runs using the K-way merging 
	mergeFiles(output_file, totalFiles);
}

// Driver program to test above 
int main(int argc, char* argv[])
{
	clock_t t1 = clock();
	constexpr int run_size = 30*1024*1024;
	/*char input_file[] = "D:\\Sample\\input";
	char output_file[] = "D:\\Sample\\output";*/
	char input_file[] = "F:\\Sample\\input";
	char output_file[] = "F:\\Sample\\output";
	if (!fileExist(input_file))
	{
		auto tmp=std::fstream(output_file, std::ios::binary | istream::out);
		tmp.close();
		return EXIT_SUCCESS;
	}
	externalSort(input_file, output_file,
		run_size);

	clock_t t2 = clock();
	cout << (t2 - t1 + .0) / CLOCKS_PER_SEC << endl;
	system("pause");
	return EXIT_SUCCESS;
}
