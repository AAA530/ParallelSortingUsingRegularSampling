#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <random>
#include <fstream>
using namespace std;
#pragma warning(disable:4996)


#define MASTER 0


void printArray(int arr[], int start, int end) {
	for (int i = start; i <= end; ++i) {
		std::cout << arr[i] << " ";
	}
	std::cout << std::endl;
}


int partition(int arr[], int start, int end)
{

	int pivot = arr[start];

	int count = 0;
	for (int i = start + 1; i <= end; i++) {
		if (arr[i] <= pivot)
			count++;
	}

	// Giving pivot element its correct position
	int pivotIndex = start + count;
	swap(arr[pivotIndex], arr[start]);

	// Sorting left and right parts of the pivot element
	int i = start, j = end;

	while (i < pivotIndex && j > pivotIndex) {

		while (arr[i] <= pivot) {
			i++;
		}

		while (arr[j] > pivot) {
			j--;
		}

		if (i < pivotIndex && j > pivotIndex) {
			swap(arr[i++], arr[j--]);
		}
	}

	return pivotIndex;
}

void quickSort(int arr[], int start, int end)
{

	// base case
	if (start >= end)
		return;

	// partitioning the array
	int p = partition(arr, start, end);

	// Sorting the left part
	quickSort(arr, start, p - 1);

	// Sorting the right part
	quickSort(arr, p + 1, end);
}


void merge_elements(int arr[], int counts[], int k) {
	int totalelements = 0;

	for (int i = 0; i < k; i++) {
		totalelements = totalelements + counts[i];
	}

	quickSort(arr, 0, totalelements - 1);
}



void multiPivotPartition(int arr[], int n, int pivots[], int numPivots, int** segments) {
	int currentSegmentStart = 0;

	for (int i = 0; i < numPivots; ++i) {
		int pivot = pivots[i];
		int segmentEnd = currentSegmentStart;

		while (segmentEnd < n && arr[segmentEnd] <= pivot) {
			++segmentEnd;
		}

		if (segmentEnd > currentSegmentStart) {
			// Copy the current segment to the result array
			segments[i][0] = currentSegmentStart;
			segments[i][1] = segmentEnd - 1;
		}
		else {
			// Empty segment
			segments[i][0] = segments[i][1] = -1;
		}

		currentSegmentStart = segmentEnd;
	}

	// Handle the last segment
	if (currentSegmentStart < n) {
		segments[numPivots][0] = currentSegmentStart;
		segments[numPivots][1] = n - 1;
	}
	else {
		// Empty segment
		segments[numPivots][0] = segments[numPivots][1] = -1;
	}
}


int read_file_txt(std::ifstream& input, int* array) {
	int i = 0;
	while (input >> array[i]) {
		i++;
	}
	return i;
}

int main(int argc, char* argv[]) {
	if (argc != 3) { std::cout << "Usage: enter input filename to perform sorting.\n"; return 1; }

	std::ifstream input(argv[1]);
	if (!input) { std::cout << "Cannot open input file.\n"; return 1; }        

	int n;
	sscanf(argv[2], "%d", &n);	// Total number of elements

	int numProcess,processId;

	double startTime, elapsedTime;


	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	printf("MPI process %d has started...\n", processId);
	//int temp;

	//int n = 100; 
	int n_per_process = n / numProcess; // Number of elements per process


	/*---------------------------------scatter n/p elements on p processes---------------------------*/
	// scatter n/p elements on p processes
	// Allocate memory for the receive buffer
	int* receive_buffer = (int*)malloc(n_per_process * sizeof(int));

	// Only one process (root) will have the data to scatter
	int* send_buffer = NULL;
	if (processId == MASTER) {
		// Allocate memory and initialize data on the root process


		send_buffer = (int*)malloc(n * sizeof(int));
		int arraySize = read_file_txt(input, send_buffer);
		//for (int i = 0; i < n; i++) {
		//	send_buffer[i] = distribution(gen);
		//	//cout << "Type a number for : " << i<<":"; // Type a number and press enter
		//	//cin >> send_buffer[i]; // Get user input from the keyboard
		//}
	}

	startTime = MPI_Wtime();
	// Scatter data from the root process to all other processes
	MPI_Scatter(send_buffer, n_per_process, MPI_INT, receive_buffer, n_per_process, MPI_INT, 0, MPI_COMM_WORLD);


	// local quick sort 
	quickSort(receive_buffer, 0, n_per_process-1);

	

	//if (processId == 1) {

	//	// Print the received data on each process
	//	printf("Process %d received: ", processId);
	//	for (int i = 0; i < n_per_process; i++) {
	//		printf("%d ", receive_buffer[i]);
	//	}
	//	printf("\n");

	//}


	



	//2 step 


	int add_to_i = n / (numProcess * numProcess);
	int* items_to_send_p0 = (int*)malloc(numProcess * sizeof(int));
	for (int i = 0, j = 0; i < n_per_process; i = i+add_to_i) {
		items_to_send_p0[j] = receive_buffer[i];
		j++;
	}


	

	/*if (processId == 1) {
		printf("Send ");
		printf("add_to_i:%d \n", add_to_i);

		for (int i = 0; i < numProcess; i = i++) {
			printf("i = %d , element=%d ", i, items_to_send[i]);
		}
	}*/

	
	printf("\n");


	int* receive_buffer_for_p_root_second_step = (int*)malloc(numProcess * numProcess * sizeof(int));

	MPI_Gather(items_to_send_p0, numProcess, MPI_INT, receive_buffer_for_p_root_second_step, numProcess, MPI_INT, 0, MPI_COMM_WORLD);


	int* items_to_send_step3 = (int*)malloc((numProcess - 1) * sizeof(int));
	if (processId == MASTER) {


		// Step 3 
		// local sort
		//printf("P*P process here \n");
		quickSort(receive_buffer_for_p_root_second_step, 0, (numProcess * numProcess) - 1);

		//for (int i = 0; i < numProcess*numProcess; i++) {
			//printf("%d ", receive_buffer_for_p_root_second_step[i]);
		//}
		//printf("\n");


		

		

		
		int add_to_i_for_step3 = numProcess + (int)floor(numProcess / 2.0) - 1;

		for (int i = add_to_i_for_step3, j = 0; i < numProcess * numProcess; i = i + numProcess) {
			//printf("i=%d ", i);
			items_to_send_step3[j] = receive_buffer_for_p_root_second_step[i];
			j++;
		}


		/*printf("Send step 3");
		printf("add_to_i_for_step3:%d \n", add_to_i_for_step3);*/

		/*for (int i = 0; i < numProcess-1; i = i++) {
			printf("element=%d ", items_to_send_step3[i]);
		}*/
		/*printf("---------------------------------\n");*/
	}

	MPI_Bcast(items_to_send_step3, numProcess - 1, MPI_INT, 0, MPI_COMM_WORLD);


	/*if (processId == 0) {
		for (int i = 0; i < numProcess - 1; i = i++) {
			printf("process id %d element=%d ,", processId, items_to_send_step3[i]);
		}
	}

	printf("\n");*/


	//Step 4
	int numPivots = numProcess - 1;

	int maxSegments = numPivots + 1;
	int** segments;
	segments = (int**)malloc(maxSegments * sizeof(int*));
	for (int i = 0; i < maxSegments; i++) {
		segments[i] = (int*)malloc(2 * sizeof(int));
	}

	multiPivotPartition(receive_buffer, n_per_process, items_to_send_step3, numPivots, segments);





	// Print the segments
	
	for (int i = 0; i < maxSegments; ++i) {
		int start = segments[i][0];
		int end = segments[i][1];

		if (start != -1 && end != -1) {
			//printf("processId :%d ", processId);
			//std::cout << "Segment " << i << ": ";
			//printArray(receive_buffer, start, end);
		}
	}


	int* partition_sizes = (int*)malloc((numProcess-1) * sizeof(int));
	for (int i = 0; i < maxSegments; ++i) {
		int start = segments[i][0];
		int end = segments[i][1];

		/*if (start != -1 && end != -1) {
			std::cout << "Segment " << i << ": ";
			printArray(receive_buffer, start, end);
		}*/

		if (start != -1 && end != -1) {
			partition_sizes[i] = end - start + 1;
		}
		else {
			partition_sizes[i] = 0;
		}
		
	}


	int* sendDisp = (int*)malloc(numProcess * sizeof(int));
	int* recvDisp = (int*)malloc(numProcess * sizeof(int));
	MPI_Barrier(MPI_COMM_WORLD);

	// sending number of elements that every process will get at last
	int* newPartitionSizes = (int*)malloc((numProcess - 1) * sizeof(int));
	MPI_Alltoall(partition_sizes, 1, MPI_INT, newPartitionSizes, 1, MPI_INT, MPI_COMM_WORLD);

	int totalSize=0;
	for (int i = 0; i < numProcess; i++) {
		//printf("procid=%d %d,",processId, newPartitionSizes[i]);
		totalSize += newPartitionSizes[i];
	}
	printf("\n");


	int* newPartitions = (int*)malloc(totalSize * sizeof(int));
	//*newPartitions = (double*)malloc(totalSize * sizeof(double));

	sendDisp[0] = 0;
	recvDisp[0] = 0; //Calculate the displacement relative to recvbuf, this displacement stores the data received from the process
	for (int i = 1; i < numProcess; i++) {
		sendDisp[i] = partition_sizes[i - 1] + sendDisp[i - 1];
		recvDisp[i] = newPartitionSizes[i - 1] + recvDisp[i - 1];
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Alltoallv(receive_buffer, partition_sizes, sendDisp, MPI_INT, newPartitions, newPartitionSizes, recvDisp, MPI_INT, MPI_COMM_WORLD);

	


	merge_elements(newPartitions, newPartitionSizes, numProcess);
	for (int i = 0; i < totalSize; i++) {
		printf("procid=%d %d,", processId, newPartitions[i]);
	}
	printf("\n");

	int totalelements = 0;

	for (int i = 0; i < numProcess; i++) {
		totalelements = totalelements + newPartitionSizes[i];
	}


	int* subListSizes = (int*)malloc(numProcess * sizeof(int));
	int* recvDisp2 = (int*)malloc(numProcess * sizeof(int));

	MPI_Gather(&totalelements, 1, MPI_INT, subListSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);


	if (processId == 0) {
		//printf("displacement ");
		recvDisp2[0] = 0;
		for (int i = 1; i < numProcess; i++) {
			recvDisp2[i] = subListSizes[i - 1] + recvDisp2[i - 1];
			
		}
		printf("\n ");

		/*for (int i = 0; i < numProcess; i++) {
			printf("i=%d %d ", i, recvDisp2[i]);
		}*/
		
	}

	MPI_Barrier(MPI_COMM_WORLD);
	int* last_array_sorted = (int*)malloc(n * sizeof(int));
	//Send each sorted sublist back to the root process
	MPI_Gatherv(newPartitions, totalelements, MPI_INT, last_array_sorted, subListSizes, recvDisp2, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (processId == 0) {

		/*printf("final result ");
		for (int i = 0; i < n; i++) {
			printf(" %d", last_array_sorted[i]);
		}

		printf("\n");*/


		FILE* file = fopen("output.txt", "w");
		if (file != NULL) {
			for (int i = 0; i < n; i++) {
				fprintf(file, "%d ", last_array_sorted[i]);
			}
			fclose(file);
		}
		else {
			fprintf(stderr, "Error opening file.\n");
		}
		//MPI_Finalize();
	}









	//if (processId == 0) {
	//	for (int i = 1; i < numProcess; i++) {
	//		MPI_Send(&arr[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	//	}

	//	temp = 0;
	//}
	//else {
	//	MPI_Recv(&temp, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//}
	//printf("process %d reiceved element %d\n", processId, temp);

	elapsedTime = MPI_Wtime() - startTime;
	if (processId == MASTER) {
		//printf("\nSequential Quicksort time is : \n");
		printf("Total Elapsed time during the process : %f\n", elapsedTime);
	}

	MPI_Finalize();
	return 0;


}