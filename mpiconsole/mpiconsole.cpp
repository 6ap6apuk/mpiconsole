#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <Windows.h>
#include <vector>
#include <iomanip>
#include <limits>
#include <numbers>
#include <mutex>

using namespace std;

void task15(int rank, int size) {
    printf("I am %d process from %d processes!\n", rank, size);
}

void task16(int rank, int size) {
    if (rank == 0) 
        printf("%d processes!\n", size);
    else
    {
        if (rank % 2 == 0)
            printf("I am %d: FIRST!\n", rank);
        else
            printf("I am %d: SECOND!\n", rank);
    }
}

void task17(int rank, int size) {
    if (rank == 0) {
        int message = 12;
        MPI_Send(&message, 1, MPI_INT, 1, 5, MPI_COMM_WORLD);
    }
    else {
        MPI_Status status;
        int message;
        MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("receive message '%d'\n", message);
    }
}

void task18(int rank, int size) {
    if (rank == 0) {
        MPI_Status status;
        int procmessage = 0;
        MPI_Ssend(&procmessage, 1, MPI_INT, 1, 5, MPI_COMM_WORLD);
        MPI_Recv(&procmessage, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("[%d]: receive message '%d'\n", rank, procmessage);
    }
    else {
        MPI_Status status;
        int procmessage;
        MPI_Recv(&procmessage, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("[%d]: receive message '%d'\n", rank, procmessage);
        procmessage++;
        if (rank == size - 1) {
            MPI_Ssend(&procmessage, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
        }
        else {
            MPI_Ssend(&procmessage, 1, MPI_INT, procmessage + 1, 5, MPI_COMM_WORLD);
        }
    }
}

void task19(int rank, int size) {
    if (rank == 0) {
        MPI_Status status;
        int procmessage;
        for(int i = 1; i < size; i++) {
            MPI_Recv(&procmessage, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("receive message '%d'\n", procmessage);
        }
    }
    else {
        MPI_Send(&rank, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
}

void task20(int rank, int size) {
    if (rank == 0) {
        MPI_Request send_request;
        MPI_Status status;
        int message = 12;
        MPI_Isend(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &send_request);
    }
    else {
        int message;
        MPI_Status status;
        MPI_Request recv_request;
        MPI_Irecv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Wait(&recv_request, &status);
        printf("receive message '%d'\n", message);
    }
}

void task21(int rank, int size) {
    MPI_Request send_request, recv_request;
    MPI_Status status;
    int message = rank;

    if (rank == size - 1)
        MPI_Isend(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);
    else 
        MPI_Isend(&message, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD, &send_request);

    if (rank == 0)
        MPI_Irecv(&message, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &recv_request);
    else
        MPI_Irecv(&message, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &recv_request);

    MPI_Wait(&recv_request, &status);
    printf("receive message '%d'\n", message);
}

void task22(int rank, int size) {
    MPI_Request send_request, recv_request;
    MPI_Status status;
    int message;
    switch (rank) {
        case 1:
            MPI_Isend(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);
            MPI_Isend(&rank, 1, MPI_INT, 2, 2, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(&message, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);
            MPI_Irecv(&message, 1, MPI_INT, 2, 1, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);
            break;
        case 2:
            MPI_Isend(&rank, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &send_request);
            MPI_Isend(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(&message, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);
            MPI_Irecv(&message, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);
            break;
        case 0:
            MPI_Isend(&rank, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &send_request);
            MPI_Isend(&rank, 1, MPI_INT, 2, 2, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);
            MPI_Irecv(&message, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            printf("[%d]; receive message '%d'\n", rank, message);  
            break;
    }
}

void task23(int rank, int size) {
    string text;

    if (rank == 0) {
        cout << "Enter text: ";
        cin >> text;
        MPI_Bcast(&text, text.length() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    int count[26] = { 0 };

    for (int i = 0; i < text.length(); i++) {
        count[text[i] - 'a']++;
    }

    int total_count[26] = { 0 };
    MPI_Reduce(count, total_count, 26, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Summarized chars: " << endl;
        for (char c = 'a'; c <= 'z'; c++) {
            if (total_count[c - 'a'] > 0) {
                cout << c << " = " << total_count[c - 'a'] << endl;
            }
        }
    }
}

float calculatePartialSum(int localStart, int localEnd) {
    float partialSum = 0;
    for (int i = localStart; i < localEnd; i++) {
        partialSum += 4 / (1 + pow((i + 0.5) * 1.0 / localEnd, 2));
    }
    return partialSum;
}

double CalculatePi(int localStart, int localEnd, unsigned int totalN) {
    double sum = 0.0;
    double h = 1.0 / totalN;
    for (int i = localStart; i < localEnd; i++) {
        double xi = (i + 0.5) * h;
        sum += 4.0 / (1.0 + xi * xi);
    }
    return sum * h;
}

void task24(int rank, int size) {
    long double massPi = 0.0;
    int Nprec = 0;

    if (rank == 0) {
        cout << "Enter precision: ";
        cin >> Nprec;
        if (Nprec <= 0) {
            cerr << "Precision must be a positive integer." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        MPI_Bcast(&Nprec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&Nprec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (Nprec < size) {
        cerr << "Precision must be greater than or equal to the number of processes." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int localNprec = Nprec / size;
    int localStart = rank * localNprec;
    int localEnd = (rank + 1) * localNprec;

    long double localPartialPiSum = CalculatePi(localStart, localEnd, Nprec);

    long double globalSum = 0.0;
    MPI_Reduce(&localPartialPiSum, &globalSum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << std::fixed << std::setprecision(15);
        cout << "pi: " << globalSum << endl;
    }
}

void task25(int rank, int size) {
    int n = 0;
    int blockSize = 0;
    vector<double> localA, localB, localC;
    vector<double> A, B, C;

    if (rank == 0) {
        cout << "Matrix size n: ";
        cin >> n;
        blockSize = n * n / size;
        A.resize(n * n);
        B.resize(n * n);
        C.resize(n * n);

        cout << "A elems: " << endl;
        for (int i = 0; i < n * n; i++) {
            cin >> A[i];
        }

        cout << "B elems :" << endl;
        for (int i = 0; i < n * n; i++) {
            cin >> B[i];
        }

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        localA.resize(n * n / size);
        localB.resize(n * n / size);

        MPI_Scatter(A.data(), blockSize, MPI_DOUBLE, localA.data(), blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(B.data(), blockSize, MPI_DOUBLE, localB.data(), blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    localC.resize(n * n / size);

    vector<double> B_T(n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B_T[j * n + i] = B[i * n + j];
        }
    }

    int localRows = blockSize / n; 
    for (int i = 0; i < localRows; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                localC[i * n + j] += localA[i * n + k] * B_T[j * n + k];
            }
        }
    }

    if (rank == 0) {
        MPI_Gather(localC.data(), n * n / size, MPI_DOUBLE, C.data(), n * n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        cout << "Result matrix C:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << C[i * n + j] << " ";
            }
            cout << endl;
        }
    }
}

void task26(int rank, int size) {
    char message[11] = "";

    MPI_Comm new_comm;
    int color = rank % 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &new_comm);

    int new_rank, new_size;
    MPI_Comm_rank(new_comm, &new_rank);
    MPI_Comm_size(new_comm, &new_size);

    if (rank == 0) {
        while (true) {
            cout << "Enter message: ";
            cin >> message;
            if (strlen(message) > 0 && strlen(message) <= 10) {
                break;
            }
        }
    }

    if (rank % 2 == 0) {
        MPI_Bcast(message, 11, MPI_CHAR, 0, new_comm);

        cout << "MPI_COMM_WORLD: " << rank << " from " << size
            << ". New comm: " << new_rank << " from " << new_size
            << ". Message = " << message << endl;
    }
    else {
        cout << "MPI_COMM_WORLD: " << rank << " from " << size
            << ". New comm: no from no. Message = no" << endl;
    }
}

void task27(int rank, int size) {

    if (rank == 0) {
        int n;
        cout << "Enter the number of processes to spawn: ";
        cin >> n;

        MPI_Comm intercomm;
        vector<int> spawned_ranks(n);

        MPI_Comm_spawn("consoleChild.exe", MPI_ARGV_NULL, n, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);

        for (int i = 0; i < n; i++) {
            spawned_ranks[i] = i; 
        }

        MPI_Comm_free(&intercomm);
    }
    else {
        cout << "I am " << rank << " process from " << size << " processes! My parent is none." << endl;
    }
}

void task28(int rank, int size) {
    double massPi = 0.0;
    int Nprec = 0;
    int localNprec = 0;
    double localPartialPiSum = 0.0;
    MPI_Request send_request, recv_request;
    MPI_Status status;

    if (rank == 0) {
        cout << "Enter precision: ";
        cin >> Nprec;
        for (int i = 1; i < size; i++) {
            MPI_Isend(&Nprec, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &send_request);
        }
    }
    else {
        MPI_Irecv(&Nprec, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Wait(&recv_request, &status);
    }

    localNprec = Nprec / size;
    int localStart = rank * localNprec;
    int localEnd = (rank + 1) * localNprec;

    localPartialPiSum = CalculatePi(localStart, localEnd, Nprec);

    if (rank == 0) {
        massPi += localPartialPiSum;
        for (int i = 1; i < size; i++) {
            MPI_Irecv(&localPartialPiSum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, &status);
            massPi += localPartialPiSum;
        }
        cout << "Final massPi= " << massPi << endl;
    }
    else {
        MPI_Isend(&localPartialPiSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &send_request);
    }
}

void countTimeForPi(int rank, int size, unsigned int Nprec) {
    long double massPi = 0.0;

    if (Nprec < size) {
        cerr << "Precision must be greater than or equal to the number of processes." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    unsigned int localNprec = Nprec / size;
    unsigned int remainder = Nprec % size; 

    if (rank < remainder) {
        localNprec += 1;
    }

    unsigned int localStart = rank * (Nprec / size) + (rank < remainder ? rank : remainder);
    unsigned int localEnd = localStart + localNprec;

    long double localPartialPiSum = CalculatePi(localStart, localEnd, Nprec);

    long double globalSum = 0.0;
    MPI_Reduce(&localPartialPiSum, &globalSum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << std::fixed << std::setprecision(15);
        cout << "pi: " << globalSum << endl;
    }
}

void task29(int rank, int size) {
    //unsigned int Nprec = 100;
    //unsigned int Nprec = 10000000;
    unsigned int Nprec = 4000000000;

    MPI_Bcast(&Nprec, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    double startTime = MPI_Wtime();
    countTimeForPi(rank, size, Nprec);
    double endTime = MPI_Wtime();

    if (rank == 0) {
        cout << "Total time taken: " << (endTime - startTime) << " seconds" << endl;
    }
}

void task30(int rank, int size) {
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        std::cout << "Process " << rank << ", thread " << omp_get_thread_num() << ": iteration " << i << std::endl;
    }
}

std::mutex cout_mutex;

void task31(int rank, int size, int n) {
    int total_threads = n * size;

    // Параллельный блок OpenMP
#pragma omp parallel num_threads(n)
    {
        int thread_id = omp_get_thread_num();
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "I am " << thread_id << " thread from " << rank
                << " process. Number of hybrid threads = " << total_threads << std::endl;
        }
    }
}

void task32(int rank, int size) {
    unsigned int Nprec;

    if (rank == 0) {
        std::cout << "Enter precision: ";
        std::cin >> Nprec;
    }

    MPI_Bcast(&Nprec, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    int localN = Nprec / size;
    int localStart = rank * localN;
    int localEnd = (rank + 1) * localN;

    if (rank == size - 1) {
        localEnd = Nprec;
    }

    long double localSum = 0.0;

    #pragma omp parallel
    {
        long double threadSum = 0.0;
        long double h = 1.0 / Nprec;
    #pragma omp for
        for (int i = localStart; i < localEnd; i++) {
            long double xi = (i + 0.5) * h; 
            threadSum += 4.0 / (1.0 + xi * xi);
        }
    #pragma omp atomic
        localSum += threadSum;
    }

    localSum *= (1.0 / Nprec);

    long double totalSum = 0.0;
    MPI_Reduce(&localSum, &totalSum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << std::setprecision(15) << "Calculated value of pi: " << totalSum << std::endl;
    }
}

int main(int argc, char* argv[]) {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc > 1) {
        int task = atoi(argv[1]);
        switch (task) {
        case 15:
            task15(rank, size);
            break;
        case 16:
            task16(rank, size);
            break;
        case 17:
            task17(rank, size);
            break;
        case 18:
            task18(rank, size);
            break;
        case 19:
            task19(rank, size);
            break;
        case 20:
            task20(rank, size);
            break;
        case 21:
            task21(rank, size);
            break;
        case 22:
            task22(rank, size);
            break;
        case 23:
            task23(rank, size);
            break;
        case 24:
            task24(rank, size);
            break;
        case 25:
            task25(rank, size);
            break;
        case 26:
            task26(rank, size);
            break;
        case 27:
            task27(rank, size);
            break;
        case 28:
            task28(rank, size);
            break;
        case 29:
            task29(rank, size);
            break;
        case 30:
            task30(rank, size);
            break;
        case 31:
            task31(rank, size, 3);
            break;
        case 32:
            task32(rank, size);
            break;
        default:
            std::cerr << "Неверный номер.\n";
            break;
        }
    }

    MPI_Finalize();
    return 0;
}