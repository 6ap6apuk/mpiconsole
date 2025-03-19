#include "mpi.h"
#include <iostream>
#include <Windows.h>
#include <vector>

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

void calculatePartialSum(int localStart, int localEnd, float& partialSum) {
    for (int i = localStart; i < localEnd; i++) {
        partialSum += (4 / (1 + pow((i + 0.5) * 1.0 / localEnd, 2))) / localEnd;
    }
}

void task24(int rank, int size) {
    float massPi;
    int Nprec = 0;

    if (rank == 0) {
        cout << "Enter precision: ";
        cin >> Nprec;
        MPI_Bcast(&Nprec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int localNprec = Nprec / size;
    float localPartialPiSum = 0;
    calculatePartialSum(rank * localNprec, (rank + 1) * localNprec, localPartialPiSum);

    float globalSum = 0;
    MPI_Reduce(&localPartialPiSum, &globalSum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "pi: " << globalSum << endl;
    }
}

void task25(int rank, int size) {
    int n = 0;
    int start = 10;
    int end = 20;
    vector<double> A, B;

    if (rank == 0) {
        cout << "Matrix size n: ";
        cin >> n;
        A.resize(n * n);
        B.resize(n * n);

        cout << "A elems: " << endl;
        for (int i = 0; i < n * n; i++) {
            cin >> A[i];
        }

        cout << "B elems :" << endl;
        for (int i = 0; i < n * n; i++) {
            cin >> B[i];
        }

        // Распределяем размеры матриц
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(A.data(), n * n / size, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(B.data(), n * n / size, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        A.resize(n * n / size);
        B.resize(n * n / size);

        MPI_Scatter(nullptr, n * n / size, MPI_DOUBLE, A.data(), n * n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(nullptr, n * n / size, MPI_DOUBLE, B.data(), n * n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Создаем массив для локального результата
    vector<double> localC(n * n / size, 0.0);
    double* B_T = new double[n * n]; // Временная матрица для транспонированной B

    // Транспонируем матрицу B
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B_T[j * n + i] = B[j * n + i]; // Транспонируем B для удобства доступа
        }
    }

    // Умножаем локальные матрицы A на транспонированную B
    for (int i = 0; i < n / size; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                localC[i * n + j] += A[i * n + k] * B_T[j * n + k];
            }
        }
    }

    // Сбор результатов в матрицу C на процессе 0
    vector<double> C(n * n);
    MPI_Gather(localC.data(), n * n / size, MPI_DOUBLE, C.data(), n * n / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Процесс 0 выводит конечный результат
    if (rank == 0) {
        cout << "Result C:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << C[i * n + j] << " ";
            }
            cout << endl;
        }
    }
}

void task26(int rank, int size) {

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
        default:
            std::cerr << "Неверный номер.\n";
            break;
        }
    }

    MPI_Finalize();
    return 0;
}