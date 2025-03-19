#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Получаем родительский процесс
    int parent_rank = 0; // У всех дочерних процессов родитель - 0
    cout << "I am " << rank << " process from " << size << " processes! My parent is " << parent_rank << "." << endl;

    MPI_Finalize();
    return 0;
}