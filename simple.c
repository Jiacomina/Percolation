#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
        int rank, size;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // reports number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &size); // reports the rank, a number between 0 and size-1 identifying the calling process
        printf("I am %d of %d\n", rank, size);
        MPI_Finalize();
        return 0;
}
