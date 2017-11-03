#include <mpi.h>
#include <stdio.h>

int fib(int n) {
        return n <= 2 ? 1 : fib(n-1)+fib(n-2);
}

 

int main(int argc, char** argv) {
        MPI_Init(NULL, NULL);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);
        printf("processor %s, rank %d/%d processors. fib(42) = %d \n",processor_name, world_rank, world_size, fib(42));

        MPI_Finalize();

}
