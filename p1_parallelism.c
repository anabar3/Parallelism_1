#include <stdio.h>
#include <math.h>
#include <mpi.h>

int MPI_FlattreeCollective(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int rank, numprocs, i;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numprocs);
    if (rank != root) { //If we are not the root process, send the result of all the calculations to the root
        MPI_Send(&sendbuf, count, datatype, root, 0, MPI_COMM_WORLD);
    }
    else {
        for (i = 1; i < numprocs; i++) { //If we are the root process, receive the results from the other processes
            void *temp;
            MPI_Recv(&temp, count, datatype, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Apply the operation to the result

        }
    }
}

int main(int argc, char *argv[]) {
    //Initial variables
    int i, n, rank, numprocs;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    while (1) {
        if (rank == 0) { //If we are the master process, get n
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d", &n);
        }

        //Broadcast n to all processes
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (n == 0) break;

        h = 1.0 / (double)n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) { //beginning at rank+1, we make the calculations with a step of numprocs
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x * x);
        }

        //Collective operation to sum all the results
        MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            pi *= h;
            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    return 0;
}