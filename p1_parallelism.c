#include <stdio.h>
#include <math.h>
#include <mpi.h>

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
        if (rank == 0) { //If we are the master process, get n and send to the others
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d", &n);
            for (i = 1; i < numprocs; i++) {
                MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        else { //If we are not the master process, receive n
            MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (n == 0) break;

        h = 1.0 / (double)n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) { //beginning at rank+1, we make the calculations with a step of numprocs
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x * x);
        }

        if (rank != 0) { //If we are not the master process, send the result of all the calculations to the master
            MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        else {
            for (i = 1; i < numprocs; i++) { //If we are the master process, receive the results from the other processes
                double temp;
                MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sum += temp;
            }
            pi = h * sum;
            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize(); //End all processes
    return 0;
}


patata ddddddddd

