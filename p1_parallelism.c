#include <stdio.h>
#include <math.h>
#include <mpi.h>

int MPI_FlattreeCollective(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int rank, numprocs, i;
    double pi;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numprocs);

    //Check all possible MPI errors
    if (sendbuf == NULL || recvbuf == NULL)     return MPI_ERR_BUFFER;
    if (count < 0)                              return MPI_ERR_COUNT;
    if (datatype != MPI_DOUBLE)                 return MPI_ERR_TYPE;
    if (op != MPI_SUM)                          return MPI_ERR_OP;
    if (root < 0 || root >= numprocs)           return MPI_ERR_ROOT;
    if (comm == MPI_COMM_NULL)                  return MPI_ERR_COMM;


    if (rank != root) { //If we are not the root process, send the result of all the calculations to the root
        MPI_Send(sendbuf, count, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
    }
    else {
        pi = *(double *)sendbuf;
        for (i = 1; i < numprocs; i++) { //If we are the root process, receive the results from the other processes
            double temp;
            MPI_Recv(&temp, count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Apply the operation to the result
            pi += temp;
        }
        *(double *)recvbuf = pi;
    }
    return MPI_SUCCESS;
}

int MPI_BinomialColective(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int rank, numprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numprocs);

    for (int i = 1; i <= log2(numprocs); i++) { //we use log2 iterations because of the tree structure
        int addFactor = pow(2, i - 1); //following the formula given to calculate the pair of the process
        if (rank < addFactor) { //If we are in the position to send
            //printf("Rank %d sending to %d\n", rank, rank + addFactor);
            MPI_Send(buffer, count, datatype, rank + addFactor, 0, comm);
        }
        else if(rank-addFactor < addFactor){ //if we don't have to send, we receive ONLY if our pair is in the position to send
            //printf("Rank %d receiving from %d\n", rank, rank - addFactor);
            MPI_Recv(buffer, count, datatype, rank - addFactor, 0, comm, MPI_STATUS_IGNORE);
        }
    }

    return MPI_SUCCESS;
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
        MPI_BinomialColective(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (n == 0) break;

        h = 1.0 / (double)n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) { //beginning at rank+1, we make the calculations with a step of numprocs
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x * x);
        }

        //Collective operation to sum all the results
        MPI_FlattreeCollective(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            pi *= h;
            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    return 0;
}