#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | main                                                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
void main(int argc, char *argv[])
{
/*------------------------------- the size of the mpi buffer in doubles */
/*                                     (do not touch if things go well) */
#ifdef PARALLEL 
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
MPI_Comm_size(MPI_COMM_WORLD, &par.nprocs);
#else
par.myrank=0;
par.nprocs=1;
#endif
if (par.myrank==0)
{
   printf("\n");
   printf("************************\n");
   printf("*  CARAT TEST VERSION  *\n");
   printf("*                      *\n");
#ifdef PARALLEL 
   printf("*   parallel version   *\n");
#else /* SEQUENTIEL */
   printf("*  sequentiel version  *\n");
#endif
   printf("*                      *\n");
   printf("*                      *\n");
   printf("************************\n");
#ifdef PARALLEL 
   printf("number of processors: %d\n",par.nprocs);
#else
   printf("-sequentiell executable-\n");
#endif
}

/*----------------------------------------------- everything is in here */
ntam(argc,argv);
/*----------------------------------------------------------------------*/

#ifdef PARALLEL
MPI_Barrier(MPI_COMM_WORLD);
printf("processor %d finished normally\n",par.myrank);
MPI_Finalize();
#else
printf("processor %d finished normally\n",par.myrank);
#endif
return;
} /* end of main */
