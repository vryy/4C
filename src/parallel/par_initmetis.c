#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "../ale3/ale3.h"
#ifdef PARALLEL 
#ifdef SUN
#include "../../../lib_sun/metis-4.0/Lib/metis.h"
#else
#include "/bau/stat33/users/statik/lib/METIS/metis.h"
#endif
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      

/*----------------------------------------------------------------------*
 |  do initial partitioning of fields                    m.gee 5/01     |
 |  the paritioning of all fields is performed on all procs,            |
 |  so at least veryone nows which piece of every field is owned by who |
 |  -> this routine lives in theMPI_COMM_WORLD space                    |
 *----------------------------------------------------------------------*/
void part_fields()
{
int  i,j,k,l,m,n;
int      counter;
int      adjcounter;
long int max,min;
int      proc;

INTRA   *actintra;
int      imyrank;
int      inprocs;
FIELD   *actfield;
NODE    *actnode;
ELEMENT *actele;
ARRAY    stack;

ARRAY    xadj[MAXFIELD];
ARRAY    adjncy[MAXFIELD];
ARRAY    vwgt[MAXFIELD];

int      ione=1;
int      options[5];
int      numflag=0;
int      edgecut;
int      wgtflag=2;
int      nparts;
ARRAY    part;
ARRAY    part_proc;
ARRAY    ele_per_proc;
#ifdef DEBUG 
dstrc_enter("part_fields");
#endif
/*--------- sequentiell version assign all elements and nodes to proc 0 */
/* the partitioning is not performed, all nodes and elements are assigned to
   proc 0, but the graphs of the fields are build in sequentiell and
   parallel version
*/
if (par.nprocs<=1)
{
   for (i=0; i<genprob.numfld; i++)
   {
      actfield = &(field[i]);
      for (j=0; j<actfield->dis[0].numele; j++) 
      actfield->dis[0].element[j].proc = 0;
      for (j=0; j<actfield->dis[0].numnp; j++)  
      actfield->dis[0].node[j].proc    = 0;
   }
}
imyrank=0;
inprocs=1;
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
/*------ check proc belonging to this intra-communicator group of procs */
   if (actintra->intra_fieldtyp==none) continue;
   imyrank  = actintra->intra_rank;
   inprocs  = actintra->intra_nprocs;
#endif
/*----------------------------------------------------------------------*/
   actfield = &(field[i]);
/*---------------------------- init the local numbering of the elements */   
/*------------------------------ numbering is c style, starts with zero */
   counter=0;
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actfield->dis[0].element[j].Id_loc = counter;
      counter++;
   }   
/*--------------------------------------init the numbering of the nodes */   
   counter=0; 
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actfield->dis[0].node[j].Id_loc = counter;
      counter++;
   }
/*------------------------------------------------- calculate the graph */   
/*--------------------------------------- size of ARRAY xadj is numnp+1 */   
   amdef("xadj",&(xadj[i]),(actfield->dis[0].numnp+1),1,"IV");
   amzero(&(xadj[i]));
/*------------------------------------- the vertex weights of the graph */
   amdef("vwgt",&(vwgt[i]),(actfield->dis[0].numnp)  ,1,"IV");
   aminit(&(vwgt[i]),&ione);
/*----------------------------- size of array adjncy has to be computed */  
   amdef("stack",&stack,1,1,"IV");
   amzero(&stack);
/*------------------------------------- estimate size of adjncy to 1000 */
   amdef("adjncy",&(adjncy[i]),1000,1,"IV");
   amzero(&(adjncy[i]));
   adjcounter=0;
/*----------------------------------------------- loop the fields nodes */
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      counter=0;
      actnode = &(actfield->dis[0].node[j]);
/*--------------------------------------------- determine size of stack */
      for (k=0; k<actnode->numele; k++)
      {
         counter += actnode->element[k]->numnp;
      }
/*------------------------------------------------------ allocate stack */
      amdel(&stack);
      amdef("stack",&stack,counter,1,"IV");
      amzero(&stack);
      counter=0;
/*-------------------------------------------------- put nodes on stack */
      for (k=0; k<actnode->numele; k++)
      {
         actele = actnode->element[k];
         for (l=0; l<actele->numnp; l++)
         {
            stack.a.iv[counter] = actele->node[l]->Id_loc;
            counter++;
         }
      }
/*--------------------------------------------- delete doubles on stack */         
      for (k=0; k<stack.fdim; k++)
      {
         counter=stack.a.iv[k];
         if (counter==-1) continue;
         for (l=k+1; l<stack.fdim; l++)
         {
            if (stack.a.iv[l]==counter) stack.a.iv[l]=-1;
         }
      }      
/*-------------------------------------- count number of nodes on stack */      
      counter=0;
      for (k=0; k<stack.fdim; k++)
      {
         if (stack.a.iv[k] != -1) counter++;
      }
/*--------------- number of vertices is from actnode to each other node */
      counter--;
/*------------------------------- put edges on stack to xadj and adjncy */
      xadj[i].a.iv[actnode->Id_loc]   = adjcounter;
      xadj[i].a.iv[actnode->Id_loc+1] = adjcounter+counter;
      for (k=0; k<stack.fdim; k++)
      {
         if (stack.a.iv[k] != -1 && stack.a.iv[k]!= actnode->Id_loc)
         {
            if (adjcounter<adjncy[i].fdim)
            {
               adjncy[i].a.iv[adjcounter] = stack.a.iv[k];
               adjcounter++;
            }
            else
            {
               amredef(&(adjncy[i]),(adjncy[i].fdim+1000),1,"IV");
               adjncy[i].a.iv[adjcounter] = stack.a.iv[k];
               adjcounter++;
            }
         } 
      }
   }  /* end of loop over nodes */
   amredef(&(adjncy[i]),(adjcounter),1,"IV");    
   amdel(&stack);
/*------------------------------------------- do not partition ale field */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*   if (actfield->fieldtyp==ale) continue;*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*-------------------------------- do not partition sequentiell version */
   if (inprocs<=1) 
   {
      amdel(&(adjncy[i]));
      amdel(&(xadj[i]));
      amdel(&(vwgt[i]));
      continue;
   }
/*----------------- now allocate the rest of the arrays needed by metis */   
   amdef("part",&part,actfield->dis[0].numnp,1,"IV");   
/*---- set the default options of metis (play with other options later) */
   options[0]=0;
   options[1]=3;
   options[2]=1;
   options[3]=1;
   options[4]=0;
/*---------------------------------------------------------- call metis */
   nparts=inprocs;
/* partitioning is not necessarily deterministic, so it is only performed
   by proc 0 and then broadcasted                                       */
   if (imyrank==0)
   {
      if (nparts < 8) /*----- better for a smaller number of partitions */
      {
#ifdef PARALLEL 
      METIS_PartGraphRecursive(
                               &(actfield->dis[0].numnp),
                               &(xadj[i].a.iv[0]),
                               &(adjncy[i].a.iv[0]),
                               &(vwgt[i].a.iv[0]),
                               NULL,
                               &wgtflag,
                               &numflag,
                               &nparts,
                               options,
                               &edgecut,
                               &(part.a.iv[0])
                              );
#endif
      }
      else /*----------------- better for a larger number of partitions */
      {
#ifdef PARALLEL 
      METIS_PartGraphKway(
                          &(actfield->dis[0].numnp),
                          &(xadj[i].a.iv[0]),
                          &(adjncy[i].a.iv[0]),
                          &(vwgt[i].a.iv[0]),
                          NULL,
                          &wgtflag,
                          &numflag, 
                          &nparts,
                          options,
                          &edgecut,
                          &(part.a.iv[0])
                         );
#endif
      }
   }
/*-------------------------------------- broadcast partitioning results */
#ifdef PARALLEL 
   MPI_Bcast(&(part.a.iv[0]),
             part.fdim,
             MPI_INT,
             0,
             actintra->MPI_INTRA_COMM);
#endif
/*------------------------------------ put the proc number to the nodes */
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actfield->dis[0].node[j].proc=part.a.iv[j];
   }
   amdel(&part);
/*---------- check the nodes of all elements and assign element to proc */   
   amdef("part",&part,nparts,1,"IV");
   amzero(&part);
   amdef("part_proc",&part_proc,nparts,1,"IV");
   amzero(&part_proc);
   amdef("ele_proc",&ele_per_proc,nparts,1,"IV");
   amzero(&ele_per_proc);
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      amzero(&part);
      amzero(&part_proc);
      for (k=0; k<actele->numnp; k++)
      {
         part.a.iv[actele->node[k]->proc]+=1;
      }
      proc = 0;
      max  = 0;
      for (k=0; k<part.fdim; k++)
      {
         if (part.a.iv[k]>=max) 
         {
            max  = part.a.iv[k];
            proc = k;
         }
      }
      counter=0;
      for (k=0; k<part.fdim; k++)
      {
         if (part.a.iv[k]==max) 
         {
            counter++;
            part_proc.a.iv[k]=1;
         }
      }
      if (counter>1)
      {
         min=VERYLARGEINT;
         for (k=0; k<part_proc.fdim; k++)
         {
            if (part_proc.a.iv[k]==1 && ele_per_proc.a.iv[k]<min)
            {
               min = ele_per_proc.a.iv[k];
               proc = k;
            }
         }
      }
      actele->proc = proc;
      ele_per_proc.a.iv[proc]+=1;
   }
   amdel(&part);
   amdel(&part_proc);
   amdel(&ele_per_proc);
   amdel(&(xadj[i]));
   amdel(&(adjncy[i]));
   amdel(&(vwgt[i]));
}/*-------------------------------------------- end of loop over fields */
/*---------------------------------------- assign procs to ale elements */
/* NOTICE: This is not ideal, as there is not an ale element to every fluid,
   so the partitioning as it is done for the fluid elements is not necesarily
   balanced for the ale field. Also this can only be done with compatible
   ale fields */  
#if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (inprocs>1)
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   if (actfield->fieldtyp!=ale) continue;
/*--------------------------------------------------- loop ale elements */   
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      actele->proc = actele->e.ale3->my_fluid->proc;
/*----------------------------------------------- loop nodes of element */
      for (k=0; k<actele->numnp; k++)
      {
         actele->node[k]->proc = actele->e.ale3->my_fluid->node[k]->proc;
      }
   }
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of part_fields */
