/*!---------------------------------------------------------------------
\file
\brief assign mesh to fields and procs

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../solver/solver.h"

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
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*! 
\addtogroup PARALLEL 
*/
/*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief one procs info about his partition

<pre>                                                         m.gee 8/00
the partition of one proc (all discretizations)
the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;                      

/*!---------------------------------------------------------------------

\brief assign fields to procs                                             

<pre>                                                        m.gee 5/01 
-allocates variable partition
-allocates variable pdis in partition
fills nodes and elements to the pdis variable in two different styles
-cuts are through elements, every node is exactly on 1 proc
-cuts are through nodes, every element is exactly on 1 proc (not used yet)
</pre>

\return void                                               

------------------------------------------------------------------------*/
void part_assignfield()
{
INT        i,j,k,kk;
INT        counter,counter2;
INT        part;
INT        proc;
INT        isbou;
FIELD     *actfield;
INTRA     *actintra;
INT        imyrank;
INT        inprocs;
NODE      *actnode;
PARTITION *actpart;
SOLVAR    *actsolv;

#ifdef DEBUG 
dstrc_enter("part_assignfield");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------- every proc alloc numfld PARTITIONs */
partition = (PARTITION*)CCACALLOC(genprob.numfld,sizeof(PARTITION));
if (partition==NULL) dserror("Allocation of PARTITION failed");
/*----------------------------------------------------- loop all fields */
for (i=0; i<genprob.numfld; i++)
{
   actsolv  = &(solv[i]);
   actfield = &(field[i]);
   /*------------------ every partition[i] allocates one discretization */
   partition[i].ndis = actfield->ndis;
   partition[i].pdis = (PARTDISCRET*)CCACALLOC(partition[i].ndis,sizeof(PARTDISCRET));
   if (!partition[i].pdis) dserror("Allocation of memory failed");
   /*-------------------------------------------------------------------*/
   actpart  = &(partition[i]);
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
#else
   actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
   if (!actintra) dserror("Allocation of INTRA failed");
   actintra->intra_fieldtyp = actfield->fieldtyp;
   actintra->intra_rank   = 0;
   actintra->intra_nprocs   = 1;
#endif
   imyrank  = actintra->intra_rank;
   inprocs  = actintra->intra_nprocs;
   actpart->fieldtyp = actfield->fieldtyp;
   /*--------------- if there is only one proc, theres is nothing to do */
   if (inprocs<=1)
   {
      for (kk=0;kk<actpart->ndis;kk++)
      {
         actpart->pdis[kk].numnp       = actfield->dis[kk].numnp;
         actpart->pdis[kk].numele      = actfield->dis[kk].numele;
         actpart->pdis[kk].bou_numnp   = 0;
         actpart->pdis[kk].bou_numele  = 0;
         actpart->pdis[kk].bou_element = NULL;
         actpart->pdis[kk].bou_node    = NULL;
         actpart->pdis[kk].element = (ELEMENT**)CCACALLOC(actpart->pdis[kk].numele,sizeof(ELEMENT*));
         actpart->pdis[kk].node    = (NODE**)CCACALLOC(actpart->pdis[kk].numnp,sizeof(NODE*));
         if (actpart->pdis[kk].element==NULL) dserror("Allocation of element pointer in PARTITION failed");
         if (actpart->pdis[kk].node==NULL)    dserror("Allocation of node pointer in PARTITION failed");
         for (j=0; j<actfield->dis[kk].numele; j++) 
         actpart->pdis[kk].element[j] = &(actfield->dis[kk].element[j]);
         for (j=0; j<actfield->dis[kk].numnp; j++)  
         actpart->pdis[kk].node[j] = &(actfield->dis[kk].node[j]);      
      }
   }
   else
   {
      /*-------------- do only for procs wqithin the intra-communicator */
      if (actintra->intra_fieldtyp==none) continue;
      /*------------------------------------ check for typ of partition */
      part=0;
      if (actsolv->parttyp == cut_elements) part = 1;
      if (actsolv->parttyp == cut_nodes)    part = 2;
      if (part==0) dserror("Typ of Partitioning unknown");


      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      /* parition field, such that cuts are through elements and each node
         is assigned exactly one proc.
         Attention: From now on information on different procs can differ!!
         severely !! (Additive Schwartz)*/
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      for (kk=0;kk<actpart->ndis;kk++) /* loop over discretisations */
      {
      if (part==1)
      {
      /*----------- loop and count elements, do pointers to my elements */
         counter=0;
         for (j=0; j<actfield->dis[kk].numele; j++)
         {
            for (k=0; k<actfield->dis[kk].element[j].numnp; k++)
            {
               if (actfield->dis[kk].element[j].node[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->pdis[kk].numele  = counter;
         actpart->pdis[kk].element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
         if (!actpart->pdis[kk].element) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[kk].numele; j++)
         {
            for (k=0; k<actfield->dis[kk].element[j].numnp; k++)
            {
               if (actfield->dis[kk].element[j].node[k]->proc == imyrank)
               {
                  actpart->pdis[kk].element[counter] = &(actfield->dis[kk].element[j]);
                  counter++;
                  break;
               }
            }
         }
      /*------------------ loop and count nodes, do pointers to my nodes */   
         counter=0;
         for (j=0; j<actfield->dis[kk].numnp; j++)
         {
            if (actfield->dis[kk].node[j].proc == imyrank) counter++;
         }
         actpart->pdis[kk].numnp = counter;
         actpart->pdis[kk].node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
         if (!actpart->pdis[kk].node) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[kk].numnp; j++)
         {
            if (actfield->dis[kk].node[j].proc == imyrank) 
            {
               actpart->pdis[kk].node[counter] = &(actfield->dis[kk].node[j]);
               counter++;
            }
         }
       /*----paritioning this way does not need inner and boundary nodes */
         actpart->pdis[kk].inner_numnp=0;
         actpart->pdis[kk].bou_numnp=0;
         actpart->pdis[kk].inner_node=NULL;
         actpart->pdis[kk].bou_node=NULL;
       /*----------------------------- now count the pure inner & boundary
                                                                elements */
       /*------------ Juhu !!! This is the first real parallel loop !!!! */
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[kk].numele; j++)
         {
            isbou=0;
            proc = actpart->pdis[kk].element[j]->node[0]->proc;
            for (k=1; k<actpart->pdis[kk].element[j]->numnp; k++)
            {
               if (actpart->pdis[kk].element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) counter++;
            else          counter2++;
         }
         actpart->pdis[kk].inner_numele = counter;
         actpart->pdis[kk].bou_numele   = counter2;
         actpart->pdis[kk].inner_element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
         actpart->pdis[kk].bou_element = (ELEMENT**)CCACALLOC(counter2,sizeof(ELEMENT*));
         if (actpart->pdis[kk].inner_element==NULL || actpart->pdis[kk].bou_element==NULL)
         dserror("Allocation of PARTITION to ELEMENT pointer failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[kk].numele; j++)
         {
            isbou=0;
            proc = actpart->pdis[kk].element[j]->node[0]->proc;
            for (k=1; k<actpart->pdis[kk].element[j]->numnp; k++)
            {
               if (actpart->pdis[kk].element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) 
            {
               actpart->pdis[kk].inner_element[counter] = actpart->pdis[kk].element[j];
               counter++;
            }
            else
            {
               actpart->pdis[kk].bou_element[counter2]  = actpart->pdis[kk].element[j];
               counter2++;
            }
         }
      } /* end of part==1 */


      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      /* partition field, such that cuts are through nodes and boundary nodes
         are on the common procs */
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      if (part==2)
      {
      /* loop and count elements, every element belongs exactly to one proc */
         counter=0;
         for (j=0; j<actfield->dis[kk].numele; j++)
         {
            if (actfield->dis[kk].element[j].proc == imyrank) counter++;
         }
         actpart->pdis[kk].numele = counter;
         actpart->pdis[kk].element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
         if (!actpart->pdis[kk].element) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[kk].numele; j++)
         {
            if (actfield->dis[kk].element[j].proc == imyrank) 
            {
               actpart->pdis[kk].element[counter] = &(actfield->dis[kk].element[j]);
               counter++;
            }
         }
      /* loop and counter nodes, every node belonging to one of the partitions
         elements belongs to this partition */
         counter=0;
         for (j=0; j<actfield->dis[kk].numnp; j++)
         {
            actnode = &(actfield->dis[kk].node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->pdis[kk].numnp = counter;
         actpart->pdis[kk].node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
         if (!actpart->pdis[kk].node) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[kk].numnp; j++)
         {
            actnode = &(actfield->dis[kk].node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  actpart->pdis[kk].node[counter] = actnode;
                  counter++;
                  break;
               }
            }
         }
      /*----------------- there is no such thing as boundary and inner
                                                               elements */
         actpart->pdis[kk].bou_numele=0;
         actpart->pdis[kk].inner_numele=0;
         actpart->pdis[kk].bou_element=NULL;
         actpart->pdis[kk].inner_element=NULL;
      /*---------------------------- count the inner and boundary nodes */      
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[kk].numnp; j++)
         {
            isbou=0;
            actnode = actpart->pdis[kk].node[j];
            proc = actnode->element[0]->proc;
            for (k=1; k<actnode->numele; k++)
            {
               if (proc != actnode->element[k]->proc)
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==1) counter2++;
            else          counter++;
         }
         actpart->pdis[kk].inner_numnp=counter;
         actpart->pdis[kk].bou_numnp  =counter2;
         actpart->pdis[kk].inner_node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
         actpart->pdis[kk].bou_node   = (NODE**)CCACALLOC(counter2,sizeof(NODE*));
         if (actpart->pdis[kk].inner_node==NULL || actpart->pdis[kk].bou_node==NULL)
         dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[kk].numnp; j++)
         {
            isbou=0;
            actnode = actpart->pdis[kk].node[j];
            proc = actnode->element[0]->proc;
            for (k=1; k<actnode->numele; k++)
            {
               if (proc != actnode->element[k]->proc)
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==1) 
            {
               actpart->pdis[kk].bou_node[counter2] = actnode;
               counter2++;
            }
            else          
            {
               actpart->pdis[kk].inner_node[counter] = actnode;
               counter++;
            }
         }
         
      } /* end of part==2 */
      } /* end loop over discretisations */


   } /* end of inprocs > 1 */
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
}/* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of part_assignfield */

/*! @} (documentation module close)*/
