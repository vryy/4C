/*!---------------------------------------------------------------------
\file
\brief assign mesh to fields and procs

---------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution.h"

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
int        i,j,k,l;
int        counter,counter2;
int        part;
int        proc;
int        isbou;
FIELD     *actfield;
INTRA     *actintra;
int        imyrank;
int        inprocs;
NODE      *actnode;
ELEMENT   *actele;
PARTITION *actpart;
SOLVAR    *actsolv;

#ifdef DEBUG 
dstrc_enter("part_assignfield");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------- every proc alloc numfld PARTITIONs */
partition = (PARTITION*)CALLOC(genprob.numfld,sizeof(PARTITION));
if (partition==NULL) dserror("Allocation of PARTITION failed");
/*----------------------------------------------------- loop all fields */
for (i=0; i<genprob.numfld; i++)
{
   /*------------------ every partition[i] allocates one discretization */
   partition[i].ndis = 1;
   partition[i].pdis = (PARTDISCRET*)CALLOC(partition[i].ndis,sizeof(PARTDISCRET));
   if (!partition[i].pdis) dserror("Allocation of memory failed");
   /*-------------------------------------------------------------------*/
   actsolv  = &(solv[i]);
   actfield = &(field[i]);
   actpart  = &(partition[i]);
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
#else
   actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
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
      actpart->pdis[0].numnp       = actfield->dis[0].numnp;
      actpart->pdis[0].numele      = actfield->dis[0].numele;
      actpart->pdis[0].bou_numnp   = 0;
      actpart->pdis[0].bou_numele  = 0;
      actpart->pdis[0].bou_element = NULL;
      actpart->pdis[0].bou_node    = NULL;
      actpart->pdis[0].element = (ELEMENT**)CALLOC(actpart->pdis[0].numele,sizeof(ELEMENT*));
      actpart->pdis[0].node    = (NODE**)CALLOC(actpart->pdis[0].numnp,sizeof(NODE*));
      if (actpart->pdis[0].element==NULL) dserror("Allocation of element pointer in PARTITION failed");
      if (actpart->pdis[0].node==NULL)    dserror("Allocation of node pointer in PARTITION failed");
      for (j=0; j<actfield->dis[0].numele; j++) 
      actpart->pdis[0].element[j] = &(actfield->dis[0].element[j]);
      for (j=0; j<actfield->dis[0].numnp; j++)  
      actpart->pdis[0].node[j] = &(actfield->dis[0].node[j]);
      
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
         Attention: From now on information on different procs can difer!!
         severely !! (Additive Schwartz)*/
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      if (part==1)
      {
      /*----------- loop and count elements, do pointers to my elements */
         counter=0;
         for (j=0; j<actfield->dis[0].numele; j++)
         {
            for (k=0; k<actfield->dis[0].element[j].numnp; k++)
            {
               if (actfield->dis[0].element[j].node[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->pdis[0].numele  = counter;
         actpart->pdis[0].element = (ELEMENT**)CALLOC(counter,sizeof(ELEMENT*));
         if (!actpart->pdis[0].element) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[0].numele; j++)
         {
            for (k=0; k<actfield->dis[0].element[j].numnp; k++)
            {
               if (actfield->dis[0].element[j].node[k]->proc == imyrank)
               {
                  actpart->pdis[0].element[counter] = &(actfield->dis[0].element[j]);
                  counter++;
                  break;
               }
            }
         }
      /*------------------ loop and count nodes, do pointers to my nodes */   
         counter=0;
         for (j=0; j<actfield->dis[0].numnp; j++)
         {
            if (actfield->dis[0].node[j].proc == imyrank) counter++;
         }
         actpart->pdis[0].numnp = counter;
         actpart->pdis[0].node = (NODE**)CALLOC(counter,sizeof(NODE*));
         if (!actpart->pdis[0].node) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[0].numnp; j++)
         {
            if (actfield->dis[0].node[j].proc == imyrank) 
            {
               actpart->pdis[0].node[counter] = &(actfield->dis[0].node[j]);
               counter++;
            }
         }
       /*----paritioning this way does not need inner and boundary nodes */
         actpart->pdis[0].inner_numnp=0;
         actpart->pdis[0].bou_numnp=0;
         actpart->pdis[0].inner_node=NULL;
         actpart->pdis[0].bou_node=NULL;
       /*----------------------------- now count the pure inner & boundary
                                                                elements */
       /*------------ Juhu !!! This is the first real parallel loop !!!! */
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[0].numele; j++)
         {
            isbou=0;
            proc = actpart->pdis[0].element[j]->node[0]->proc;
            for (k=1; k<actpart->pdis[0].element[j]->numnp; k++)
            {
               if (actpart->pdis[0].element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) counter++;
            else          counter2++;
         }
         actpart->pdis[0].inner_numele = counter;
         actpart->pdis[0].bou_numele   = counter2;
         actpart->pdis[0].inner_element = (ELEMENT**)CALLOC(counter,sizeof(ELEMENT*));
         actpart->pdis[0].bou_element = (ELEMENT**)CALLOC(counter2,sizeof(ELEMENT*));
         if (actpart->pdis[0].inner_element==NULL || actpart->pdis[0].bou_element==NULL)
         dserror("Allocation of PARTITION to ELEMENT pointer failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[0].numele; j++)
         {
            isbou=0;
            proc = actpart->pdis[0].element[j]->node[0]->proc;
            for (k=1; k<actpart->pdis[0].element[j]->numnp; k++)
            {
               if (actpart->pdis[0].element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) 
            {
               actpart->pdis[0].inner_element[counter] = actpart->pdis[0].element[j];
               counter++;
            }
            else
            {
               actpart->pdis[0].bou_element[counter2]  = actpart->pdis[0].element[j];
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
         for (j=0; j<actfield->dis[0].numele; j++)
         {
            if (actfield->dis[0].element[j].proc == imyrank) counter++;
         }
         actpart->pdis[0].numele = counter;
         actpart->pdis[0].element = (ELEMENT**)CALLOC(counter,sizeof(ELEMENT*));
         if (!actpart->pdis[0].element) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[0].numele; j++)
         {
            if (actfield->dis[0].element[j].proc == imyrank) 
            {
               actpart->pdis[0].element[counter] = &(actfield->dis[0].element[j]);
               counter++;
            }
         }
      /* loop and counter nodes, every node belonging to one of the partitions
         elements belongs to this partition */
         counter=0;
         for (j=0; j<actfield->dis[0].numnp; j++)
         {
            actnode = &(actfield->dis[0].node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->pdis[0].numnp = counter;
         actpart->pdis[0].node = (NODE**)CALLOC(counter,sizeof(NODE*));
         if (!actpart->pdis[0].node) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->dis[0].numnp; j++)
         {
            actnode = &(actfield->dis[0].node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  actpart->pdis[0].node[counter] = actnode;
                  counter++;
                  break;
               }
            }
         }
      /*----------------- there is no such thing as boundary and inner
                                                               elements */
         actpart->pdis[0].bou_numele=0;
         actpart->pdis[0].inner_numele=0;
         actpart->pdis[0].bou_element=NULL;
         actpart->pdis[0].inner_element=NULL;
      /*---------------------------- count the inner and boundary nodes */      
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[0].numnp; j++)
         {
            isbou=0;
            actnode = actpart->pdis[0].node[j];
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
         actpart->pdis[0].inner_numnp=counter;
         actpart->pdis[0].bou_numnp  =counter2;
         actpart->pdis[0].inner_node = (NODE**)CALLOC(counter,sizeof(NODE*));
         actpart->pdis[0].bou_node   = (NODE**)CALLOC(counter2,sizeof(NODE*));
         if (actpart->pdis[0].inner_node==NULL || actpart->pdis[0].bou_node==NULL)
         dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->pdis[0].numnp; j++)
         {
            isbou=0;
            actnode = actpart->pdis[0].node[j];
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
               actpart->pdis[0].bou_node[counter2] = actnode;
               counter2++;
            }
            else          
            {
               actpart->pdis[0].inner_node[counter] = actnode;
               counter++;
            }
         }
         
      } /* end of part==2 */



   } /* end of inprocs > 1 */
#ifndef PARALLEL 
FREE(actintra);
#endif
}/* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of part_assignfield */

/*! @} (documentation module close)*/
