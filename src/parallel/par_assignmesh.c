#include "../headers/standardtypes.h"
#include "../headers/solution.h"

/*----------------------------------------------------------------------*
 |  do initial partitioning of fields                    m.gee 5/01     |
 *----------------------------------------------------------------------*/
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
partition = (PARTITION*)calloc(genprob.numfld,sizeof(PARTITION));
if (partition==NULL) dserror("Allocation of PARTITION failed");
/*----------------------------------------------------- loop all fields */
for (i=0; i<genprob.numfld; i++)
{
   actsolv  = &(solv[i]);
   actfield = &(field[i]);
   actpart  = &(partition[i]);
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
#else
   actintra    = (INTRA*)calloc(1,sizeof(INTRA));
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
      actpart->numnp       = actfield->numnp;
      actpart->numele      = actfield->numele;
      actpart->bou_numnp   = 0;
      actpart->bou_numele  = 0;
      actpart->bou_element = NULL;
      actpart->bou_node    = NULL;
      actpart->element = (ELEMENT**)calloc(actpart->numele,sizeof(ELEMENT*));
      actpart->node    = (NODE**)calloc(actpart->numnp,sizeof(NODE*));
      if (actpart->element==NULL) dserror("Allocation of element pointer in PARTITION failed");
      if (actpart->node==NULL)    dserror("Allocation of node pointer in PARTITION failed");
      for (j=0; j<actfield->numele; j++) actpart->element[j] = &(actfield->element[j]);
      for (j=0; j<actfield->numnp; j++)  actpart->node[j] = &(actfield->node[j]);
      
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
         for (j=0; j<actfield->numele; j++)
         {
            for (k=0; k<actfield->element[j].numnp; k++)
            {
               if (actfield->element[j].node[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->numele  = counter;
         actpart->element = (ELEMENT**)calloc(counter,sizeof(ELEMENT*));
         if (actpart->element == NULL) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->numele; j++)
         {
            for (k=0; k<actfield->element[j].numnp; k++)
            {
               if (actfield->element[j].node[k]->proc == imyrank)
               {
                  actpart->element[counter] = &(actfield->element[j]);
                  counter++;
                  break;
               }
            }
         }
      /*------------------ loop and count nodes, do pointers to my nodes */   
         counter=0;
         for (j=0; j<actfield->numnp; j++)
         {
            if (actfield->node[j].proc == imyrank) counter++;
         }
         actpart->numnp = counter;
         actpart->node = (NODE**)calloc(counter,sizeof(NODE*));
         if (actpart->node==NULL) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->numnp; j++)
         {
            if (actfield->node[j].proc == imyrank) 
            {
               actpart->node[counter] = &(actfield->node[j]);
               counter++;
            }
         }
       /*----paritioning this way does not need inner and boundary nodes */
         actpart->inner_numnp=0;
         actpart->bou_numnp=0;
         actpart->inner_node=NULL;
         actpart->bou_node=NULL;
       /*----------------------------- now count the pure inner & boundary
                                                                elements */
       /*------------ Juhu !!! This is the first real parallel loop !!!! */
         counter=0;
         counter2=0;
         for (j=0; j<actpart->numele; j++)
         {
            isbou=0;
            proc = actpart->element[j]->node[0]->proc;
            for (k=1; k<actpart->element[j]->numnp; k++)
            {
               if (actpart->element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) counter++;
            else          counter2++;
         }
         actpart->inner_numele = counter;
         actpart->bou_numele   = counter2;
         actpart->inner_element = (ELEMENT**)calloc(counter,sizeof(ELEMENT*));
         actpart->bou_element = (ELEMENT**)calloc(counter2,sizeof(ELEMENT*));
         if (actpart->inner_element==NULL || actpart->bou_element==NULL)
         dserror("Allocation of PARTITION to ELEMENT pointer failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->numele; j++)
         {
            isbou=0;
            proc = actpart->element[j]->node[0]->proc;
            for (k=1; k<actpart->element[j]->numnp; k++)
            {
               if (actpart->element[j]->node[k]->proc != proc) 
               {
                  isbou=1;
                  break;
               }
            }
            if (isbou==0) 
            {
               actpart->inner_element[counter] = actpart->element[j];
               counter++;
            }
            else
            {
               actpart->bou_element[counter2]  = actpart->element[j];
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
         for (j=0; j<actfield->numele; j++)
         {
            if (actfield->element[j].proc == imyrank) counter++;
         }
         actpart->numele = counter;
         actpart->element = (ELEMENT**)calloc(counter,sizeof(ELEMENT*));
         if (actpart->element==NULL) dserror("Allocation of ELEMENT ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->numele; j++)
         {
            if (actfield->element[j].proc == imyrank) 
            {
               actpart->element[counter] = &(actfield->element[j]);
               counter++;
            }
         }
      /* loop and counter nodes, every node belonging to one of the partitions
         elements belongs to this partition */
         counter=0;
         for (j=0; j<actfield->numnp; j++)
         {
            actnode = &(actfield->node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  counter++;
                  break;
               }
            }
         }
         actpart->numnp = counter;
         actpart->node = (NODE**)calloc(counter,sizeof(NODE*));
         if (actpart->node==NULL) dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         for (j=0; j<actfield->numnp; j++)
         {
            actnode = &(actfield->node[j]);
            for (k=0; k<actnode->numele; k++)
            {
               if (actnode->element[k]->proc == imyrank)
               {
                  actpart->node[counter] = actnode;
                  counter++;
                  break;
               }
            }
         }
      /*----------------- there is no such thing as boundary and inner
                                                               elements */
         actpart->bou_numele=0;
         actpart->inner_numele=0;
         actpart->bou_element=NULL;
         actpart->inner_element=NULL;
      /*---------------------------- count the inner and boundary nodes */      
         counter=0;
         counter2=0;
         for (j=0; j<actpart->numnp; j++)
         {
            isbou=0;
            actnode = actpart->node[j];
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
         actpart->inner_numnp=counter;
         actpart->bou_numnp  =counter2;
         actpart->inner_node = (NODE**)calloc(counter,sizeof(NODE*));
         actpart->bou_node   = (NODE**)calloc(counter2,sizeof(NODE*));
         if (actpart->inner_node==NULL || actpart->bou_node==NULL)
         dserror("Allocation of NODE ptr in PARTITION failed");
         counter=0;
         counter2=0;
         for (j=0; j<actpart->numnp; j++)
         {
            isbou=0;
            actnode = actpart->node[j];
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
               actpart->bou_node[counter2] = actnode;
               counter2++;
            }
            else          
            {
               actpart->inner_node[counter] = actnode;
               counter++;
            }
         }
         
      } /* end of part==2 */



   } /* end of inprocs > 1 */
#ifndef PARALLEL 
free(actintra);
#endif
}/* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of part_assignfield */
