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

  INT        i,j,k,disnum;
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


  /* every proc alloc numfld PARTITIONs */
  partition = (PARTITION*)CCACALLOC(genprob.numfld,sizeof(PARTITION));
  if (partition==NULL) dserror("Allocation of PARTITION failed");


  /* loop all fields */
  for (i=0; i<genprob.numfld; i++)
  {
    actsolv  = &(solv[i]);
    actfield = &(field[i]);
    actfield->part = &(partition[i]);

    /* every partition[i] allocates discretizations */
    partition[i].ndis = actfield->ndis;
    partition[i].pdis = (PARTDISCRET*)CCACALLOC(partition[i].ndis,sizeof(PARTDISCRET));
    if (!partition[i].pdis) dserror("Allocation of memory failed");


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


    /* if there is only one proc, there is nothing to do */
    if (inprocs<=1)
    {
      for (disnum=0;disnum<actpart->ndis;disnum++)
      {
        PARTDISCRET     *actpdis;
        DISCRET         *actdis;

        actpdis = &(actpart->pdis[disnum]);
        actdis  = &(actfield->dis[disnum]);

        actpdis->numnp       = actdis->numnp;
        actpdis->numele      = actdis->numele;
        actpdis->numlele     = actdis->numele;
        actpdis->bou_numnp   = 0;
        actpdis->bou_numele  = 0;
        actpdis->bou_element = NULL;
        actpdis->bou_node    = NULL;
        actpdis->element = (ELEMENT**)CCACALLOC(actpdis->numele,sizeof(ELEMENT*));
        actpdis->node    = (NODE**)CCACALLOC(actpdis->numnp,sizeof(NODE*));
        if (actpdis->element==NULL) dserror("Allocation of element pointer in PARTITION failed");
        if (actpdis->node==NULL)    dserror("Allocation of node pointer in PARTITION failed");
        for (j=0; j<actdis->numele; j++)
          actpdis->element[j] = &(actdis->element[j]);
        for (j=0; j<actdis->numnp; j++)
          actpdis->node[j] = &(actdis->node[j]);
      }
    }
    else
    {

      /* do only for procs within the intra-communicator */
      if (actintra->intra_fieldtyp==none) continue;

      /* check for typ of partition */
      part=0;
      if (actsolv->parttyp == cut_elements) part = 1;
      if (actsolv->parttyp == cut_nodes)    part = 2;
      if (part==0) dserror("Typ of Partitioning unknown");


      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
       * parition field, such that cuts are through elements and each node
       * is assigned exactly one proc.
       *
       * Attention: From now on information on different procs can differ!!
       *            severely !! (Additive Schwartz)
       *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


      for (disnum=0;disnum<actpart->ndis;disnum++) /* loop over discretisations */
      {
        PARTDISCRET     *actpdis;
        DISCRET         *actdis;

        actpdis = &(actpart->pdis[disnum]);
        actdis  = &(actfield->dis[disnum]);

        if (part==1)
        {
          /* count local elements */
          /* these are needed for discontinuous pressure calculations */
          counter=0;
          for (j=0; j<actdis->numele; j++)
          {
            if (actdis->element[j].proc == imyrank)
            {
              counter++;
            }
          }
          actpdis->numlele  = counter;

          /* loop and count elements, do pointers to my elements */
          counter=0;
          for (j=0; j<actdis->numele; j++)
          {
            for (k=0; k<actdis->element[j].numnp; k++)
            {
              if (actdis->element[j].node[k]->proc == imyrank)
              {
                counter++;
                break;
              }
            }
          }
          actpdis->numele  = counter;
          actpdis->element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
          if (!actpdis->element)
            dserror("Allocation of ELEMENT ptr in PARTITION failed");

          counter=0;
          for (j=0; j<actdis->numele; j++)
          {
            for (k=0; k<actdis->element[j].numnp; k++)
            {
              if (actdis->element[j].node[k]->proc == imyrank)
              {
                actpdis->element[counter] = &(actdis->element[j]);
                counter++;
                break;
              }
            }
          }


          /* loop and count nodes, do pointers to my nodes */
          counter=0;
          for (j=0; j<actdis->numnp; j++)
          {
            if (actdis->node[j].proc == imyrank) counter++;
          }
          actpdis->numnp = counter;
          actpdis->node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
          if (!actpdis->node) dserror("Allocation of NODE ptr in PARTITION failed");
          counter=0;
          for (j=0; j<actdis->numnp; j++)
          {
            if (actdis->node[j].proc == imyrank)
            {
              actpdis->node[counter] = &(actdis->node[j]);
              counter++;
            }
          }

          /* partitioning this way does not need inner and boundary nodes */
          actpdis->inner_numnp=0;
          actpdis->bou_numnp=0;
          actpdis->inner_node=NULL;
          actpdis->bou_node=NULL;


          /* now count the pure inner & boundary elements */
          counter=0;
          counter2=0;
          for (j=0; j<actpdis->numele; j++)
          {
            isbou=0;
            proc = actpdis->element[j]->node[0]->proc;
            for (k=1; k<actpdis->element[j]->numnp; k++)
            {
              if (actpdis->element[j]->node[k]->proc != proc)
              {
                isbou=1;
                break;
              }
            }
            if (isbou==0) counter++;
            else          counter2++;
          }
          actpdis->inner_numele = counter;
          actpdis->bou_numele   = counter2;
          actpdis->inner_element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
          actpdis->bou_element = (ELEMENT**)CCACALLOC(counter2,sizeof(ELEMENT*));

          if (actpdis->inner_element==NULL || actpdis->bou_element==NULL)
            dserror("Allocation of PARTITION to ELEMENT pointer failed");

          counter=0;
          counter2=0;
          for (j=0; j<actpdis->numele; j++)
          {
            isbou=0;
            proc = actpdis->element[j]->node[0]->proc;
            for (k=1; k<actpdis->element[j]->numnp; k++)
            {
              if (actpdis->element[j]->node[k]->proc != proc)
              {
                isbou=1;
                break;
              }
            }
            if (isbou==0)
            {
              actpdis->inner_element[counter] = actpdis->element[j];
              counter++;
            }
            else
            {
              actpdis->bou_element[counter2]  = actpdis->element[j];
              counter2++;
            }
          }
        }  /* if (part==1) */


        /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
         * partition field, such that cuts are through nodes and boundary nodes
         * are on the common procs *
         *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

        if (part==2)
        {

          /* loop and count elements, every element belongs exactly to one proc */
          counter=0;
          for (j=0; j<actdis->numele; j++)
          {
            if (actdis->element[j].proc == imyrank) counter++;
          }
          actpdis->numele = counter;
          actpdis->element = (ELEMENT**)CCACALLOC(counter,sizeof(ELEMENT*));
          if (!actpdis->element)
            dserror("Allocation of ELEMENT ptr in PARTITION failed");

          counter=0;
          for (j=0; j<actdis->numele; j++)
          {
            if (actdis->element[j].proc == imyrank)
            {
              actpdis->element[counter] = &(actdis->element[j]);
              counter++;
            }
          }


          /* loop and counter nodes, every node belonging to one of the partitions
             elements belongs to this partition */
          counter=0;
          for (j=0; j<actdis->numnp; j++)
          {
            actnode = &(actdis->node[j]);
            for (k=0; k<actnode->numele; k++)
            {
              if (actnode->element[k]->proc == imyrank)
              {
                counter++;
                break;
              }
            }
          }

          actpdis->numnp = counter;
          actpdis->node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
          if (!actpdis->node) dserror("Allocation of NODE ptr in PARTITION failed");
          counter=0;
          for (j=0; j<actdis->numnp; j++)
          {
            actnode = &(actdis->node[j]);
            for (k=0; k<actnode->numele; k++)
            {
              if (actnode->element[k]->proc == imyrank)
              {
                actpdis->node[counter] = actnode;
                counter++;
                break;
              }
            }
          }


          /* there is no such thing as boundary and inner elements */
          actpdis->bou_numele=0;
          actpdis->inner_numele=0;
          actpdis->bou_element=NULL;
          actpdis->inner_element=NULL;


          /* count the inner and boundary nodes */
          counter=0;
          counter2=0;
          for (j=0; j<actpdis->numnp; j++)
          {
            isbou=0;
            actnode = actpdis->node[j];
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


          actpdis->inner_numnp=counter;
          actpdis->bou_numnp  =counter2;
          actpdis->inner_node = (NODE**)CCACALLOC(counter,sizeof(NODE*));
          actpdis->bou_node   = (NODE**)CCACALLOC(counter2,sizeof(NODE*));
          if (actpdis->inner_node==NULL || actpdis->bou_node==NULL)
            dserror("Allocation of NODE ptr in PARTITION failed");
          counter=0;
          counter2=0;


          for (j=0; j<actpdis->numnp; j++)
          {
            isbou=0;
            actnode = actpdis->node[j];
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
              actpdis->bou_node[counter2] = actnode;
              counter2++;
            }
            else
            {
              actpdis->inner_node[counter] = actnode;
              counter++;
            }
          }

        }  /* if (part==2) */

      }  /* for (disnum=0;disnum<actpart->ndis;disnum++) */


    }  /* end of inprocs > 1 */


#ifndef PARALLEL
    CCAFREE(actintra);
#endif


  }  /* end of loop over fields */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of part_assignfield */

/*! @} (documentation module close)*/
