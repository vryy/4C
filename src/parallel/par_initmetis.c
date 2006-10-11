/*!---------------------------------------------------------------------
\file
\brief domain decomposition and metis routines

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"

#ifdef D_ALE
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
#endif


#ifdef PARALLEL

#ifdef HPUX10
#include "metis/metis.h"
#endif

#ifdef HPUX11
#include "metis/metis.h"
#endif

#ifdef HPUXITA
#include "metis/metis.h"
#endif

#ifdef AZUSA
#include "../../../../lib_ita1/metis/metis.h"
#endif

#ifdef SUN
#include "../../../lib_sun/metis-4.0/Lib/metis.h"
#endif

#ifdef LINUX_MUENCH
#include "metis.h"
#endif

#endif



static INT MAXNODPERGLINE=3;
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

/*!
\addtogroup PARALLEL
*/

/*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

/*!---------------------------------------------------------------------
\brief initial partitioning of fields

<pre>                                                        m.gee 5/01
the partitioning of all fields is performed on all procs,
so at least everyone nows which piece of every field is owned by who
this routine lives in the MPI_COMM_WORLD space
it uses the sequentiell metis lib to do graph partitioning
every field is partitioned separately, except for ale, that inherits it's
partition from the corresponding fluid elements. This feature is switched
off by malte neumann at the moment (Aug/2002)
</pre>

\return void
\sa part_assignfield()

------------------------------------------------------------------------*/
void part_fields()
{

  INT      i,j,k,l,disnum;
  INT      counter;
  INT      adjcounter;
  long int max,min;
  INT      proc, proc2;
  INT      sameproc;
  INT      ngline;

  INT      imyrank;
  INT      inprocs;
  NODE    *actnode;
  FIELD   *actfield = NULL;
  DISCRET *actdis   = NULL;
  GNODE   *actgnode;
  GLINE   *actgline;
  ELEMENT *actele;

#ifdef D_FLUID2_PRO
  ELEMENT *actvele, *actpele;
#endif

  ELEMENT *actele2;

  ARRAY    stack;

  ARRAY    xadj[MAXFIELD];
  ARRAY    adjncy[MAXFIELD];
  ARRAY    vwgt[MAXFIELD];

  INT      ione=1;
  INT      options[5];
  INT      nparts;
  ARRAY    part;
  ARRAY    part_proc;
  ARRAY    ele_per_proc;
  ARRAY    gl_per_proc_a;
  INT     *gl_per_proc;
  ARRAY    lineproc_a;
  INT     *lineproc;

  INT      upbound;

#ifdef PARALLEL
  INTRA   *actintra;
  INT      numflag=0;
  INT      edgecut;
  INT      wgtflag=2;
#endif


#ifdef DEBUG
  dstrc_enter("part_fields");
#endif


  /* sequentiell version:
   * --------------------
   * assign all elements and nodes to proc 0;
   * the partitioning is not performed, all nodes and elements are assigned to
   * proc 0, but the graphs of the fields are build in sequentiell and
   * parallel version
   */
  if (par.nprocs<=1)
  {
    for (i=0; i<genprob.numfld; i++)
    {
      actfield = &(field[i]);
      for(disnum=0;disnum<actfield->ndis;disnum++)
      {
        for (j=0; j<actfield->dis[disnum].numele; j++)
          actfield->dis[disnum].element[j].proc = 0;
        for (j=0; j<actfield->dis[disnum].numnp; j++)
          actfield->dis[disnum].node[j].proc    = 0;
      }
    }
  }
  imyrank=0;
  inprocs=1;




  for (i=0; i<genprob.numfld; i++)
  {
    actfield = &(field[i]);

#ifdef PARALLEL
    actintra = &(par.intra[i]);
    /*------ check proc belonging to this intra-communicator group of procs */
    if (actintra->intra_fieldtyp==none) continue;
    imyrank  = actintra->intra_rank;
    inprocs  = actintra->intra_nprocs;
#endif

      for(disnum=0;disnum<actfield->ndis;disnum++)
      {

        actdis = &(actfield->dis[disnum]);


        /* init the local numbering of the elements */
        /* numbering is c style, starts with zero */
        counter=0;
        for (j=0; j<actdis->numele; j++)
        {
          actdis->element[j].Id_loc = counter;
          counter++;
        }


        /* init the numbering of the nodes */
        counter=0;
        for (j=0; j<actdis->numnp; j++)
        {
          actdis->node[j].Id_loc = counter;
          counter++;
        }




        /* calculate the graph */
        /*=====================*/


        /* size of ARRAY xadj is numnp+1 */
        amdef("xadj",&(xadj[i]),(actdis->numnp+1),1,"IV");
        amzero(&(xadj[i]));

        /* the vertex weights of the graph */
        amdef("vwgt",&(vwgt[i]),(actdis->numnp)  ,1,"IV");
        aminit(&(vwgt[i]),&ione);

        /* size of array adjncy has to be computed */
        amdef("stack",&stack,1,1,"IV");
        amzero(&stack);



        /* Upper Bound for adjncy:
         * numele * ele.numnp * (ele.numnp-1) */
        upbound = actdis->numele *
          actdis->element[0].numnp *
          (actdis->element[0].numnp - 1);

        if (upbound < 10)
          upbound = 10;

#if 0
        printf("Upper Bound for adjncy: %d\n",upbound);
        fflush(stdout);
#endif


        amdef("adjncy",&(adjncy[i]),0.6*upbound,1,"IV");
        amzero(&(adjncy[i]));
        adjcounter=0;


        /* loop the nodes of this discretization */
        for (j=0; j<actdis->numnp; j++)
        {

          counter=0;
          actnode = &(actdis->node[j]);

          /* determine size of stack */
          /* use all nodes of the neighbour elements */
          for (k=0; k<actnode->numele; k++)
          {
            counter += actnode->element[k]->numnp;
          }


          /* allocate stack */
          amdel(&stack);
          amdef("stack",&stack,counter,1,"IV");
          amzero(&stack);
          counter=0;


          /* put nodes on stack */
          for (k=0; k<actnode->numele; k++)
          {
            actele = actnode->element[k];
            for (l=0; l<actele->numnp; l++)
            {
              stack.a.iv[counter] = actele->node[l]->Id_loc;
              counter++;
            }
          }


          /* delete doubles on stack */
          for (k=0; k<stack.fdim; k++)
          {
            counter=stack.a.iv[k];
            if (counter==-1) continue;
            for (l=k+1; l<stack.fdim; l++)
            {
              if (stack.a.iv[l]==counter) stack.a.iv[l]=-1;
            }
          }


          /* count number of nodes on stack */
          counter=0;
          for (k=0; k<stack.fdim; k++)
          {
            if (stack.a.iv[k] != -1) counter++;
          }


          /* number of vertices is from actnode to each other node */
          counter--;



          /* put edges on stack to xadj and adjncy */
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
                amredef(&(adjncy[i]),(adjncy[i].fdim+0.1*upbound),1,"IV");
#if 0
                printf("REDEFINE: New size: %d\n",adjncy[i].fdim);
                fflush(stdout);
#endif
                adjncy[i].a.iv[adjcounter] = stack.a.iv[k];
                adjcounter++;
              }
            }  /* if (stack.a.iv[k] != -1 && stack.a.iv[k]!= actnode->Id_loc) */

          }  /* for (k=0; k<stack.fdim; k++) */

        }  /* for (j=0; j<actdis->numnp; j++) */

#if 0
        printf("Check size of adjncy!! Current size: %d  --- Requiered size: %d\n",
            adjncy[i].fdim,adjcounter);
        fflush(stdout);
#endif

        amredef(&(adjncy[i]),(adjcounter),1,"IV");
        amdel(&stack);



        /* do NOT partition sequentiell version */
        if (inprocs<=1)
        {
          amdel(&(adjncy[i]));
          amdel(&(xadj[i]));
          amdel(&(vwgt[i]));
          continue;
        }



        /* do NOT partition some classes of discretizations */
        if (actdis->disclass == dc_created_ale ||
            actdis->disclass == dc_subdiv_io_created_ale ||
            actdis->disclass == dc_created_tu  ||
            actdis->disclass == dc_created_f2p)
        {
          goto after_part;
        }


        /* now allocate the rest of the arrays needed by metis */
        amdef("part",&part,actdis->numnp,1,"IV");


        /* set the default options of metis (play with other options later) */
        options[0]=0;
        options[1]=3;
        options[2]=1;
        options[3]=1;
        options[4]=0;


        /* call metis */
        nparts=inprocs;
        /* partitioning is not necessarily deterministic, so it is only performed
           by proc 0 and then broadcasted                                       */
        if (imyrank==0)
        {
          if (nparts < 8) /*----- better for a smaller number of partitions */
          {
#ifdef PARALLEL
            METIS_PartGraphRecursive(
                &(actdis->numnp),
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
                &(actdis->numnp),
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

          }  /* if (nparts < 8) */
        }  /* if (imyrank==0) */




#ifdef PARALLEL
        /* broadcast partitioning results */
        MPI_Bcast(&(part.a.iv[0]), part.fdim, MPI_INT, 0, actintra->MPI_INTRA_COMM);
#endif


        /* put the proc number to the nodes */
        for (j=0; j<actdis->numnp; j++)
        {
          actdis->node[j].proc=part.a.iv[j];
        }
        amdel(&part);



        /* check the nodes of all elements and assign element to proc */
        amdef("part",&part,nparts,1,"IV");
        amzero(&part);
        amdef("part_proc",&part_proc,nparts,1,"IV");
        amzero(&part_proc);
        amdef("ele_proc",&ele_per_proc,nparts,1,"IV");
        amzero(&ele_per_proc);

        for (j=0; j<actdis->numele; j++)
        {
          actele = &(actdis->element[j]);

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

        }  /* for (j=0; j<actdis->numele; j++) */

        amdel(&part);
        amdel(&part_proc);
        amdel(&ele_per_proc);
        amdel(&(xadj[i]));
        amdel(&(adjncy[i]));
        amdel(&(vwgt[i]));



after_part:

        /* Special cases:
         * ==============
         * Now assign the nodes to the processors for the special cases.
         */



#ifdef D_FLUID2_PRO
        /*
           For the FLUID2_PRO element we have two discretisations
           (velocity=0 and pressure=1):
           discretisation 0 was partitioned by METIS and this is no copied to
           discretisation 1:
           */
        if (actdis->disclass==dc_created_f2p)
        {
          for(j=0; j<actfield->dis[0].numele; j++)
          {
            actvele = &(actfield->dis[0].element[j]);
            actpele = &(actfield->dis[1].element[j]);
            dsassert(actvele->eltyp==el_fluid2_pro,
                "actfield=fluid & ndis>1 but eltyp!=el_fluid2_pro\n");
            dsassert(actpele->eltyp==el_fluid2_pro,
                "actfield=fluid & ndis>1 but eltyp!=el_fluid2_pro\n");
            actpele->proc = actvele->proc;
            for (k=0; k<actpele->numnp; k++)
              actpele->node[k]->proc = actvele->node[k]->proc;
          }
        }
#endif


        /*
           for solving the turbulence-models we need a 2nd discretisation.
           discretisation 0 was partitioned by METIS and this is no copied to
           discretisation 1:
           */
        if (actdis->disclass==dc_created_tu)
        {
          /* loop fluid elements for TURBULENCE*/
          for (j=0; j<actfield->dis[1].numele; j++)
          {
            actele  = &(actfield->dis[1].element[j]);
            actele2 = &(actfield->dis[0].element[j]);
            actele->proc = actele2->proc;

            /* loop nodes of element */
            for (k=0; k<actele->numnp; k++)
              actele->node[k]->proc = actele2->node[k]->proc;
          }
        }


#ifdef D_ALE
        /* assign procs to created ale elements:
         * -------------------------------------
         */
        if (actdis->disclass == dc_created_ale ||
            actdis->disclass == dc_subdiv_io_created_ale)
        {
          /* loop ale elements */
          for (j=0; j<actdis->numele; j++)
          {
            actele = &(actfield->dis[disnum].element[j]);

            switch (actele->eltyp)
            {
              case el_ale3:
                actele->proc = actele->e.ale3->fluid_ele->proc;

                /* loop nodes of element */
                for (k=0; k<actele->numnp; k++)
                  actele->node[k]->proc = actele->e.ale3->fluid_ele->node[k]->proc;

                break;


              case el_ale2:
                actele->proc = actele->e.ale2->fluid_ele->proc;

                /* loop nodes of element */
                for (k=0; k<actele->numnp; k++)
                  actele->node[k]->proc = actele->e.ale2->fluid_ele->node[k]->proc;

                break;

              default:
                dserror("Unknown element type in ale discretization!!");
            }
          }
        }
#endif


#if 0
#ifdef SUBDIV
        /* assign procs to subdivided elements:
         * ------------------------------------
         */
        if (actdis->disclass == dc_subdiv_calc)
        {
          /* loop elements */
          for (j=0; j<actdis->numele; j++)
          {
            actele = &(actfield->dis[disnum].element[j]);

            actele->proc = actele->master_ele->proc;

            /* loop nodes of element */
            for (k=0; k<actele->numnp; k++)
            {
              if (actele->node[k]->master_node != NULL)
                actele->node[k]->proc = actele->node[k]->master_node->proc;
              else
                actele->node[k]->proc = actele->master_ele->proc;
            }
          }
        }
#endif
#endif





        /* assign gline to proc */
        gl_per_proc = amdef("gl_per_proc",&gl_per_proc_a,nparts,1,"IV");
        amzero(&gl_per_proc_a);
        lineproc = amdef("lineproc",&lineproc_a,MAXNODPERGLINE,1,"IV");
        ngline=actdis->ngline;
        for (j=0;j<ngline;j++)
        {
          actgline = &(actdis->gline[j]);
          actgline->proc = -1;
        }


        /* loop all glines and assign proc to these glines where all nodes
           belong to the same domain                                          */
        for (j=0;j<ngline;j++)
        {
          actgline = &(actdis->gline[j]);
          dsassert(actgline->ngnode<=3,
              "number of nodes at a gline > 3 not possible!\n");

          counter=0;
          for (k=0;k<actgline->ngnode;k++)
          {
            actgnode = actgline->gnode[k];
            lineproc[k]=actgnode->node->proc;
          }

          /* check if gline belongs only to one proc */
          proc = lineproc[0];
          sameproc = 0;
          for (k=1;k<actgline->ngnode;k++)
            if (proc != lineproc[k]) sameproc=1;
          if (sameproc==0)
          {
            actgline->proc=proc;
            gl_per_proc[proc]+=1;
          }
        }


        /* now loop again all glines and assign the procs to the rest */
        for (j=0;j<ngline;j++)
        {
          actgline = &(actdis->gline[j]);
          if (actgline->proc>-1) continue;

          actgnode = actgline->gnode[0];
          proc = actgnode->node->proc;

          for (k=1;k<actgline->ngnode;k++)
          {
            actgnode=actgline->gnode[k];
            proc2 = actgnode->node->proc;
            if (gl_per_proc[proc2]<gl_per_proc[proc]) proc = proc2;
          }

          actgline->proc=proc;
          gl_per_proc[proc]+=1;
        }

        amdel(&gl_per_proc_a);
        amdel(&lineproc_a);


      }  /* for(disnum=0;disnum<actfield->ndis;disnum++) */

  }  /* for (i=0; i<genprob.numfld; i++) */





#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of part_fields */



/*! @} (documentation module close)*/






