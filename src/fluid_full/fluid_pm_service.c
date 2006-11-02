/*!
\file
\brief service functions for projection algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*!
\addtogroup FLUID_PM
*//*! @{ (documentation module open)*/

#ifdef D_FLUID_PM

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../solver/solver_trilinos_service.H"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"

#ifdef D_FLUID2_PRO
#include "../fluid2_pro/fluid2pro.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#endif

#ifdef D_FLUID3_PRO
#include "../fluid3_pro/fluid3pro.h"
#include "../fluid3_pro/fluid3pro_prototypes.h"
#endif


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;



/*----------------------------------------------------------------------*/
/*!
  \brief build the sparse mask of the C^T*ML^-1*C matrix

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param numpdof            (i) number of pressure dofs
  \param pmat               (o) pressure sparse mask

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_build_pmat_sparse_mask(FIELD* actfield,
                               PARTITION *actpart,
                               INT disnum,
                               INTRA* actintra,
                               INT numpdof,
                               PARALLEL_SPARSE* pmat)
{
  INT i;
  SPARSE_MASK_LIST sml;

#ifdef DEBUG
  dstrc_enter("pm_build_pmat_sparse_mask");
#endif

  sparse_mask_list_init(&sml,
                        numpdof*actpart->pdis[disnum].numlele,
                        numpdof*actfield->dis[disnum].numele,
                        100);

  /*
   * The sparse pattern is enlarged here. Each velocity dof is coupled
   * with the velocities of all neighbouring elements. And that is why
   * the pressures of neighbouring elements are coupled, too.
   *
   * There are at least two ways to find the mask. One can go to each
   * element, loop the neighbouring ones and mask the matrix. (Needs
   * different behaviour depending on the space dimension.) That is
   * what we do here. This way each position that gets an entry is
   * marked more that once. (The mask list can handle that.)
   *
   * An alternative approach is to start from the know mask of the
   * gradient matrix and test all entries in the new mask. This
   * results in a huge loop (all entries including the zero ones) and
   * requires communication. But it seems conceptually simpler. */

  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT j;
    ELEMENT* actele = actpart->pdis[disnum].element[i];
    if (actele->proc == actintra->intra_rank)
    {
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:

        /*
         * The strategy is to visit all nodes of the original element
         * and go to all elements that are connected to these. This loop
         * can be improved depending on the element type, but right now
         * I'm interested in the general case. */
        for (j=0; j<actele->numnp; ++j)
        {
          INT e;
          NODE* actnode;
          actnode = actele->node[j];
          for (e=0; e<actnode->numele; ++e)
          {
            INT k;
            ELEMENT* otherele;
            otherele = actnode->element[e];

            for (k=0; k<numpdof; ++k)
            {
              INT l;
              for (l=0; l<numpdof; ++l)
              {
                dsassert(otherele->eltyp==el_fluid2_pro,
                         "Projection method requires projection elements");

                /*
                 * Note the asymmetry: the local dof number from actele
                 * becomes the row index and the global dof of the
                 * connected element becomes the column index. That is
                 * because we slice pmat horizontally in the parallel
                 * case. actele belongs to the local pdis, but otherele
                 * might be a foreign element that is not calculated
                 * here. */
                sparse_mask_list_mark(&sml,
                                      actele->e.f2pro->ldof[k],
                                      otherele->e.f2pro->dof[l]);
              }
            }
          }
        }
        break;
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:

        /*
         * The strategy is to visit all nodes of the original element
         * and go to all elements that are connected to these. This loop
         * can be improved depending on the element type, but right now
         * I'm interested in the general case. */
        for (j=0; j<actele->numnp; ++j)
        {
          INT e;
          NODE* actnode;
          actnode = actele->node[j];
          for (e=0; e<actnode->numele; ++e)
          {
            INT k;
            ELEMENT* otherele;
            otherele = actnode->element[e];

            for (k=0; k<numpdof; ++k)
            {
              INT l;
              for (l=0; l<numpdof; ++l)
              {
                dsassert(otherele->eltyp==el_fluid3_pro,
                         "Projection method requires projection elements");

                /*
                 * Note the asymmetry: the local dof number from actele
                 * becomes the row index and the global dof of the
                 * connected element becomes the column index. That is
                 * because we slice pmat horizontally in the parallel
                 * case. actele belongs to the local pdis, but otherele
                 * might be a foreign element that is not calculated
                 * here. */
                sparse_mask_list_mark(&sml,
                                      actele->e.f3pro->ldof[k],
                                      otherele->e.f3pro->dof[l]);
              }
            }
          }
        }
        break;
#endif
      default:
        dserror("unsupported element type %d", actele->eltyp);
      }
    }
  }

  sparse_fix_mask(&(pmat->slice), &sml);
  sparse_mask_list_destroy(&sml);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief assign dof numbers to the discontinuous pressure dofs

  The dofs of the discontinuous pressure discretization belong to the
  elements. Thus we loop the elements here and assign dof numbers to
  the pressure dofs. These numbers start from zero. They cannot be
  confused with dofs that life in nodes.

  Dirichlet boundaries on the pressure are currently not supported.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
INT pm_assign_press_dof(FIELD *actfield,
                        PARTITION *actpart,
                        INT disnum,
                        INTRA* actintra)
{
  INT counter;
  INT i, numpdof;

#ifdef DEBUG
  dstrc_enter("assign_press_dof");
#endif

  counter = 0;
  numpdof = 0;

  for (i=0; i<actfield->dis[disnum].numele; ++i)
  {
    ELEMENT* actele;
    actele = &(actfield->dis[disnum].element[i]);

    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      FLUID2_PRO* f2pro;
      f2pro = actele->e.f2pro;
      switch (f2pro->dm)
      {
      case dm_q2pm1:
	if ((numpdof != 0) && (numpdof != 3))
	{
	  dserror("just one type of element allowed in one mesh");
	}
	numpdof = 3;

        /* the array have been allocated on reading the element
         * (inpele) */

#if 1
        /* assign numbers as we go */
        /* might result in a very bad sparse pattern */
        f2pro->dof[0] = counter++;
        f2pro->dof[1] = counter++;
        f2pro->dof[2] = counter++;
#else
        /* follow the element numbers */
        /* this relies on the mesh to contain nothing but f2pro
         * elements */
        f2pro->dof[0] = 3*actele->Id;
        f2pro->dof[1] = 3*actele->Id+1;
        f2pro->dof[2] = 3*actele->Id+2;
#endif
        f2pro->ldof[0] = -1;
        f2pro->ldof[1] = -1;
        f2pro->ldof[2] = -1;
        break;
      case dm_q1p0:
	if ((numpdof != 0) && (numpdof != 1))
	{
	  dserror("just one type of element allowed in one mesh");
	}
	numpdof = 1;

        /* the array have been allocated on reading the element
         * (inpele) */

#if 1
        /* assign numbers as we go */
        /* might result in a very bad sparse pattern */
        f2pro->dof[0] = counter++;
#else
        /* follow the element numbers */
        /* this relies on the mesh to contain nothing but f2pro
         * elements */
        f2pro->dof[0] = 3*actele->Id;
#endif
        f2pro->ldof[0] = -1;
        break;
      case dm_q2q1:
      default:
        dserror("discretization mode %d currently unsupported", f2pro->dm);
      }
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      FLUID3_PRO* f3pro;
      f3pro = actele->e.f3pro;
      switch (f3pro->dm)
      {
      case dm_q2pm1:
	if ((numpdof != 0) && (numpdof != 4))
	{
	  dserror("just one type of element allowed in one mesh");
	}
	numpdof = 4;

        /* the array have been allocated on reading the element
         * (inpele) */

#if 1
        /* assign numbers as we go */
        /* might result in a very bad sparse pattern */
        f3pro->dof[0] = counter++;
        f3pro->dof[1] = counter++;
        f3pro->dof[2] = counter++;
        f3pro->dof[3] = counter++;
#else
        /* follow the element numbers */
        /* this relies on the mesh to contain nothing but f3pro
         * elements */
        f3pro->dof[0] = 3*actele->Id;
        f3pro->dof[1] = 3*actele->Id+1;
        f3pro->dof[2] = 3*actele->Id+2;
#endif
        f3pro->ldof[0] = -1;
        f3pro->ldof[1] = -1;
        f3pro->ldof[2] = -1;
        f3pro->ldof[3] = -1;
        break;
      case dm_q1p0:
	if ((numpdof != 0) && (numpdof != 1))
	{
	  dserror("just one type of element allowed in one mesh");
	}
	numpdof = 1;

        /* the array have been allocated on reading the element
         * (inpele) */

#if 1
        /* assign numbers as we go */
        /* might result in a very bad sparse pattern */
        f3pro->dof[0] = counter++;
#else
        /* follow the element numbers */
        /* this relies on the mesh to contain nothing but f3pro
         * elements */
        f3pro->dof[0] = 3*actele->Id;
#endif
        f3pro->ldof[0] = -1;
        break;
      default:
        dserror("discretization mode %d currently unsupported", f3pro->dm);
      }
      break;
    }
#endif
    default:
      dserror("unsupported element type for projection method: %d", actele->eltyp);
    }
  }

  /* We need processor local dof numbers, too. */

  counter = 0;
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    ELEMENT* actele;
    actele = actpart->pdis[disnum].element[i];
    if (actele->proc == actintra->intra_rank)
    {
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
      {
        FLUID2_PRO* f2pro;
        f2pro = actele->e.f2pro;
        switch (f2pro->dm)
        {
        case dm_q2pm1:
          f2pro->ldof[0] = counter++;
          f2pro->ldof[1] = counter++;
          f2pro->ldof[2] = counter++;
          break;
        case dm_q1p0:
          f2pro->ldof[0] = counter++;
          break;
        case dm_q2q1:
        default:
          dserror("discretization mode %d currently unsupported", f2pro->dm);
        }
        break;
      }
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
      {
        FLUID3_PRO* f3pro;
        f3pro = actele->e.f3pro;
        switch (f3pro->dm)
        {
        case dm_q2pm1:
          f3pro->ldof[0] = counter++;
          f3pro->ldof[1] = counter++;
          f3pro->ldof[2] = counter++;
          f3pro->ldof[3] = counter++;
          break;
        case dm_q1p0:
          f3pro->ldof[0] = counter++;
          break;
        default:
          dserror("discretization mode %d currently unsupported", f3pro->dm);
        }
        break;
      }
#endif
      default:
        dserror("unsupported element type for projection method: %d", actele->eltyp);
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return(numpdof);
}


#if defined(PARALLEL) || defined(PM_TRILINOS)

/*----------------------------------------------------------------------*/
/*!
  \brief create the mapping from local to global dof numbers

  In parallel execution the sparse matrix's update array must be
  provided by the user. Here we do just that.

  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param numpdof            (i) number of pressure dofs
  \param numldof            (i) number of local pressure dofs
  \param update             (o) update array to be filled

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_fill_gradient_update(PARTITION *actpart,
                             INT disnum,
                             INTRA* actintra,
                             INT numpdof,
			     INT numldof,
                             INT* update)
{
  INT i;
#ifdef DEBUG
  dstrc_enter("pm_fill_gradient_update");
#endif

  /* fill the column id array */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT j;
    ELEMENT* actele = actpart->pdis[disnum].element[i];
    if (actele->proc == actintra->intra_rank)
    {
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        for (j=0; j<numpdof; ++j)
        {
          /*grad->update[j+numpdof*i] = actele->e.f2pro->dof[j];*/
          dsassert((actele->e.f2pro->ldof[j] >= 0) &&
                   (actele->e.f2pro->ldof[j] < numldof),
                   "local dof out of range");
          update[actele->e.f2pro->ldof[j]] = actele->e.f2pro->dof[j];
        }
        break;
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
        for (j=0; j<numpdof; ++j)
        {
          /*grad->update[j+numpdof*i] = actele->e.f3pro->dof[j];*/
          dsassert((actele->e.f3pro->ldof[j] >= 0) &&
                   (actele->e.f3pro->ldof[j] < numldof),
                   "local dof out of range");
          update[actele->e.f3pro->ldof[j]] = actele->e.f3pro->dof[j];
        }
        break;
#endif
      default:
        dserror("unsupported element type for projection method: %d", actele->eltyp);
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief create gradient G sparse mask

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param numpdof            (i) number of pressure dofs
  \param grad             (i/o) sparse pressure gradient

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_gradient_mask_mat(FIELD *actfield,
                          PARTITION *actpart,
                          INT disnum,
                          INTRA* actintra,
                          INT numpdof,
                          PARALLEL_SPARSE *grad)
{
  SPARSE_MASK_LIST sml;
  INT i;

#ifdef DEBUG
  dstrc_enter("pm_gradient_mask_mat");
#endif

  sparse_mask_list_init(&sml,
                        actfield->dis[disnum].numeq,
                        numpdof*actpart->pdis[disnum].numlele,
                        100);

  /* build sparse mask */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT j;
    ELEMENT* actele = actpart->pdis[disnum].element[i];
    if (actele->proc == actintra->intra_rank)
    {
      for (j=0; j<actele->numnp; ++j)
      {
        INT dof;
        NODE* actnode = actele->node[j];

        for (dof=0; dof<genprob.ndim; ++dof)
        {
          if (actnode->dof[dof] < actfield->dis[disnum].numeq)
          {
            INT k;
            for (k=0; k<numpdof; ++k)
            {
              /* Now we have the entry (actnode->dof[dof],
               * actele->ldof[k]) to be stored. */
              /*
               * The gradient matrix is vertically sliced, so the column
               * number matches the pressure's local dof number. */
              switch (actele->eltyp)
              {
#ifdef D_FLUID2_PRO
              case el_fluid2_pro:
                sparse_mask_list_mark(&sml,
                                      actnode->dof[dof],
                                      actele->e.f2pro->ldof[k]);
                break;
#endif
#ifdef D_FLUID3_PRO
              case el_fluid3_pro:
                sparse_mask_list_mark(&sml,
                                      actnode->dof[dof],
                                      actele->e.f3pro->ldof[k]);
                break;
#endif
              default:
                dserror("element type %d not suitable for projection method", actele->eltyp);
              }
            }
          }
        }
      }
    }
  }

  /* Mask done. Make it fix. */
  sparse_fix_mask(&(grad->slice), &sml);
  sparse_mask_list_destroy(&sml);

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create gradient and mass matrix

  All the element call and element matrix assembling is in here. We do
  not use global_calelm because it is such a mess and our requirements
  are specific. We have to build the inverted diagonalized mass
  matrix, too.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actsolv            (i) solver
  \param sysarray           (i) solver type
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param grad             (i/o) sparse pressure gradient
  \param lmass_vec          (o) inverted diagonal mass matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calelm(FIELD *actfield,
               PARTITION *actpart,
               INT disnum,
               SOLVAR *actsolv,
               INT sysarray,
               INTRA *actintra,
               ARRAY_POSITION *ipos,
               INT numpdof,
#ifdef PM_TRILINOS
	       TRILINOSMATRIX* grad,
	       TRILINOSMATRIX* lmass
#else
               PARALLEL_SPARSE* grad,
               DOUBLE* lmass_vec
#endif
  )
{
  INT i;
#ifdef PM_TRILINOS
#else
  DOUBLE* lmass;
#endif

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY emass_global;
  extern struct _ARRAY lmass_global;
  extern struct _ARRAY gradopr_global;

#ifdef DEBUG
  dstrc_enter("pm_calelm");
#endif

#ifdef PM_TRILINOS
#else
#ifdef PARALLEL
  lmass = (DOUBLE*)CCAMALLOC(actfield->dis[disnum].numeq*sizeof(DOUBLE));
#else
  lmass = lmass_vec;
#endif

  memset(lmass, 0, actfield->dis[disnum].numeq*sizeof(DOUBLE));
#endif

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT j;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    /* We need to assemble the global mass matrix. To do this in
     * parallel we cannot calculate the local elements only but have
     * to calculate any neighbouring elements. This in turn demands
     * that these elements know their pressure values. That is a
     * truely data parallel approach for discontinuous pressure
     * demands for element based communication (in contrast to the
     * node based communication we get by with in the continuous
     * case.) */

    /* Calculate gradient and mass matrix */
    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      f2pro_calgradp(actele, ipos);
      break;
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
      f3pro_calgradp(actele, ipos);
      break;
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }

    /* Assemble */

    /* At first lets do the global mass matrix */
    /*
     * We do this just once, thus there is no need to have this inside
     * the fluid element... */
    assemble(sysarray,
             &emass_global,
             -1,
             NULL,
             actpart,
             actsolv,
             actintra,
             actele,
             assemble_one_matrix,
             NULL);

    /* And now the lumped mass matrix (in vector form) and the
     * gradient matrix. */
    for (j=0; j<actele->numnp; ++j)
    {
      INT dof;
      NODE* actnode = actele->node[j];

      for (dof=0; dof<actnode->numdf; ++dof)
      {
	/* non-dirichlet dofs */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	if ((actnode->gnode->dirich==NULL) ||
	    (actnode->gnode->dirich->dirich_onoff.a.iv[dof]==0))
#else
        if (actnode->dof[dof] < actfield->dis[disnum].numeq)
#endif
        {
          INT k;

          /* But of course we assemble just to those nodes that belong to us. */
          if (actnode->proc == actintra->intra_rank)
          {
#ifdef PM_TRILINOS
	    add_trilinos_value(lmass,
			       lmass_global.a.dv[j*actnode->numdf+dof],
			       actnode->dof[dof],
			       actnode->dof[dof]);
#else
            lmass[actnode->dof[dof]] += lmass_global.a.dv[j*actnode->numdf+dof];
#endif
          }

          if (actele->proc == actintra->intra_rank)
          {
            for (k=0; k<numpdof; ++k)
            {
              /* Now we have the entry (actnode->dof[dof],
               * actele->ldof[k]) to be stored. */
              /*
               * Again: We access the columns via the local dof numbers. */
              switch (actele->eltyp)
              {
#ifdef D_FLUID2_PRO
              case el_fluid2_pro:
#ifdef PM_TRILINOS
		/*
		 * The orientation of the gradient matrix differs from
		 * the non-trilinos version! */
		add_trilinos_value(grad,
				   gradopr_global.a.da[actnode->numdf*j+dof][k],
				   /*
				    * I suppose trilinos expects the
				    * global dof number... */
				   /*actele->e.f2pro->ldof[k],*/
				   actele->e.f2pro->dof[k],
				   actnode->dof[dof]);
#else
                *sparse_entry(&(grad->slice),
                              actnode->dof[dof],
                              actele->e.f2pro->ldof[k]) +=
                  gradopr_global.a.da[actnode->numdf*j+dof][k];
#endif
                break;
#endif
#ifdef D_FLUID3_PRO
              case el_fluid3_pro:
#ifdef PM_TRILINOS
		/*
		 * The orientation of the gradient matrix differs from
		 * the non-trilinos version! */
		add_trilinos_value(grad,
				   gradopr_global.a.da[actnode->numdf*j+dof][k],
				   /*actele->e.f3pro->ldof[k],*/
				   actele->e.f3pro->dof[k],
				   actnode->dof[dof]);
#else
                *sparse_entry(&(grad->slice),
                              actnode->dof[dof],
                              actele->e.f3pro->ldof[k]) +=
                  gradopr_global.a.da[actnode->numdf*j+dof][k];
#endif
                break;
#endif
              default:
                dserror("element type %d unsupported", actele->eltyp);
              }
            }
          }
        }
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	/* dirichlet dofs */
	else
	{
          if (actnode->proc == actintra->intra_rank)
          {
#ifdef PM_TRILINOS
	    add_trilinos_value(lmass,
			       1,
			       actnode->dof[dof],
			       actnode->dof[dof]);
#else
            lmass[actnode->dof[dof]] = 1;
#endif
          }
	}
#endif
      }
    }
  }

#ifdef PM_TRILINOS

  close_nonquad_trilinos_matrix(grad,lmass);
  close_trilinos_matrix(lmass);
  invert_trilinos_diagonal_matrix(lmass);

#else /* PM_TRILINOS */

#ifdef PARALLEL
  /* We need the lumped masses globally. */
  MPI_Allreduce(lmass, lmass_vec, actfield->dis[disnum].numeq, MPI_DOUBLE,
                MPI_SUM, actintra->MPI_INTRA_COMM);
  CCAFREE(lmass);
  lmass = lmass_vec;
#endif

  /* Invert the lumed masses. It's such a pleasure. :)*/
  for (i=0; i<actfield->dis[disnum].numeq; ++i)
  {
    lmass[i] = 1./lmass[i];
  }

#endif /* PM_TRILINOS */

  /* Close the full mass matrix to enable the trilinos matrix-vector
   * product. */

  if (genprob.usetrilinosalgebra)
  {
    if (actsolv->sysarray_typ[sysarray]==trilinos)
    {
#ifdef TRILINOS_PACKAGE
      close_trilinos_matrix(actsolv->sysarray[sysarray].trilinos);
#else
      dserror("trilinos missing");
#endif
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#ifdef PM_TRILINOS

/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create gradient and mass matrix

  All the element call and element matrix assembling is in here. We do
  not use global_calelm because it is such a mess and our requirements
  are specific. We have to build the inverted diagonalized mass
  matrix, too.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actsolv            (i) solver
  \param sysarray           (i) solver type
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param grad             (i/o) sparse pressure gradient
  \param lmass_vec          (o) inverted diagonal mass matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calelm_cont(FIELD *actfield,
		    PARTITION *actpart,
		    INT vel_dis,
		    INT press_dis,
		    SOLVAR *actsolv,
		    INT sysarray,
		    INTRA *actintra,
		    ARRAY_POSITION *ipos,
		    TRILINOSMATRIX* grad,
		    TRILINOSMATRIX* lmass
  )
{
  INT i;
#ifdef PM_TRILINOS
#else
  DOUBLE* lmass;
#endif

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY emass_global;
  extern struct _ARRAY lmass_global;
  extern struct _ARRAY gradopr_global;

#ifdef DEBUG
  dstrc_enter("pm_calelm_cont");
#endif

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[vel_dis].numele; ++i)
  {
    INT j;
    ELEMENT* vele;
    ELEMENT* pele;
    vele = actpart->pdis[vel_dis].element[i];
    pele = actpart->pdis[press_dis].element[i];

    /* We need to assemble the global mass matrix. To do this in
     * parallel we cannot calculate the local elements only but have
     * to calculate any neighbouring elements. This in turn demands
     * that these elements know their pressure values. That is a
     * truely data parallel approach for discontinuous pressure
     * demands for element based communication (in contrast to the
     * node based communication we get by with in the continuous
     * case.) */

    /* Calculate gradient and mass matrix */
    switch (vele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      f2pro_calgradp(vele, ipos);
      break;
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
      f3pro_calgradp(vele, ipos);
      break;
#endif
    default:
      dserror("element type %d unsupported", vele->eltyp);
    }

    /* Assemble */

    /* At first lets do the global mass matrix */
    /*
     * We do this just once, thus there is no need to have this inside
     * the fluid element... */
    assemble(sysarray,
             &emass_global,
             -1,
             NULL,
             actpart,
             actsolv,
             actintra,
             vele,
             assemble_one_matrix,
             NULL);

    /* And now the lumped mass matrix (in vector form) and the
     * gradient matrix. */
    for (j=0; j<vele->numnp; ++j)
    {
      INT dof;
      NODE* vnode = vele->node[j];

      for (dof=0; dof<vnode->numdf; ++dof)
      {
	/* non-dirichlet dofs */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	if ((vnode->gnode->dirich==NULL) ||
	    (vnode->gnode->dirich->dirich_onoff.a.iv[dof]==0))
#else
        if (vnode->dof[dof] < actfield->dis[vel_dis].numeq)
#endif
        {
          INT l;

          /* But of course we assemble just to those nodes that belong to us. */
          if (vnode->proc == actintra->intra_rank)
          {
	    add_trilinos_value(lmass,
			       lmass_global.a.dv[j*vnode->numdf+dof],
			       vnode->dof[dof],
			       vnode->dof[dof]);
          }

	  for (l=0; l<pele->numnp; ++l)
	  {
	    NODE* pnode = pele->node[l];
	    if (pnode->proc == actintra->intra_rank)
	    {
	      /* Now we have the entry (actnode->dof[dof],
	       * actele->ldof[k]) to be stored. */
	      /*
	       * Again: We access the columns via the local dof numbers. */
	      switch (vele->eltyp)
	      {
#ifdef D_FLUID2_PRO
	      case el_fluid2_pro:
		/*
		 * The orientation of the gradient matrix differs from
		 * the non-trilinos version! */
		add_trilinos_value(grad,
				   gradopr_global.a.da[vnode->numdf*j+dof][l],
				   /* trilinos expects the global dof
				    * number... */
				   pnode->dof[0],
				   vnode->dof[dof]);
		break;
#endif
#ifdef D_FLUID3_PRO
	      case el_fluid3_pro:
		/*
		 * The orientation of the gradient matrix differs from
		 * the non-trilinos version! */
		add_trilinos_value(grad,
				   gradopr_global.a.da[vnode->numdf*j+dof][l],
				   pnode->dof[0],
				   vnode->dof[dof]);
		break;
#endif
	      default:
		dserror("element type %d unsupported", vele->eltyp);
	      }
	    }
	  }
        }
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	/* dirichlet dofs */
	else
	{
          if (vnode->proc == actintra->intra_rank)
          {
	    add_trilinos_value(lmass,
			       1,
			       vnode->dof[dof],
			       vnode->dof[dof]);
          }
	}
#endif
      }
    }
  }

  close_nonquad_trilinos_matrix(grad,lmass);
  close_trilinos_matrix(lmass);
  invert_trilinos_diagonal_matrix(lmass);

  /* Close the full mass matrix to enable the trilinos matrix-vector
   * product. */

  if (actsolv->sysarray_typ[sysarray]==trilinos)
  {
    close_trilinos_matrix(actsolv->sysarray[sysarray].trilinos);
  }
  else
  {
    dserror("trilinos matrix expected");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create gradient and mass matrix

  All the element call and element matrix assembling is in here. We do
  not use global_calelm because it is such a mess and our requirements
  are specific. We have to build the inverted diagonalized mass
  matrix, too.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actsolv            (i) solver
  \param sysarray           (i) solver type
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param grad             (i/o) sparse pressure gradient
  \param lmass_vec          (o) inverted diagonal mass matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calelm_laplace(FIELD *actfield,
		       PARTITION *actpart,
		       INT vdisnum,
		       INT pdisnum,
		       SOLVAR *actsolv,
		       INT press_array,
		       INT mass_array,
		       INTRA *actintra,
		       ARRAY_POSITION *ipos
  )
{
  INT i;
  TRILINOSMATRIX *press_mat;

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY emass_global;
  extern struct _ARRAY estif_global;

#ifdef DEBUG
  dstrc_enter("pm_calelm_laplace");
#endif

  if (genprob.usetrilinosalgebra)
  {
    if (actsolv->sysarray_typ[press_array]==trilinos)
    {
#ifdef TRILINOS_PACKAGE
      press_mat = actsolv->sysarray[press_array].trilinos;
#else
      dserror("trilinos missing");
#endif
    }
  }

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[vdisnum].numele; ++i)
  {
    INT j;
    ELEMENT* other = NULL;
    ELEMENT* actele = actpart->pdis[vdisnum].element[i];

    /* We need to assemble the global mass matrix. To do this in
     * parallel we cannot calculate the local elements only but have
     * to calculate any neighbouring elements. This in turn demands
     * that these elements know their pressure values. That is a
     * truely data parallel approach for discontinuous pressure
     * demands for element based communication (in contrast to the
     * node based communication we get by with in the continuous
     * case.) */

    /* Calculate gradient and mass matrix */
    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      other = actele->e.f2pro->other;
      f2pro_calpress(actele, &estif_global, &emass_global, ipos);
      break;
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
      other = actele->e.f3pro->other;
      f3pro_calpress(actele, &estif_global, &emass_global, ipos);
      break;
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }

    /* Assemble */

    /* At first lets do the global mass matrix */
    /*
     * We do this just once, thus there is no need to have this inside
     * the fluid element... */
#if 1
    /*
     * We need a different set of boundary conditions here, so we do
     * the assembling by hand. */
    assemble(mass_array,
             &emass_global,
             -1,
             NULL,
             actpart,
             actsolv,
             actintra,
             actele,
             assemble_one_matrix,
             NULL);
#endif

#if 0
    for (j=0; j<actele->numnp; ++j)
    {
      NODE* actnode = actele->node[j];

      if (actele->proc == actintra->intra_rank)
      {
	INT dj;
	for (dj=0; dj<actnode->numdf; ++dj)
	{
	  add_trilinos_value(actsolv->sysarray[mass_array].trilinos,
			     1,
			     actnode->dof[dj],
			     actnode->dof[dj]);
	}
      }
    }
#endif

#if 0
    for (j=0; j<actele->numnp; ++j)
    {
      NODE* actnode = actele->node[j];

      if (actele->proc == actintra->intra_rank)
      {
	INT dj;
	for (dj=0; dj<actnode->numdf; ++dj)
	{
	  INT k;

	  /* only normal dirichlet hack! */
	  if ((dj==1) &&
	      (actnode->gnode->dirich!=NULL) &&
	      (actnode->gnode->dirich->dirich_onoff.a.iv[0]!=0) &&
	      (actnode->gnode->dirich->dirich_onoff.a.iv[1]!=0) &&
	      (actnode->gnode->dirich->dirich_val.a.dv[0]==0) &&
	      (actnode->gnode->dirich->dirich_val.a.dv[1]==0))
	  {
	    add_trilinos_value(actsolv->sysarray[mass_array].trilinos,
			       1,
			       actnode->dof[dj],
			       actnode->dof[dj]);
	  }
	  else if ((dj==0) &&
	      (actnode->gnode->dirich!=NULL) &&
	      (actnode->gnode->dirich->dirich_onoff.a.iv[0]!=0) &&
	      (actnode->gnode->dirich->dirich_onoff.a.iv[1]!=0) &&
	      (actnode->gnode->dirich->dirich_val.a.dv[0]==1) &&
	      (actnode->gnode->dirich->dirich_val.a.dv[1]==0))
	  {
	    add_trilinos_value(actsolv->sysarray[mass_array].trilinos,
			       1,
			       actnode->dof[dj],
			       actnode->dof[dj]);
	  }
	  else
	  {
	    for (k=0; k<actele->numnp; ++k)
	    {
	      NODE* actnode2 = actele->node[k];
	      INT dk;
	      for (dk=0; dk<actnode->numdf; ++dk)
	      {
		if ((dk==1) &&
		    (actnode2->gnode->dirich!=NULL) &&
		    (actnode2->gnode->dirich->dirich_onoff.a.iv[0]!=0) &&
		    (actnode2->gnode->dirich->dirich_onoff.a.iv[1]!=0) &&
		    (actnode2->gnode->dirich->dirich_val.a.dv[0]==0) &&
		    (actnode2->gnode->dirich->dirich_val.a.dv[1]==0))
		{
		}
		else if ((dk==0) &&
		    (actnode2->gnode->dirich!=NULL) &&
		    (actnode2->gnode->dirich->dirich_onoff.a.iv[0]!=0) &&
		    (actnode2->gnode->dirich->dirich_onoff.a.iv[1]!=0) &&
		    (actnode2->gnode->dirich->dirich_val.a.dv[0]==1) &&
		    (actnode2->gnode->dirich->dirich_val.a.dv[1]==0))
		{
		}
		else
		{
		  add_trilinos_value(actsolv->sysarray[mass_array].trilinos,
				     emass_global.a.da[j*actnode->numdf+dj][k*actnode->numdf+dk],
				     actnode->dof[dj],
				     actnode2->dof[dk]);
		}
	      }
	    }
	  }
	}
      }
    }
#endif

    /* laplace pressure matrix */
    for (j=0; j<other->numnp; ++j)
    {
      NODE* actnode = other->node[j];

      if (actele->proc == actintra->intra_rank)
      {
	/* non-dirichlet dofs */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	if ((actnode->gnode->dirich==NULL) ||
	    (actnode->gnode->dirich->dirich_onoff.a.iv[genprob.ndim]==0))
#else
	  dserror("not implemented");
        if (actnode->dof[0] < actfield->dis[pdisnum].numeq)
#endif
	{
	  INT k;
	  for (k=0; k<other->numnp; ++k)
	  {
	    NODE* node2 = other->node[k];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
	    if ((node2->gnode->dirich==NULL) ||
		(node2->gnode->dirich->dirich_onoff.a.iv[genprob.ndim]==0))
#else
	      dserror("not implemented");
	    if (node2->dof[0] < actfield->dis[pdisnum].numeq)
#endif
	    {
	      add_trilinos_value(press_mat,
				 estif_global.a.da[j][k],
				 actnode->dof[0],
				 node2->dof[0]);
	    }
	  }
	}
	else
	{
	  add_trilinos_value(press_mat,
			     1,
			     actnode->dof[0],
			     actnode->dof[0]);
	}
      }
    }
  }

  if (genprob.usetrilinosalgebra)
  {
    if (actsolv->sysarray_typ[mass_array]==trilinos)
    {
#ifdef TRILINOS_PACKAGE
      close_trilinos_matrix(actsolv->sysarray[mass_array].trilinos);
      close_trilinos_matrix(actsolv->sysarray[press_array].trilinos);
#else
      dserror("trilinos missing");
#endif
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create gradient

  Calculate pressure discrete gradient. For simplicity we do it
  elementwise here. Maybe a global matrix-vector multiplication would
  be faster. Could be improved.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param rhs                (o) result vector

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calprhs(FIELD *actfield,
                PARTITION *actpart,
                INT disnum,
                INTRA *actintra,
                ARRAY_POSITION *ipos,
                INT numpdof,
                DIST_VECTOR* rhs)
{
  INT i;

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY eforce_global;

#ifdef DEBUG
  dstrc_enter("pm_calprhs");
#endif

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT k;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    if (actele->proc == actintra->intra_rank)
    {
      /* Calculate gradient and mass matrix */
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        f2pro_calprhs(actele, ipos);

        /* Assemble */
        /* We have discontinuous pressure here. No need to loop the
         * nodes. */
        for (k=0; k<numpdof; ++k)
        {
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */
          dsassert((actele->e.f2pro->ldof[k] >= 0) &&
                   (actele->e.f2pro->ldof[k] < rhs->numeq),
                   "local dof number out of range");
          rhs->vec.a.dv[actele->e.f2pro->ldof[k]] += eforce_global.a.dv[k];
        }
        break;
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
        f3pro_calprhs(actele, ipos);

        /* Assemble */
        /* We have discontinuous pressure here. No need to loop the
         * nodes. */
        for (k=0; k<numpdof; ++k)
        {
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */
          dsassert((actele->e.f3pro->ldof[k] >= 0) &&
                   (actele->e.f3pro->ldof[k] < rhs->numeq),
                   "local dof number out of range");
          rhs->vec.a.dv[actele->e.f3pro->ldof[k]] += eforce_global.a.dv[k];
        }
        break;
#endif
      default:
        dserror("element type %d unsupported", actele->eltyp);
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create gradient

  Calculate pressure discrete gradient. For simplicity we do it
  elementwise here. Maybe a global matrix-vector multiplication would
  be faster. Could be improved.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param rhs                (o) result vector

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calprhs_cont(FIELD *actfield,
		     PARTITION *actpart,
		     INT disnum,
		     INTRA *actintra,
		     ARRAY_POSITION *ipos,
		     DIST_VECTOR* rhs,
		     DOUBLE* full_rhs)
{
  INT i;

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY eforce_global;

#ifdef DEBUG
  dstrc_enter("pm_calprhs_cont");
#endif

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT k;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    /* Calculate gradient and mass matrix */
    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      ELEMENT* other;
      other = actele->e.f2pro->other;
      f2pro_calprhs(actele, ipos);

      /* Assemble */
      /* We have discontinuous pressure here. No need to loop the
       * nodes. */
      for (k=0; k<other->numnp; ++k)
      {
	NODE* node = other->node[k];

	if (node->proc == actintra->intra_rank)
	{
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */
          dsassert((node->dof[0] >= 0) &&
                   (node->dof[0] < rhs->numeq),
                   "local dof number out of range");
          full_rhs[node->dof[0]] += eforce_global.a.dv[k];
	}
      }
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      ELEMENT* other;
      other = actele->e.f3pro->other;
      f3pro_calprhs(actele, ipos);

      /* Assemble */
      /* We have discontinuous pressure here. No need to loop the
       * nodes. */
      for (k=0; k<other->numnp; ++k)
      {
	NODE* node = other->node[k];

	if (node->proc == actintra->intra_rank)
	{
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */
          dsassert((node->dof[0] >= 0) &&
                   (node->dof[0] < rhs->numeq),
                   "local dof number out of range");
          full_rhs[node->dof[0]] += eforce_global.a.dv[k];
	}
      }
      break;
    }
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief call elements to create rhs for velocity correction

  This is used if we do not lump the mass matrix and have to solve for
  the velocity correction...

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param rhs                (o) result vector

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_calvrhs(FIELD *actfield,
		PARTITION *actpart,
		INT disnum,
		INTRA *actintra,
		ARRAY_POSITION *ipos,
		DIST_VECTOR* rhs,
		DOUBLE* full_rhs)
{
  INT i;

  /* Variables from global_calelm. These are filled by the element
   * routines. */
  extern struct _ARRAY eforce_global;

#ifdef DEBUG
  dstrc_enter("pm_calvrhs");
#endif

  /* calculate matrix values */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT k;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    /* Calculate gradient and mass matrix */
    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      f2pro_calvrhs(actele, ipos);
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      f3pro_calvrhs(actele, ipos);
      break;
    }
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }

    /* Assemble */
    /* We have discontinuous pressure here. No need to loop the
     * nodes. */
    for (k=0; k<actele->numnp; ++k)
    {
      NODE* node = actele->node[k];

      if (node->proc == actintra->intra_rank)
      {
	INT j;
	for (j=0; j<node->numdf; ++j)
	{
	  dsassert((node->dof[j] >= 0) &&
		   (node->dof[j] < rhs->numeq),
		   "dof number out of range");
	  full_rhs[node->dof[j]] += eforce_global.a.dv[node->numdf*k+j];
	}
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief update pressure values after the pressure equation has been
  solved

  Unfortunately we need the new pressure values in more that just the
  local elements (to be able to calculate the pressure gradient in the
  momentum equation). So we follow the ccarat way, introduce total
  vectors and allreduce the pressure increments.

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param numpdof            (i) number of pressure dofs
  \param sol                (i) pressure solution to be distributed
  \param dta                (i) current time step size

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_press_update(FIELD *actfield,
                     PARTITION *actpart,
                     INT disnum,
                     INTRA *actintra,
                     ARRAY_POSITION *ipos,
                     INT numpdof,
                     DIST_VECTOR* sol,
                     DOUBLE dta)
{
  INT i;

  DOUBLE* p;
  DOUBLE* press1;
#ifdef PARALLEL
  DOUBLE* press2;
#endif

#ifdef DEBUG
  dstrc_enter("pm_press_update");
#endif

  press1 = (DOUBLE*)CCACALLOC(numpdof*actfield->dis[disnum].numele, sizeof(DOUBLE));
#ifdef PARALLEL
  press2 = (DOUBLE*)CCACALLOC(numpdof*actfield->dis[disnum].numele, sizeof(DOUBLE));
#endif

  /* gather pressure */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT k;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    if (actele->proc == actintra->intra_rank)
    {
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
      {
        FLUID2_PRO* f2pro;
        f2pro = actele->e.f2pro;
        for (k=0; k<numpdof; ++k)
        {
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */

          /* Remember the pressure increment, we need it to update the
           * velocity, too. */
          dsassert((f2pro->ldof[k] >= 0) &&
                   (f2pro->ldof[k] < sol->numeq),
                   "local dof number out of range");
          dsassert((f2pro->dof[k] >= 0) &&
                   (f2pro->dof[k] < numpdof*actfield->dis[disnum].numele),
                   "global dof number out of range");
          press1[f2pro->dof[k]] = sol->vec.a.dv[f2pro->ldof[k]];
        }
        break;
      }
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
      {
        FLUID3_PRO* f3pro;
        f3pro = actele->e.f3pro;
        for (k=0; k<numpdof; ++k)
        {
          /* there are no dirichlet conditions on the pressure dofs
           * allowed... currently. */

          /* Remember the pressure increment, we need it to update the
           * velocity, too. */
          dsassert((f3pro->ldof[k] >= 0) &&
                   (f3pro->ldof[k] < sol->numeq),
                   "local dof number out of range");
          dsassert((f3pro->dof[k] >= 0) &&
                   (f3pro->dof[k] < numpdof*actfield->dis[disnum].numele),
                   "global dof number out of range");
          press1[f3pro->dof[k]] = sol->vec.a.dv[f3pro->ldof[k]];
        }
        break;
      }
#endif
      default:
        dserror("element type %d unsupported", actele->eltyp);
      }
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(press1, press2, numpdof*actfield->dis[disnum].numele, MPI_DOUBLE,
                MPI_SUM, actintra->MPI_INTRA_COMM);
  p = press2;
#else
  p = press1;
#endif

  /* Update the pressure on all elements, including foreign ones. */
  for (i=0; i<actfield->dis[disnum].numele; ++i)
  {
    INT k;
    ELEMENT* actele = &(actfield->dis[disnum].element[i]);

    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      FLUID2_PRO* f2pro;
      f2pro = actele->e.f2pro;
      for (k=0; k<numpdof; ++k)
      {
        /* there are no dirichlet conditions on the pressure dofs
         * allowed... currently. */

        /* Remember the pressure increment, we need it to update the
         * velocity, too. */
        dsassert((f2pro->dof[k] >= 0) &&
                 (f2pro->dof[k] < numpdof*actfield->dis[disnum].numele),
                 "global dof number out of range");

#ifdef DEBUG
        /* yet again paranoia tests */
        if (actele->proc==actintra->intra_rank)
        {
          dsassert(p[f2pro->dof[k]]==sol->vec.a.dv[f2pro->ldof[k]], "allreduce failed");
        }
#endif

        f2pro->phi[k] = p[f2pro->dof[k]];
        f2pro->press[k] += 2/dta*f2pro->phi[k];
      }
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      FLUID3_PRO* f3pro;
      f3pro = actele->e.f3pro;
      for (k=0; k<numpdof; ++k)
      {
        /* there are no dirichlet conditions on the pressure dofs
         * allowed... currently. */

        /* Remember the pressure increment, we need it to update the
         * velocity, too. */
        dsassert((f3pro->dof[k] >= 0) &&
                 (f3pro->dof[k] < numpdof*actfield->dis[disnum].numele),
                 "global dof number out of range");

#ifdef DEBUG
        /* yet again paranoia tests */
        if (actele->proc==actintra->intra_rank)
        {
          dsassert(p[f3pro->dof[k]]==sol->vec.a.dv[f3pro->ldof[k]], "allreduce failed");
        }
#endif

        f3pro->phi[k] = p[f3pro->dof[k]];
        f3pro->press[k] += 2/dta*f3pro->phi[k];
      }
      break;
    }
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

  CCAFREE(press1);
#ifdef PARALLEL
  CCAFREE(press2);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief update velocity after the pressure increment is known

  \param actfield           (i) actual field
  \param actpart            (i) actual partition
  \param disnum             (i) number of discretization
  \param actintra           (i) communicator
  \param ipos               (i) node positions
  \param lmass              (i) global inverted lumped mass matrix
  \param rhs1               (i) temp array (velocity dofs)
  \param rhs2               (i) temp array (velocity dofs)

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_vel_update(FIELD *actfield,
                   PARTITION *actpart,
                   INT disnum,
                   INTRA *actintra,
                   ARRAY_POSITION *ipos,
#ifdef PM_TRILINOS
		   TRILINOSMATRIX* lmass,
		   SOLVAR *actsolv,
		   INT actsysarray,
#else
                   DOUBLE* lmass,
#endif
                   DOUBLE* rhs1,
                   DOUBLE* rhs2)
{
  INT i;
  INT velnp;
#ifndef PM_TRILINOS
  DOUBLE* gradip;
#endif

#ifdef DEBUG
  dstrc_enter("pm_vel_update");
#endif

  /* Due to missing processor local dof numbers we can only assemble
   * to global vectors. In order to allreduce it a second one is
   * needed. */
  memset(rhs1, 0, actfield->dis[disnum].numeq*sizeof(DOUBLE));
#if defined(PARALLEL) || defined(PM_TRILINOS)
  memset(rhs2, 0, actfield->dis[disnum].numeq*sizeof(DOUBLE));
#endif

  /* We need to gather the term G*phi in an distributed vector */
  for (i=0; i<actpart->pdis[disnum].numele; ++i)
  {
    INT j;
    ELEMENT* actele = actpart->pdis[disnum].element[i];

    if (actele->proc == actintra->intra_rank)
    {
      /* Calculate pressure increment gradient */
      switch (actele->eltyp)
      {
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        f2pro_calvelupdate(actele, ipos);
        break;
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
        f3pro_calvelupdate(actele, ipos);
        break;
#endif
      default:
        dserror("element type %d unsupported", actele->eltyp);
      }

      /* Assemble */
      for (j=0; j<actele->numnp; ++j)
      {
        INT dof;
        NODE* actnode = actele->node[j];

        for (dof=0; dof<actnode->numdf; ++dof)
        {
          INT gdof;
          gdof = actnode->dof[dof];
          if (gdof < actfield->dis[disnum].numeq)
          {
            /* eforce_global filled by element call above */
            extern ARRAY eforce_global;
            rhs1[gdof] += eforce_global.a.dv[actnode->numdf*j+dof];
          }
        }
      }
    }
  }

#ifdef PM_TRILINOS

  /* multiply the new rhs with lmass and allreduce the result */

  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_vec(&(actsolv->rhs[1]));

  assemble_vec(actintra,
	       &(actsolv->sysarray_typ[actsysarray]),
	       &(actsolv->sysarray[actsysarray]),
	       &(actsolv->rhs[0]),
	       rhs1,
	       1.0
    );

  matvec_trilinos(&(actsolv->rhs[1]),
		  &(actsolv->rhs[0]),
		  lmass);

  solserv_reddistvec(&(actsolv->rhs[1]),
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]),
		     rhs2,
		     actsolv->rhs[1].numeq_total,
		     actintra);

#else /* PM_TRILINOS */

#ifdef PARALLEL
  MPI_Allreduce(rhs1, rhs2, actfield->dis[disnum].numeq, MPI_DOUBLE, MPI_SUM,
                actintra->MPI_INTRA_COMM);
  gradip = rhs2;
#else
  gradip = rhs1;
#endif

#endif /* PM_TRILINOS */

  velnp = ipos->velnp;

  /* Now the velocity dofs can be updated. */
  for (i=0; i<actfield->dis[disnum].numnp; ++i)
  {
    INT dof;
    NODE* actnode = &(actfield->dis[disnum].node[i]);

    for (dof=0; dof<actnode->numdf; ++dof)
    {
      INT gdof;
      gdof = actnode->dof[dof];

      /*
       * We do not want to change the prescribed boundaries. Do we
       * have to worry because the calculated decrement would change
       * them? */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
      if ((actnode->gnode->dirich==NULL) ||
	  (actnode->gnode->dirich->dirich_onoff.a.iv[dof]==0))
#else
      if (gdof < actfield->dis[disnum].numeq)
#endif
      {
#ifdef PM_TRILINOS

        actnode->sol_increment.a.da[velnp][dof] -= rhs2[gdof];

#else /* PM_TRILINOS */

        actnode->sol_increment.a.da[velnp][dof] -= lmass[gdof]*gradip[gdof];

#endif /* PM_TRILINOS */
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief output to the screen

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void pm_out_screen_header(INT numeq,
                          INT numeq_total,
                          INTRA *actintra,
                          FILE *out,
                          FLUID_DYNAMIC *fdyn)
{
  INT i;
#ifdef DEBUG
  dstrc_enter("pm_screen_header");
#endif

  /*------------------------------------------------ output to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0;i<par.nprocs;i++)
    if (par.myrank==i)
      printf("PROC  %3d | FIELD FLUID     | number of equations      : %10d \n",
             par.myrank,numeq);
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  if (par.myrank==0)
    printf("          | FIELD FLUID     | total number of equations: %10d \n",numeq_total);
  if (par.myrank==0) printf("\n\n");

  /* write general data to .out */
  if (par.myrank==0)
  {
    fprintf(out,"max. values:\n");
    fprintf(out,"============\n");


    /* table head */
    fprintf(out," time |            |fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
    case fncc_Linf: /* infinity norm */
      fprintf(out,"inf-norm");
      break;
    case fncc_L1: /* L_1 norm */
      fprintf(out,"L_1-norm");
      break;
    case fncc_L2: /* L_2 norm */
      fprintf(out,"L_2-norm");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
    }

    fprintf(out," | steady state in ");

    switch(fdyn->stnorm)
    {
    case fnst_Linf: /* infinity norm */
      fprintf(out,"inf-norm|");
      break;
    case fnst_L1: /* L_1 norm */
      fprintf(out,"L_1-norm|");
      break;
    case fnst_L2: /* L_2 norm */
      fprintf(out,"L_2-norm|");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
    }
    fprintf(out,"    total   |\n");

    fprintf(out," step |  sim. time | ite |     vel.   |     pre.   |     vel.   |     pre.   | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------\n");

    fprintf(out,"%5d | %10.3E | %3d |        %10.3E       |        %10.3E       |            |\n",
            fdyn->nstep,fdyn->maxtime,fdyn->itemax,fdyn->ittol,fdyn->sttol);
    fprintf(out,"-------------------------------------------------------------------------------------------\n");

    fprintf(out,"\n\ntimeloop:  ");

    switch(fdyn->iop)
    {
    case 1:
      fprintf(out,"Generalised Alpha\n");
      break;
    case 4:
      fprintf(out,"One-Step-Theta\n");
      break;
    case 7:
      fprintf(out,"BDF2\n");
      break;
    default:
      dserror("parameter out of range: IOP\n");
    }

    fprintf(out,"=========\n");

    /* table head */
    fprintf(out," time |            |fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
    case fncc_Linf: /* infinity norm */
      fprintf(out,"inf-norm");
      break;
    case fncc_L1: /* L_1 norm */
      fprintf(out,"L_1-norm");
      break;
    case fncc_L2: /* L_2 norm */
      fprintf(out,"L_2-norm");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
    }

    fprintf(out," | steady state in ");

    switch(fdyn->stnorm)
    {
    case fnst_Linf: /* infinity norm */
      fprintf(out,"inf-norm|");
      break;
    case fnst_L1: /* L_1 norm */
      fprintf(out,"L_1-norm|");
      break;
    case fnst_L2: /* L_2 norm */
      fprintf(out,"L_2-norm|");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
    }
    fprintf(out,"    total   |\n");

    fprintf(out," step |  sim. time | ite |     vel.   |     pre.   |     vel.   |     pre.   | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------\n");

    fflush(out);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
