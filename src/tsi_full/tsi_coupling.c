/*----------------------------------------------------------------------*/
/*!
\file
\brief Coupling of TSI problems

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/

/*!
\addtogroup TSI
*//*! @{ (documentation module open)*/
/*----------------------------------------------------------------------*/


#ifdef D_TSI


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"  /* this includes enums.h */
#include "../solver/solver.h"
#include "tsi_prototypes.h"
#ifdef D_WALL1
#include "../wall1/wall1.h"
#endif
#ifdef D_BRICK1
#include "../brick1/brick1.h"
#endif
#ifdef D_THERM2
#include "../therm2/therm2.h"
#endif
#ifdef D_THERM3
#include "../therm3/therm3.h"
#endif


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable genprob (of struct _GENPROB) is defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Pointer to dynamic variables

Dynamic controls of the differennt fields, defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern ALLDYNA *alldyn;


/*----------------------------------------------------------------------*/
/*!
\brief Pointer to design entities

Design entities, such as DNODE, DLINE, DSURF, DVOL; 
defined in global_control.c 

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*----------------------------------------------------------------------*/
/*!
\brief Initialise TSI coupling conditions

Create multifield solution history and set pointers to the corresponding
nodes of the other fields

\param *structfield   FIELD         (i)      structure field
\param  disnum_s      INT           (i)      structure discretisation
\param *thermfield    FIELD         (i)      thermal field
\param  disnum_t      INT           (i)      thermal discretisation

\return void

\author bborn
\date 03/06
*/
/*-----------------------------------------------------------------------*/
void tsi_coupling(FIELD *structfield,
		    INT disnum_s,
		    FIELD *thermfield,
		    INT disnum_t)

{

  INT numfld;  /* number of fields */
  INT numsnp;  /* number of structure nodes */
  INT numtnp;  /* number of thermal nodes */
  INT numdf;  /* number of dofs */
  INT numc = 0;  /* number of columns in mf */
  INT numsf;  /* number of structural fields */
  INT numtf;  /* number of thermal fields */
  INT e, i, ii, j, k;  /* simply some counters */
  INT ierr;  /* flag */
  INT dim;  /* dimension of problem */
  INT tfound;  /* flag for found thermal node */
  INT partner;  /* flag for checking thermal partnership */

  DOUBLE tol = EPS4;  /* critcal tolerance for node distance

  /* declare and nullify actual nodes */
  NODE *actsnode = NULL;  /* actual structure node */
  NODE *acttnode = NULL;  /* actual thermal node */
  GNODE *actsgnode = NULL;  /* actual structure geometry node */
  GNODE *acttgnode = NULL;  /* actual thermal geometry node */

  /* actual element pointer */
  ELEMENT *actele;
  ELEMENT *actsele;  /* for structure */
  ELEMENT *acttele;  /* for thermal */

  /* The array holds flags for each thermal node.
   * The array avoids redundant repetitions. */
  ARRAY tindex_a;
  INT *tindex;

  /* TSI dynamic configuration */
  TSI_DYNAMIC *tsidyn;  /* pointer to TSI dynamic control */

  INT numsel;  /* number od structure elements */
  INT nen;  /* number of element nodes */
  INT maxadjele;
  INT actId;  /* current element ID */

  /* Array of element IDs of each element node */
  ARRAY adjele_a;
  INT **adjele;


#ifdef DEBUG
  dstrc_enter("tsi_initcoupling");
#endif

  /* total number of fields */
  numfld = genprob.numfld;


  /* set pointer of TSI dynamic control */
  /* This should be the 3rd entry, ie C-index is 2:
   * numfld=2, ie alldyn[2].tsidyn */
  tsidyn = alldyn[numfld].tsidyn;


  /* find number of nodes in different fields */
  numsnp = structfield->dis[disnum_s].numnp;
  numtnp = thermfield->dis[disnum_t].numnp;
  numsf = genprob.numsf;
  numtf = genprob.numtf;
  dim = genprob.ndim;

#ifdef PERF
  perf_begin(53);
#endif

  /*--------------------------------------------------------------------*/
  /* allocate space for multifield solutions
   * 
   * requests:
   * (1) structure needs temperature
   * (2) thermal needs structural displacements
   *
   * offers:
   * (1) structure offers displacements
   * (2) thermal offers temperature
   */ 

  /*--------------------------------------------------------------------*/
  /* structure */

  /* multifield solution history of structural nodes: numc=6 entries:
   * actsnode->sol_mf.a.da[0][i]: actual displacements transfered to thermal
   * actsnode->sol_mf.a.da[1][i]: displacements of old iteration step
   * actsnode->sol_mf.a.da[2][i]: displacements of old time step
   * actsnode->sol_mf.a.da[3][i]: dispi
   * actsnode->sol_mf.a.da[4][i]: couplingforces at the end of time step
   * actsnode->sol_mf.a.da[5][i]: coupling forces at the beginning of 
   *                              time step */
  numc = 6;

  /* loop struct nodes */
  for (i=0; i<numsnp; i++)
  {
    actsnode  = &(structfield->dis[disnum_s].node[i]);
    actsgnode = actsnode->gnode;
    numdf = actsnode->numdf;  /* number of DOFs at node */

    /* Allocate space for multifield solution, i.e. the global
     * displacement solution can be copied/disassembled to
     * each structural node. 
     * The thermal field element can access the structure 
     * displacements by getting hold of the displacements of the
     * conforming element and its nodes */
    /* at node, actsnode->sol_mf is not a pointer */
    amdef("sol_mf", &actsnode->sol_mf, numc, numdf, "DA");
    amzero(&(actsnode->sol_mf));

    /* at geometry node */  /* why dimension 3 ?! */
    actsgnode->mfcpnode = (NODE**) CCACALLOC(numfld, sizeof(NODE*));
    if (actsgnode->mfcpnode == NULL)
    {
      dserror("Allocation of coupling node pointers failed");
    }
    for (j=0; j<numfld; j++) 
    {
      actsgnode->mfcpnode[j] = NULL;
    }

  } /* end of loop over struct nodes */


  /*--------------------------------------------------------------------*/
  /* thermal */

  /* multifield solution history of thermal nodes: numc=3 entries:
   * acttnode->sol_mf.a.da[0][i]: actual temperature transfered to structure
   * acttnode->sol_mf.a.da[1][i]: temperature of old iteration step
   * acttnode->sol_mf.a.da[2][i]: temperature of old time step */
  numc = 3;
  
  /* loop thermal nodes */
  for (i=0; i<numtnp; i++)
  {
    acttnode  = &(thermfield->dis[disnum_t].node[i]);
    acttgnode = acttnode->gnode;
    numdf = acttnode->numdf;

    amdef("sol_mf", &acttnode->sol_mf, numc, numdf, "DA");
    amzero(&(acttnode->sol_mf));

    acttgnode->mfcpnode = (NODE**) CCACALLOC(numfld, sizeof(NODE*));
    if (acttgnode->mfcpnode == NULL)
    {
      dserror("Allocation of coupling node pointers failed");
    }
    for (j=0; j<numfld; j++)
    {
      acttgnode->mfcpnode[j] = NULL;
    }

  } /* end of loop over thermal nodes */


#ifdef PERF
  perf_end(53);
#endif


#ifdef PERF
  perf_begin(54);
#endif


  /*--------------------------------------------------------------------*/
  /* find TSI coupled nodes
   * and set pointers to the corresponding nodes of the other fields */

  /* create and initialise index arrays */
  tindex = amdef("tindex", &tindex_a, numtnp, 1, "IV");
  for (i=0; i<numtnp; i++)  /* loop over all thermal nodes */
  {
    tindex[i] = 1;
  }



  /* loop structure nodes  */
  for (i=0; i<numsnp; i++)
  {
    actsnode = &(structfield->dis[disnum_s].node[i]);
    actsgnode = actsnode->gnode;

    /* initialise partner counter */
    partner = 0;
    /* initialise adjacent element counter */
    j = 0;
    /* loop all adjacent elements of node */
    while ((j < actsnode->numele) && (partner <= 0))
    {

      actele = actsnode->element[j];

      switch (actele->eltyp)
      {
#ifdef D_WALL1
        case el_wall1:
          if (actele->e.w1->tsi_couptyp == tsi_coup_thermconf) 
          {
            partner++;
          }
          break;
#endif
#ifdef D_BRICK1
        case el_brick1:
          if (actele->e.c1->tsi_couptyp == tsi_coup_thermconf)
          {
            partner++;
          }
          break;
#endif
        default:
          dserror("eltyp unknown\n");
      }  /* switch (actele->eltyp) */

      /* increment element counter */
      j++;

    }  /* end of for (j=0;j<actfnode->numele;j++) */


#ifdef PERF
    perf_begin(56);
#endif


#ifdef PERF
    perf_end(56);
#endif


#ifdef PERF
    perf_begin(57);
#endif

    /* initialise flag marking a found corresponding thermal node */
    tfound = 0;

    /* check wether structure node has a thermal partner */
    if (partner > 0)
    {
      /* tfound enlightens if coincident node is found
       * tfound is initialised above */
      /* initialise thermal node counter */
      j = 0;
      /* loop all thermal nodes and find corresponding */
      while ((j < numtnp) && (tfound == 0))
      {
        /* node j was not processed previously */
        if (tindex[j] != 0)
        {
          /* get pointer on actual thermal node */
          acttnode = &(thermfield->dis[disnum_t].node[j]);
          /* check distance of coords */
          cheque_distance(&(actsnode->x[0]), &(acttnode->x[0]), 
                          tol, &ierr);
          /* strike! actsnode and acttnode have identical coords */
          if (ierr == 1)
          {
            /* increment found coincident thermal nodes */
            tfound++;
            /* get actual thermal geometry node */
            acttgnode = acttnode->gnode;
            /* jth thermal node processed */
            tindex[j] = 0;
          }
        }
        /* increment nodal counter */
        j++;
      }  /* end of while ((j < numtnp) && (tfound > 0)) */
      /* The loop stops if a coincident thermal node is found,
       * thus, the variables acttnode and acttgnode
       * point at this coincident thermal node! */
    }  /* end of if (partner > 0) */


#ifdef PERF
    perf_end(57);
#endif


    /* set pointer of geometry node to node (always possible) */
    actsgnode->mfcpnode[numsf] = actsnode;
    /* set pointers to corresponding nodes */
    if (tfound > 0)
    {
      actsgnode->mfcpnode[numtf] = acttnode;
      acttgnode->mfcpnode[numsf] = actsnode;
      acttgnode->mfcpnode[numtf] = acttnode;
    }
    else
    {
      dserror("Corresponding thermal node was not found!");
    }

  }  /* end of for (i=0;i<numsnp;i++) */


#ifdef PERF
  perf_end(54);
#endif

  /* clean up temporary index */
  amdel(&tindex_a);

  /*--------------------------------------------------------------------*/
  /* CHECK THIS --- POSSIBLY UNNECESSARY */
  if (genprob.visual > 0) goto end;


#ifdef PERF
  perf_begin(55);
#endif


  /*--------------------------------------------------------------------*/
  /* plausibility checks */

  /* loop all structure nodes */
  for (i=0; i<numsnp; i++)
  {
    /* get pointers at current structure node */
    actsnode  = &(structfield->dis[disnum_s].node[i]);
    actsgnode = actsnode->gnode;

    /* if (actfgnode->fsicouple==NULL) continue; */

    /* get corresponding current thermal node */
    acttnode = actsgnode->mfcpnode[numtf];
    acttgnode = acttnode->gnode;


    if (acttgnode == NULL)
    {
      dserror("No thermal node for structure node: #%d", actsnode->Id);
    }


    /* not sure if we need the latter */
    /* check locsys */
    if (actsnode->locsysId != 0)
    {
      dserror("No locsys at TSI coupling node: #%d", actsnode->Id);
    }
    if (acttnode->locsysId != 0)
    {
      dserror("No locsys at TSI coupling node: #%d", acttnode->Id);
    }

  }  /* for (i=0;i<numanp;i++) */



#ifdef PERF
  perf_end(55);
#endif


  /*--------------------------------------------------------------------*/
  /* identify conforming elements */

  /* number of structure elements */
  numsel = structfield->dis[disnum_s].numele;
  
  /* loop all structure elements */
  for (e=0; e<numsel; e++)
  {
    /* set pointer to current element */
    actsele = &(structfield->dis[disnum_s].element[e]);

    /* number of element nodes */
    nen = actsele->numnp;

    /* maximum number adjacent thermal elements of one element node */
    maxadjele = 0;
    for (k=0; k<nen; k++)
    {
      actsnode = actsele->node[k];
      actsgnode = actsnode->gnode;
      acttnode = actsgnode->mfcpnode[numtf];
      maxadjele = IMAX(maxadjele, acttnode->numele);
    }

    /* allocate arrays for storing element IDs */
    /* for every element node an array row is of maxadjele length
     * plus one additional row */
    adjele = amdef("adjele", &(adjele_a), nen, maxadjele, "IA");
    j = -1;
    aminit(&adjele_a, &(j));  /* initialise with -1 */
      
    /* grab thermal element IDs of first conforming thermal node */
    actsnode = actsele->node[0];
    actsgnode = actsnode->gnode;
    acttnode = actsgnode->mfcpnode[numtf];
    for (i=0; i<acttnode->numele; i++)
    {
      adjele[0][i] = acttnode->element[i]->Id_loc;
    }

    /* grab thermal element IDs of j-th conforming thermal node
     * add only those which also occur in (j-1)th node */
    for (k=1; k<nen; k++)
    {
      actsnode = actsele->node[k];
      actsgnode = actsnode->gnode;
      acttnode = actsgnode->mfcpnode[numtf];
      ii = 0;
      for (i=0; i<acttnode->numele; i++)  /* loop of all adj elements */
      {
        actId = acttnode->element[i]->Id_loc;
        j = 0;
        while ( (adjele[k-1][j] != -1) && (j < maxadjele) )
        {
          if (actId == adjele[k-1][j])
          {
            adjele[k][ii] = actId;
            ii = ii + 1;
            break;
          }  /* end of if */
          j = j + 1;
        }  /* end of while */
      }  /* end of for */
    }  /* end of for */

    /* set ID of conforming thermal element */
    actId = adjele[nen-1][0];
    if (actId == -1)
    {
      dserror("Conforming thermal element was not found!");
    }
    else
    {
      /* printf("The ID: %d\n", actId); */
    }
    acttele = &(thermfield->dis[disnum_s].element[actId]);

    /* set conforming element pointers */

    /* structure element gets conforming thermal element pointer */
    switch (actsele->eltyp)
    {
#ifdef D_WALL1
        case el_wall1:
          actsele->e.w1->therm_ele = acttele;
          break;
#endif
#ifdef D_BRICK1
        case el_brick1:
          actsele->e.c1->therm_ele = acttele;
          break;
#endif
        default:
          dserror("eltyp unknown\n");
    }  /* end of switch */

    /* thermal element gets conforming structure element (pointer) */
    switch (acttele->eltyp)
    {
#ifdef D_THERM2
        case el_therm2:
          acttele->e.th2->struct_ele = actsele;
          break;
#endif
#ifdef D_THERM3
        case el_therm3:
          acttele->e.th3->struct_ele = actsele;
          break;
#endif
        default:
          dserror("eltyp unknown\n");
    }  /* end of switch */
    
    /* deallocate array of adjacent element IDs */
    amdel(&adjele_a);
    
  }  /* end of for */



  /*--------------------------------------------------------------------*/
  /* print out coupling */  
  if (genprob.visual == 0)
  {
    out_tsi(structfield);
  }

  /* this is the end, my friend */
  end:


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

} /* end of tsi_coupling */



#endif  /* end of #ifdef D_TSI */

/*! @} (documentation module close)*/

