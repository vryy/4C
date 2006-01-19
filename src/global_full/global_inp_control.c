/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/


#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fsi_full/fsi_prototypes.h"

#ifdef D_SSI
#include "../ssi_full/ssi_prototypes.h"
#endif



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;


#ifdef DEBUG
/*!----------------------------------------------------------------------
  \brief the tracing variable

  <pre>                                                         m.gee 8/00
  defined in pss_ds.c, declared in tracing.h
  </pre>
 *----------------------------------------------------------------------*/
extern struct _CCA_TRACE         trace;
#endif


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


#ifdef D_AXISHELL
/*!----------------------------------------------------------------------
  \brief prototypes callable only in this file
 *-----------------------------------------------------------------------*/
void interpolate_axishell_conds(DISCRET  *actdis);
#endif /* D_AXISHELL */


/*----------------------------------------------------------------------*
  | input of control, element and load information         m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntainp()
{

  INT i,j;


  /* the input of the tracing option has not been done yet, so
     we have to make the dstrc_enter 'by hand'
     */
#ifdef DEBUG
  trace.actroutine = trace.actroutine->next;
  trace.actroutine->name = "ntainp";
  trace.actroutine->dsroutcontrol=dsin;
  trace.deepness++;
#endif


  /* input of not mesh or time based problem data  */
#ifdef PERF
  perf_begin(3);
#endif

  inpctr();

#ifdef PERF
  perf_end(3);
#endif



#ifdef D_OPTIM
  /* input of optimization */
  if (genprob.probtyp == prb_opt) inpctropt();
#endif



  /* input of design data:
   * ---------------------
   * read description of design volumes, surfaces, lines and nodes, allocate
   * and fill the structures DVOL, DSURF, DLINE and DNODE
   */
#ifdef PERF
  perf_begin(4);
#endif

  inpdesign();

#ifdef PERF
  perf_end(4);
#endif



  /* build the topology among the design elements, allocate and create the
   * pointers connecting DVOL, DSURF, DLINE and DNODE
   */
#ifdef PERF
  perf_begin(5);
#endif

  inpdesign_topology_design();

#ifdef PERF
  perf_end(5);
#endif



  /*----------------------------------------------------------------------*/
  /* NOTE: the materials have to be read before the input of the elements */
  /*       as these informations are needed for shell9 element            */
  /*                                                             sh 10/02 */
  /*                                                                      */
  /*-------------------------------------------------- input of materials */
#ifdef PERF
  perf_begin(6);
#endif

  inp_material();

  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();

#ifdef PERF
  perf_end(6);
#endif



  /* input of meshes:
   * ----------------
   * read the fe-nodes in NODE, the fe-elements in ELEMENT and build the
   * topology among ELEMENTs and NODEs. Assign the disctretization to
   * the different fields in a multifield problem dependent on type
   * of element
   */
#ifdef PERF
  perf_begin(7);
#endif

  inpfield();

#ifdef PERF
  perf_end(7);
#endif



  /* built the detailed topology of the discretizations:
   * ---------------------------------------------------
   * for each existing field and discretization build a detailed topology
   * of GNODEs, GLINEs, GSURFs and GVOLs and connect them to the topology
   * of NODEs and ELEMENTs
   */
#ifdef PERF
  perf_begin(8);
#endif

  for (i=0; i<genprob.numfld; i++)
    for (j=0; j<field[i].ndis; j++)
      if (field[i].dis[j].disclass != dc_subdiv_calc)
        inp_detailed_topology(&(field[i].dis[j]));

#ifdef PERF
  perf_end(8);
#endif



  /* design-fe topology:
   * -------------------
   * Read which node is on which design object from file.
   * For each field and discretization build the pointers among design
   * and FE-objects DVOL<->GVOL,DSURF<->GSURF,DLINE<->GLINE,DNODE<->GNODE.
   */
#ifdef PERF
  perf_begin(9);
#endif

  for (i=0; i<genprob.numfld; i++)
    for (j=0; j<field[i].ndis; j++)
      if (field[i].dis[j].disclass != dc_subdiv_calc)
        inpdesign_topology_fe(&(field[i].dis[j]));



#ifdef PERF
  perf_end(9);
#endif



  /* generate the second discretization for subdivision of elements */
  /*================================================================*/
#ifdef SUBDIV

    if (genprob.create_dis != 1 && genprob.create_ale != 1)
    {
      design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));
      design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;
      design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;
      design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;

      for (i=0; i<design->ndnode; i++)
      {
        design->dnode_fenode2[i]= (INT*)CCAMALLOC(sizeof(INT));
      }

      for (i=0; i<design->ndline; i++)
      {
        design->dline_fenode2[i]= (INT*)CCAMALLOC(sizeof(INT));
      }

      for (i=0; i<design->ndsurf; i++)
      {
        design->dsurf_fenode2[i]= (INT*)CCAMALLOC(sizeof(INT));
      }

      for (i=0; i<design->ndvol; i++)
      {
        design->dvol_fenode2[i] = (INT*)CCAMALLOC(sizeof(INT));
      }
    }



  for (i=0; i<genprob.numfld; i++)
    if (field[i].subdivide > 0)
    {
      global_subdivide(&(field[i]));

    }  /* if (field[i].subdivide > 0) */


  /* write warnings */
  dswarning(2,0);


#endif /* ifdef SUBDIV */


  /* tidy up */
  if (design->dnode_fenode2 != NULL)
  if (design->dnode_fenode2[0] != NULL)
  {
    for (i=0; i<design->ndnode; i++) CCAFREE(design->dnode_fenode2[i]);
    CCAFREE(design->dnode_fenode2);
    for (i=0; i<design->ndline; i++) CCAFREE(design->dline_fenode2[i]);
    CCAFREE(design->dline_fenode2);
    for (i=0; i<design->ndsurf; i++) CCAFREE(design->dsurf_fenode2[i]);
    CCAFREE(design->dsurf_fenode2);
    for (i=0; i<design->ndvol; i++) CCAFREE(design->dvol_fenode2[i]);
    CCAFREE(design->dvol_fenode2);
  }

  if (design->dnode_fenode != NULL)
  if (design->dnode_fenode[0] != NULL)
  {
    for (i=0; i<design->ndnode; i++) CCAFREE(design->dnode_fenode[i]);
    CCAFREE(design->dnode_fenode);
    for (i=0; i<design->ndline; i++) CCAFREE(design->dline_fenode[i]);
    CCAFREE(design->dline_fenode);
    for (i=0; i<design->ndsurf; i++) CCAFREE(design->dsurf_fenode[i]);
    CCAFREE(design->dsurf_fenode);
    for (i=0; i<design->ndvol; i++) CCAFREE(design->dvol_fenode[i]);
    CCAFREE(design->dvol_fenode);
  }

  if (design->ndnode_fenode != NULL)
  {
    CCAFREE(design->ndnode_fenode);
    CCAFREE(design->ndline_fenode);
    CCAFREE(design->ndsurf_fenode);
    CCAFREE(design->ndvol_fenode);
  }




  /*--------------------------------------- input of general dynamic data */
  if (genprob.timetyp==time_dynamic) inpctrdyn();
  /*---------------------------------------- input of general static data */
  else inpctrstat();
  /*----------------------------------------- input of eigensolution data */
  inpctreig();
  /*------------------------------------------------- input of conditions */
  /* dirichlet/coupling/neumann conditions are read from file to the
     DVOLS/DSURFS/DLINES/DNODES
     */
#ifdef PERF
  perf_begin(10);
#endif
  inp_conditions();
#ifdef PERF
  perf_end(10);
#endif
  /*----------- inherit the fsi coupling conditions inside the design
    condition is transformed in a dirichlet condition for fluid- and
    ale-field and into a neumann condition for structure-field           */
#ifdef D_FSI
  if (genprob.probtyp==prb_fsi) fsi_createfsicoup();
#endif
#ifdef D_SSI
  if (genprob.probtyp==prb_ssi) ssi_createssicoup();
#endif
#ifdef D_FLUID
  /*----------------- inherit the freesurface condition inside the design
    condition is transformed into a dirichlet condition for ale fiedl    */
  if ((genprob.probtyp==prb_fluid && genprob.numfld>=2) || genprob.probtyp==prb_fsi)
    fluid_createfreesurf();

  if (genprob.probtyp==prb_fluid && genprob.numfld>=2)
    fluid_createpseudofsi();

  /*------ inherit stabilisation condition from design to the elements ---*/
  for (i=0; i<genprob.numfld; i++)
  {
    if (field[i].fieldtyp == fluid)
      for (j=0; j<field[i].ndis; j++)
        inherit_design_ele(&(field[i].dis[j]));
  }
#endif
  /*----- inherit the dirichlet and coupling conditions inside the design */
  /* conditions are inherited 'downwards': DVOL->DSURF->DLINE->DNODE      */
  /* BUT: if a 'lower object already has its own condition, it does NOT   */
  /* inherit from above, it does inherit its own condition further down   */
#ifdef PERF
  perf_begin(11);
#endif


  inherit_dirich_coup_indesign();


  /* set pointers in the discretization to the design dirichlet conditions*/
  for (i=0; i<genprob.numfld; i++)
    for (j=0; j<field[i].ndis; j++)
      inherit_design_dis_dirichlet(&field[i], &(field[i].dis[j]));


  /*--- set pointers in the discretization to the design couple conditions*/
  for (i=0; i<genprob.numfld; i++)
    for (j=0; j<field[i].ndis; j++)
      inherit_design_dis_couple(&(field[i].dis[j]));


#ifdef D_FSI
  /*---- set pointers in the discretisations to the design fsi conditions */
  if (genprob.probtyp==prb_fluid || genprob.probtyp==prb_fsi)
  {
    for (i=0; i<genprob.numfld; i++)
      for (j=0; j<field[i].ndis; j++)
        inherit_design_dis_fsicouple(&(field[i].dis[j]));
  }
#endif


#ifdef D_SSI
  /*---- set pointers in the discretisations to the design ssi conditions */
  if (genprob.probtyp==prb_ssi)
  {
    for (i=0; i<genprob.numfld; i++)
      for (j=0; j<field[i].ndis; j++)
        inherit_design_dis_ssicouple(&(field[i].dis[j]));
  }
#endif


#ifdef D_FLUID
  /*------ set pointers in the discr. to the design freesurface condition */
  if (genprob.probtyp==prb_fluid || genprob.probtyp==prb_fsi)
  {
    for (i=0; i<genprob.numfld; i++)
      for (j=0; j<field[i].ndis; j++)
      {
        inherit_design_dis_freesurf(&(field[i].dis[j]));
        inherit_design_dis_slipdirich(&(field[i].dis[j]));
      }
  }
#endif


  /*------------------------ interpolate axishell conditions to the nodes */
#ifdef D_AXISHELL
  if (genprob.probtyp==prb_structure)
  {
    for (j=0; j<field[genprob.numsf].ndis; j++)
      interpolate_axishell_conds(&(field[genprob.numsf].dis[j]));
  }
#endif
  /*---------------------------- input of submesh for material-multiscale */
#ifdef D_MLSTRUCT
  if (genprob.multisc_struct == 1)  inp_submesh();
#endif /* D_MLSTRUCT */

#ifdef PERF
  perf_end(11);
#endif

  /*-------------------------------------------- input of monitoring data */
  inp_monitor();

#ifdef RESULTTEST
  /*---------------------------------------- input of result descriptions */
  inp_resultdescr();
#endif

  /*--------------------------------------------- all reading is over here*/
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of ntainp */
