/*!----------------------------------------------------------------------
\file
\brief input of mesh data

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/


#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_pro/fluid3pro_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
#include "../axishell/axishell.h"
#include "../beam3/beam3.h"
#include "../interf/interf.h"

#ifdef D_WALLGE
#include "../wallge/wallge.h"
#endif

#ifdef D_THERM2
#include "../therm2/therm2.h"
#endif

#ifdef D_THERM3
#include "../therm3/therm3.h"
#endif



/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | vector of numfld FIELDs, defined in global_control.c               |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*
  | Global variables for this file                        m.gee 11/00  |
 *----------------------------------------------------------------------*/
static ARRAY tmpnodes1;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;


/*----------------------------------------------------------------------*
  | prototypes for this file                                           |
 *----------------------------------------------------------------------*/
#ifdef D_SSI
void inp_struct_field_ssi(FIELD *masterfield, FIELD *slavefield);
#endif


/*-----------------------------------------------------------------------*/
/*!
  \brief read the fields

  Depending on the problem types the elemnts for one or several fields
  are read from the input file and crated

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inpfield()
{
  INT  i,j;
  INT  numnd;
  INT  nnode_total = 0;


#ifdef DEBUG
  dstrc_enter("inpfield");
#endif


  genprob.create_dis = 0;
  genprob.create_ale = 0;
  genprob.maxnode    = 0;
  genprob.nodeshift  = genprob.nnode;


  /* read node coordinates from file */
  /* nodes coordinates are temporarily stored in global variable
   * tmpnodes1, which used later on */
  inpnodes();


  /* read elements for all fields from field  */


  /* FSI typ of problem:
   * -------------------
   */
  if (genprob.probtyp == prb_fsi)
  {
    if (genprob.numfld!=3) dserror("numfld != 3 for FSI");

    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    inp_struct_field(&(field[genprob.numsf]));

    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));
    inp_fluid_field (&(field[genprob.numff]));

    if (genprob.create_ale == 0)
    {
      field[genprob.numaf].fieldtyp = ale;
      inpdis(&(field[genprob.numaf]));
      inp_ale_field  (&(field[genprob.numaf]));
    }

  }


  /* SSI type of problem:
   * -------------------
   */
#ifdef D_SSI
  if (genprob.probtyp == prb_ssi)
  {
    if (genprob.numfld!=2) dserror("numfld != 2 for SSI");

    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[0].fieldtyp = structure;
    field[1].fieldtyp = structure;
    field[0].ndis=1;
    field[1].ndis=1;
    inp_struct_field_ssi(&(field[0]),&(field[1]));
  }
#endif


  /* structure type of problem:
   * --------------------------
   */
  if (genprob.probtyp==prb_structure)
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    inp_struct_field(&(field[genprob.numsf]));
  }


  /* Optimisation type of problem:
   * -----------------------------
   */
  if (genprob.probtyp == prb_opt)
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for optimization problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[0].fieldtyp = structure;
    inpdis(&(field[0]));
    inp_struct_field(&(field[0]));
  }


  /* fluid type of problem:
   * ----------------------
   */
  if (genprob.probtyp==prb_fluid)
  {
    if (genprob.numfld==1) /* single field fluid problem  */
    {
      field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

      field[genprob.numff].fieldtyp = fluid;
      inpdis(&(field[genprob.numff]));
      inp_fluid_field (&(field[genprob.numff]));
    }
    else if (genprob.numfld==2) /* two field fluid problem (fluid+ale) */
    {
      field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

      field[genprob.numff].fieldtyp = fluid;
      inpdis(&(field[genprob.numff]));
      inp_fluid_field (&(field[genprob.numff]));

      field[genprob.numaf].fieldtyp = ale;
      inpdis(&(field[genprob.numaf]));
      inp_ale_field  (&(field[genprob.numaf]));
    }
    else dserror("NUMFLD>2 not allowed for Problemtype FLUID\n");
  }


  /* fluid type of problem:
   * ----------------------
   */
  if (genprob.probtyp==prb_fluid_pm)
  {
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));
    inp_fluid_field (&(field[genprob.numff]));
  }


  /* ale type of problem:
   * --------------------
   */
  if (genprob.probtyp==prb_ale)
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));
    inp_ale_field(&(field[genprob.numaf]));
  }


#ifdef D_TSI
  /* TSI typ of problem:
   * -------------------
   */
  if (genprob.probtyp == prb_tsi)
  {
    if (genprob.numfld != 2)
    {
      dserror("numfld != 2 for TSI");
    }

    /* allocate global field variable */
    field = (FIELD*) CCACALLOC(genprob.numfld, sizeof(FIELD));

    /* structure field */
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    inp_struct_field(&(field[genprob.numsf]));

    /* thermal field */
    field[genprob.numtf].fieldtyp = thermal;
    inpdis(&(field[genprob.numtf]));
    inp_therm_field(&(field[genprob.numtf]));

    /* automatic creation of thermal mesh */
    /* ==> needs to be accomplished */
  }
#endif



  /* copy nodes for the second discretisation  or the ale field */
  if (genprob.create_dis == 1 || genprob.create_ale == 1)
  {

    numnd = genprob.nodeshift;
    genprob.nnode = 2*genprob.nnode;
    genprob.maxnode = genprob.nnode;
    amredef(&tmpnodes1,(genprob.nnode),3,"DA");

    for (i=numnd; i<genprob.nnode; i++)
    {
      tmpnodes1.a.da[i][0] = tmpnodes1.a.da[i-numnd][0];
      tmpnodes1.a.da[i][1] = tmpnodes1.a.da[i-numnd][1];
      tmpnodes1.a.da[i][2] = tmpnodes1.a.da[i-numnd][2];
    }
  }


  /* assign the nodes to the fields */
  for (i=0; i<genprob.numfld; i++)
  {
    for (j=0;j<field[i].ndis;j++)
    {
      if (field[i].dis[j].disclass != dc_subdiv_calc)
        inp_assign_nodes(&(field[i].dis[j]));

      nnode_total += field[i].dis[j].numnp;
    }
  }

  amdel(&tmpnodes1);



  /* make element-node-element topology */
  genprob.nnode = nnode_total;

  for (i=0; i<genprob.numfld; i++)
    for (j=0;j<field[i].ndis;j++)
      if (field[i].dis[j].disclass != dc_subdiv_calc)
        inp_topology(&(field[i].dis[j]));


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inpfield */





/*-----------------------------------------------------------------------*/
/*!
  \brief sort nodes to the fields

  Count the number of nodes in the current discretization, allocate the
  necessary array and create the nodes, i.e. assign the ids and coordinates

  \param actdis      *DISCRET  (i) the current discretization

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inp_assign_nodes(
    DISCRET       *actdis)
{
  INT        i,j;
  INT        node_Id;
  INT        counter;
  INT        minusone = -1;
  ARRAY      nodeflag;
  ELEMENT   *actele;


#ifdef DEBUG
  dstrc_enter("inp_assign_nodes");
#endif


  /* allocate a temporary flag array */
  amdef("nodeflag",&nodeflag,genprob.nnode,1,"IV");
  aminit(&nodeflag,&minusone);


  /* set a flag to the node_id for each node in the field */
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    for (j=0; j<actele->numnp; j++)
    {
      node_Id = actele->lm[j];
      nodeflag.a.iv[node_Id]=node_Id;
    }
  }


  /* count the flags */
  counter=0;
  for (i=0; i<genprob.nnode; i++)
  {
    if (nodeflag.a.iv[i]!=-1) counter++;
  }
  actdis->numnp=counter;
  dsassert(actdis->numnp>0, "No nodes in discretization");


  /* Allocate the nodes to the field */
  actdis->node = (NODE*)CCACALLOC(counter,sizeof(NODE));


  /* assign the node Ids and coords to the NODE structure */
  counter=0;
  for (i=0; i<genprob.nnode; i++)
  {
    if (nodeflag.a.iv[i]!=-1)
    {
      actdis->node[counter].Id = nodeflag.a.iv[i];
      for (j=0; j<3; j++)
      {
        actdis->node[counter].x[j] = tmpnodes1.a.da[i][j];
      }

      if (genprob.maxnode < (actdis->node[counter].Id+1))
        genprob.maxnode = actdis->node[counter].Id + 1;

      counter++;
    }
  }

  amdel(&nodeflag);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_assign_nodes */





/*-----------------------------------------------------------------------*/
/*!
  \brief read discretization data

  Read how many discretization the current field will have

  \param actfield    *FIELD  (i) the currenct field

  \return void

  \author mn
  \date   08/04

 */
/*-----------------------------------------------------------------------*/
void inpdis(
    FIELD       *actfield)
{

  INT  ierr=0;


#ifdef DEBUG
  dstrc_enter("inpdis");
#endif


  /* set default values */
  actfield->ndis      = 1;
  actfield->subdivide = 0;



  /* read discretisation */
  if (frfind("--DISCRETISATION")==0) goto end_dis;
  frread();
  switch(actfield->fieldtyp)
  {
    case fluid:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("NUMFLUIDDIS", &(actfield->ndis),&ierr);
        frread();
      }
      break;


    case structure:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("NUMSTRUCDIS", &(actfield->ndis),&ierr);
        frread();
      }
      break;


    case ale:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("NUMALEDIS", &(actfield->ndis),&ierr);
        frread();
      }
      break;


#ifdef D_TSI
    case thermal:
      while (strncmp(allfiles.actplace, "------", 6) != 0)
      {
        frint("NUMTHERMDIS", &(actfield->ndis), &ierr);
        frread();
      }
      break;
#endif


    default:
      dserror("Unknown fieldtype");
      break;
  }

  frrewind();

  if (actfield->ndis > MAXDIS)
    dserror("Too many discretizations: Increase the value of MAXDIS!!");


end_dis:



  /* read subdivide */
  if (frfind("--SUBDIVIDE")==0) goto end_sub;
  frread();
  switch(actfield->fieldtyp)
  {
    case fluid:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("FLUID_SUBDIVIDE", &(actfield->subdivide),&ierr);
        frread();
      }
      break;


    case structure:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("STRUCT_SUBDIVIDE", &(actfield->subdivide),&ierr);
        frread();
      }
      break;


    case ale:
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        frint("ALE_SUBDIVIDE", &(actfield->subdivide),&ierr);
        frread();
      }
      break;


#ifdef D_TSI
    case thermal:
      dserror("SUBDIVIDE is not available for thermal field!");
      break;
#endif


    default:
      dserror("Unknown fieldtype");
      break;
  }
  frrewind();


  /* read output_dis */
  if (frfind("--SUBDIVIDE")==0) goto end_sub;
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frint("OUTPUT_DIS", &(ioflags.output_dis),&ierr);
    frread();
  }


  frrewind();



end_sub:

#ifndef SUBDIV
  if (actfield->subdivide > 0 )
    dserror("Element subdivision chosen but not compiled!!");
#endif

  if (actfield->subdivide > 0 && actfield->ndis != 2 )
    dserror("Two discretizations needed for element subdivision!!");


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inpdis */





/*-----------------------------------------------------------------------*/
/*!
  \brief input of node coords

  Read the coordinates of all nodes and store them in a temporary array

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inpnodes()
{

  INT  ierr=0;
  INT  counter;
  INT  tmp_nodeid;


#ifdef DEBUG
  dstrc_enter("inpnodes");
#endif


  /* allocate temporary array for coordinates of all nodes */
  amdef("tempnod1",&tmpnodes1,(genprob.nnode),3,"DA");


  /* read nodes */
  if (frfind("--NODE COORDS")==0) dserror("frfind: NODE COORDS is not in input file");
  frread();
  counter=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frint("NODE",&(tmp_nodeid),&ierr);
    if (ierr!=1)
      dserror("reading of nodes failed");
    if (tmp_nodeid-1 != counter)
      dserror("Reading of nodes failed: Nodes must be numbered consecutive!!");

    frdouble_n("COORD",&(tmpnodes1.a.da[counter][0]),3,&ierr);
    if (ierr!=1) dserror("reading of nodes failed");

    counter++;

    frread();
  }
  frrewind();


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of inpnodes */





/*-----------------------------------------------------------------------*/
/*!
  \brief input of structure field

  Create the structure field: allocate the discretizations, the required
  number of elements and then read and create the elements

  \param structfield    FIELD  (i) pointer to the structure field

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inp_struct_field(
    FIELD       *structfield)
{

  INT         ierr;
  INT         counter=0;
  INT         elenumber;
  char       *colpointer;


#ifdef DEBUG
  dstrc_enter("inp_struct_field");
#endif


  /* allocate discretizations */
  structfield->dis = (DISCRET*)CCACALLOC(structfield->ndis,sizeof(DISCRET));


  if (structfield->subdivide > 0)
  {
    structfield->dis[0].disclass = dc_subdiv_io;
    structfield->dis[1].disclass = dc_subdiv_calc;

    /* initialize array positions with -1 */
    memset(&structfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
    memset(&structfield->dis[1].ipos, 0xff, sizeof(ARRAY_POSITION));
  }
  else
  {
    structfield->dis[0].disclass = dc_normal;

    /* initialize array positions with -1 */
    memset(&structfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
  }


  /* count number of elements */
  if (frfind("--STRUCTURE ELEMENTS")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      counter++;
      frread();
    }
  }
  structfield->dis[0].numele = counter;


  /* allocate elements */
  structfield->dis[0].element =
    (ELEMENT*)CCACALLOC(structfield->dis[0].numele,sizeof(ELEMENT));


  /* read elements */
  if (frfind("--STRUCTURE ELEMENTS")==0) goto end;
  frread();
  counter=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    colpointer = allfiles.actplace;
    elenumber  = strtol(colpointer,&colpointer,10);
    structfield->dis[0].element[counter].Id = --elenumber;

    /* read the typ of element and call element readning function */


    /*   elementtyp is SHELL8:
     *   ---------------------
     */
    frchk("SHELL8",&ierr);
    if (ierr==1)
    {
#ifndef D_SHELL8
      dserror("SHELL8 needed but not defined in Makefile");
#endif
    }
#ifdef D_SHELL8
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_shell8;
      s8inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is SHELL9:
     *   ---------------------
     */
    frchk("SHELL9",&ierr);
    if (ierr==1)
    {
#ifndef D_SHELL9
      dserror("SHELL9 needed but not defined in Makefile");
#endif
    }
#ifdef D_SHELL9
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_shell9;
      s9inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is BRICK1:
     *   ---------------------
     */
    frchk("BRICK1",&ierr);
    if (ierr==1)
    {
#ifndef D_BRICK1
      dserror("BRICK1 needed but not defined in Makefile");
#endif
    }
#ifdef D_BRICK1
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_brick1;
      c1inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is WALL:
     *   -------------------
     */
    frchk("WALL",&ierr);
    if (ierr==1)
    {
#ifndef D_WALL1
      dserror("WALL1 needed but not defined in Makefile");
#endif
    }
#ifdef D_WALL1
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_wall1;
      w1inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is BEAM3:
     *   --------------------
     */
    frchk("BEAM3",&ierr);
    if (ierr==1)
    {
#ifndef D_BEAM3
      dserror("BEAM3 needed but not defined in Makefile");
#endif
    }
#ifdef D_BEAM3
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_beam3;
      b3inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is SAXI:
     *   -------------------
     */
    frchk("SAXI",&ierr);
    if (ierr==1)
    {
#ifndef D_AXISHELL
      dserror("AXISHELL needed but not defined in Makefile");
#endif
    }
#ifdef D_AXISHELL
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_axishell;
      saxi_inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is Interf:
     *   ---------------------
     */
    frchk("IF",&ierr);
    if (ierr==1)
    {
#ifndef D_INTERF
      dserror("INTERF needed but not defined in Makefile");
#endif
    }
#ifdef D_INTERF
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_interf;
      interf_inp(&(structfield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is Wallge:
     *   ---------------------
     */
    frchk("WGE",&ierr);
    if (ierr==1)
    {
#ifndef D_WALLGE
      dserror("WALLGE needed but not defined in Makefile");
#endif
    }
#ifdef D_WALLGE
    if (ierr==1)
    {
      structfield->dis[0].element[counter].eltyp=el_wallge;
      wge_inp(&(structfield->dis[0].element[counter]));
    }
#endif


    counter++;
    frread();

  }  /* while(strncmp(allfiles.actplace,"------",6)!=0) */

  frrewind();


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_struct_field */





#ifdef D_SSI
/*-----------------------------------------------------------------------*/
/*!
  \brief input of structure field for SSI

  Create the structure field for SSI: allocate the discretizations, the
  required number of elements and then read and create the elements

  \param structfield    FIELD  (i) pointer to the structure field

  \return void

  \author chfoe
  \date   07/04

 */
/*-----------------------------------------------------------------------*/
void inp_struct_field_ssi(
    FIELD       *masterfield,
    FIELD       *slavefield)
{

  INT        ierr,ierr_s,ierr_m;
  INT        slavecounter=0;
  INT        mastercounter=0;
  INT        counter=0;
  INT        elenumber;
  char      *colpointer;
  INT       *flag;
  ARRAY      flag_a;


#ifdef DEBUG
  dstrc_enter("inp_struct_field_ssi");
#endif


  /* allocate discretizations */
  if (masterfield->ndis>1)
    dserror("different discretisations not implemented yet for structural elements\n");
  if (slavefield->ndis>1)
    dserror("different discretisations not implemented yet for structural elements\n");
  masterfield->dis = (DISCRET*)CCACALLOC(masterfield->ndis,sizeof(DISCRET));
  slavefield->dis  = (DISCRET*)CCACALLOC(slavefield->ndis,sizeof(DISCRET));

  /* initialize array positions with -1 */
  memset(&masterfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
  memset(&slavefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));

  slavefield->dis[0].disclass  = dc_normal;
  masterfield->dis[0].disclass = dc_normal;


#ifndef D_WALL1
  dserror("WALL1 needed but not defined in Makefile");
#endif


  /* count number of elements */
  if (frfind("--STRUCTURE ELEMENTS")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      counter++;
      frread();
    }
  }

  flag = amdef("flag",&flag_a,counter,1,"IV");
  counter=0;

  if (frfind("--STRUCTURE ELEMENTS")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frchk("WALL",&ierr);
      if (ierr!=1)
        dserror("SSI only possible with wall elements!");
      frchk("Master",&ierr_m);
      frchk("Slave",&ierr_s);
      if (ierr_s==1)
      {
        slavecounter++;
        flag[counter]=1;
      }
      else if (ierr_m==1)
      {
        mastercounter++;
        flag[counter]=0;
      }
      else dserror("SSI_COUPTYP not possible for wall element!");
      counter++;
      frread();
    }
  }
  masterfield->dis[0].numele = mastercounter;
  slavefield->dis[0].numele  = slavecounter;


  /* allocate elements */
  masterfield->dis[0].element =
    (ELEMENT*)CCACALLOC(masterfield->dis[0].numele,sizeof(ELEMENT));
  slavefield->dis[0].element =
    (ELEMENT*)CCACALLOC(slavefield->dis[0].numele,sizeof(ELEMENT));


  /* read elements */
  if (frfind("--STRUCTURE ELEMENTS")==0) goto end;
  frread();
  counter=0;
  mastercounter=0;
  slavecounter=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    switch (flag[counter])
    {

      case 0: /*master field */
        colpointer = allfiles.actplace;
        elenumber  = strtol(colpointer,&colpointer,10);
        masterfield->dis[0].element[mastercounter].Id = --elenumber;


        /* read the typ of element and call element readning function */

        /*   elementtyp is WALL:
         *   -------------------
         */
        frchk("WALL",&ierr);
        if (ierr==1)
        {
#ifndef D_WALL1
          dserror("WALL1 needed but not defined in Makefile");
#endif
        }
#ifdef D_WALL1
        if (ierr==1)
        {
          masterfield->dis[0].element[mastercounter].eltyp=el_wall1;
          w1inp(&(masterfield->dis[0].element[mastercounter]));
        }
#endif
        mastercounter++;
        break;



      case 1:
        colpointer = allfiles.actplace;
        elenumber  = strtol(colpointer,&colpointer,10);
        slavefield->dis[0].element[slavecounter].Id = --elenumber;


        /* read the typ of element and call element readning function */

        /*   elementtyp is WALL:
         *   -------------------
         */
        frchk("WALL",&ierr);
        if (ierr==1)
        {
#ifndef D_WALL1
          dserror("WALL1 needed but not defined in Makefile");
#endif
        }
#ifdef D_WALL1
        if (ierr==1)
        {
          slavefield->dis[0].element[slavecounter].eltyp=el_wall1;
          w1inp(&(slavefield->dis[0].element[slavecounter]));
        }
#endif
        slavecounter++;
        break;


      default:
        dserror("flag out of range!\n");
    }

    counter++;
    frread();

  }  /* while(strncmp(allfiles.actplace,"------",6)!=0) */

  frrewind();


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_struct_field_ssi */

#endif /* #ifdef D_SSI */




/*-----------------------------------------------------------------------*/
/*!
  \brief input of fluid field

  Create the fluid field: allocate the discretizations, the required
  number of elements and then read and create the elements

  \param fluidfield    FIELD  (i) pointer to the fluid field

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inp_fluid_field(
    FIELD       *fluidfield)
{

  INT        ierr;
  INT        counter = 0;
  INT        ale_counter;
  INT        ale_element;
  INT        elenumber;
  char      *colpointer;
  FIELD     *alefield;
  char      *actplace_save;
  INT        actrow_save;

#if defined(D_FLUID2TU) || defined(D_FLUID2_PRO) || defined(D_FLUID3_PRO)
  INT        cpro = 0;
#endif


#ifdef DEBUG
  dstrc_enter("inp_fluid_field");
#endif

  alefield = &(field[genprob.numaf]);

  /* allocate discretizations */
  fluidfield->dis = (DISCRET*)CCACALLOC(fluidfield->ndis,sizeof(DISCRET));

  /* remarks about different discretisations:
   * ----------------------------------------
   * we asume to read in one "global" discretisation from the input file.
   * from this discretisation all the other ones can be directly derived!!!
   */

  if (fluidfield->subdivide > 0)
  {
    fluidfield->dis[0].disclass = dc_subdiv_io;
    fluidfield->dis[1].disclass = dc_subdiv_calc;

    /* initialize array positions with -1 */
    memset(&fluidfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
    memset(&fluidfield->dis[1].ipos, 0xff, sizeof(ARRAY_POSITION));
  }
  else
  {
    fluidfield->dis[0].disclass = dc_normal;

    /* initialize array positions with -1 */
    memset(&fluidfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
  }

  /* count number of elements */
  ale_counter = 0;
  if (frfind("--FLUID ELEMENTS")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      INT create_ale = 0;
      if (!frreadyes("CA", &create_ale))
        create_ale = 0;
      if (create_ale != 0)
        ale_counter++;

      counter++;
      frread();
    }
  }
  fluidfield->dis[0].numele = counter;


  /* allocate elements */
  fluidfield->dis[0].element =
    (ELEMENT*)CCACALLOC(fluidfield->dis[0].numele,sizeof(ELEMENT));


  /* read elements */
  if (frfind("--FLUID ELEMENTS")==0) goto end;
  frread();
  counter=0;
  ale_element=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    colpointer = allfiles.actplace;
    elenumber  = strtol(colpointer,&colpointer,10);
    fluidfield->dis[0].element[counter].Id = --elenumber;


    /* read the typ of element and call element readning function */


    /*   elementtyp is FLUID2_PRO
     *   ------------------------
     */
    frchk("FLUID2_PRO",&ierr);
    if (ierr==1)
    {
#ifndef D_FLUID2_PRO
      dserror("FLUID2_PRO needed but not defined in Makefile");
#else
      fluidfield->dis[0].element[counter].eltyp = el_fluid2_pro;
      f2pro_inp(&(fluidfield->dis[0].element[counter]));

      /*
       * We always create a second discretization, even if we do not
       * need it for this projection method. The particular kind of
       * method is not yet known... */
      if (cpro==0)
      {
	dsassert(fluidfield->subdivide==0,"no subdiv with pm at the moment");
	dsassert(fluidfield->ndis==1,"one discretization expected");
	fluidfield->ndis += 1;
	fluidfield->dis = (DISCRET*)CCAREALLOC(fluidfield->dis,fluidfield->ndis*sizeof(DISCRET));

	fluidfield->dis[1].numele  = fluidfield->dis[0].numele;
	fluidfield->dis[1].element = (ELEMENT*)CCACALLOC(fluidfield->dis[1].numele,sizeof(ELEMENT));
	cpro++;
	genprob.create_dis = 1;
	fluidfield->dis[1].disclass = dc_created_tu;
      }
      fluidfield->dis[1].element[counter].eltyp=el_fluid2_pro;
      f2pro_dis(&(fluidfield->dis[0].element[counter]),
		&(fluidfield->dis[1].element[counter]),
		genprob.nele,genprob.nodeshift);
#endif
    }


    /*   elementtyp is FLUID3_PRO
     *   ------------------------
     */
    frchk("FLUID3_PRO",&ierr);
    if (ierr==1)
    {
#ifndef D_FLUID3_PRO
      dserror("FLUID3_PRO needed but not defined in Makefile");
#else
      fluidfield->dis[0].element[counter].eltyp = el_fluid3_pro;
      f3pro_inp(&(fluidfield->dis[0].element[counter]));

      /*
       * We always create a second discretization, even if we do not
       * need it for this projection method. The particular kind of
       * method is not yet known... */
      if (cpro==0)
      {
	dsassert(fluidfield->subdivide==0,"no subdiv with pm at the moment");
	dsassert(fluidfield->ndis==1,"one discretization expected");
	fluidfield->ndis += 1;
	fluidfield->dis = (DISCRET*)CCAREALLOC(fluidfield->dis,fluidfield->ndis*sizeof(DISCRET));

	fluidfield->dis[1].numele  = fluidfield->dis[0].numele;
	fluidfield->dis[1].element = (ELEMENT*)CCACALLOC(fluidfield->dis[1].numele,sizeof(ELEMENT));
	cpro++;
	genprob.create_dis = 1;
	fluidfield->dis[1].disclass = dc_created_tu;
      }
      fluidfield->dis[1].element[counter].eltyp=el_fluid3_pro;
      f3pro_dis(&(fluidfield->dis[0].element[counter]),
		&(fluidfield->dis[1].element[counter]),
		genprob.nele,genprob.nodeshift);
#endif
    }


    /*   elementtyp is FLUID3
     *   --------------------
     */
    frchk("FLUID3 ",&ierr);
    if (ierr==1)
    {
#ifndef D_FLUID3
      dserror("FLUID3 needed but not defined in Makefile");
#endif
    }
#ifdef D_FLUID3
    if (ierr==1)
    {
      fluidfield->dis[0].element[counter].eltyp=el_fluid3;
      f3inp(&(fluidfield->dis[0].element[counter]),counter);

      if (fluidfield->dis[0].element[counter].e.f3->create_ale > 0)
      {
        /* allocate discretization, if necessary */
        if( alefield->dis == NULL)
        {
          dsassert(ale_counter>0, "implicit ale expected");
          alefield->fieldtyp = ale;
          actplace_save = allfiles.actplace;
          actrow_save   = allfiles.actrow;
          inpdis(alefield);
          allfiles.actplace = actplace_save;
          allfiles.actrow   = actrow_save;

          alefield->dis = (DISCRET*)CCACALLOC(alefield->ndis,sizeof(DISCRET));

          if (alefield->subdivide > 0)
          {
            alefield->dis[0].disclass = dc_subdiv_io_created_ale;
            alefield->dis[1].disclass = dc_subdiv_calc;

	    /* initialize array positions with -1 */
	    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
	    memset(&alefield->dis[1].ipos, 0xff, sizeof(ARRAY_POSITION));
          }
          else
	  {
            alefield->dis[0].disclass = dc_created_ale;

	    /* initialize array positions with -1 */
	    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
	  }
        }


        /* allocate elements */
        if( alefield->dis[0].element == NULL)
        {
          alefield->dis[0].numele = ale_counter;
          alefield->dis[0].element =
            (ELEMENT*)CCACALLOC(ale_counter,sizeof(ELEMENT));
        }


        /* create the corresponding ale element */
        f3_createale(&(fluidfield->dis[0].element[counter]),
            &(alefield->dis[0].element[ale_element]),genprob.nele,genprob.nodeshift);
        ale_element++;
      }  /* if (genprob.create_ale == 1) */
    }
#endif



    /*   elementtyp is FLUID3_FAST
     *   -------------------------
     */
    frchk("FLUID3_FAST",&ierr);
    if (ierr==1)
    {
#ifndef D_FLUID3_F
      dserror("FLUID3_F needed but not defined in Makefile");
#else
      fluidfield->dis[0].element[counter].eltyp=el_fluid3_fast;
      f3inp(&(fluidfield->dis[0].element[counter]),counter);

      if (genprob.create_ale == 1)
      {
        /* allocate discretization, if necessary */
        if( alefield->dis == NULL)
        {
          alefield->fieldtyp = ale;
          actplace_save = allfiles.actplace;
          actrow_save   = allfiles.actrow;
          inpdis(alefield);
          allfiles.actplace = actplace_save;
          allfiles.actrow   = actrow_save;

          alefield->dis = (DISCRET*)CCACALLOC(alefield->ndis,sizeof(DISCRET));

          if (alefield->subdivide > 0)
          {
            alefield->dis[0].disclass = dc_subdiv_io_created_ale;
            alefield->dis[1].disclass = dc_subdiv_calc;

	    /* initialize array positions with -1 */
	    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
	    memset(&alefield->dis[1].ipos, 0xff, sizeof(ARRAY_POSITION));
          }
          else
	  {
            alefield->dis[0].disclass = dc_created_ale;

	    /* initialize array positions with -1 */
	    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
	  }
        }


        /* allocate elements */
        if( alefield->dis[0].element == NULL)
        {
          alefield->dis[0].numele = fluidfield->dis[0].numele;
          alefield->dis[0].element =
            (ELEMENT*)CCACALLOC(alefield->dis[0].numele,sizeof(ELEMENT));
        }


        /* create the corresponding ale element */
        f3_createale(&(fluidfield->dis[0].element[counter]),
            &(alefield->dis[0].element[ale_element]),genprob.nele,genprob.nodeshift);
        ale_element++;
      }  /* if (genprob.create_ale == 1) */

#endif
    }



    /*   elementtyp is FLUID2
     *   --------------------
     */
    frchk("FLUID2 ",&ierr);
    if (ierr==1)
    {
#ifndef D_FLUID2
      dserror("FLUID2 needed but not defined in Makefile");
#endif
    }
#ifdef D_FLUID2
    if (ierr==1)
    {
      fluidfield->dis[0].element[counter].eltyp=el_fluid2;
      f2_inp(&(fluidfield->dis[0].element[counter]),counter);


      if (genprob.create_ale == 1)
      {
        /* allocate discretization, if necessary */
        if( alefield->dis == NULL)
        {
          alefield->fieldtyp = ale;
          actplace_save = allfiles.actplace;
          actrow_save   = allfiles.actrow;
          inpdis(alefield);
          allfiles.actplace = actplace_save;
          allfiles.actrow   = actrow_save;

          if (alefield->ndis>1)
            dserror("different discretisations not implemented yet for structural elements\n");
          alefield->dis = (DISCRET*)CCACALLOC(alefield->ndis,sizeof(DISCRET));
          alefield->dis[0].disclass = dc_created_ale;

	  /* initialize array positions with -1 */
	  memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
	}


        /* allocate elements */
        if( alefield->dis[0].element == NULL)
        {
          alefield->dis[0].numele = fluidfield->dis[0].numele;
          alefield->dis[0].element =
            (ELEMENT*)CCACALLOC(alefield->dis[0].numele,sizeof(ELEMENT));
        }


        /* create the corresponding ale element */
        f2_createale(&(fluidfield->dis[0].element[counter]),
            &(alefield->dis[0].element[counter]),genprob.nele,genprob.nodeshift);
      }  /* if (genprob.create_ale == 1) */



#ifdef D_FLUID2TU
      if (fluidfield->dis[0].element[counter].e.f2->turbu == 2 ||
          fluidfield->dis[0].element[counter].e.f2->turbu == 3 )
      {
        if(cpro==0)
        {
          fluidfield->dis[1].numele  = fluidfield->dis[0].numele;
          fluidfield->dis[1].element =
            (ELEMENT*)CCACALLOC(fluidfield->dis[1].numele,sizeof(ELEMENT));
          cpro++;
          genprob.create_dis = 1;
          fluidfield->dis[1].disclass = dc_created_tu;
        } /* endif (cpro==0) */

        fluidfield->dis[1].element[counter].eltyp=el_fluid2_tu;
        f2tu_dis(&(fluidfield->dis[0].element[counter]),
            &(fluidfield->dis[1].element[counter]),genprob.nele,genprob.nodeshift);
      } /* endif (e.f2->turbu == 2 || 3) */
#endif

    }
#endif


    counter++;
    frread();

  }  /* while(strncmp(allfiles.actplace,"------",6)!=0) */
  dsassert(ale_counter==ale_element, "ale element count mismatch");

  frrewind();

end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_fluid_field */






/*-----------------------------------------------------------------------*/
/*!
  \brief input of ale field

  Create the ale field: allocate the discretizations, the required
  number of elements and then read and create the elements

  \param alefield    FIELD  (i) pointer to the ale field

  \return void

  \author m.gee
  \date   04/01

 */
/*-----------------------------------------------------------------------*/
void inp_ale_field(
    FIELD       *alefield)
{

  INT  ierr;
  INT  counter=0;
  INT  elenumber;
  char *colpointer;


#ifdef DEBUG
  dstrc_enter("inp_ale_field");
#endif


  /* allocate discretizations */
  alefield->dis = (DISCRET*)CCACALLOC(alefield->ndis,sizeof(DISCRET));

  /* remarks about different discretisations:
   * ----------------------------------------
   * we asume to read in one "global" discretisation from the input file.
   * from this discretisation all the other ones can be directly derived!!!
   */


  if (alefield->subdivide > 0)
  {
    alefield->dis[0].disclass = dc_subdiv_io;
    alefield->dis[1].disclass = dc_subdiv_calc;

    /* initialize array positions with -1 */
    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
    memset(&alefield->dis[1].ipos, 0xff, sizeof(ARRAY_POSITION));
  }
  else
  {
    alefield->dis[0].disclass = dc_normal;

    /* initialize array positions with -1 */
    memset(&alefield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));
  }

  /* count number of elements */
  if (frfind("--ALE ELEMENTS")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      counter++;
      frread();
    }
  }
  frrewind();
  alefield->dis[0].numele = counter;


  /* allocate elements */
  alefield->dis[0].element=(ELEMENT*)CCACALLOC(alefield->dis[0].numele,sizeof(ELEMENT));


  /* read elements */
  frrewind();
  if (frfind("--ALE ELEMENTS")==0) goto end;
  frread();
  counter=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    colpointer = allfiles.actplace;
    elenumber  = strtol(colpointer,&colpointer,10);
    alefield->dis[0].element[counter].Id = --elenumber;


    /* read the typ of element and call element reading function */


    /*   elementtyp is ALE3:
     *   -------------------
     */
    frchk("ALE3",&ierr);
    if (ierr==1)
    {
#ifndef D_ALE
      dserror("ALE3 needed but not defined in Makefile");
#endif
    }
#ifdef D_ALE
    if (ierr==1)
    {
      alefield->dis[0].element[counter].eltyp=el_ale3;
      ale3inp(&(alefield->dis[0].element[counter]));
    }
#endif



    /*   elementtyp is ALE2:
     *   -------------------
     */
    frchk("ALE2",&ierr);
    if (ierr==1)
    {
#ifndef D_ALE
      dserror("ALE2 needed but not defined in Makefile");
#endif
    }
#ifdef D_ALE
    if (ierr==1)
    {
      alefield->dis[0].element[counter].eltyp=el_ale2;
      ale2inp(&(alefield->dis[0].element[counter]));
    }
#endif


    counter++;
    frread();

  }  /* while(strncmp(allfiles.actplace,"------",6)!=0) */

  frrewind();


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_ale_field */




/*-----------------------------------------------------------------------*/
/*!
  \brief  Input of thermal field

  Create the thermal field: Allocate the discretisations, the required
  number of elements and then read and create the elements

  \param  thermfield    FIELD  (i) pointer to the thermal field
  \return  void

  \author  bborn
  \date  03/06

 */
/*-----------------------------------------------------------------------*/
#ifdef D_TSI
void inp_therm_field(
  FIELD *thermfield)
{
  INT idis;  /* index of "current" discretisation */
  INT ierr;
  INT elecounter;  /* all element counter */
  INT elenumber;  /* index/ID of currently read element */
  char *colpointer;


#ifdef DEBUG
  dstrc_enter("inp_therm_field");
#endif


  /* allocate discretizations */
  thermfield->dis
    = (DISCRET*) CCACALLOC(thermfield->ndis, sizeof(DISCRET));

  /* several thermal discretisations */
  if (thermfield->ndis <= 0)
  {
    dserror("Thermal field does not have any discretisations!");
  }
  else if (thermfield->ndis > 1)
  {
    dserror("Thermal field can only deal with a single discretisation!");
  }
  else
  {
    idis = 0;
  }

  /* subdivison */
  if (thermfield->subdivide > 0)
  {
    dserror("Subdivision is not available for thermal discretisations!");
    /*
    structfield->dis[0].disclass = dc_subdiv_io;
    structfield->dis[1].disclass = dc_subdiv_calc;
    */
  }
  else
  {
    thermfield->dis[idis].disclass = dc_normal;
  }


  /* count number of elements */
  elecounter = 0;  /* initialise element counter */
  if (frfind("--THERMAL ELEMENTS") == 1)
  {
    frread();
    while (strncmp(allfiles.actplace,"------",6) != 0)
    {
      elecounter++;
      frread();
    }
  }
  thermfield->dis[idis].numele = elecounter;


  /* allocate elements */
  thermfield->dis[idis].element
    = (ELEMENT*) CCACALLOC(thermfield->dis[idis].numele, sizeof(ELEMENT));


  /* read elements */
  if (frfind("--THERMAL ELEMENTS")==0) goto end;
  frread();
  elecounter = 0;  /* initialise element counter */
  while (strncmp(allfiles.actplace,"------",6) != 0)
  {
    colpointer = allfiles.actplace;
    elenumber  = strtol(colpointer, &colpointer, 10);  /* ? */
    /* set ID */
    thermfield->dis[idis].element[elecounter].Id = --elenumber;

    /*------------------------------------------------------------------*/
    /* read the typ of element and call element readning function */


    /*------------------------------------------------------------------*/
    /* element type : THERM2 */
    frchk("THERM2", &ierr);
    if (ierr == 1)
    {
#ifndef D_THERM2
      dserror("THERM2 needed but not defined in Makefile");
#else
      thermfield->dis[idis].element[elecounter].eltyp = el_therm2;
      th2_inp(&(thermfield->dis[idis].element[elecounter]));
#endif
    }

    /*------------------------------------------------------------------*/
    /*  element type : THERM3 */
    frchk("THERM3", &ierr);
    if (ierr == 1)
    {
#ifndef D_THERM3
      dserror("THERM3 needed but not defined in Makefile");
#else
      thermfield->dis[idis].element[elecounter].eltyp = el_therm3;
      th3_inp(&(thermfield->dis[idis].element[elecounter]));
#endif
    }


    elecounter++;  /* increment element counter */
    frread();

  }  /* while(strncmp(allfiles.actplace,"------",6)!=0) */

  frrewind();


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of inp_therm_field */
#endif  /* end of #ifdef D_TSI */
