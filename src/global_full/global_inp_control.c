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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief the tracing variable

<pre>                                                         m.gee 8/00
defined in pss_ds.c, declared in tracing.h
</pre>
*----------------------------------------------------------------------*/
#ifdef DEBUG
extern struct _TRACE         trace;
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/* global variable: flag for the creation of a second discretisation */
INT      create_dis;


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
INT i,j,k,id;
INT  ngnode;
INT  ngline;
INT  ngsurf;
INT  ngvol;
INT  counter1, counter2, counter3;

/* 
   the input of the tracing option has not been done yet, so
   we have to make the dstrc_enter 'by hand' 
*/
#ifdef DEBUG 
trace.actroutine = trace.actroutine->next;
strncpy(trace.actroutine->name,"ntainp",49);
trace.actroutine->dsroutcontrol=dsin;
trace.deepness++;
#endif

/*--------------------------------------------- input of tracing option */
#ifdef DEBUG 
inptrace();
#endif
/*=========================== tracing is active and working from now on */
/*----------------------- input of not mesh or time based problem data  */
inpctr();
#ifdef D_OPTIM
/*------------------------------------------------input of optimization */
if (genprob.probtyp == prb_opt) inpctropt();
#endif
/*------------------------------------------------ input of design data */
/* read description of design volumes, surfaces, lines and nodes, allocate
   and fill the structures DVOL, DSURF, DLINE and DNODE 
*/
inpdesign();
/* build the topology among the design elements, allocate and create the
   pointers connecting DVOL, DSURF, DLINE and DNODE 
*/
inpdesign_topology_design();
/*----------------------------------------------------------------------*/
/* NOTE: the materials have to be read before the input of the elements */
/*       as these informations are needed for shell9 element            */
/*                                                             sh 10/02 */
/*                                                                      */
/*-------------------------------------------------- input of materials */
inp_material();
/*---------------- sh 10/02 --- input of multilayer materials -> shell9 */
inp_multimat();
/*------------------------------------------------------input of meshes */
/* read the fe-nodes in NODE, the fe-elements in ELEMENT and build the 
   topology among ELEMENTs and NODEs. Assign the disctretization to
   the different fields in a multifield problem dependent on type
   of element
*/    
inpfield();
/*------------------ built the detailed topology of the discretizations */
/* for each existing field and discretization build a detailed topology 
   of GNODEs, GLINEs, GSURFs and GVOLs and connect them to the topology
   of NODEs and ELEMENTs
*/   
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inp_detailed_topology(&(field[i].dis[j]));
/*---------------------------------------------------design-fe topology */
/* Read which node is on which design object from file.
   For each field and discretization build the pointers among design
   and FE-objects DVOL<->GVOL,DSURF<->GSURF,DLINE<->GLINE,DNODE<->GNODE.
*/   
  /* count number of gnodes, glines, gsurfs and gvols in whole problem */
  ngnode =0;
  ngline =0;
  ngsurf =0;
  ngvol  =0;
  for (i=0; i<genprob.numfld; i++)
  {
    for (j=0;j<field[i].ndis;j++)
    {
      ngnode += field[i].dis[j].ngnode;
      ngline += field[i].dis[j].ngline;
      ngsurf += field[i].dis[j].ngsurf;
      ngvol  += field[i].dis[j].ngvol;
    }
  }
  genprob.ngnode = ngnode;
  genprob.gnodes = (GNODE**)CCACALLOC(genprob.maxnode,sizeof(GNODE*));
  genprob.ngline = ngline;
  genprob.glines = (GLINE**)CCACALLOC(ngline,sizeof(GLINE*));
  genprob.ngsurf = ngsurf;
  genprob.gsurfs = (GSURF**)CCACALLOC(ngsurf,sizeof(GSURF*));
  genprob.ngvol  = ngvol;
  genprob.gvols  = (GVOL**)CCACALLOC(ngvol,sizeof(GVOL*));

  counter1 = 0;
  counter2 = 0;
  counter3 = 0;
  for (i=0; i<genprob.numfld; i++)
  {
    for (j=0;j<field[i].ndis;j++)
    {
      /* make pointers to all gnodes in genprob.gnodes */
      for (k=0; k<field[i].dis[j].ngnode; k++)
      {
        id = field[i].dis[j].gnode[k].node->Id;
        dsassert(id <= genprob.maxnode,"Zu wenig KNOTEN");
        genprob.gnodes[id] = &(field[i].dis[j].gnode[k]);
      }
      for (k=0; k<field[i].dis[j].ngline; k++)
        genprob.glines[counter1++] = &(field[i].dis[j].gline[k]);
      for (k=0; k<field[i].dis[j].ngsurf; k++)
        genprob.gsurfs[counter2++] = &(field[i].dis[j].gsurf[k]);
      for (k=0; k<field[i].dis[j].ngvol; k++)
        genprob.gvols[counter3++] = &(field[i].dis[j].gvol[k]);
    }
  }

  inpdesign_topology_fe();

  CCAFREE(genprob.nodes);
  CCAFREE(genprob.gnodes);
  CCAFREE(genprob.glines);
  CCAFREE(genprob.gsurfs);
  CCAFREE(genprob.gvols);

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
inp_conditions();
/*----------- inherit the fsi coupling conditions inside the design 
   condition is transformed in a dirichlet condition for fluid- and
   ale-field and into a neumann condition for structure-field           */
#ifdef D_FSI
if (genprob.probtyp==prb_fsi) fsi_createfsicoup();  
#endif
/*----------------- inherit the freesurface condition inside the design 
   condition is transformed into a dirichlet condition for ale fiedl    */
#ifdef D_FLUID
if (genprob.numff>=0 && genprob.numfld>1) fluid_createfreesurf();
#endif   
/*----- inherit the dirichlet and coupling conditions inside the design */
/* conditions are inherited 'downwards': DVOL->DSURF->DLINE->DNODE      */
/* BUT: if a 'lower object already has its own condition, it does NOT   */
/* inherit from above, it does inherit its own condition further down   */
inherit_dirich_coup_indesign();
/* set pointers in the discretization to the design dirichlet conditions*/
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inherit_design_dis_dirichlet(&(field[i].dis[j]));
/*--- set pointers in the discretization to the design couple conditions*/
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inherit_design_dis_couple(&(field[i].dis[j]));
#ifdef D_FSI
/*--------------------------------------------- do we really need this? 
  I don't think so - check it!!!                                        */
if (genprob.probtyp==prb_fsi) 
{   
   for (i=0; i<genprob.numfld; i++)
   for (j=0; j<field[i].ndis; j++)
   inherit_design_dis_fsicouple(&(field[i].dis[j]));  
}
#endif
/*------ set pointers n the discretisation to the freesurface condition */
#ifdef D_FLUID
if (genprob.probtyp==prb_fluid || genprob.probtyp==prb_fsi)
{
  for (i=0; i<genprob.numfld; i++)
    for (j=0; j<field[i].ndis; j++)
      inherit_design_dis_freesurf(&(field[i].dis[j]));
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
/*-------------------------------------------- input of monitoring data */
inp_monitor();

/*--------------------------------------------- all reading is over here*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ntainp */
