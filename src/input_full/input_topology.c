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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 | prototypes of routines which are only to be called inside this file  |
 |                                                        m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inptop_findadjele(
    ELEMENT       *centerele,
    ELEMENT       *elepatch[400],
    INT            nelepatch,
    ELEMENT       *adjele[400],
    INT            adjelelinenum[400],
    INT           *nadjele,
    INT            linenodes[3]);

static void inptop_makelinestoele(
    ELEMENT       *actele,
    INT            linenodes[12][3]);

static void inptop_findotherele(
    ELEMENT       *firstele,
    ELEMENT      **otherele,
    INT           *facenumber,
    ELEMENT       *elepatch[400],
    INT            npatch,
    INT            surfnodes[4]);

static void inptop_makepatch(
    ELEMENT       *centerele,
    ELEMENT       *elepatch[400],
    INT           *nelepatch);

static void inptop_makesurfnodes(
    ELEMENT       *actele,
    INT            surfnodes[6][4]);




/*----------------------------------------------------------------------*
 | create the node-element topology for this field        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_topology(
    DISCRET       *actdis)
{

  INT  i,j,k;
  INT  node_id;
  ELEMENT *actele;
  NODE    *actnode;

#ifdef DEBUG
  dstrc_enter("inp_topology");
#endif
  /* create pointer from elements to nodes */

  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    actele->node = (NODE**)CCACALLOC(actele->numnp,sizeof(NODE*));

    for (j=0; j<actele->numnp; j++)
    {
      node_id = actele->lm[j];
      actele->node[j] = genprob.nodes[node_id];
    } /* end of loop over elements nodes */
  }/* end of loop over elements */

/*------------------------------ create pointers from nodes to elements */
for (i=0; i<actdis->numnp; i++) actdis->node[i].numele=0;
/*---------------------------- count the number of elements to one node */
for (i=0; i<actdis->numele; i++)
{
   actele = &(actdis->element[i]);
   for (j=0; j<actele->numnp; j++)
   {
      (actele->node[j]->numele)++;
   }
}
/*------------------------- allocate space for element pointers in NODE */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   actnode->element = (ELEMENT**)CCACALLOC(actnode->numele,sizeof(ELEMENT*));
   for (j=0; j<actnode->numele; j++) actnode->element[j]=NULL;
}
/*---------------- loop elements and point from their nodes to themself */
for (i=0; i<actdis->numele; i++)
{
   actele = &(actdis->element[i]);
   for (j=0; j<actele->numnp; j++)
   {
      actnode = actele->node[j];
      for (k=0; k<actnode->numele; k++)
      {
         if (actnode->element[k]==NULL) break;
      }
      actnode->element[k]=actele;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_topology */






/*-----------------------------------------------------------------------*
 | create the node-line-surface-volume topology            m.gee 3/02    |
 | This routine is rather complex, so here is the order in which         |
 | the whole topology is built:                                          |
 |                                                                       |
 | [[1]]   Allocate GNODEs                                               |
 |         make pointer GNODE <-> NODE                                   |
 |                                                                       |
 | [[2]]   Count number of GVOLS                                         |
 |         Allocate GVOLS                                                |
 |                                                                       |
 | [[3]]   make pointer GVOL <-> ELEMENT                                 |
 |         make GVOL.ngsurf                                              |
 |         Allocate GVOL.gsurf                                           |
 |                                                                       |
 | [[4]]   count number of 2D elements                                   |
 |         make upper estimate for DISCRET.ngsurf                        |
 |         Allocate GSURF                                                |
 |                                                                       |
 | [[5]]   make patches around GVOL                                      |
 |         make GSURF.ngvol                                              |
 |         allocate GSURF.gvol                                           |
 |         make GVOL <-> GSURF                                           |
 |                                                                       |
 | [[6]]   make GSURF <-> ELEMENT for 2D elements                        |
 |                                                                       |
 | [[7]]   realloc GSURFs, shift pointers to GSURFs if necessary         |
 |                                                                       |
 | [[8]]   make GSURF.ngline for 2D elements                             |
 |         make GVOL.ngline for 3D elements                              |
 |         allocate GSURF.gline for 2D elements                          |
 |         allocate GVOL.gline for 3D elements                           |
 | [[9]]   counter upper estimate for DISCRET.nglines                    |
 |         allocate GLINES                                               |
 |                                                                       |
 | [[9a]]  count number of 1D elements                          mn 05/03 |
 |                                                                       |
 | [[10]]  make GSURF -> GLINE for 2D elements                           |
 |         make GVOL -> GLINE for 3D elements                            |
 |                                                                       |
 | [[10a]] make GLINE <-> ELEMENT for 1D elements               mn 05/03 |
 |                                                                       |
 | [[11]]  realloc GLINEs , move pointers if necessary                   |
 |                                                                       |
 | [[12]]  built GSURF -> GLINE for 3D elements                          |
 |                                                                       |
 | [[13]]  count GLINE.ngsurf                                            |
 |         allocate GLINE.gsurf                                          |
 |         make GLINE -> GSURF                                           |
 |                                                                       |
 | [[14]]  count GLINE.ngnode                                            |
 |         allocate GLINE.gnode                                          |
 |                                                                       |
 | [[15]]  make GLINE -> GNODE                                           |
 |         make GNODE -> GLINE                                           |
 |         make GSURF -> GNODE                                           |
 *-----------------------------------------------------------------------*/
void inp_detailed_topology(DISCRET   *actdis)
{
  PTRSIZE    ptrdistance;
  INT        i,j,k,counter;
  INT        nsurfelement;
  INT        nsurftovol;
  INT        nlineelement;
  INT        nlinetosurf;
  INT        nline = 0;
  INT        nnodeperline = 0;
  INT        surfnodes[6][4];
  INT        linenodes[12][3];
  ELEMENT   *actele = NULL;
  ELEMENT   *otherele;
  GVOL      *actgvol = NULL;
  GSURF     *actgsurf = NULL;
  GSURF     *tmpgsurf;
  GLINE     *tmpgline;
  GLINE    **gline = NULL;
  GLINE     *actgline = NULL;
  GNODE     *actgnode;
  ELEMENT   *elepatch[400];
  INT        nelepatch;
  ELEMENT   *adjelepatch[400];
  INT        adjelelinenum[400];
  INT        nadjelepatch;

  INT        facenumber;
  INT        isgline;
  INT        isgsurf;
  INT        isgvol;
#ifdef DEBUG
  dstrc_enter("inp_detailed_topology");
#endif
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  /*                                                          make gnodes */
  /*----------------------------------------------------------------------*/
  /*---------------------------- set number of gnodes and allocate memory */
  /*
     [[1]]
     */
  actdis->ngnode = actdis->numnp;
  actdis->gnode = (GNODE*)CCACALLOC(actdis->ngnode,sizeof(GNODE));
#ifdef DEBUG
  for (i=0; i<actdis->ngnode; i++) actdis->gnode[i].Id = i;
#endif
  /*----------------------------------------- set pointers gnode <-> node */
  for (i=0; i<actdis->numnp; i++)
  {
    actdis->node[i].gnode = &(actdis->gnode[i]);
    actdis->gnode[i].node = &(actdis->node[i]);
  }

  /*----------------------------------------------------------------------*/
  /*                                                           make gvols */
  /*----------------------------------------------------------------------*/
  /*------------- count number of volumetric elements and allocate memory */
  /*
     [[2]]
     */
  counter=0;
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case hex8:
      case hex20:
      case hex27:
      case tet4:
      case tet10:
        counter++;
        break;
      default:
        break;
    }
  }
  actdis->ngvol=counter;
  actdis->gvol = (GVOL*)CCACALLOC(actdis->ngvol,sizeof(GVOL));
#ifdef DEBUG
  for (i=0; i<actdis->ngvol; i++) actdis->gvol[i].Id = i;
#endif
  /*----------------------------------------------------------------------*/
  /*                                                          make gsurfs */
  /*----------------------------------------------------------------------*/
  /* surfaces are made in this steps:
     - count number of surface to volumetric elements and allocate pointers
     - count number of surfaces to surface elements
     - allocate and upper guess number of GURFs
     - make topology GVOLs <-> GSURFs
     - make topology ELEMENTs <-> GSURFs for 2D elements
     - realloc the total number of surfaces
     */
  /*------------------------- set pointers gvol <-> (volumetric) elements */
  /*                                  and allocate pointers for gsurfaces */
  /*
     [[3]]
     */
  counter=0;
  nsurftovol=0;
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case hex8:
      case hex20:
      case hex27:
        actdis->gvol[counter].element = actele;
        actele->g.gvol                = &(actdis->gvol[counter]);
        actdis->gvol[counter].ngsurf  = 6;
        actdis->gvol[counter].gsurf   = (GSURF**)CCACALLOC(actdis->gvol[counter].ngsurf,sizeof(GSURF*));
        nsurftovol += actdis->gvol[counter].ngsurf;
        counter++;
        break;
      case tet4:
      case tet10:
        actdis->gvol[counter].element = actele;
        actele->g.gvol                = &(actdis->gvol[counter]);
        actdis->gvol[counter].ngsurf  = 4;
        actdis->gvol[counter].gsurf   = (GSURF**)CCACALLOC(actdis->gvol[counter].ngsurf,sizeof(GSURF*));
        nsurftovol += actdis->gvol[counter].ngsurf;
        counter++;
        break;
      default:
        break;
    }
  }

  /*------------------------------------ count number of surface elements */
  /*
     [[4]]
     */
  nsurfelement=0;
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case quad4:
      case quad8:
      case quad9:
      case tri3:
      case tri6:
        nsurfelement++;
        break;
      default:
        break;
    }
  }
  /* an upper estimate to the number of surfaces is nsurfelement+nsurftovol */
  /* allocate (temporarily) some surfaces */
  actdis->ngsurf = nsurftovol+nsurfelement;
  actdis->gsurf = (GSURF*)CCACALLOC(actdis->ngsurf,sizeof(GSURF));
#ifdef DEBUG
  for (i=0; i<actdis->ngsurf; i++) actdis->gsurf[i].Id = i;
#endif
  /*------------------------------------ now loop all volumetric elements */
  /*
     [[5]]
     */
  counter=0;
  for (i=0; i<actdis->ngvol; i++)
  {
    /*------------------------------ set actual element and actual gvol */
    actgvol = &(actdis->gvol[i]);
    actele = actgvol->element;
    /*------------------- built a patch of elements around this element */
    /* elepatch holds pointers of on all elements on patch and nelepatch
       is the number of elements on patch. Note that actele is also included
       in the patch
       */
    inptop_makepatch(actele,elepatch,&nelepatch);
    /* make the array surfnodes containing the nodes of each surface to actele */
    /* the array surfnodes[6][4] hold the maximal 4 nodes to each of maximal
       6 surfaces for the element actele
       */
    inptop_makesurfnodes(actele,surfnodes);
    /*----------------------------------------- loop surfaces of actele */
    for (j=0; j<actgvol->ngsurf; j++)
    {
      /*----------- this surface of element has been recognized before */
      if (actgvol->gsurf[j] != NULL) continue;
      /* find the element, which shares the actual surface (if it exists) */
      /* function inptop_findotherele returns elementpointer otherele if there
         is an element 'on the other side of gsurf[j], it returns otherele=NULL
         if there is no other element to this surface.
         The function returns facenumber=the local surface number of the otherele
         */
      inptop_findotherele(actele,&otherele,&facenumber,elepatch,nelepatch,surfnodes[j]);
      /*------------------ if there is another element to this surface */
      if (otherele)
      {
        /*-------------------------- this is a surface to the actgvol */
        /*--------------------- set the pointer in the actual actgvol */
        actgvol->gsurf[j] = &(actdis->gsurf[counter]);
        /*-------------------------------------------- check otherele */
        dsassert(otherele->g.gvol->gsurf[facenumber]==NULL,"surfaces of volumetric elements got mixed up");
        dsassert(otherele->g.gvol->ngsurf > facenumber,"surfaces of volumetric elements got mixed up");
        /*----------------- this is a surface to the otherele as well */
        /*--------------------------- set the pointer in the otherele */
        otherele->g.gvol->gsurf[facenumber] = &(actdis->gsurf[counter]);
        /*-------------- so the gsurf has 2 volumetric elements to it */
        /* set number of volumes in the surface to 2, allocate 2 pointers
           and point to actgvol and otherele.g.gvol */
        dsassert(actdis->gsurf[counter].ngvol==0,"surfaces of volumetric elements got mixed up");
        actdis->gsurf[counter].ngvol = 2;
        actdis->gsurf[counter].gvol = (GVOL**)CCACALLOC(actdis->gsurf[counter].ngvol,sizeof(GVOL*));
        actdis->gsurf[counter].gvol[0] = actgvol;
        actdis->gsurf[counter].gvol[1] = otherele->g.gvol;
        /* increment gsurf, check whether we still have space in the initial guess */
        counter++;
        dsassert(counter<actdis->ngsurf,"Initial guess of ngsurf too small");
      }
      /*---------------- there is no otherele to this surface gsurf[j] */
      else
      {
        /*-------------------------- this is a surface to the actgvol */
        /*--------------------- set the pointer in the actual actgvol */
        actgvol->gsurf[j] = &(actdis->gsurf[counter]);
        /*--------------- so the gsurf has 1 volumetric element to it */
        /* set number of volumes in the surface to 1, allocate 1 pointer
           and point to actgvol */
        dsassert(actdis->gsurf[counter].ngvol==0,"surfaces of volumetric elements got mixed up");
        actdis->gsurf[counter].ngvol = 1;
        actdis->gsurf[counter].gvol = (GVOL**)CCACALLOC(actdis->gsurf[counter].ngvol,sizeof(GVOL*));
        actdis->gsurf[counter].gvol[0] = actgvol;
        /* increment gsurf, check whether we still have space in the initial guess */
        dsassert(counter<actdis->ngsurf,"Initial guess of ngsurf too small");
        counter++;
      }
    }/* end loop surfaces of volumetric elements */
  } /* end loop over volumetric elements */
  nsurftovol=counter;
  /*-------- now loop the 2D elements and make pointers gsurf <-> element */
  /*
     [[6]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case quad4:
      case quad8:
      case quad9:
      case tri3:
      case tri6:
        dsassert(actdis->gsurf[counter].ngvol==0,"surfaces of 2D elements got mixed up with 3D elements");
        actdis->gsurf[counter].element = actele;
        actele->g.gsurf = &(actdis->gsurf[counter]);
        counter++;
        dsassert(counter<=actdis->ngsurf,"Initial guess of ngsurf too small");
        break;
      default:
        break;
    }
  }
  dsassert(counter-nsurftovol==nsurfelement,"surfaces of 2D elements got mixed up with 3D elements");
  /*--------- surface topology is complete here, redefine the gsurf array */
  /* since the new array is smaller than the initial guess, I expect to keep
     the same pointer in an REALLOC, but this is checked here and in case it is not
     the case all pointers to gsurfs are moved
     */
  /*
     [[7]
     */
  if (actdis->ngsurf != 0 || counter != 0)
  {
    actdis->ngsurf = counter;
    tmpgsurf       = actdis->gsurf;
    actdis->gsurf  = (GSURF*)CCAREALLOC(actdis->gsurf,(actdis->ngsurf)*sizeof(GSURF));
    /*---- if the realloc moved the gsurf to another place, adjust pointers */
    if (actdis->gsurf != tmpgsurf)
    {
      /*--------------- find the pointer offset between new and old gsurfs */
      ptrdistance = (PTRSIZE)(actdis->gsurf) - (PTRSIZE)(tmpgsurf);
      /*---------------------------- loop elements and adjust the pointers */
      for (i=0; i<actdis->numele; i++)
      {
        actele = &(actdis->element[i]);
        switch (actele->distyp)
        {
          case line2:
          case line3:
            break;
          case quad4:
          case quad8:
          case quad9:
          case tri3:
          case tri6:
            ShiftPointer((void**)&(actele->g.gsurf),ptrdistance);
            break;
          case hex8:
          case hex20:
          case hex27:
          case tet4:
          case tet10:
            actgvol = actele->g.gvol;
            for (j=0; j<actgvol->ngsurf; j++)
              ShiftPointer((void**)&(actgvol->gsurf[j]),ptrdistance);
            break;
          default:
            dserror("Unknown type of discretization 1");
            break;
        }
      }/* end i loop over elements */
    }/* end of if (actdis->gsurf != tmpgsurf) */
  }

  /*----------------------------------------------------------------------*/
  /*                                                          make glines */
  /*----------------------------------------------------------------------*/
  /*----------------------------------now start making the lines topology */
  /*------------------------------------ prepare the surfaces and volumes */
  /*
     set number of lines in the gsurfs
     set number of lines in the gvols
     allocate pointers to lines in the gsurfs and in the gvols
     */
  /*
     [[8]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case line2:
      case line3:
        break;
      case quad4:
      case quad8:
      case quad9:
        actgsurf         = actele->g.gsurf;
        actgsurf->ngline = 4;
        actgsurf->gline = (GLINE**)CCACALLOC(actgsurf->ngline,sizeof(GLINE*));
        break;
      case tri3:
      case tri6:
        actgsurf         = actele->g.gsurf;
        actgsurf->ngline = 3;
        actgsurf->gline = (GLINE**)CCACALLOC(actgsurf->ngline,sizeof(GLINE*));
        break;
      case hex8:
      case hex20:
      case hex27:
        actgvol  = actele->g.gvol;
        actgvol->ngline = 12;
        actgvol->gline = (GLINE**)CCACALLOC(actgvol->ngline,sizeof(GLINE*));
        for (j=0; j<actgvol->ngsurf; j++)
        {
          actgvol->gsurf[j]->ngline = 4;
          if (!actgvol->gsurf[j]->gline)
            actgvol->gsurf[j]->gline = (GLINE**)CCACALLOC(actgvol->gsurf[j]->ngline,sizeof(GLINE*));
        }
        break;
      case tet4:
      case tet10:
        actgvol  = actele->g.gvol;
        actgvol->ngline = 6;
        actgvol->gline = (GLINE**)CCACALLOC(actgvol->ngline,sizeof(GLINE*));
        for (j=0; j<actgvol->ngsurf; j++)
        {
          actgvol->gsurf[j]->ngline = 3;
          if (!actgvol->gsurf[j]->gline)
            actgvol->gsurf[j]->gline = (GLINE**)CCACALLOC(actgvol->gsurf[j]->ngline,sizeof(GLINE*));
        }
        break;
      default:
        dserror("Unknown type of discretization 2: %s",actele->distyp);
        break;
    }
  }
  /*------ loop all gsurfs and count how many lines there could be maximal */
  /*
     [[9]]
     */
  counter=0;
  nlinetosurf = 0;
  for (i=0; i<actdis->ngsurf; i++)
  {
    nlinetosurf += actdis->gsurf[i].ngline;
    dsassert(actdis->gsurf[i].ngline!=0,"lines got mixed up ");
  }
  /*------------------------------------ count number of line elements */
  /*
     [[9a]]
     */
  nlineelement=0;
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case line2:
      case line3:
        nlineelement++;
        break;
      default:
        break;
    }
  }
  /* -an upper estimate to the number of lines is nlineelement+nlinetosurf */
  /*-------------------------------------------- allocate vector of glines */
  actdis->ngline = nlinetosurf+nlineelement;
  actdis->gline = (GLINE*)CCACALLOC(actdis->ngline,sizeof(GLINE));
#ifdef DEBUG
  for (i=0; i<actdis->ngline; i++) actdis->gline[i].Id = i;
#endif
  /*---------------------------------------------------- loop all elements */
  /*
     [[10]]
     */
  counter=0;
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    /*------------------------------------------ make patch around actele */
    /* elepatch holds pointers of on all elements on patch and nelepatch
       is the number of elements on patch. Note that actele is also included
       in the patch
       */
    inptop_makepatch(actele,elepatch,&nelepatch);
    /*--------------------------------- build list of lines to an element */
    /* linenodes[12][3] holds the global node ids of the nodes on each line
       to the element actele.
       Dependent on the element type there are linenodes[numlines][numpointsonline]
       On 2-noded lines, the nodes are in linenodes[line][0 and 2]
       On 3-noded lines, the nodes are in linenodes[line][0 and 1 and 2]
       */
    inptop_makelinestoele(actele,linenodes);
    /*------------------------------------------------------ check distyp */
    /*-------------- get number of lines to element and its gline pointer */
    switch (actele->distyp)
    {
      case line2:
      case line3:
        nline = 0;
        /*gline = actele->g.gline;*/
        break;
      case quad4:
      case quad8:
      case quad9:
      case tri3:
      case tri6:
        nline = actele->g.gsurf->ngline;
        gline = actele->g.gsurf->gline;
        break;
      case hex8:
      case hex20:
      case hex27:
      case tet4:
      case tet10:
        nline = actele->g.gvol->ngline;
        gline = actele->g.gvol->gline;
        break;
      default:
        dserror("Unknown type of discretization 3");
        break;
    }
    /*----------------------------------------- loop all lines to element */
    for (j=0; j<nline; j++)
    {
      /*--------- if the actual gline has been done before then continue */
      if (gline[j]!=NULL) continue;
      /* find all adjacent elements and which line numbers of them to this line */
      /*
         actele is the centerelement of the elementpatch
         elepatch is the patch of element around actele
         nelepatch is the number of element on the patch
         (elepatch includes actele)
         */
      inptop_findadjele(actele,elepatch,nelepatch,adjelepatch,adjelelinenum,
          &nadjelepatch,linenodes[j]);
      /*----------------------- set pointer to the actual line in gsurf of gvol */
      gline[j] = &(actdis->gline[counter]);
      /*----------- set pointer to this line in all other elements found */
      /*---------------------------- loop elements adjacent to this line */
      for (k=0; k<nadjelepatch; k++)
      {
        switch (adjelepatch[k]->distyp)
        {
          case line2:
          case line3:
            break;
          case quad4:
          case quad8:
          case quad9:
          case tri3:
          case tri6:
            dsassert(adjelepatch[k]->g.gsurf->gline[adjelelinenum[k]]==NULL,"Lines got mixed up");
            adjelepatch[k]->g.gsurf->gline[adjelelinenum[k]] = &(actdis->gline[counter]);
            break;
          case hex8:
          case hex20:
          case hex27:
          case tet4:
          case tet10:
            dsassert(adjelepatch[k]->g.gvol->gline[adjelelinenum[k]]==NULL,"Lines got mixed up");
            adjelepatch[k]->g.gvol->gline[adjelelinenum[k]] = &(actdis->gline[counter]);
            break;
          default: dserror("Unknown type of discretization 4"); break;
        }
      }
      /*--------------------------------------- increment the line number */
      counter++;
    }/* end j loop over lines of actele */
  } /* end i loop over all elements */
  nlinetosurf = counter;
  /*-------- now loop the 1D elements and make pointers gline <-> element */
  /*
     [[10a]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    actele = &(actdis->element[i]);
    switch (actele->distyp)
    {
      case line2:
      case line3:
        actdis->gline[counter].element = actele;
        actele->g.gline = &(actdis->gline[counter]);
        counter++;
        dsassert(counter<=actdis->ngline,"Initial guess of ngline too small");
        break;
      default:
        break;
    }
  }
  dsassert(counter-nlinetosurf==nlineelement,"lines of 1D elements got mixed up with 2D elements");
  /*------------------------------------- reallocate the vector of glines */
  /*
     [[11]]
     */
  actdis->ngline = counter;
  tmpgline = actdis->gline;
  actdis->gline = (GLINE*)CCAREALLOC(actdis->gline,(actdis->ngline)*sizeof(GLINE));
  /* if the realloc moved the lines to another place, adjust the pointers */
  if (tmpgline != actdis->gline)
  {
    /*--------------- find the pointer offset between new and old glines */
    ptrdistance = (PTRSIZE)(actdis->gline) - (PTRSIZE)(tmpgline);
    /*--------------------- shift pointers gsurf -> gline in 2D elements */
    /*----------------------- shift pointers gvol->glines in 3D elements */
    /*------------------- all other pointers to glines are not yet built */
    for (i=0; i<actdis->numele; i++)
    {
      actele = &(actdis->element[i]);
      switch (actele->distyp)
      {
        case line2:
        case line3:
          ShiftPointer((void**)&(actele->g.gline),ptrdistance);
          break;
        case quad4:
        case quad8:
        case quad9:
        case tri3:
        case tri6:
          nline = actele->g.gsurf->ngline;
          gline = actele->g.gsurf->gline;
          break;
        case hex8:
        case hex20:
        case hex27:
        case tet4:
        case tet10:
          nline = actele->g.gvol->ngline;
          gline = actele->g.gvol->gline;
          break;
        default: dserror("Unknown type of discretization 5"); break;
      }
      for (j=0; j<nline; j++)
      {
        ShiftPointer((void**)&(gline[j]),ptrdistance);
      }
    }
  }/* end of (tmpgline != actdis->gline) */
  /* in the volumetric elements gvol we have to build the gsurf to gline topology */
  /*
     [[12]]
     */
  for (i=0; i<actdis->ngvol; i++)
  {
    actgvol = &(actdis->gvol[i]);
    actele  = actgvol->element;
    /*------------------------------------------- switch type of element */
    switch (actele->distyp)
    {
      case hex8:
      case hex20:
      case hex27:
        /* surface 0 has lines 0 1 2 3 */
        dsassert(actgvol->gsurf[0]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[0]->gline[0] = actgvol->gline[0];
        actgvol->gsurf[0]->gline[1] = actgvol->gline[1];
        actgvol->gsurf[0]->gline[2] = actgvol->gline[2];
        actgvol->gsurf[0]->gline[3] = actgvol->gline[3];
        /* surface 1 has lines 0 5 8 4 */
        dsassert(actgvol->gsurf[1]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[1]->gline[0] = actgvol->gline[0];
        actgvol->gsurf[1]->gline[1] = actgvol->gline[5];
        actgvol->gsurf[1]->gline[2] = actgvol->gline[8];
        actgvol->gsurf[1]->gline[3] = actgvol->gline[4];
        /* surface 2 has lines 1 6 9 5 */
        dsassert(actgvol->gsurf[2]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[2]->gline[0] = actgvol->gline[1];
        actgvol->gsurf[2]->gline[1] = actgvol->gline[6];
        actgvol->gsurf[2]->gline[2] = actgvol->gline[9];
        actgvol->gsurf[2]->gline[3] = actgvol->gline[5];
        /* surface 3 has lines 2 7 10 6 */
        dsassert(actgvol->gsurf[3]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[3]->gline[0] = actgvol->gline[2];
        actgvol->gsurf[3]->gline[1] = actgvol->gline[7];
        actgvol->gsurf[3]->gline[2] = actgvol->gline[10];
        actgvol->gsurf[3]->gline[3] = actgvol->gline[6];
        /* surface 4 has lines 3 4 11 7 */
        dsassert(actgvol->gsurf[4]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[4]->gline[0] = actgvol->gline[3];
        actgvol->gsurf[4]->gline[1] = actgvol->gline[4];
        actgvol->gsurf[4]->gline[2] = actgvol->gline[11];
        actgvol->gsurf[4]->gline[3] = actgvol->gline[7];
        /* surface 5 has lines 8 9 10 11 */
        dsassert(actgvol->gsurf[5]->ngline==4,"wrong number of lines in surface");
        actgvol->gsurf[5]->gline[0] = actgvol->gline[8];
        actgvol->gsurf[5]->gline[1] = actgvol->gline[9];
        actgvol->gsurf[5]->gline[2] = actgvol->gline[10];
        actgvol->gsurf[5]->gline[3] = actgvol->gline[11];
        break;
      case tet4:
      case tet10:
        /* surface 0 has lines 0 1 2 */
        dsassert(actgvol->gsurf[0]->ngline==3,"wrong number of lines in surface");
        actgvol->gsurf[0]->gline[0] = actgvol->gline[0];
        actgvol->gsurf[0]->gline[1] = actgvol->gline[1];
        actgvol->gsurf[0]->gline[2] = actgvol->gline[2];
        /* surface 1 has lines 0 3 4 */
        dsassert(actgvol->gsurf[1]->ngline==3,"wrong number of lines in surface");
        actgvol->gsurf[1]->gline[0] = actgvol->gline[0];
        actgvol->gsurf[1]->gline[1] = actgvol->gline[3];
        actgvol->gsurf[1]->gline[2] = actgvol->gline[4];
        /* surface 2 has lines 2 5 3 */
        dsassert(actgvol->gsurf[2]->ngline==3,"wrong number of lines in surface");
        actgvol->gsurf[2]->gline[0] = actgvol->gline[2];
        actgvol->gsurf[2]->gline[1] = actgvol->gline[5];
        actgvol->gsurf[2]->gline[2] = actgvol->gline[3];
        /* surface 3 has lines 1 4 5 */
        dsassert(actgvol->gsurf[3]->ngline==3,"wrong number of lines in surface");
        actgvol->gsurf[3]->gline[0] = actgvol->gline[1];
        actgvol->gsurf[3]->gline[1] = actgvol->gline[4];
        actgvol->gsurf[3]->gline[2] = actgvol->gline[5];
        break;
      default: dserror("Unknown type of discretization 6"); break;
    }
  }/* end loop i over gvols */

  /*------------------------------------ build topology gsurfs <-> glines */
  /*----------------------------- count the number of gsurf to each gline */
  /*
     [[13]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    /*----------------------------------------------- set active element */
    actele = &(actdis->element[i]);
    /*------------------------------------------------------ check distyp */
    switch (actele->distyp)
    {
      case line2:
      case line3:
        break;
      case quad4:
      case quad8:
      case quad9:
      case tri3:
      case tri6:
        actgsurf = actele->g.gsurf; isgvol=0;
        for (j=0; j<actgsurf->ngline; j++)
        {
          dsassert(actgsurf->gline[j]!=NULL,"missing line to surface");
          actgsurf->gline[j]->ngsurf++;
        }
        break;
      case hex8:
      case hex20:
      case hex27:
      case tet4:
      case tet10:
        actgvol = actele->g.gvol;   isgvol=1;
        for (k=0; k<actgvol->ngsurf; k++)
          for (j=0; j<actgvol->gsurf[k]->ngline; j++)
          {
            dsassert(actgvol->gsurf[k]->gline[j]!=NULL,"missing line to surface");
            actgvol->gsurf[k]->gline[j]->ngsurf++;
          }
        break;
      default: dserror("Unknown type of discretization 7"); break;
    }
  }/* end i loop over all elements */
  /*------------ loop all glines and allocate pointers in gline to gsurfs */
  for (i=0; i<actdis->ngline; i++)
  {
    actgline = &(actdis->gline[i]);
    /*dsassert(actgline->ngsurf!=0,"Error in gline<->gsurf connectivity");*/
    actgline->gsurf = (GSURF**)CCACALLOC(actgline->ngsurf,sizeof(GSURF*));
  }
  /*----------------- loop all surfaces and make pointer glines -> gsurfs */
  for (i=0; i<actdis->ngsurf; i++)
  {
    actgsurf = &(actdis->gsurf[i]);
    for (j=0; j<actgsurf->ngline; j++)
    {
      actgline = actgsurf->gline[j];
      for (k=0; k<actgline->ngsurf; k++)
      {
        if (actgline->gsurf[k]==NULL) break;
      }
      dsassert(k!=actgline->ngsurf,"No empty pointer in gline to gsurf found");
      actgline->gsurf[k] = actgsurf;
    }
  }
  /*------------------------------------ count number of nodes on 1 gline */
  /*------------------------------------ allocate pointers gline to gnode */
  /*
     [[14]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    isgline=0; isgsurf=0;isgvol=0;
    /*----------------------------------------------- set active element */
    actele = &(actdis->element[i]);
    /*------------------------------------------------------ check distyp */
    switch (actele->distyp)
    {
      case line2: actgline = actele->g.gline; isgline=1; nnodeperline=2; break;
      case line3: actgline = actele->g.gline; isgline=1; nnodeperline=3; break;
      case quad4: actgsurf = actele->g.gsurf; isgsurf=1; nnodeperline=2; break;
      case quad8: actgsurf = actele->g.gsurf; isgsurf=1; nnodeperline=3; break;
      case quad9: actgsurf = actele->g.gsurf; isgsurf=1; nnodeperline=3; break;
      case tri3:  actgsurf = actele->g.gsurf; isgsurf=1; nnodeperline=2; break;
      case tri6:  actgsurf = actele->g.gsurf; isgsurf=1; nnodeperline=3; break;
      case hex8:  actgvol = actele->g.gvol;   isgvol=1; nnodeperline=2; break;
      case hex20: actgvol = actele->g.gvol;   isgvol=1; nnodeperline=3; break;
      case hex27: actgvol = actele->g.gvol;   isgvol=1; nnodeperline=3; break;
      case tet4:  actgvol = actele->g.gvol;   isgvol=1; nnodeperline=2; break;
      case tet10: actgvol = actele->g.gvol;   isgvol=1; nnodeperline=3; break;
      default: dserror("Unknown type of discretization 8"); break;
    }
    if (isgvol)
    {
      for (j=0; j<actgvol->ngline; j++)
      {
        dsassert(actgvol->gline[j]->ngnode==0 || actgvol->gline[j]->ngnode==nnodeperline,
            "Nonconforming grid detected, this wil lead to severe topology mix-up");
        actgvol->gline[j]->ngnode =  nnodeperline;
        if (!actgvol->gline[j]->gnode)
          actgvol->gline[j]->gnode = (GNODE**)CCACALLOC(nnodeperline,sizeof(GNODE*));
      }
    }
    if (isgsurf)
    {
      for (j=0; j<actgsurf->ngline; j++)
      {
        dsassert(actgsurf->gline[j]->ngnode==0 || actgsurf->gline[j]->ngnode==nnodeperline,
            "Nonconforming grid detected, this wil lead to severe topology mix-up");
        actgsurf->gline[j]->ngnode =  nnodeperline;
        if (!actgsurf->gline[j]->gnode)
          actgsurf->gline[j]->gnode = (GNODE**)CCACALLOC(nnodeperline,sizeof(GNODE*));
      }
    }
    if (isgline)
    {
      dsassert(actgline->ngnode==0 || actgline->ngnode==nnodeperline,
          "Nonconforming grid detected, this wil lead to severe topology mix-up");
      actgline->ngnode =  nnodeperline;
      if (!actgline->gnode)
        actgline->gnode = (GNODE**)CCACALLOC(nnodeperline,sizeof(GNODE*));
    }


  } /* end i loop over elements */
  /*----------------------------------------- make pointer gline -> gnode */
  /*
     [[15]]
     */
  for (i=0; i<actdis->numele; i++)
  {
    /*----------------------------------------------- set active element */
    actele = &(actdis->element[i]);
    /*------------------------------------------------------ check distyp */
    switch (actele->distyp)
    {
      case line2:
        actgline = actele->g.gline;
        /* link and count new gline at gnond only, if not yet pointered */
        /* 1st gnode */
        if(actgline->gnode[0] == NULL)
        {
          actgline->gnode[0] = actgline->element->node[0]->gnode;
          actgline->gnode[0]->ngline++;
        }
        /* 2nd gnode */
        if(actgline->gnode[1] == NULL)
        {
          actgline->gnode[1] = actgline->element->node[1]->gnode;
          actgline->gnode[1]->ngline++;
        }
        break;
      case line3:
        actgline = actele->g.gline;
        /* link and count new gline at gnond only, if not yet pointered */
        /* 1st gnode */
        if(actgline->gnode[0] == NULL)
        {
          actgline->gnode[0] = actgline->element->node[0]->gnode;
          actgline->gnode[0]->ngline++;
        }
        /* 2nd gnode */
        if(actgline->gnode[1] == NULL)
        {
          actgline->gnode[1] = actgline->element->node[1]->gnode;
          actgline->gnode[1]->ngline++;
        }
        /* 3rd gnode */
        if(actgline->gnode[2] == NULL)
        {
          actgline->gnode[2] = actgline->element->node[2]->gnode;
          actgline->gnode[2]->ngline++;
        }
        break;
      case quad4:
        actgsurf = actele->g.gsurf;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgsurf->gline[0]->gnode[0] == NULL)
        {
          actgsurf->gline[0]->gnode[0] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[0]->gnode[0]->ngline++;
        }
        /*         2nd gnode */
        if(actgsurf->gline[0]->gnode[1] == NULL)
        {
          actgsurf->gline[0]->gnode[1] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[0]->gnode[1]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgsurf->gline[1]->gnode[0] == NULL)
        {
          actgsurf->gline[1]->gnode[0] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[1]->gnode[0]->ngline++;
        }
        /*         2nd gnode */
        if(actgsurf->gline[1]->gnode[1] == NULL)
        {
          actgsurf->gline[1]->gnode[1] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[1]->gnode[1]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgsurf->gline[2]->gnode[0] == NULL)
        {
          actgsurf->gline[2]->gnode[0] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[2]->gnode[0]->ngline++;
        }
        /*         2nd gnode */
        if(actgsurf->gline[2]->gnode[1] == NULL)
        {
          actgsurf->gline[2]->gnode[1] = actgsurf->element->node[3]->gnode;
          actgsurf->gline[2]->gnode[1]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgsurf->gline[3]->gnode[0] == NULL)
        {
          actgsurf->gline[3]->gnode[0] = actgsurf->element->node[3]->gnode;
          actgsurf->gline[3]->gnode[0]->ngline++;
        }
        /*         2nd gnode */
        if(actgsurf->gline[3]->gnode[1] == NULL)
        {
          actgsurf->gline[3]->gnode[1] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[3]->gnode[1]->ngline++;
        }
        break;
      case quad8:
      case quad9:
        actgsurf = actele->g.gsurf;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgsurf->gline[0]->gnode[0] == NULL)
        {
          actgsurf->gline[0]->gnode[0] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgsurf->gline[0]->gnode[1] == NULL)
        {
          actgsurf->gline[0]->gnode[1] = actgsurf->element->node[4]->gnode;
          actgsurf->gline[0]->gnode[1]->ngline++;
        }
        /* line 0: 3rd gnode */
        if(actgsurf->gline[0]->gnode[2] == NULL)
        {
          actgsurf->gline[0]->gnode[2] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[0]->gnode[2]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgsurf->gline[1]->gnode[0] == NULL)
        {
          actgsurf->gline[1]->gnode[0] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgsurf->gline[1]->gnode[1] == NULL)
        {
          actgsurf->gline[1]->gnode[1] = actgsurf->element->node[5]->gnode;
          actgsurf->gline[1]->gnode[1]->ngline++;
        }
        /* line 1: 3rd gnode */
        if(actgsurf->gline[1]->gnode[2] == NULL)
        {
          actgsurf->gline[1]->gnode[2] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[1]->gnode[2]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgsurf->gline[2]->gnode[0] == NULL)
        {
          actgsurf->gline[2]->gnode[0] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgsurf->gline[2]->gnode[1] == NULL)
        {
          actgsurf->gline[2]->gnode[1] = actgsurf->element->node[6]->gnode;
          actgsurf->gline[2]->gnode[1]->ngline++;
        }
        /* line 2: 3rd gnode */
        if(actgsurf->gline[2]->gnode[2] == NULL)
        {
          actgsurf->gline[2]->gnode[2] = actgsurf->element->node[3]->gnode;
          actgsurf->gline[2]->gnode[2]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgsurf->gline[3]->gnode[0] == NULL)
        {
          actgsurf->gline[3]->gnode[0] = actgsurf->element->node[3]->gnode;
          actgsurf->gline[3]->gnode[0]->ngline++;
        }
        /* line 3: 2nd gnode */
        if(actgsurf->gline[3]->gnode[1] == NULL)
        {
          actgsurf->gline[3]->gnode[1] = actgsurf->element->node[7]->gnode;
          actgsurf->gline[3]->gnode[1]->ngline++;
        }
        /* line 3: 3rd gnode */
        if(actgsurf->gline[3]->gnode[2] == NULL)
        {
          actgsurf->gline[3]->gnode[2] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[3]->gnode[2]->ngline++;
        }
        break;
      case tri3:
        actgsurf = actele->g.gsurf;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgsurf->gline[0]->gnode[0] == NULL)
        {
          actgsurf->gline[0]->gnode[0] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgsurf->gline[0]->gnode[1] == NULL)
        {
          actgsurf->gline[0]->gnode[1] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[0]->gnode[1]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgsurf->gline[1]->gnode[0] == NULL)
        {
          actgsurf->gline[1]->gnode[0] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgsurf->gline[1]->gnode[1] == NULL)
        {
          actgsurf->gline[1]->gnode[1] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[1]->gnode[1]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgsurf->gline[2]->gnode[0] == NULL)
        {
        actgsurf->gline[2]->gnode[0] = actgsurf->element->node[2]->gnode;
        actgsurf->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgsurf->gline[2]->gnode[1] == NULL)
        {
          actgsurf->gline[2]->gnode[1] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[2]->gnode[1]->ngline++;
        }
        break;
      case tri6:
        actgsurf = actele->g.gsurf;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgsurf->gline[0]->gnode[0] == NULL)
        {
          actgsurf->gline[0]->gnode[0] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgsurf->gline[0]->gnode[1] == NULL)
        {
          actgsurf->gline[0]->gnode[1] = actgsurf->element->node[3]->gnode;
          actgsurf->gline[0]->gnode[1]->ngline++;
        }
        /* line 0: 3rd gnode */
        if(actgsurf->gline[0]->gnode[2] == NULL)
        {
          actgsurf->gline[0]->gnode[2] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[0]->gnode[2]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgsurf->gline[1]->gnode[0] == NULL)
        {
          actgsurf->gline[1]->gnode[0] = actgsurf->element->node[1]->gnode;
          actgsurf->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgsurf->gline[1]->gnode[1] == NULL)
        {
          actgsurf->gline[1]->gnode[1] = actgsurf->element->node[4]->gnode;
          actgsurf->gline[1]->gnode[1]->ngline++;
        }
        /* line 1: 3rd gnode */
        if(actgsurf->gline[1]->gnode[2] == NULL)
        {
          actgsurf->gline[1]->gnode[2] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[1]->gnode[2]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgsurf->gline[2]->gnode[0] == NULL)
        {
          actgsurf->gline[2]->gnode[0] = actgsurf->element->node[2]->gnode;
          actgsurf->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgsurf->gline[2]->gnode[1] == NULL)
        {
          actgsurf->gline[2]->gnode[1] = actgsurf->element->node[5]->gnode;
          actgsurf->gline[2]->gnode[1]->ngline++;
        }
        /* line 2: 3rd gnode */
        if(actgsurf->gline[2]->gnode[2] == NULL)
        {
          actgsurf->gline[2]->gnode[2] = actgsurf->element->node[0]->gnode;
          actgsurf->gline[2]->gnode[2]->ngline++;
        }
        break;
      case hex8:
        actgvol = actele->g.gvol;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgvol->gline[0]->gnode[0] == NULL)
        {
          actgvol->gline[0]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgvol->gline[0]->gnode[1] == NULL)
        {
          actgvol->gline[0]->gnode[1] = actgvol->element->node[1]->gnode;
          actgvol->gline[0]->gnode[1]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgvol->gline[1]->gnode[0] == NULL)
        {
          actgvol->gline[1]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgvol->gline[1]->gnode[1] == NULL)
        {
          actgvol->gline[1]->gnode[1] = actgvol->element->node[2]->gnode;
          actgvol->gline[1]->gnode[1]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgvol->gline[2]->gnode[0] == NULL)
        {
          actgvol->gline[2]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgvol->gline[2]->gnode[1] == NULL)
        {
          actgvol->gline[2]->gnode[1] = actgvol->element->node[3]->gnode;
          actgvol->gline[2]->gnode[1]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgvol->gline[3]->gnode[0] == NULL)
        {
          actgvol->gline[3]->gnode[0] = actgvol->element->node[3]->gnode;
          actgvol->gline[3]->gnode[0]->ngline++;
        }
        /* line 3: 2nd gnode */
        if(actgvol->gline[3]->gnode[1] == NULL)
        {
          actgvol->gline[3]->gnode[1] = actgvol->element->node[0]->gnode;
          actgvol->gline[3]->gnode[1]->ngline++;
        }
        /* line 4: 1st gnode */
        if(actgvol->gline[4]->gnode[0] == NULL)
        {
          actgvol->gline[4]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[4]->gnode[0]->ngline++;
        }
        /* line 4: 2nd gnode */
        if(actgvol->gline[4]->gnode[1] == NULL)
        {
          actgvol->gline[4]->gnode[1] = actgvol->element->node[4]->gnode;
          actgvol->gline[4]->gnode[1]->ngline++;
        }
        /* line 5: 1st gnode */
        if(actgvol->gline[5]->gnode[0] == NULL)
        {
          actgvol->gline[5]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[5]->gnode[0]->ngline++;
        }
        /* line 5: 2nd gnode */
        if(actgvol->gline[5]->gnode[1] == NULL)
        {
          actgvol->gline[5]->gnode[1] = actgvol->element->node[5]->gnode;
          actgvol->gline[5]->gnode[1]->ngline++;
        }
        /* line 6: 1st gnode */
        if(actgvol->gline[6]->gnode[0] == NULL)
        {
          actgvol->gline[6]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[6]->gnode[0]->ngline++;
        }
        /* line 6: 2nd gnode */
        if(actgvol->gline[6]->gnode[1] == NULL)
        {
          actgvol->gline[6]->gnode[1] = actgvol->element->node[6]->gnode;
          actgvol->gline[6]->gnode[1]->ngline++;
        }
        /* line 7: 1st gnode */
        if(actgvol->gline[7]->gnode[0] == NULL)
        {
          actgvol->gline[7]->gnode[0] = actgvol->element->node[3]->gnode;
          actgvol->gline[7]->gnode[0]->ngline++;
        }
        /* line 7: 2nd gnode */
        if(actgvol->gline[7]->gnode[1] == NULL)
        {
          actgvol->gline[7]->gnode[1] = actgvol->element->node[7]->gnode;
          actgvol->gline[7]->gnode[1]->ngline++;
        }
        /* line 8: 1st gnode */
        if(actgvol->gline[8]->gnode[0] == NULL)
        {
          actgvol->gline[8]->gnode[0] = actgvol->element->node[4]->gnode;
          actgvol->gline[8]->gnode[0]->ngline++;
        }
        /* line 8: 2nd gnode */
        if(actgvol->gline[8]->gnode[1] == NULL)
        {
          actgvol->gline[8]->gnode[1] = actgvol->element->node[5]->gnode;
          actgvol->gline[8]->gnode[1]->ngline++;
        }
        /* line 9: 1st gnode */
        if(actgvol->gline[9]->gnode[0] == NULL)
        {
          actgvol->gline[9]->gnode[0] = actgvol->element->node[5]->gnode;
          actgvol->gline[9]->gnode[0]->ngline++;
        }
        /* line 9: 2nd gnode */
        if(actgvol->gline[9]->gnode[1] == NULL)
        {
          actgvol->gline[9]->gnode[1] = actgvol->element->node[6]->gnode;
          actgvol->gline[9]->gnode[1]->ngline++;
        }
        /* line 10: 1st gnode */
        if(actgvol->gline[10]->gnode[0] == NULL)
        {
          actgvol->gline[10]->gnode[0] = actgvol->element->node[6]->gnode;
          actgvol->gline[10]->gnode[0]->ngline++;
        }
        /* line 10: 2nd gnode */
        if(actgvol->gline[10]->gnode[1] == NULL)
        {
          actgvol->gline[10]->gnode[1] = actgvol->element->node[7]->gnode;
          actgvol->gline[10]->gnode[1]->ngline++;
        }
        /* line 11: 1st gnode */
        if(actgvol->gline[11]->gnode[0] == NULL)
        {
          actgvol->gline[11]->gnode[0] = actgvol->element->node[7]->gnode;
          actgvol->gline[11]->gnode[0]->ngline++;
        }
        /* line 11: 2nd gnode */
        if(actgvol->gline[11]->gnode[1] == NULL)
        {
          actgvol->gline[11]->gnode[1] = actgvol->element->node[4]->gnode;
          actgvol->gline[11]->gnode[1]->ngline++;
        }
        break;
      case hex20:
      case hex27:
        actgvol = actele->g.gvol;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgvol->gline[0]->gnode[0] == NULL)
        {
          actgvol->gline[0]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgvol->gline[0]->gnode[1] == NULL)
        {
          actgvol->gline[0]->gnode[1] = actgvol->element->node[8]->gnode;
          actgvol->gline[0]->gnode[1]->ngline++;
        }
        /* line 0: 3rd gnode */
        if(actgvol->gline[0]->gnode[2] == NULL)
        {
          actgvol->gline[0]->gnode[2] = actgvol->element->node[1]->gnode;
          actgvol->gline[0]->gnode[2]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgvol->gline[1]->gnode[0] == NULL)
        {
          actgvol->gline[1]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgvol->gline[1]->gnode[1] == NULL)
        {
          actgvol->gline[1]->gnode[1] = actgvol->element->node[9]->gnode;
          actgvol->gline[1]->gnode[1]->ngline++;
        }
        /* line 1: 3rd gnode */
        if(actgvol->gline[1]->gnode[2] == NULL)
        {
          actgvol->gline[1]->gnode[2] = actgvol->element->node[2]->gnode;
          actgvol->gline[1]->gnode[2]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgvol->gline[2]->gnode[0] == NULL)
        {
          actgvol->gline[2]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgvol->gline[2]->gnode[1] == NULL)
        {
          actgvol->gline[2]->gnode[1] = actgvol->element->node[10]->gnode;
          actgvol->gline[2]->gnode[1]->ngline++;
        }
        /* line 2: 3rd gnode */
        if(actgvol->gline[2]->gnode[2] == NULL)
        {
          actgvol->gline[2]->gnode[2] = actgvol->element->node[3]->gnode;
          actgvol->gline[2]->gnode[2]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgvol->gline[3]->gnode[0] == NULL)
        {
          actgvol->gline[3]->gnode[0] = actgvol->element->node[3]->gnode;
          actgvol->gline[3]->gnode[0]->ngline++;
        }
        /* line 3: 2nd gnode */
        if(actgvol->gline[3]->gnode[1] == NULL)
        {
          actgvol->gline[3]->gnode[1] = actgvol->element->node[11]->gnode;
          actgvol->gline[3]->gnode[1]->ngline++;
        }
        /* line 3: 3rd gnode */
        if(actgvol->gline[3]->gnode[2] == NULL)
        {
          actgvol->gline[3]->gnode[2] = actgvol->element->node[0]->gnode;
          actgvol->gline[3]->gnode[2]->ngline++;
        }
        /* line 4: 1st gnode */
        if(actgvol->gline[4]->gnode[0] == NULL)
        {
          actgvol->gline[4]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[4]->gnode[0]->ngline++;
        }
        /* line 4: 2nd gnode */
        if(actgvol->gline[4]->gnode[1] == NULL)
        {
          actgvol->gline[4]->gnode[1] = actgvol->element->node[12]->gnode;
          actgvol->gline[4]->gnode[1]->ngline++;
        }
        /* line 4: 3rd gnode */
        if(actgvol->gline[4]->gnode[2] == NULL)
        {
          actgvol->gline[4]->gnode[2] = actgvol->element->node[4]->gnode;
          actgvol->gline[4]->gnode[2]->ngline++;
        }
        /* line 5: 1st gnode */
        if(actgvol->gline[5]->gnode[0] == NULL)
        {
          actgvol->gline[5]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[5]->gnode[0]->ngline++;
        }
        /* line 5: 2nd gnode */
        if(actgvol->gline[5]->gnode[1] == NULL)
        {
          actgvol->gline[5]->gnode[1] = actgvol->element->node[13]->gnode;
          actgvol->gline[5]->gnode[1]->ngline++;
        }
        /* line 5: 3rd gnode */
        if(actgvol->gline[5]->gnode[2] == NULL)
        {
          actgvol->gline[5]->gnode[2] = actgvol->element->node[5]->gnode;
          actgvol->gline[5]->gnode[2]->ngline++;
        }
        /* line 6: 1st gnode */
        if(actgvol->gline[6]->gnode[0] == NULL)
        {
          actgvol->gline[6]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[6]->gnode[0]->ngline++;
        }
        /* line 6: 2nd gnode */
        if(actgvol->gline[6]->gnode[1] == NULL)
        {
          actgvol->gline[6]->gnode[1] = actgvol->element->node[14]->gnode;
          actgvol->gline[6]->gnode[1]->ngline++;
        }
        /* line 6: 3rd gnode */
        if(actgvol->gline[6]->gnode[2] == NULL)
        {
          actgvol->gline[6]->gnode[2] = actgvol->element->node[6]->gnode;
          actgvol->gline[6]->gnode[2]->ngline++;
        }
        /* line 7: 1st gnode */
        if(actgvol->gline[7]->gnode[0] == NULL)
        {
          actgvol->gline[7]->gnode[0] = actgvol->element->node[3]->gnode;
          actgvol->gline[7]->gnode[0]->ngline++;
        }
        /* line 7: 2nd gnode */
        if(actgvol->gline[7]->gnode[1] == NULL)
        {
          actgvol->gline[7]->gnode[1] = actgvol->element->node[15]->gnode;
          actgvol->gline[7]->gnode[1]->ngline++;
        }
        /* line 7: 3rd gnode */
        if(actgvol->gline[7]->gnode[2] == NULL)
        {
          actgvol->gline[7]->gnode[2] = actgvol->element->node[7]->gnode;
          actgvol->gline[7]->gnode[2]->ngline++;
        }
        /* line 8: 1st gnode */
        if(actgvol->gline[8]->gnode[0] == NULL)
        {
          actgvol->gline[8]->gnode[0] = actgvol->element->node[4]->gnode;
          actgvol->gline[8]->gnode[0]->ngline++;
        }
        /* line 8: 2nd gnode */
        if(actgvol->gline[8]->gnode[1] == NULL)
        {
          actgvol->gline[8]->gnode[1] = actgvol->element->node[16]->gnode;
          actgvol->gline[8]->gnode[1]->ngline++;
        }
        /* line 8: 3rd gnode */
        if(actgvol->gline[8]->gnode[2] == NULL)
        {
          actgvol->gline[8]->gnode[2] = actgvol->element->node[5]->gnode;
          actgvol->gline[8]->gnode[2]->ngline++;
        }
        /* line 9: 1st gnode */
        if(actgvol->gline[9]->gnode[0] == NULL)
        {
          actgvol->gline[9]->gnode[0] = actgvol->element->node[5]->gnode;
          actgvol->gline[9]->gnode[0]->ngline++;
        }
        /* line 9: 2nd gnode */
        if(actgvol->gline[9]->gnode[1] == NULL)
        {
          actgvol->gline[9]->gnode[1] = actgvol->element->node[17]->gnode;
          actgvol->gline[9]->gnode[1]->ngline++;
        }
        /* line 9: 3rd gnode */
        if(actgvol->gline[9]->gnode[2] == NULL)
        {
          actgvol->gline[9]->gnode[2] = actgvol->element->node[6]->gnode;
          actgvol->gline[9]->gnode[2]->ngline++;
        }
        /* line 10: 1st gnode */
        if(actgvol->gline[10]->gnode[0] == NULL)
        {
          actgvol->gline[10]->gnode[0] = actgvol->element->node[6]->gnode;
          actgvol->gline[10]->gnode[0]->ngline++;
        }
        /* line 10: 2nd gnode */
        if(actgvol->gline[10]->gnode[1] == NULL)
        {
          actgvol->gline[10]->gnode[1] = actgvol->element->node[18]->gnode;
          actgvol->gline[10]->gnode[1]->ngline++;
        }
        /* line 10: 3rd gnode */
        if(actgvol->gline[10]->gnode[2] == NULL)
        {
          actgvol->gline[10]->gnode[2] = actgvol->element->node[7]->gnode;
          actgvol->gline[10]->gnode[2]->ngline++;
        }
        /* line 11: 1st gnode */
        if(actgvol->gline[11]->gnode[0] == NULL)
        {
          actgvol->gline[11]->gnode[0] = actgvol->element->node[7]->gnode;
          actgvol->gline[11]->gnode[0]->ngline++;
        }
        /* line 11: 2nd gnode */
        if(actgvol->gline[11]->gnode[1] == NULL)
        {
          actgvol->gline[11]->gnode[1] = actgvol->element->node[19]->gnode;
          actgvol->gline[11]->gnode[1]->ngline++;
        }
        /* line 11: 3rd gnode */
        if(actgvol->gline[11]->gnode[2] == NULL)
        {
          actgvol->gline[11]->gnode[2] = actgvol->element->node[4]->gnode;
          actgvol->gline[11]->gnode[2]->ngline++;
        }
        break;
      case tet4:
        actgvol = actele->g.gvol;
        /* line 0: 1st gnode */
        if(actgvol->gline[0]->gnode[0] == NULL)
        {
          actgvol->gline[0]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgvol->gline[0]->gnode[1] == NULL)
        {
          actgvol->gline[0]->gnode[1] = actgvol->element->node[1]->gnode;
          actgvol->gline[0]->gnode[1]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgvol->gline[1]->gnode[0] == NULL)
        {
          actgvol->gline[1]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgvol->gline[1]->gnode[1] == NULL)
        {
          actgvol->gline[1]->gnode[1] = actgvol->element->node[2]->gnode;
          actgvol->gline[1]->gnode[1]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgvol->gline[2]->gnode[0] == NULL)
        {
          actgvol->gline[2]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgvol->gline[2]->gnode[1] == NULL)
        {
          actgvol->gline[2]->gnode[1] = actgvol->element->node[0]->gnode;
          actgvol->gline[2]->gnode[1]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgvol->gline[3]->gnode[0] == NULL)
        {
          actgvol->gline[3]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[3]->gnode[0]->ngline++;
        }
        /* line 3: 2nd gnode */
        if(actgvol->gline[3]->gnode[1] == NULL)
        {
          actgvol->gline[3]->gnode[1] = actgvol->element->node[3]->gnode;
          actgvol->gline[3]->gnode[1]->ngline++;
        }
        /* line 4: 1st gnode */
        if(actgvol->gline[4]->gnode[0] == NULL)
        {
          actgvol->gline[4]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[4]->gnode[0]->ngline++;
        }
        /* line 4: 2nd gnode */
        if(actgvol->gline[4]->gnode[1] == NULL)
        {
          actgvol->gline[4]->gnode[1] = actgvol->element->node[3]->gnode;
          actgvol->gline[4]->gnode[1]->ngline++;
        }
        /* line 5: 1st gnode */
        if(actgvol->gline[5]->gnode[0] == NULL)
        {
          actgvol->gline[5]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[5]->gnode[0]->ngline++;
        }
        /* line 5: 2nd gnode */
        if(actgvol->gline[5]->gnode[1] == NULL)
        {
          actgvol->gline[5]->gnode[1] = actgvol->element->node[3]->gnode;
          actgvol->gline[5]->gnode[1]->ngline++;
        }
        break;
      case tet10:
        actgvol = actele->g.gvol;
        /* link and count new gline at gnond only, if not yet pointered */
        /* line 0: 1st gnode */
        if(actgvol->gline[0]->gnode[0] == NULL)
        {
          actgvol->gline[0]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[0]->gnode[0]->ngline++;
        }
        /* line 0: 2nd gnode */
        if(actgvol->gline[0]->gnode[1] == NULL)
        {
          actgvol->gline[0]->gnode[1] = actgvol->element->node[4]->gnode;
          actgvol->gline[0]->gnode[1]->ngline++;
        }
        /* line 0: 3rd gnode */
        if(actgvol->gline[0]->gnode[2] == NULL)
        {
          actgvol->gline[0]->gnode[2] = actgvol->element->node[1]->gnode;
          actgvol->gline[0]->gnode[2]->ngline++;
        }
        /* line 1: 1st gnode */
        if(actgvol->gline[1]->gnode[0] == NULL)
        {
          actgvol->gline[1]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[1]->gnode[0]->ngline++;
        }
        /* line 1: 2nd gnode */
        if(actgvol->gline[1]->gnode[1] == NULL)
        {
          actgvol->gline[1]->gnode[1] = actgvol->element->node[5]->gnode;
          actgvol->gline[1]->gnode[1]->ngline++;
        }
        /* line 1: 3rd gnode */
        if(actgvol->gline[1]->gnode[2] == NULL)
        {
          actgvol->gline[1]->gnode[2] = actgvol->element->node[2]->gnode;
          actgvol->gline[1]->gnode[2]->ngline++;
        }
        /* line 2: 1st gnode */
        if(actgvol->gline[2]->gnode[0] == NULL)
        {
          actgvol->gline[2]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[2]->gnode[0]->ngline++;
        }
        /* line 2: 2nd gnode */
        if(actgvol->gline[2]->gnode[1] == NULL)
        {
          actgvol->gline[2]->gnode[1] = actgvol->element->node[6]->gnode;
          actgvol->gline[2]->gnode[1]->ngline++;
        }
        /* line 2: 3rd gnode */
        if(actgvol->gline[2]->gnode[2] == NULL)
        {
          actgvol->gline[2]->gnode[2] = actgvol->element->node[0]->gnode;
          actgvol->gline[2]->gnode[2]->ngline++;
        }
        /* line 3: 1st gnode */
        if(actgvol->gline[3]->gnode[0] == NULL)
        {
          actgvol->gline[3]->gnode[0] = actgvol->element->node[0]->gnode;
          actgvol->gline[3]->gnode[0]->ngline++;
        }
        /* line 3: 2nd gnode */
        if(actgvol->gline[3]->gnode[1] == NULL)
        {
          actgvol->gline[3]->gnode[1] = actgvol->element->node[7]->gnode;
          actgvol->gline[3]->gnode[1]->ngline++;
        }
        /* line 3: 3rd gnode */
        if(actgvol->gline[3]->gnode[2] == NULL)
        {
          actgvol->gline[3]->gnode[2] = actgvol->element->node[3]->gnode;
          actgvol->gline[3]->gnode[2]->ngline++;
        }
        /* line 4: 1st gnode */
        if(actgvol->gline[4]->gnode[0] == NULL)
        {
          actgvol->gline[4]->gnode[0] = actgvol->element->node[1]->gnode;
          actgvol->gline[4]->gnode[0]->ngline++;
        }
        /* line 4: 2nd gnode */
        if(actgvol->gline[4]->gnode[1] == NULL)
        {
          actgvol->gline[4]->gnode[1] = actgvol->element->node[8]->gnode;
          actgvol->gline[4]->gnode[1]->ngline++;
        }
        /* line 4: 3rd gnode */
        if(actgvol->gline[4]->gnode[2] == NULL)
        {
          actgvol->gline[4]->gnode[2] = actgvol->element->node[3]->gnode;
          actgvol->gline[4]->gnode[2]->ngline++;
        }
        /* line 5: 1st gnode */
        if(actgvol->gline[5]->gnode[0] == NULL)
        {
          actgvol->gline[5]->gnode[0] = actgvol->element->node[2]->gnode;
          actgvol->gline[5]->gnode[0]->ngline++;
        }
        /* line 5: 2nd gnode */
        if(actgvol->gline[5]->gnode[1] == NULL)
        {
          actgvol->gline[5]->gnode[1] = actgvol->element->node[9]->gnode;
          actgvol->gline[5]->gnode[1]->ngline++;
        }
        /* line 5: 3rd gnode */
        if(actgvol->gline[5]->gnode[2] == NULL)
        {
          actgvol->gline[5]->gnode[2] = actgvol->element->node[3]->gnode;
          actgvol->gline[5]->gnode[2]->ngline++;
        }
        break;
      default: dserror("Unknown type of discretization 9"); break;
    }
  } /* end i loop over elements */

  /*----------------------------------------- make pointer gnode -> gline */
  for (i=0; i<actdis->ngline; i++)
  {
    actgline = &(actdis->gline[i]);
    for (j=0; j<actgline->ngnode; j++)
    {
      actgnode = actgline->gnode[j];
      dsassert(actgnode->ngline!=0,"gnode <-> gline topology mixed up");
      if (!actgnode->gline)
        actgnode->gline = (GLINE**)CCACALLOC(actgnode->ngline,sizeof(GLINE*));
      for (k=0; k<actgnode->ngline; k++)
        if (actgnode->gline[k]==NULL) break;
      dsassert(k!=actgnode->ngline,"gnode <-> gline topology mixed up");
      actgnode->gline[k] = actgline;
    }
  }
  /*----------------------------------------- make pointers GSURF ->GNODE */
  for (i=0; i<actdis->numele; i++)
  {
    /*----------------------------------------------- set active element */
    actele = &(actdis->element[i]);
    /*------------------------------------------------------ check distyp */
    switch (actele->distyp)
    {
      case line2:
      case line3:
        break;
      case quad4:
        actgsurf = actele->g.gsurf;
        actgsurf->ngnode = 4;
        actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
        for (j=0; j<actele->numnp; j++)
          actgsurf->gnode[j] = actele->node[j]->gnode;
        break;
      case quad8:
        actgsurf = actele->g.gsurf;
        actgsurf->ngnode = 8;
        actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
        for (j=0; j<actele->numnp; j++)
          actgsurf->gnode[j] = actele->node[j]->gnode;
        break;
      case quad9:
        actgsurf = actele->g.gsurf;
        actgsurf->ngnode = 9;
        actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
        for (j=0; j<actele->numnp; j++)
          actgsurf->gnode[j] = actele->node[j]->gnode;
        break;
      case tri3:
        actgsurf = actele->g.gsurf;
        actgsurf->ngnode = 3;
        actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
        for (j=0; j<actele->numnp; j++)
          actgsurf->gnode[j] = actele->node[j]->gnode;
        break;
      case tri6:
        actgsurf = actele->g.gsurf;
        actgsurf->ngnode = 6;
        actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
        for (j=0; j<actele->numnp; j++)
          actgsurf->gnode[j] = actele->node[j]->gnode;
        break;
      case hex8:
        actgvol = actele->g.gvol;
        if (actgvol->gsurf[0]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[0];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[2]->gnode;
          actgsurf->gnode[3] = actele->node[3]->gnode;
        }
        if (actgvol->gsurf[1]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[1];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[5]->gnode;
          actgsurf->gnode[3] = actele->node[4]->gnode;
        }
        if (actgvol->gsurf[2]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[2];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[1]->gnode;
          actgsurf->gnode[1] = actele->node[2]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[5]->gnode;
        }
        if (actgvol->gsurf[3]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[3];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[2]->gnode;
          actgsurf->gnode[1] = actele->node[3]->gnode;
          actgsurf->gnode[2] = actele->node[7]->gnode;
          actgsurf->gnode[3] = actele->node[6]->gnode;
        }
        if (actgvol->gsurf[4]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[4];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[3]->gnode;
          actgsurf->gnode[1] = actele->node[0]->gnode;
          actgsurf->gnode[2] = actele->node[4]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
        }
        if (actgvol->gsurf[5]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[5];
          actgsurf->ngnode = 4;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[4]->gnode;
          actgsurf->gnode[1] = actele->node[5]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
        }
        break;
      case hex20:
        actgvol = actele->g.gvol;
        if (actgvol->gsurf[0]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[0];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[2]->gnode;
          actgsurf->gnode[3] = actele->node[3]->gnode;
          actgsurf->gnode[4] = actele->node[8]->gnode;
          actgsurf->gnode[5] = actele->node[9]->gnode;
          actgsurf->gnode[6] = actele->node[10]->gnode;
          actgsurf->gnode[7] = actele->node[11]->gnode;
        }
        if (actgvol->gsurf[1]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[1];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[5]->gnode;
          actgsurf->gnode[3] = actele->node[4]->gnode;
          actgsurf->gnode[4] = actele->node[8]->gnode;
          actgsurf->gnode[5] = actele->node[13]->gnode;
          actgsurf->gnode[6] = actele->node[15]->gnode;
          actgsurf->gnode[7] = actele->node[12]->gnode;
        }
        if (actgvol->gsurf[2]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[2];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[1]->gnode;
          actgsurf->gnode[1] = actele->node[2]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[5]->gnode;
          actgsurf->gnode[4] = actele->node[9]->gnode;
          actgsurf->gnode[5] = actele->node[14]->gnode;
          actgsurf->gnode[6] = actele->node[17]->gnode;
          actgsurf->gnode[7] = actele->node[13]->gnode;
        }
        if (actgvol->gsurf[3]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[3];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[2]->gnode;
          actgsurf->gnode[1] = actele->node[3]->gnode;
          actgsurf->gnode[2] = actele->node[7]->gnode;
          actgsurf->gnode[3] = actele->node[6]->gnode;
          actgsurf->gnode[4] = actele->node[10]->gnode;
          actgsurf->gnode[5] = actele->node[15]->gnode;
          actgsurf->gnode[6] = actele->node[18]->gnode;
          actgsurf->gnode[7] = actele->node[14]->gnode;
        }
        if (actgvol->gsurf[4]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[4];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[3]->gnode;
          actgsurf->gnode[1] = actele->node[0]->gnode;
          actgsurf->gnode[2] = actele->node[4]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
          actgsurf->gnode[4] = actele->node[11]->gnode;
          actgsurf->gnode[5] = actele->node[12]->gnode;
          actgsurf->gnode[6] = actele->node[19]->gnode;
          actgsurf->gnode[7] = actele->node[15]->gnode;
        }
        if (actgvol->gsurf[5]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[5];
          actgsurf->ngnode = 8;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[4]->gnode;
          actgsurf->gnode[1] = actele->node[5]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
          actgsurf->gnode[4] = actele->node[16]->gnode;
          actgsurf->gnode[5] = actele->node[17]->gnode;
          actgsurf->gnode[6] = actele->node[18]->gnode;
          actgsurf->gnode[7] = actele->node[19]->gnode;
        }
        break;
      case hex27:
        actgvol = actele->g.gvol;
        if (actgvol->gsurf[0]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[0];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[2]->gnode;
          actgsurf->gnode[3] = actele->node[3]->gnode;
          actgsurf->gnode[4] = actele->node[8]->gnode;
          actgsurf->gnode[5] = actele->node[9]->gnode;
          actgsurf->gnode[6] = actele->node[10]->gnode;
          actgsurf->gnode[7] = actele->node[11]->gnode;
          actgsurf->gnode[8] = actele->node[20]->gnode;
        }
        if (actgvol->gsurf[1]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[1];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[5]->gnode;
          actgsurf->gnode[3] = actele->node[4]->gnode;
          actgsurf->gnode[4] = actele->node[8]->gnode;
          actgsurf->gnode[5] = actele->node[13]->gnode;
          actgsurf->gnode[6] = actele->node[15]->gnode;
          actgsurf->gnode[7] = actele->node[12]->gnode;
          actgsurf->gnode[8] = actele->node[21]->gnode;
        }
        if (actgvol->gsurf[2]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[2];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[1]->gnode;
          actgsurf->gnode[1] = actele->node[2]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[5]->gnode;
          actgsurf->gnode[4] = actele->node[9]->gnode;
          actgsurf->gnode[5] = actele->node[14]->gnode;
          actgsurf->gnode[6] = actele->node[17]->gnode;
          actgsurf->gnode[7] = actele->node[13]->gnode;
          actgsurf->gnode[8] = actele->node[22]->gnode;
        }
        if (actgvol->gsurf[3]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[3];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[2]->gnode;
          actgsurf->gnode[1] = actele->node[3]->gnode;
          actgsurf->gnode[2] = actele->node[7]->gnode;
          actgsurf->gnode[3] = actele->node[6]->gnode;
          actgsurf->gnode[4] = actele->node[10]->gnode;
          actgsurf->gnode[5] = actele->node[15]->gnode;
          actgsurf->gnode[6] = actele->node[18]->gnode;
          actgsurf->gnode[7] = actele->node[14]->gnode;
          actgsurf->gnode[8] = actele->node[23]->gnode;
        }
        if (actgvol->gsurf[4]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[4];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[3]->gnode;
          actgsurf->gnode[1] = actele->node[0]->gnode;
          actgsurf->gnode[2] = actele->node[4]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
          actgsurf->gnode[4] = actele->node[11]->gnode;
          actgsurf->gnode[5] = actele->node[12]->gnode;
          actgsurf->gnode[6] = actele->node[19]->gnode;
          actgsurf->gnode[7] = actele->node[15]->gnode;
          actgsurf->gnode[8] = actele->node[24]->gnode;
        }
        if (actgvol->gsurf[5]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[5];
          actgsurf->ngnode = 9;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[4]->gnode;
          actgsurf->gnode[1] = actele->node[5]->gnode;
          actgsurf->gnode[2] = actele->node[6]->gnode;
          actgsurf->gnode[3] = actele->node[7]->gnode;
          actgsurf->gnode[4] = actele->node[16]->gnode;
          actgsurf->gnode[5] = actele->node[17]->gnode;
          actgsurf->gnode[6] = actele->node[18]->gnode;
          actgsurf->gnode[7] = actele->node[19]->gnode;
          actgsurf->gnode[8] = actele->node[25]->gnode;
        }
        break;
      case tet4:
        actgvol = actele->g.gvol;
        if (actgvol->gsurf[0]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[0];
          actgsurf->ngnode = 3;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[2]->gnode;
        }
        if (actgvol->gsurf[1]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[1];
          actgsurf->ngnode = 3;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
        }
        if (actgvol->gsurf[2]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[2];
          actgsurf->ngnode = 3;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[2]->gnode;
          actgsurf->gnode[1] = actele->node[0]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
        }
        if (actgvol->gsurf[3]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[3];
          actgsurf->ngnode = 3;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[1]->gnode;
          actgsurf->gnode[1] = actele->node[2]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
        }
        break;
      case tet10:
        actgvol = actele->g.gvol;
        if (actgvol->gsurf[0]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[0];
          actgsurf->ngnode = 6;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[2]->gnode;
          actgsurf->gnode[3] = actele->node[4]->gnode;
          actgsurf->gnode[4] = actele->node[5]->gnode;
          actgsurf->gnode[5] = actele->node[6]->gnode;
        }
        if (actgvol->gsurf[1]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[1];
          actgsurf->ngnode = 6;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[0]->gnode;
          actgsurf->gnode[1] = actele->node[1]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
          actgsurf->gnode[3] = actele->node[4]->gnode;
          actgsurf->gnode[4] = actele->node[8]->gnode;
          actgsurf->gnode[5] = actele->node[7]->gnode;
        }
        if (actgvol->gsurf[2]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[2];
          actgsurf->ngnode = 6;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[2]->gnode;
          actgsurf->gnode[1] = actele->node[0]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
          actgsurf->gnode[3] = actele->node[6]->gnode;
          actgsurf->gnode[4] = actele->node[7]->gnode;
          actgsurf->gnode[5] = actele->node[9]->gnode;
        }
        if (actgvol->gsurf[3]->ngnode == 0)
        {
          actgsurf = actgvol->gsurf[3];
          actgsurf->ngnode = 6;
          actgsurf->gnode = (GNODE**)CCACALLOC(actgsurf->ngnode,sizeof(GNODE*));
          actgsurf->gnode[0] = actele->node[1]->gnode;
          actgsurf->gnode[1] = actele->node[2]->gnode;
          actgsurf->gnode[2] = actele->node[3]->gnode;
          actgsurf->gnode[3] = actele->node[5]->gnode;
          actgsurf->gnode[4] = actele->node[9]->gnode;
          actgsurf->gnode[5] = actele->node[8]->gnode;
        }
        break;
      default: dserror("Unknown type of discretization 10"); break;
    }
  } /* end i loop over elements */




  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of inp_detailed_topology */


/*----------------------------------------------------------------------*
 | to a given element and a given line find the           m.gee 3/02    |
 | elements on the patch conected to this line                          |
 *----------------------------------------------------------------------*/
static void inptop_findadjele(ELEMENT *centerele, ELEMENT *elepatch[400], INT nelepatch,
                              ELEMENT *adjele[400], INT adjelelinenum[400], INT *nadjele,
                              INT linenodes[3])
{
INT        i,j,counter;
ELEMENT   *actele;
INT        firstnode,scndnode,thirdnode;
INT        matchfirst,matchscnd,matchthird;
INT        searchlnodes[12][3];
INT        ngline = 0;
#ifdef DEBUG
dstrc_enter("inptop_findadjele");
#endif
/*----------------------------------------------------------------------*/
firstnode = linenodes[0];
scndnode  = linenodes[1];
thirdnode = linenodes[2];
/*---------------------------------------- loop elements on input patch */
counter=0;
for (i=0; i<nelepatch; i++)
{
   actele = elepatch[i];
   if (actele==centerele) continue;
   matchfirst=matchscnd=matchthird=0;
   /*-------------------------------------- make lines array for actele */
   inptop_makelinestoele(actele,searchlnodes);
   /*---------------------------------------- check for type of element */
   switch (actele->distyp)
   {
   case line2: ngline = 1;break;
   case line3: ngline = 1;break;
   case quad4: ngline = actele->g.gsurf->ngline;break;
   case quad8: ngline = actele->g.gsurf->ngline;break;
   case quad9: ngline = actele->g.gsurf->ngline;break;
   case tri3:  ngline = actele->g.gsurf->ngline;break;
   case tri6:  ngline = actele->g.gsurf->ngline;break;
   case hex8:  ngline = actele->g.gvol->ngline; break;
   case hex20: ngline = actele->g.gvol->ngline; break;
   case hex27: ngline = actele->g.gvol->ngline; break;
   case tet4:  ngline = actele->g.gvol->ngline; break;
   case tet10: ngline = actele->g.gvol->ngline; break;
   default: dserror("Unknown type of discretization 11"); break;
   }
   for (j=0; j<ngline; j++)
   {
      if (searchlnodes[j][0]==firstnode) matchfirst=1;
      if (searchlnodes[j][2]==firstnode) matchfirst=1;
      if (searchlnodes[j][2]==thirdnode) matchthird=1;
      if (searchlnodes[j][0]==thirdnode) matchthird=1;
      if (matchfirst && matchthird)
      {
         adjele[counter] = actele;
         adjelelinenum[counter]=j;
         counter++;
         break;
      }
      else
      {
         matchfirst=matchscnd=matchthird=0;
      }
   }/* end j loop over lines */
}/* end i loop over input patch */
*nadjele=counter;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inptop_findadjele */





/*----------------------------------------------------------------------*
 | fill an array with all nodes of the lines to an element    m.gee 3/02|
 *----------------------------------------------------------------------*/
static void inptop_makelinestoele(ELEMENT *actele, INT linenodes[12][3])
{

#ifdef DEBUG
dstrc_enter("inptop_makelinestoele");
#endif
/*----------------------------------------------------------------------*/
switch (actele->distyp)
{
case line2:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][2] = actele->node[1]->Id;
break;
case line3:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[1]->Id;
   linenodes[0][2] = actele->node[2]->Id;
break;
case quad4:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][2] = actele->node[0]->Id;
break;
case quad8:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[4]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[5]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[6]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][1] = actele->node[7]->Id;
   linenodes[3][2] = actele->node[0]->Id;
break;
case quad9:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[4]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[5]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[6]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][1] = actele->node[7]->Id;
   linenodes[3][2] = actele->node[0]->Id;
break;
case tri3:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][2] = actele->node[0]->Id;
break;
case tri6:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[3]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[4]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[5]->Id;
   linenodes[2][2] = actele->node[0]->Id;
break;
case hex8:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][2] = actele->node[0]->Id;
   /* line 4 */
   linenodes[4][0] = actele->node[0]->Id;
   linenodes[4][2] = actele->node[4]->Id;
   /* line 5 */
   linenodes[5][0] = actele->node[1]->Id;
   linenodes[5][2] = actele->node[5]->Id;
   /* line 6 */
   linenodes[6][0] = actele->node[2]->Id;
   linenodes[6][2] = actele->node[6]->Id;
   /* line 7 */
   linenodes[7][0] = actele->node[3]->Id;
   linenodes[7][2] = actele->node[7]->Id;
   /* line 8 */
   linenodes[8][0] = actele->node[4]->Id;
   linenodes[8][2] = actele->node[5]->Id;
   /* line 9 */
   linenodes[9][0] = actele->node[5]->Id;
   linenodes[9][2] = actele->node[6]->Id;
   /* line 10 */
   linenodes[10][0] = actele->node[6]->Id;
   linenodes[10][2] = actele->node[7]->Id;
   /* line 11 */
   linenodes[11][0] = actele->node[7]->Id;
   linenodes[11][2] = actele->node[4]->Id;
break;
case hex20:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[8]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[9]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[10]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][1] = actele->node[11]->Id;
   linenodes[3][2] = actele->node[0]->Id;
   /* line 4 */
   linenodes[4][0] = actele->node[0]->Id;
   linenodes[4][1] = actele->node[12]->Id;
   linenodes[4][2] = actele->node[4]->Id;
   /* line 5 */
   linenodes[5][0] = actele->node[1]->Id;
   linenodes[5][1] = actele->node[13]->Id;
   linenodes[5][2] = actele->node[5]->Id;
   /* line 6 */
   linenodes[6][0] = actele->node[2]->Id;
   linenodes[6][1] = actele->node[14]->Id;
   linenodes[6][2] = actele->node[6]->Id;
   /* line 7 */
   linenodes[7][0] = actele->node[3]->Id;
   linenodes[7][1] = actele->node[15]->Id;
   linenodes[7][2] = actele->node[7]->Id;
   /* line 8 */
   linenodes[8][0] = actele->node[4]->Id;
   linenodes[8][1] = actele->node[16]->Id;
   linenodes[8][2] = actele->node[5]->Id;
   /* line 9 */
   linenodes[9][0] = actele->node[5]->Id;
   linenodes[9][1] = actele->node[17]->Id;
   linenodes[9][2] = actele->node[6]->Id;
   /* line 10 */
   linenodes[10][0] = actele->node[6]->Id;
   linenodes[10][1] = actele->node[18]->Id;
   linenodes[10][2] = actele->node[7]->Id;
   /* line 11 */
   linenodes[11][0] = actele->node[7]->Id;
   linenodes[11][1] = actele->node[19]->Id;
   linenodes[11][2] = actele->node[4]->Id;
break;
case hex27:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[8]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[9]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[10]->Id;
   linenodes[2][2] = actele->node[3]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[3]->Id;
   linenodes[3][1] = actele->node[11]->Id;
   linenodes[3][2] = actele->node[0]->Id;
   /* line 4 */
   linenodes[4][0] = actele->node[0]->Id;
   linenodes[4][1] = actele->node[12]->Id;
   linenodes[4][2] = actele->node[4]->Id;
   /* line 5 */
   linenodes[5][0] = actele->node[1]->Id;
   linenodes[5][1] = actele->node[13]->Id;
   linenodes[5][2] = actele->node[5]->Id;
   /* line 6 */
   linenodes[6][0] = actele->node[2]->Id;
   linenodes[6][1] = actele->node[14]->Id;
   linenodes[6][2] = actele->node[6]->Id;
   /* line 7 */
   linenodes[7][0] = actele->node[3]->Id;
   linenodes[7][1] = actele->node[15]->Id;
   linenodes[7][2] = actele->node[7]->Id;
   /* line 8 */
   linenodes[8][0] = actele->node[4]->Id;
   linenodes[8][1] = actele->node[16]->Id;
   linenodes[8][2] = actele->node[5]->Id;
   /* line 9 */
   linenodes[9][0] = actele->node[5]->Id;
   linenodes[9][1] = actele->node[17]->Id;
   linenodes[9][2] = actele->node[6]->Id;
   /* line 10 */
   linenodes[10][0] = actele->node[6]->Id;
   linenodes[10][1] = actele->node[18]->Id;
   linenodes[10][2] = actele->node[7]->Id;
   /* line 11 */
   linenodes[11][0] = actele->node[7]->Id;
   linenodes[11][1] = actele->node[19]->Id;
   linenodes[11][2] = actele->node[4]->Id;
break;
case tet4:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][2] = actele->node[0]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[0]->Id;
   linenodes[3][2] = actele->node[3]->Id;
   /* line 4 */
   linenodes[4][0] = actele->node[1]->Id;
   linenodes[4][2] = actele->node[3]->Id;
   /* line 5 */
   linenodes[5][0] = actele->node[2]->Id;
   linenodes[5][2] = actele->node[3]->Id;
break;
case tet10:
   /* line 0 */
   linenodes[0][0] = actele->node[0]->Id;
   linenodes[0][1] = actele->node[4]->Id;
   linenodes[0][2] = actele->node[1]->Id;
   /* line 1 */
   linenodes[1][0] = actele->node[1]->Id;
   linenodes[1][1] = actele->node[5]->Id;
   linenodes[1][2] = actele->node[2]->Id;
   /* line 2 */
   linenodes[2][0] = actele->node[2]->Id;
   linenodes[2][1] = actele->node[6]->Id;
   linenodes[2][2] = actele->node[0]->Id;
   /* line 3 */
   linenodes[3][0] = actele->node[0]->Id;
   linenodes[3][1] = actele->node[7]->Id;
   linenodes[3][2] = actele->node[3]->Id;
   /* line 4 */
   linenodes[4][0] = actele->node[1]->Id;
   linenodes[4][1] = actele->node[8]->Id;
   linenodes[4][2] = actele->node[3]->Id;
   /* line 5 */
   linenodes[5][0] = actele->node[2]->Id;
   linenodes[5][1] = actele->node[9]->Id;
   linenodes[5][2] = actele->node[3]->Id;
break;
default: dserror("Unknown type of discretization 12");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inptop_makelinestoele */



/*----------------------------------------------------------------------*
 | to a given element and a given surface find the        m.gee 3/02    |
 | volumetric element on the patch which is on the other side           |
 *----------------------------------------------------------------------*/
static void inptop_findotherele(ELEMENT *firstele, ELEMENT **otherele, INT *facenumber,
                                ELEMENT *elepatch[400], INT npatch, INT surfnodes[4])
{
INT        i,j,k;
ELEMENT   *actele = NULL;
GVOL      *actgvol;
INT        surfnodesother[6][4];
INT        matchfirst=0,matchscnd=0,matchthird=0,foundsurface=0;
INT        firstnode,scndnode,thirdnode;
INT        nodespersurf=0;
#ifdef DEBUG
dstrc_enter("inptop_findotherele");
#endif
/*----------------------------------------------------------------------*/
firstnode = surfnodes[0];
scndnode  = surfnodes[1];
thirdnode = surfnodes[2];
/*------------------------------------------ loop the patch of elements */
for (i=0; i<npatch; i++)
{
   foundsurface=-1;
   actele  = elepatch[i];
   actgvol = actele->g.gvol;
   /*---------------------------- we do not want to find firstele again */
   if (actele == firstele) continue;
   /*------------------------------- check number of surfaces to actele */
   /*----------- if actele has 6 surfaces it's a brick, else it's a tet */
   if (actgvol->ngsurf==6)      nodespersurf=4;
   else if (actgvol->ngsurf==4) nodespersurf=3;
   else dserror("Unknown type of element 13");
   /*---------------------------- make the surfacenodes array of actele */
   inptop_makesurfnodes(actele,surfnodesother);
   /* loop nodes in surfnodes, we need a match of three nodes to detect
                                                       a common surface */
   foundsurface=-1;
   for (j=0; j<actgvol->ngsurf; j++)/* loop over surfaces of actgvol */
   {
      for (k=0; k<nodespersurf; k++)/* loop over nodes of surface of actgvol */
      {
         if (surfnodesother[j][k]==firstnode)/* check match first node */
         {
            matchfirst=1;
            break;
         }
      }
      if (matchfirst)
      for (k=0; k<nodespersurf; k++)/* loop over nodes of surface of actgvol */
      {
         if (surfnodesother[j][k]==scndnode)/* check match scnd node */
         {
            matchscnd=1;
            break;
         }
      }
      if (matchscnd)
      for (k=0; k<nodespersurf; k++)/* loop over nodes of surface of actgvol */
      {
         if (surfnodesother[j][k]==thirdnode)/* check match third node */
         {
            matchthird=1;
            break;
         }
      }
      if (matchfirst==1 && matchscnd==1 && matchthird==1)
      {
         foundsurface = j;
         goto exit;
      }
      else
      {
         matchfirst = 0;
         matchscnd  = 0;
         matchthird = 0;
      }
   } /* end loop over surfaces */
} /* end loop over elements on patch */
exit:
if (foundsurface!=-1) /* found an element and a surface */
{
   *otherele = actele;
   *facenumber = foundsurface;
}
else
{
   *otherele=NULL;
   *facenumber=-1;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inptop_findotherele */




/*----------------------------------------------------------------------*
 | make element patch around a given element              m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inptop_makepatch(ELEMENT *centerele, ELEMENT *elepatch[400], INT *nelepatch)
{
INT        i,j,counter;
ELEMENT   *actele;
ELEMENT   *patch[400];
#ifdef DEBUG
dstrc_enter("inptop_makepatch");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------- put all surrounding elements on patch */
counter = 0;
for (i=0; i<centerele->numnp; i++)
{
   for (j=0; j<centerele->node[i]->numele; j++)
   {
      patch[counter] = centerele->node[i]->element[j];
      counter++;
   }
}
/*------------------------------------------ delete doubles on elepatch */
for (i=0; i<counter; i++)
{
   if (!patch[i]) continue;
   actele = patch[i];
   for (j=i+1; j<counter; j++)
   {
      if (actele == patch[j]) patch[j]=NULL;
   }
}
/*-------------------------------------- move all elements to the front */
*nelepatch = 0;
for (i=0; i<counter; i++)
{
   if (patch[i])
   {
      elepatch[*nelepatch] = patch[i];
      (*nelepatch)++;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inptop_makepatch */



/*----------------------------------------------------------------------*
 | fill an array with all nodes of the surfaces               m.gee 3/02|
 *----------------------------------------------------------------------*/
static void inptop_makesurfnodes(ELEMENT *actele, INT surfnodes[6][4])
{
INT        nsurf;
#ifdef DEBUG
dstrc_enter("inptop_makesurfnodes");
#endif
/*----------------------------------------------------------------------*/
switch(actele->distyp)
{
  case hex8:  nsurf=6;
     /* surface 0 */
     surfnodes[0][0] = actele->node[0]->Id;
     surfnodes[0][1] = actele->node[1]->Id;
     surfnodes[0][2] = actele->node[2]->Id;
     surfnodes[0][3] = actele->node[3]->Id;
     /* surface 1 */
     surfnodes[1][0] = actele->node[0]->Id;
     surfnodes[1][1] = actele->node[1]->Id;
     surfnodes[1][2] = actele->node[5]->Id;
     surfnodes[1][3] = actele->node[4]->Id;
     /* surface 2 */
     surfnodes[2][0] = actele->node[1]->Id;
     surfnodes[2][1] = actele->node[2]->Id;
     surfnodes[2][2] = actele->node[6]->Id;
     surfnodes[2][3] = actele->node[5]->Id;
     /* surface 3 */
     surfnodes[3][0] = actele->node[2]->Id;
     surfnodes[3][1] = actele->node[3]->Id;
     surfnodes[3][2] = actele->node[7]->Id;
     surfnodes[3][3] = actele->node[6]->Id;
     /* surface 4 */
     surfnodes[4][0] = actele->node[3]->Id;
     surfnodes[4][1] = actele->node[0]->Id;
     surfnodes[4][2] = actele->node[4]->Id;
     surfnodes[4][3] = actele->node[7]->Id;
     /* surface 5 */
     surfnodes[5][0] = actele->node[4]->Id;
     surfnodes[5][1] = actele->node[5]->Id;
     surfnodes[5][2] = actele->node[6]->Id;
     surfnodes[5][3] = actele->node[7]->Id;
  break;
  case hex20: nsurf=6;
     /* surface 0 */
     surfnodes[0][0] = actele->node[0]->Id;
     surfnodes[0][1] = actele->node[1]->Id;
     surfnodes[0][2] = actele->node[2]->Id;
     surfnodes[0][3] = actele->node[3]->Id;
     /* surface 1 */
     surfnodes[1][0] = actele->node[0]->Id;
     surfnodes[1][1] = actele->node[1]->Id;
     surfnodes[1][2] = actele->node[5]->Id;
     surfnodes[1][3] = actele->node[4]->Id;
     /* surface 2 */
     surfnodes[2][0] = actele->node[1]->Id;
     surfnodes[2][1] = actele->node[2]->Id;
     surfnodes[2][2] = actele->node[6]->Id;
     surfnodes[2][3] = actele->node[5]->Id;
     /* surface 3 */
     surfnodes[3][0] = actele->node[2]->Id;
     surfnodes[3][1] = actele->node[3]->Id;
     surfnodes[3][2] = actele->node[7]->Id;
     surfnodes[3][3] = actele->node[6]->Id;
     /* surface 4 */
     surfnodes[4][0] = actele->node[3]->Id;
     surfnodes[4][1] = actele->node[0]->Id;
     surfnodes[4][2] = actele->node[4]->Id;
     surfnodes[4][3] = actele->node[7]->Id;
     /* surface 5 */
     surfnodes[5][0] = actele->node[4]->Id;
     surfnodes[5][1] = actele->node[5]->Id;
     surfnodes[5][2] = actele->node[6]->Id;
     surfnodes[5][3] = actele->node[7]->Id;
  break;
  case hex27: nsurf=6;
     /* surface 0 */
     surfnodes[0][0] = actele->node[0]->Id;
     surfnodes[0][1] = actele->node[1]->Id;
     surfnodes[0][2] = actele->node[2]->Id;
     surfnodes[0][3] = actele->node[3]->Id;
     /* surface 1 */
     surfnodes[1][0] = actele->node[0]->Id;
     surfnodes[1][1] = actele->node[1]->Id;
     surfnodes[1][2] = actele->node[5]->Id;
     surfnodes[1][3] = actele->node[4]->Id;
     /* surface 2 */
     surfnodes[2][0] = actele->node[1]->Id;
     surfnodes[2][1] = actele->node[2]->Id;
     surfnodes[2][2] = actele->node[6]->Id;
     surfnodes[2][3] = actele->node[5]->Id;
     /* surface 3 */
     surfnodes[3][0] = actele->node[2]->Id;
     surfnodes[3][1] = actele->node[3]->Id;
     surfnodes[3][2] = actele->node[7]->Id;
     surfnodes[3][3] = actele->node[6]->Id;
     /* surface 4 */
     surfnodes[4][0] = actele->node[3]->Id;
     surfnodes[4][1] = actele->node[0]->Id;
     surfnodes[4][2] = actele->node[4]->Id;
     surfnodes[4][3] = actele->node[7]->Id;
     /* surface 5 */
     surfnodes[5][0] = actele->node[4]->Id;
     surfnodes[5][1] = actele->node[5]->Id;
     surfnodes[5][2] = actele->node[6]->Id;
     surfnodes[5][3] = actele->node[7]->Id;
  break;
  case tet4:  nsurf=4;
     /* surface 0 */
     surfnodes[0][0] = actele->node[0]->Id;
     surfnodes[0][1] = actele->node[1]->Id;
     surfnodes[0][2] = actele->node[2]->Id;
     /* surface 1 */
     surfnodes[1][0] = actele->node[0]->Id;
     surfnodes[1][1] = actele->node[1]->Id;
     surfnodes[1][2] = actele->node[3]->Id;
     /* surface 2 */
     surfnodes[2][0] = actele->node[0]->Id;
     surfnodes[2][1] = actele->node[2]->Id;
     surfnodes[2][2] = actele->node[3]->Id;
     /* surface 3 */
     surfnodes[3][0] = actele->node[1]->Id;
     surfnodes[3][1] = actele->node[2]->Id;
     surfnodes[3][2] = actele->node[3]->Id;
  break;
  case tet10: nsurf=4;
     /* surface 0 */
     surfnodes[0][0] = actele->node[0]->Id;
     surfnodes[0][1] = actele->node[1]->Id;
     surfnodes[0][2] = actele->node[2]->Id;
     /* surface 1 */
     surfnodes[1][0] = actele->node[0]->Id;
     surfnodes[1][1] = actele->node[1]->Id;
     surfnodes[1][2] = actele->node[3]->Id;
     /* surface 2 */
     surfnodes[2][0] = actele->node[0]->Id;
     surfnodes[2][1] = actele->node[2]->Id;
     surfnodes[2][2] = actele->node[3]->Id;
     /* surface 3 */
     surfnodes[3][0] = actele->node[1]->Id;
     surfnodes[3][1] = actele->node[2]->Id;
     surfnodes[3][2] = actele->node[3]->Id;
  break;
  default:
     dserror("Unknown type of discretization 14");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inptop_makesurfnodes */
