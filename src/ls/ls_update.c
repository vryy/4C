/*!----------------------------------------------------------------------
\file
\brief ls_fluid.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../xfem/xfem_prototypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



struct _LS_UPDATE             lsupdt;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern 	       ALLDYNA       *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD         *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB        genprob;



static  FIELD                *actfield;               /* pointer to levelset field */
static  LS_DYNAMIC           *lsdyn;                  
static  INT                   completionflag;         /* flag for completion of level set update */
const static DOUBLE           NODETHRESHOLD = 0.001;  /* node threshold (used in case anyone of
                                                       * the nodes lie on the interface) */


static  ARRAY                 lset00_a;
static 	DOUBLE               *lset00;                 /* nodal level set values at time step (n) */
static 	ARRAY                 lset01_a;
static 	DOUBLE               *lset01;                 /* nodal level set values at time step (n+1) */
static 	ARRAY                 xyze_a;
static 	DOUBLE              **xyze;                   /* coordinate array for 2D ls2 element */
static 	ARRAY                 node_to_node_a;
static 	INT                 **node_to_node;           /* node to node connectivity */
static 	ARRAY                 edge_to_edge_a;
static 	INT                 **edge_to_edge;           /* edge to edge connectivity */
static 	ARRAY                 edge_to_node_tri_a;
static 	INT                 **edge_to_node_tri;       /* edge to node connectivity (for triangle) */
static 	ARRAY                 node_to_edge_tri_a;
static 	INT                 **node_to_edge_tri;       /* node to edge connectivity (for triangle) */
static 	ARRAY                 node_to_node_tri_a;
static 	INT                 **node_to_node_tri;       /* node to node connectivity (for triangle) */
static 	ARRAY                 edge_to_triangle_a;
static 	INT                  *edge_to_triangle;       /* edge to triangle connectivity */



static  FILE                 *file01;                 /* file pointers for output */
static  FILE	             *file02;
static  FILE                 *file03;
static  FILE	             *file04;
static  FILE	             *file04_;
static  FILE	             *file05;
static  FILE                 *file06;
static  FILE                 *file07;
static  FILE                 *file08;
static  FILE                 *file09;
static  FILE                 *file10;
static  FILE                 *file11;
static  FILE                 *file12;



/*!----------------------------------------------------------------------
\brief control subroutine for main level set operations

<pre>                                                            irhan 05/04
control subroutine for main level set operations
</pre>

*----------------------------------------------------------------------*/
void ls_main(
  LSFLAG lsflag
  )      
{
  FRONTLSFLAG      frontlsflag;
  XFEMPOLYFLAG     xfempolyflag;

#ifdef DEBUG 
  dstrc_enter("ls_main");
#endif
/*----------------------------------------------------------------------*/
  
  switch (lsflag)
  {
/*---------------------------------------------------------- initialize */
      case ls_initphase:
        /* open */
        xfempolyflag = xfem_poly_open;
        xfem_polygon(xfempolyflag,NULL);
        /* initialize level set update */
        frontlsflag = front_ls_set;
        ls_update(frontlsflag);
        /* construct initial front and write it to output files */
        frontlsflag = front_ls_init;
        ls_update(frontlsflag);
        /* initialize material data */
        frontlsflag = front_ls_initmat;
        ls_update(frontlsflag);
        /* activate the nodes and elements */
        frontlsflag = front_ls_activate;
        ls_update(frontlsflag);
        /* activate layers for localization */
        frontlsflag = front_ls_localize;
        ls_update(frontlsflag);
        /* polygonize the elements cut by the interface */
        frontlsflag = front_ls_polygonize;
        ls_update(frontlsflag);
        /* write initial front into file */
        frontlsflag = front_ls_write;
        ls_update(frontlsflag);
        break;
/*-------------------------------------------------------------- update */
      case ls_updtphase:
        /* update zero level set */
        frontlsflag = front_ls_updt;
        ls_update(frontlsflag);
        /* activate the nodes and elements */
        frontlsflag = front_ls_activate;
        ls_update(frontlsflag);
        /* activate layers for localization */
        frontlsflag = front_ls_localize;
        ls_update(frontlsflag);
        /* update material data */
        frontlsflag = front_ls_updtmat;
        ls_update(frontlsflag);
        /* polygonize the elements cut by the interface */
        frontlsflag = front_ls_polygonize;
        ls_update(frontlsflag);
        if (lsdyn->lsdata->print_on_off==1)
        {
          /* write front into file */
          frontlsflag = front_ls_write;
          ls_update(frontlsflag);
          /* turn off print flag */
          lsdyn->lsdata->print_on_off = 0;
        }
        break;
/*------------------------------------------------------------ finalize */
      case ls_finaphase:
        frontlsflag = front_ls_finalize;
        ls_update(frontlsflag);
        xfempolyflag = xfem_poly_close;
        xfem_polygon(xfempolyflag,NULL);
        break;        
/*------------------------------------------------------------- default */
      default:
        dserror("action unknown\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_main */



/*!----------------------------------------------------------------------
\brief perform operations dictated by execution flag (frontlsflag)

<pre>                                                            irhan 05/04
perform operations dictated by execution flag (frontlsflag)
</pre>

*----------------------------------------------------------------------*/
void ls_update(
  FRONTLSFLAG frontlsflag
  )      
{
#ifdef DEBUG 
  dstrc_enter("ls_update");
#endif
/*----------------------------------------------------------------------*/
  
  switch (frontlsflag)
  {
/*--------------------------------------------- set the data structures */
      case front_ls_set:
        ls_setdata();
        ls_to_matlab();
        break;
/*------------------------------------------------ initialize the front */
      case front_ls_init:
        ls_reset(); /* reset the front */
        ls_check_profile();
        ls_initialize();
        break;
/*-------------------------------------------- initialize material data */
      case front_ls_initmat:
        ls_init_material();
        break;        
/*---------------------------------------------------- update the front */
      case front_ls_updt:
        ls_reset(); /* reset the front */
        ls_check_profile();
        ls_updt();
        break;
/*------------------------------------------------ update material data */
      case front_ls_updtmat:
        ls_updt_material();
        break;        
/*------------------------------------- activate the nodes and elements */
      case front_ls_activate: 
        ls_activate();
        break;
/*------------------------------------------------ perform localization */
      case front_ls_localize: 
        ls_localize();
        break;
/*------------------------ polygonize the elements cut by the interface */
      case front_ls_polygonize:
        ls_polygonize();
        break;
/*------------------------------------------------- write front to file */        
      case front_ls_write:
        ls_write();
        break;
/*------------------------------------------------------------ finalize */
      case front_ls_finalize: 
        ls_finalize();
        break;
/*------------------------------------------------------------- default */
      default:
        dserror("action unknown\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_updt */



/*!----------------------------------------------------------------------
\brief initialize front (called only at the beginning of time
history analysis)

<pre>                                                            irhan 05/04
In this subroutine initial front is constructed. A brute force search is
performed to find the first element (first_in_list) in the interface.
After the first element is found, the interface is constructed by using
the information provided by edge connectivity. Execution is terminated
when the first_in_list is reached once again. Interfaces which are not
closed (those reaching to the boundary of computational doamin) can also
be constructed depending on the flag 'lsdyn->lsdata->boundary_on_off'. The
number of interfaces may be more than one. In this case other interfaces are
able to be constructed by recursive call.
</pre>

*----------------------------------------------------------------------*/
void ls_initialize()
{
  INT              i,j;
  INT              numint;
  INT              intcnt;
  INT              ctrl;
  ELEMENT         *actele;
  LS_INT_DATA     *intdata;
  
#ifdef DEBUG 
  dstrc_enter("ls_initialize");
#endif
/*----------------------------------------------------------------------*/

  /* loop over the elements */
  for (i=0; i<actfield->dis[0].numele; i++)
  {	
    /* access to the element */
    actele = &(actfield->dis[0].element[i]);
    /*
      check whether element is already included
      in one of the previously generated interfaces
    */
    if (actele->e.ls2->is_elcut==1) continue;
    /* get information from the element */
    ls2_calset(actele,xyze,lset00,lset01);
    /* check whether element is cut by the front */
    ls_updatelement(actele,lset01);
    /* check if element became active */
    if (actele->e.ls2->is_elcut==1)
    {
      /* access to the interface number */
      intcnt = actele->e.ls2->intcnt;
      /*
        find the element which is neighbor to eedge of actele and
        set its sedge and coordinates of the start point
      */
      ls_updateneighbor(actele,actele->e.ls2->intdata[intcnt].edge[1],0);
      /* check */
      if (completionflag==1)
      {
        /* not possible! */
        fprintf(file12,"\n\n**ERROR** in ls_initialize()");
        fprintf(file12,"\n**ERROR** completionflag=1");
        fprintf(file12,"\nsomething is wrong level set initialization!");
        ls_printelinfo_to_file(actele);
#ifdef PARALLEL 
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif			
      }
      if (completionflag==2)
      {
        if (lsdyn->lsdata->boundary_on_off==0)
        {
          /*
           * boundary off the computational domain is not allowed
           * to be cut by the interface!
           */
          fprintf(file12,"\n\n**ERROR** in ls initialize()");
          fprintf(file12,"\nelement first_in_list is on boundary!");
          ls_printelinfo_to_file(actele);			
#ifdef PARALLEL 
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
          exit(EXIT_FAILURE);
#endif
        }
        else if (lsdyn->lsdata->boundary_on_off==1)
        {
          /**/
        }
        else
        {
          dserror("lsdyn->lsdata->boundary_on_off set improperly!");
        }
      }
      /* access to the interface number */
      numint = lsupdt.numinterface;
      /* add elements actele to active element list */
      lsupdt.first_in_list[numint] = *actele;
      goto end;
    }
  }
  return;
 end:
  /* construct the front */
  if (completionflag==0)
  {
    actele->e.ls2->intcnt++;
    /* increment the number of active elements (nael++) */
    lsupdt.nael[numint]++;
    ls_construct(actele->e.ls2->intdata[intcnt].nbor);
  }
  else if (completionflag==2)
  {
/*******************************BE CAREFUL*******************************/
/*******************************BE CAREFUL*******************************/
/*******************************BE CAREFUL*******************************/
    /* deactivate the element */
    actele->e.ls2->is_elcut = 0;
    actele->e.ls2->is_sedge_set = 0;
    actele->e.ls2->nlayer = 0;
    /* reset interface data */
    intdata = &(actele->e.ls2->intdata[intcnt]); 
    ls_resetintdata(intdata);
    /* construct the front */
    ls_construct(actele);
/*******************************BE CAREFUL*******************************/
/*******************************BE CAREFUL*******************************/
/*******************************BE CAREFUL*******************************/
  }
  lsupdt.numinterface++;
  /* recursive call for other possible interfaces */
  ls_initialize();
  /* print to screen */  
  printf("\n**WARNING** NUMINTERFACE =%2d",lsupdt.numinterface);  
  /* print to file */  
  fprintf(file12,"\n**WARNING** NUMINTERFACE =%2d",lsupdt.numinterface);  
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_initialize */



/*!----------------------------------------------------------------------
\brief update front (at an intermediate time step)

<pre>                                                            irhan 05/04
The idea is similar to the one used in ls_initialize(). The difference here
is that a confined search is performed instead of brute force search to
find first element (first_in_list) in the interface. To confine the search
region the elements which were previously active (cut by the interface) are
marked. It is done by turning 'ele->e.ls2->is_elsearch' flag on.
</pre>

*----------------------------------------------------------------------*/
void ls_updt()
{
  INT              i,j;
  INT              numint;
  INT              intcnt;
  INT              ctrl;
  ELEMENT         *actele;
  LS_INT_DATA     *intdata;
  
#ifdef DEBUG 
  dstrc_enter("ls_updt");
#endif
/*----------------------------------------------------------------------*/

  /* loop over the elements */
  for (i=0; i<actfield->dis[0].numele; i++)
  {	
    /* access to the element */
    actele = &(actfield->dis[0].element[i]);
    /*
      check whether element is already included
      in one of the previously generated interfaces or
      outisde the active search region
    */	
    if (actele->e.ls2->is_elcut==1 ||
        actele->e.ls2->is_elsearch==0) continue;
    /* get information from the element */
    ls2_calset(actele,xyze,lset00,lset01);
    /* check whether element is cut by the front */
    ls_updatelement(actele,lset01);
    /* check if element became active */
    if (actele->e.ls2->is_elcut==1)
    {
      /* access to the interface number */
      intcnt = actele->e.ls2->intcnt;
      /*
        find the element which is neighbor to eedge of actele and
        set its sedge and coordinates of the start point
      */
      ls_updateneighbor(actele,actele->e.ls2->intdata[intcnt].edge[1],0);
      /* check */
      if (completionflag==1)
      {
        fprintf(file12,"\n\n**ERROR** in ls_update()");
        fprintf(file12,"\ncompletionflag=1");
        ls_printelinfo_to_file(actele);
#ifdef PARALLEL 
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif			
      }
      if (completionflag==2)
      {
        if (lsdyn->lsdata->boundary_on_off==0)
        {
          /*
           * boundary off the computational domain is not allowed
           * to be cut by the interface!
           */
          fprintf(file12,"\n\n**ERROR** in ls_update()");
          fprintf(file12,"\nelement first_in_list is on boundary!");
          ls_printelinfo_to_file(actele);			
#ifdef PARALLEL 
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
          exit(EXIT_FAILURE);
#endif
        }
        else if (lsdyn->lsdata->boundary_on_off==1)
        {
          /**/
        }
        else
        {
          dserror("lsdyn->lsdata->boundary_on_off set improperly!");
        }
      }
      /* access to the interface number */
      numint = lsupdt.numinterface;
      /* add elements actele to active element list */
      lsupdt.first_in_list[numint] = *actele;
      goto end;
    }
  }
  goto end1;
 end:
  /* construct the front */
  if (completionflag==0)
  {
    actele->e.ls2->intcnt++;
    /* increment the number of active elements (nael++) */
    lsupdt.nael[numint]++;
    ls_construct(actele->e.ls2->intdata[intcnt].nbor);
  }
  else if (completionflag==2)
  {
    /* deactivate the element */
    actele->e.ls2->is_elcut = 0;
    actele->e.ls2->is_sedge_set = 0;
    actele->e.ls2->nlayer = 0;
    /* reset interface data */
    intdata = &(actele->e.ls2->intdata[intcnt]); 
    ls_resetintdata(intdata);
    /* construct the front */
    ls_construct(actele);
  }
  lsupdt.numinterface++;
  /* recursive call for other possible interfaces */
  ls_updt();
  /* print to screen */  
  printf("\n**WARNING** NUMINTERFACE =%2d",lsupdt.numinterface);  
  /* print to file */  
  fprintf(file12,"\n**WARNING** NUMINTERFACE =%2d",lsupdt.numinterface);  
 end1:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_updt */



/*!----------------------------------------------------------------------
\brief construct front (called by ls_initialize(), or ls_update())

<pre>                                                            irhan 05/04
In this subroutine interface is constructed starting from the element ele.
Note that ele is not always to be the first_in_list. It is actually the case
when the subroutine is called for reconstruction phase.
</pre>

*----------------------------------------------------------------------*/
void ls_construct(
  ELEMENT *actele
  )
{
  INT           i;
  INT           ctrl;
  INT           numint;
  INT           intcnt;
  INT           Id_first_old;
  
#ifdef DEBUG 
  dstrc_enter("ls_construct");
#endif
/*----------------------------------------------------------------------*/

  numint = lsupdt.numinterface;
  /* check */
  if (actele==NULL) dserror("\n**ERROR** in level set construction!\n");
  
  ctrl = 0;
  while (ctrl==0)
  {
    /* get information from the element */
    ls2_calset(actele,xyze,lset00,lset01);
    /* check whether element is cut by the front */
    ls_updatelement(actele,lset01);
    /*
     * NOTE =>
     * since we are in ls_construct element "must" have been cut
     * by the interface
     */
    if (actele->e.ls2->is_elcut!=1)
    {
      fprintf(file12,"\n\n**ERROR** in ls_construct()");
      fprintf(file12,"\nelement must have been cut by the interface");
      ls_printelinfo_to_file(actele);
#ifdef PARALLEL 
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
      exit(EXIT_FAILURE);
#endif			
    }
    else
    {
      /* access to the interface number */
      intcnt = actele->e.ls2->intcnt;
      /* increment the number of active elements (nael++) */
      lsupdt.nael[numint]++;		
      /*
       * find the element which is neighbor to eedge of actele
       * and set its sedge and coordinates of the start point
       */
      ls_updateneighbor(actele,actele->e.ls2->intdata[intcnt].edge[1],0);
      /* increment the number of interfaces cut the element */
      actele->e.ls2->intcnt++;
      /* check completion of interface construction */      
      if (completionflag==1) goto end;
      else if (completionflag==2) 
      {
        if (lsdyn->lsdata->boundary_on_off==0)
        {
          /*
           * boundary off the computational domain is not allowed
           * to be cut by the interface!
           */
          fprintf(file12,"\n\n**ERROR** in ls_construct()");
          fprintf(file12,"\nelement actele is on boundary!");
          ls_printelinfo_to_file(actele);			
#ifdef PARALLEL 
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
          exit(EXIT_FAILURE);
#endif
        }
        else if (lsdyn->lsdata->boundary_on_off==1)
        {
          if (lsdyn->lsdata->reconstruct_on_off==0)
          {
            /* set type of inerface to open ( 0 => closed; 1 => open ) */
            lsupdt.typeinterface[numint] = 1;
            /* print to screen */
            printf("\n**WARNING** LEVELSET TOUCHES BOUNDARY => RECONSTRUCTION ACTIVE");
            /* print to file */
            fprintf(file12,"\n**WARNING** LEVELSET TOUCHES BOUNDARY => RECONSTRUCTION ACTIVE");            
            /* turn on lsdyn->lsdata->reconstruct_on_off flag */
            lsdyn->lsdata->reconstruct_on_off = 1;
            /* RECONSTRUCTION! ( => use nbor_s ) */
            /* last element to be treated is to be old first_in_list */
            Id_first_old = lsupdt.first_in_list[numint].Id;
            /* set new first_in_list which is to be located on the boundary */
            lsupdt.first_in_list[numint] = *actele;
            /* set ele->e.ls2->recontruct flag for the elements */
            for (i=0; i<lsupdt.nael[numint]; i++)
            {
              intcnt = actele->e.ls2->intcnt-1;
              actele->e.ls2->intdata[intcnt].reconstruct = 1;
              if (i==lsupdt.nael[numint]-1) break;
              actele = actele->e.ls2->intdata[intcnt].nbor_s;
            }
            /* check */
            if (actele->Id!=Id_first_old)
            {
              fprintf(file12,"\n\n**ERROR** in ls_construct()");
              fprintf(file12,"\nactele->Id!=Id_first_old");
              ls_printelinfo_to_file(actele);
#ifdef PARALLEL 
              MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
              exit(EXIT_FAILURE);
#endif			
            }
            /*
             * find the element which is neighbor to eedge of actele and
             * set its sedge and coordinates of the start point
             */
            ls_updateneighbor(actele,actele->e.ls2->intdata[intcnt].edge[0],1);
            /* check */
            if (actele->e.ls2->intdata[intcnt].nbor_s!=NULL)
            {
              /* there are still elements to be processed */
              /* recursive call for ls_construct() */
              ls_construct(actele->e.ls2->intdata[intcnt].nbor_s);
            }
            else
            {
              /* there are no further elements to be processed */
              /* turn off lsdyn->lsdata->reconstruct_on_off flag */
              lsdyn->lsdata->reconstruct_on_off = 0;
            }
            goto end;            
          }
          else if (lsdyn->lsdata->reconstruct_on_off==1)
          {
            lsdyn->lsdata->reconstruct_on_off = 0;
            goto end;
          }
          else
          {
            dserror("lsdyn->lsdata->reconstruct_on_off set improperly!");
          }
        }
        else
        {
          dserror("lsdyn->lsdata->boundary_on_off set improperly!");
        }
      }
      /* set new actele */
      actele = actele->e.ls2->intdata[intcnt].nbor;
    }
  }
 end:
  /* print to screen */
  printf("\n**WARNING** LEVELSET UPDATE COMPLETED!");  
  /* print to screen */
  fprintf(file12,"\n**WARNING** LEVELSET UPDATE COMPLETED!");  
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_construct */



/*!----------------------------------------------------------------------
\brief update the element

<pre>                                                            irhan 05/04
In this subroutine element interface data is set if it is involved in anyone
of the interfaces.
</pre>

*----------------------------------------------------------------------*/
void ls_updatelement(
  ELEMENT* ele,
  DOUBLE* lset
  )      
{
  INT         i,j;
  DOUBLE      u[3];
  DIS_TYP     typ;
  INT         nd;
  INT         is_cut=0;
  INT         sedge,eedge,snode;
  INT         sedge_;
  INT         ntri;
  INT         counter;
  INT         intcnt;
  INT         temp;
  INT         ind[2];
  
#ifdef DEBUG 
  dstrc_enter("ls_updatelement");
#endif
/*----------------------------------------------------------------------*/

  /* access to the interface number */
  intcnt = ele->e.ls2->intcnt;  
  typ = ele->distyp;
  counter = 0;
  switch (typ)
  {
      case tri3:   /* => 3 noded triangular element */
        /*
         * check whether the sedge of element is set by previously
         * activated element
         * NOTE =>
         * it is NOT only the case for the element first_in_list
         */
        sedge_ = 0;
        if (ele->e.ls2->is_sedge_set==1)
        {
          /* access to the start edge */
          /*
           * NOTE =>
           * note that this information is also provided by
           * the element which has already been activated
           */
          sedge_ = ele->e.ls2->intdata[intcnt].edge[0];
          if (sedge_<=0) dserror("start edge not set properly!\n");
        }
        /* loop over the nodes */
        for (j=0; j<3; j++)
        {
          u[j] = lset[j];
        }
        /* check whether subtriangle is cut by zero level set */
        is_tricut(u,&is_cut,&sedge,&eedge,&snode);
        if (!is_cut) 
        {
          /* element interior not cut by interface! */
          return;
        }
        /* switch edges if necessary */         
        if (sedge_!=0 && sedge_!=sedge)
        {
          temp  = sedge;
          sedge = eedge;
          eedge = temp;
        }
        /* set edges */
        ele->e.ls2->intdata[intcnt].edge[0]  = sedge;
        ele->e.ls2->intdata[intcnt].edge[1]  = eedge;
        /* compute intersection with the edges */
        ls_compintersection(ele->e.ls2->intdata[intcnt].p[0],u,sedge,-1);
        ls_compintersection(ele->e.ls2->intdata[intcnt].p[1],u,eedge,-1);
        /* construct polygon connectivity to be used in refined integration (partially) */
        ls_polygon_con(ele->e.ls2->intdata[intcnt].polycon[0],snode,-1);
        /* set is_diagcut */
        ele->e.ls2->intdata[intcnt].is_diagcut = 0;
        if (ele->e.ls2->is_elcut!=1)
        {
          /* activate the element */          
          ele->e.ls2->is_elcut = 1;
          /* add element to the first layer */
          ele->e.ls2->nlayer = 1;
        }
        break;
      case quad4:  /* => 4 noded quadrilateral element */
        ind[0] = 1;
        ind[1] = 2;
        /* check whether the element is set by previously activated element */
        /*
         * NOTE =>
         * it is NOT only the case for the element first_in_list
         */
        sedge_ = 0;
        if (ele->e.ls2->is_sedge_set==1)
        {
          /* access to the start edge */
          /*
           * NOTE =>
           * note that this information is also provided by
           * the element which has already been activated
           */
          sedge_ = ele->e.ls2->intdata[intcnt].edge[0];
          if (sedge_<=0) dserror("start edge not set properly!\n");
          /* check to which triangle this edge belongs to */
          ntri = edge_to_triangle[sedge_-1];
          /* switch the order of triangles if necessary */
          if (ntri==2)
          {
            ind[0] = 2;
            ind[1] = 1;
          }
        }
        /* loop over the triangles */
        for (i=0; i<2; i++)
        {
          /* loop over the nodes */
          for (j=0; j<3; j++)
          {
            nd = node_to_node[ind[i]-1][j]-1;
            u[j] = lset[nd];
          }
          /* check whether subtriangle is cut by zero level set */
          is_tricut(u,&is_cut,&sedge,&eedge,&snode);
          if (!is_cut) 
          {
            counter++;
            if (counter==2)
            {
              /* element interior not cut by interface! */
              goto end1;
            }
            continue;
          }
          if (edge_to_edge[ind[i]-1][eedge-1]!=-1) goto end;
          /* set edge */
          ele->e.ls2->intdata[intcnt].edge[i]  = edge_to_edge[ind[i]-1][sedge-1];
          /* compute intersection with the edges */
          ls_compintersection(ele->e.ls2->intdata[intcnt].p[i],u,sedge,ind[i]);
          /* construct polygon connectivity to be used in refined integration (partially) */
          ls_polygon_con(ele->e.ls2->intdata[intcnt].polycon[ind[i]-1],snode,ind[i]);
          /* compute intersection with the diagonal */
          if (i==1)
          {
            ele->e.ls2->intdata[intcnt].is_diagcut = 1;
            ls_compintersection(ele->e.ls2->intdata[intcnt].pd,u,3,ind[i]);
          }
        }
        if (ele->e.ls2->is_elcut!=1)
        {
          /*  activate the element */          
          ele->e.ls2->is_elcut = 1;
          /*  add element to the first layer */
          ele->e.ls2->nlayer = 1;
        }
        goto end1;
  end:
        /* switch edges if necessary */         
        if (sedge_!=0 && sedge_!=edge_to_edge[ind[i]-1][sedge-1])
        {
          temp  = sedge;
          sedge = eedge;
          eedge = temp;
        }
        /* set edge and solitary node */
        ele->e.ls2->intdata[intcnt].edge[0] = edge_to_edge[ind[i]-1][sedge-1];
        ele->e.ls2->intdata[intcnt].edge[1] = edge_to_edge[ind[i]-1][eedge-1];
        /* compute intersection with the edges */
        ls_compintersection(ele->e.ls2->intdata[intcnt].p[0],u,sedge,ind[i]);
        ls_compintersection(ele->e.ls2->intdata[intcnt].p[1],u,eedge,ind[i]);
        /* construct polygon connectivity to be used in refined integration (partially) */
        ls_polygon_con(ele->e.ls2->intdata[intcnt].polycon[ind[i]-1],snode,ind[i]);
        if (ele->e.ls2->is_elcut!=1)
        {
          /* activate the element */          
          ele->e.ls2->is_elcut = 1;
          /* add element to the first layer */
          ele->e.ls2->nlayer = 1;
        }        
        break;
      default:
        dserror("typ unknown!");
  } /* end switch(typ) */
  
/*----------------------------------------------------------------------*/
 end1:
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_updatelement */



/*!----------------------------------------------------------------------
\brief is triangle (or subtriangle in case of quadrilateral element) cut
by the interface

<pre>                                                            irhan 05/04
In this subroutine it is checked whether triangular domain is cut by the
interface or not. Nodal level set values are used for this purpose. If the
element is cut by the interface, the local edge numbers having the start and
end points of the interface as well as the local node number corresponding to
solitary node are returned. Solitary node is the node which is located opposite
to the edge which is not cut by the interface.
</pre>

*----------------------------------------------------------------------*/
void is_tricut(
  DOUBLE* u,
  INT* is_cut,
  INT* sedge,
  INT* eedge,
  INT* snode
  )
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("ls_tricut");
#endif
/*----------------------------------------------------------------------*/
  
  if (u[0]>0.0)
  {
    for (i=0; i<3; i++)
    {
      u[i] *= -1;
    }
  }
  *is_cut = 1;
  *sedge  = 2;
  *eedge  = 3;
  *snode  = 3;
  if (u[1]<0.0 && u[2]<0.0)
  {
    *is_cut = 0;
    *sedge  = 0;
    *eedge  = 0;		
    *snode  = 0;
  }
  if (u[1]>0.0 && u[2]>0.0)
  {
    *sedge = 1;
    *eedge = 3;
    *snode = 1; 
    return;
  }
  if (u[1]>0.0)
  {
    *sedge = 1;
    *eedge = 2;
    *snode = 2;	
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_istricut */



/*!----------------------------------------------------------------------
\brief compute coordinates of the intersection point between the interface
and local edge 'edge' of the triangular domain.

<pre>                                                            irhan 05/04
compute coordinates of the intersection point between the interface
and local edge 'edge' of the triangular domain.
</pre>

*----------------------------------------------------------------------*/
void ls_compintersection(
  DOUBLE* p,
  DOUBLE* u,
  INT edge,
  INT numtri
  )
{
  INT     nd1,nd2;
  
#ifdef DEBUG 
  dstrc_enter("ls_compintersection");
#endif
/*----------------------------------------------------------------------*/

  /* get the local node numbers corresponding to edge */
  nd1 = edge_to_node_tri[0][edge-1];
  nd2 = edge_to_node_tri[1][edge-1];
  /* compute intersection point */
  ls_comppoint(p,u,nd1,nd2,numtri);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_compintersection */



/*!----------------------------------------------------------------------
\brief compute coordinates of the intersection point between the interface
and local edge defined by nodes 'nd1' and 'nd2'

<pre>                                                            irhan 05/04
compute coordinates of the intersection point between the interface
and local edge defined by nodes 'nd1' and 'nd2'
</pre>

*----------------------------------------------------------------------*/
void ls_comppoint(
  DOUBLE* p,
  DOUBLE* u,
  INT nd1,
  INT nd2,
  INT numtri
  )
{
  DOUBLE     ksi;
  INT        ndg1,ndg2;
  
#ifdef DEBUG 
  dstrc_enter("ls_comppoint");
#endif
/*----------------------------------------------------------------------*/

  /* get the global node numbers corresponding to edge */
  if (numtri==-1)            /* discretization is of triangles */
  {
    ndg1 = nd1;
    ndg2 = nd2;    
  }
  else                       /* discretization is of quadrilaterals */
  {
    ndg1 = node_to_node[numtri-1][nd1-1];
    ndg2 = node_to_node[numtri-1][nd2-1];
  }
  /* compute natural coordinate of the node */
  ksi = (u[nd1-1]+u[nd2-1])/(u[nd1-1]-u[nd2-1]);
  p[0] = 0.5*(1-ksi)*xyze[0][ndg1-1] + 0.5*(1+ksi)*xyze[0][ndg2-1];
  p[1] = 0.5*(1-ksi)*xyze[1][ndg1-1] + 0.5*(1+ksi)*xyze[1][ndg2-1];
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_comppoint */



/*!----------------------------------------------------------------------
\brief update the neighbor

<pre>                                                            irhan 05/04
In this subroutine element neighbor to the end edge of element ele is
found. Then it is appended to element ele and its start edge is set.
Element ele is also appended to neighbor. If the neighbor element is the
first element (first_in_list) in the list, completion flag is set to 1
indicating successfull completion of the interface construction.
</pre>

*----------------------------------------------------------------------*/
void ls_updateneighbor(
  ELEMENT* ele,
  INT eedge,
  INT junction
  )
{
  INT          i,j;
  INT          nelepatch;
  INT          numint;
  INT          intcntele;
  INT          intcntnbor;
  GSURF       *actgsurf,*nborgsurf;
  GLINE       *actgline,*nborgline;
  ELEMENT     *elepatch[50],*nbor;
  
#ifdef DEBUG 
  dstrc_enter("ls_updateneighbor");
#endif
/*----------------------------------------------------------------------*/

  /* access to interface number crossing the element */  
  intcntele = ele->e.ls2->intcnt;
  /* set completion flag */
  completionflag = 0;
  /* access to the interface number */
  numint = lsupdt.numinterface;
  /* access to geometry objects of ele */
  actgsurf = ele->g.gsurf;
  actgline = actgsurf->gline[eedge-1];
  /* construct an element patch around a given element */
  ls_makepatch(ele,elepatch,&nelepatch);
  /* loop */
  for (i=0; i<nelepatch; i++)
  {
    nbor = elepatch[i];
    nborgsurf = nbor->g.gsurf;
    for (j=0; j<nborgsurf->ngline; j++)
    {
      nborgline = nborgsurf->gline[j];
      if (nborgline!=actgline || nbor==ele)
      {
        continue;
      }
      else
      {
        /* access to interface number crossing the element */
        intcntnbor = nbor->e.ls2->intcnt;      
        /* set nbor_s of nbor */
        nbor->e.ls2->intdata[intcntnbor].nbor_s = ele;
        /* check whether nbor is first_in_list */
        if (nbor->Id==lsupdt.first_in_list[numint].Id)
        {
          /* interface construction completed succesfully */
          completionflag = 1;
          goto end;
        }
        else
        {
          /* set start edge */
          nbor->e.ls2->intdata[intcntnbor].edge[0] = j+1;
          /* turn is_sedge_set flag on */
          nbor->e.ls2->is_sedge_set = 1;
/*************************BE CAREFUL*************************************/
/*************************BE CAREFUL*************************************/
/*************************BE CAREFUL*************************************/
          /* set nbor of actele */
          if (junction==0)
          {
            ele->e.ls2->intdata[intcntele].nbor = nbor;
          }
          else
          {
            ele->e.ls2->intdata[intcntele-1].nbor_s = nbor;
          }
/*************************BE CAREFUL*************************************/
/*************************BE CAREFUL*************************************/
/*************************BE CAREFUL*************************************/
          if (nbor->e.ls2->is_elcut==1) /* nbor is already cut */
          {
            /* print to file */
            fprintf(file12,"\nin ls_updateneighbor()");
            fprintf(file12,"\nneighbor element is cut second times by the interface!");
            fprintf(file12,"\njunction = %2d",junction);            
            fprintf(file12,"\nactive element =>");
            ls_printelinfo_to_file(ele);
            fprintf(file12,"\nneighbor to active element =>");            
            ls_printelinfo_to_file(nbor);
          }
        }
      }
      goto end;;
    }
  }
  /* there is no neighbor element */
  nbor = NULL;
  completionflag = 2;
 end:
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_updateneighbor */



/*!----------------------------------------------------------------------
\brief make a patch of elements around centerele.

<pre>                                                            irhan 05/04
In this subroutine an element patch is constructed including the elements
surrounding centerele. Centerele is also included in the patch.
</pre>

*----------------------------------------------------------------------*/
void ls_makepatch(
  ELEMENT* centerele, 
  ELEMENT *elepatch[50], 
  INT *nelepatch
  )
{
  INT          i,j;
  INT          counter;
  ELEMENT     *actele;
  ELEMENT     *patch[50];
  
#ifdef DEBUG 
  dstrc_enter("ls_makepatch");
#endif
/*----------------------------------------------------------------------*/

  /* put all surrounding elements on patch */
  counter = 0;
  for (i=0; i<centerele->numnp; i++)
  {
    for (j=0; j<centerele->node[i]->numele; j++)
    {
      patch[counter] = centerele->node[i]->element[j];
      counter++;
    }
  }
  /* delete doubles on elepatch */
  for (i=0; i<counter; i++)
  {
    if (!patch[i]) continue;
    actele = patch[i];
    for (j=i+1; j<counter; j++)
    {
      if (actele == patch[j]) patch[j]=NULL;
    }
  }
  /* move all elements to the front */
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
} /* end of ls_makepatch */



/*!----------------------------------------------------------------------
\brief activate

<pre>                                                            irhan 05/04
In this subroutine the nodes and elements belonging to first and second
layers surrounding the interface are activated. 
</pre>

*----------------------------------------------------------------------*/
void ls_activate()
{
  INT          ii,i,j;
  INT         *flagactive,*flaginside;
  INT          dummy;
  INT         *nenode;
  DOUBLE       val;
  NODE        *actnode;
  ELEMENT     *actele,*myf2;
  DOUBLE       lset01[4];
  
#ifdef DEBUG 
  dstrc_enter("ls_activate");
#endif
/*----------------------------------------------------------------------*/

  /* activate the nodes of firstlayer */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==1)
    {
      for (j=0; j<actele->numnp; j++)
      {
        actnode =  actele->node[j];
        flagactive = &(actnode->gnode->is_node_active);
        flaginside = &(actnode->gnode->is_node_inside);
        val = actnode->sol_increment.a.da[1][0];
        if (*flagactive==0)
        {
          /* activate the node */
          *flagactive = 1;
          /* set whether node is inside or outside */
          if (val<0.0)
            *flaginside = -1;
          else
            *flaginside = 1;
          /* increment the nand */
          lsupdt.nand++;
        }
      }
    }
    else
    {
      /* do nothing for a while */
    }
  }
  
  /* activate the elements of second layer */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    /* set nenode and enode */
    /*
     * NOTE => these data structures are used
     * in co-operation with X-FEM
     */
    nenode = &(actele->e.ls2->nenode);
    *nenode = 0;
    
    if (actele->e.ls2->is_elcut==1)
    {
      for (j=0; j<actele->numnp; j++)
      {
        actele->e.ls2->enode[*nenode] = j+1;
        (*nenode)++;
      }
      continue;
    }
    /* loop */
    for (j=0; j<actele->numnp; j++)
    {
      actnode =  actele->node[j];
      flagactive = &(actnode->gnode->is_node_active);
      if (*flagactive==0) continue;
      /* store local number of enriched node */
      actele->e.ls2->enode[*nenode] = j+1;
      (*nenode)++;
      /* add element to second layer */
      actele->e.ls2->is_elcut = 2;
      actele->e.ls2->nlayer = 2;      
    }
  }
  
  /* activate the nodes of second layer */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==2)
    {
      for (j=0; j<actele->numnp; j++)
      {
        actnode =  actele->node[j];
        flagactive = &(actnode->gnode->is_node_active);
        if (*flagactive == 1)
        {
          flaginside = &(actnode->gnode->is_node_inside);
          break;
        }
      }
      for (j=0; j<actele->numnp; j++)
      {
        actnode =  actele->node[j];
        flagactive = &(actnode->gnode->is_node_active);
        if (*flagactive==0)
        {
          /* activate the node */
          *flagactive = 2;
          /* set whether node is inside or outside */
          actnode->gnode->is_node_inside = *flaginside;          
          /* increment the nand */
          lsupdt.nand++;
        }
      }
    }
  }
  
  /*
     reset sol_increment of enriched directions corresponding
     to previously active nodes
  */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode =  &(actfield->dis[0].node[i]);
    flagactive = &(actnode->gnode->is_node_active);
    if (*flagactive>2)
    {
      /* reset solution at enriched directions */
      myf2->node[j]->sol_increment.a.da[1][3] = 0.0;
      myf2->node[j]->sol_increment.a.da[1][4] = 0.0;
      myf2->node[j]->sol_increment.a.da[3][3] = 0.0;
      myf2->node[j]->sol_increment.a.da[3][4] = 0.0;
    }
  }
  /* print to screen */
  printf("\n**WARNING** LEVELSET ACTIVATE COMPLETED!");  
  /* print to screen */
  fprintf(file12,"\n**WARNING** LEVELSET ACTIVATE COMPLETED!");  
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_activate */



/*!----------------------------------------------------------------------
\brief write interface related data to files

<pre>                                                            irhan 05/04
In this subroutine interface related data is written into files. These data
files are sent to MATLAB and processed to produce some graphical outputs.
</pre>

*----------------------------------------------------------------------*/
void ls_write()
{
  INT          i,j;
  INT          prncnt;
  INT          counter,counter1;
  INT          flagactive,flaginside;
  NODE        *actnode;
  ELEMENT     *actele;
  
#ifdef DEBUG 
  dstrc_enter("ls_write");
#endif
/*----------------------------------------------------------------------*/

  /*  write number of interfaces to file */	
  fprintf(file05,"%5d\n",lsupdt.numinterface);
  
  for (i=0; i<lsupdt.numinterface; i++)
  {
    /* print out the front */
    actele = &(lsupdt.first_in_list[i]);
    /* write number of active elements to file */	
    fprintf(file08,"%5d\n",lsupdt.nael[i]);

    prncnt = actele->e.ls2->prncnt;
    if (lsupdt.typeinterface[i]==1)
    {
      fprintf(file01,"%10.5E %10.5E\n",actele->e.ls2->intdata[prncnt].p[1][0],
              actele->e.ls2->intdata[prncnt].p[1][1]);
    }
    /* loop */
    for (j=0; j<lsupdt.nael[i]; j++)
    {
      prncnt = actele->e.ls2->prncnt;
      if (actele->e.ls2->intdata[prncnt].reconstruct==0)
      {
        fprintf(file01,"%10.5E %10.5E\n",actele->e.ls2->intdata[prncnt].p[1][0],
                actele->e.ls2->intdata[prncnt].p[1][1]);
      }
      else if (actele->e.ls2->intdata[prncnt].reconstruct==1)
      {
        fprintf(file01,"%10.5E %10.5E\n",actele->e.ls2->intdata[prncnt].p[0][0],
                actele->e.ls2->intdata[prncnt].p[0][1]);
      }
      else
      {
        dserror("lsdyn->lsdata->reconstruct_on_off set improperly!");
      }
      
      actele->e.ls2->prncnt++;

      if (actele->e.ls2->intdata[prncnt].reconstruct==0)
      {
        actele = actele->e.ls2->intdata[prncnt].nbor;
      }
      else if(actele->e.ls2->intdata[prncnt].reconstruct==1)
      {
        actele = actele->e.ls2->intdata[prncnt].nbor_s;        
      }
      else
      {
        dserror("lsdyn->lsdata->reconstruct_on_off set improperly!");
      }
    }
    
    if (lsupdt.typeinterface[i]==0)
    {
      /* close the polygon */
      actele = &(lsupdt.first_in_list[i]);
      prncnt = actele->e.ls2->prncnt-1;    
      fprintf(file01,"%10.5E %10.5E\n\n",actele->e.ls2->intdata[prncnt].p[1][0],
              actele->e.ls2->intdata[prncnt].p[1][1]);
    }
  }
  fprintf(file01,"\n");													  

  /* print out the level set profile */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    fprintf(file02,"%10.5E\n",actnode->sol_increment.a.da[1][0]);
  }
  fprintf(file02,"\n");

  /* print out nodes activated */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    flagactive = actnode->gnode->is_node_active;
    flaginside = actnode->gnode->is_node_inside;
    if (flagactive!=0) fprintf(file07,"%5d %5d %5d\n",actnode->Id+1,flagactive,flaginside);
  }
  fprintf(file07,"\n");    

  fprintf(file06,"%5d\n",lsupdt.nand);

  /* print out lael */
  counter = 0;
  counter1 = 0;
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==1) counter1++;
    if (actele->e.ls2->is_elcut==1 ||
        actele->e.ls2->is_elcut==2)
    {
      counter++;
      fprintf(file03,"%5d %5d\n",actele->Id+1,actele->e.ls2->is_elcut);          
    }
  }
  fprintf(file03,"\n");

  /* print out nael */
  fprintf(file04,"%5d\n",counter);
  fprintf(file04_,"%5d\n",counter1);
  
  /* print out element layers */
  counter = 0;
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->nlayer==0)
    {
      continue;
    }
    else
    {
      counter++;
      fprintf(file10,"%5d %5d\n",actele->Id+1,actele->e.ls2->nlayer);                  
    }
  }
  fprintf(file10,"\n");
  /* print out numlayer */
  fprintf(file11,"%5d\n",counter);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_write */



/*!----------------------------------------------------------------------
\brief reset interface

<pre>                                                            irhan 05/04
In this subroutine interface data is reset, i.e., corresponding flags and
data structures are initialized for the next time step.
</pre>

*----------------------------------------------------------------------*/
void ls_reset()
{
  INT              i,j;
  ELEMENT         *actele;
  NODE            *actnode;
  LS_INT_DATA     *intdata;
  LS_POLY_DATA    *polydata;
  
#ifdef DEBUG 
  dstrc_enter("ls_reset");
#endif
/*----------------------------------------------------------------------*/

  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->nlayer!=0)
    {
      actele->e.ls2->is_elsearch = 1;
    }
    else
    {
      actele->e.ls2->is_elsearch = 0;
    }
    /* reset flags and counters */
    actele->e.ls2->is_elcut = 0;
    actele->e.ls2->is_sedge_set = 0;
    actele->e.ls2->nlayer = 0;
    actele->e.ls2->nenode = 0;
    actele->e.ls2->prncnt = 0;
    actele->e.ls2->rstcnt = 0;
    /* reset enode */
    for (j=0; j<actele->numnp; j++)
    {
      actele->e.ls2->enode[j] = 0;
    }
    /* reset interface data */
    for (j=0; j<actele->e.ls2->intcnt; j++)
    {
      intdata = &(actele->e.ls2->intdata[j]);
      polydata = &(actele->e.ls2->polydata[j]);
      ls_resetintdata(intdata);
      ls_resetpolydata(polydata);      
      intdata->nbor = NULL;
      intdata->nbor_s = NULL;
    }
    actele->e.ls2->intcnt = 0;    
  }
  
  /* reset flag for node activation */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    actnode->gnode->is_node_active = 0;
    actnode->gnode->is_node_inside = 0;
  }
  /* reset nand */
  lsupdt.nand = 0;

  /* reset LS_UPDATE struct */
  lsupdt.first_in_list = NULL;
  for (i=0; i<lsupdt.numinterface; i++)
  {
    lsupdt.nael[i] = 0;
    lsupdt.typeinterface[i] = 0;    
  }
  lsupdt.numinterface = 0;
  lsupdt.first_in_list = (ELEMENT*)CCACALLOC(50,sizeof(ELEMENT));
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_reset */



/*!----------------------------------------------------------------------
\brief reset interface data

<pre>                                                            irhan 05/04
reset interface data
</pre>

*----------------------------------------------------------------------*/
void ls_resetintdata(
  LS_INT_DATA *data
  )
{
  INT     i,j,k;
  
#ifdef DEBUG 
  dstrc_enter("ls_resetintdata");
#endif
/*----------------------------------------------------------------------*/
  
  data->edge[0] = 0;
  data->edge[1] = 0;
  data->p[0][0] = 0.0;
  data->p[0][1] = 0.0;
  data->p[1][0] = 0.0;
  data->p[1][1] = 0.0;
  data->pd[0] = 0.0;
  data->pd[1] = 0.0;
  data->polycon[0][0] = 0;
  data->polycon[1][0] = 0;
  data->polycon[0][1] = 0;
  data->polycon[1][1] = 0;
  data->polycon[0][2] = 0;
  data->polycon[1][2] = 0;
  data->is_diagcut = 0;
  data->reconstruct = 0;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_resetintdata */



/*!----------------------------------------------------------------------
\brief reset polygon data

<pre>                                                            irhan 05/04
reset polygon data
</pre>

*----------------------------------------------------------------------*/
void ls_resetpolydata(
  LS_POLY_DATA *data
  )
{
  INT     i,j,k;
  
#ifdef DEBUG 
  dstrc_enter("ls_resetpolydata");
#endif
/*----------------------------------------------------------------------*/
  
  data->ind[0] = 0;
  data->ind[1] = 0;

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      for (k=0; k<7; k++)
        data->polygonGP[i][j][k] = 0.0;

  for (i=0; i<2; i++)
    for (j=0; j<7; j++)
      data->polygonwgt[i][j] = 0.0;

  for (i=0; i<2; i++)
    for (j=0; j<7; j++)
      data->polygonmat[i][j] = 0;  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_resetpolydata */



/*!----------------------------------------------------------------------
\brief finalize

<pre>                                                            irhan 05/04
finalize
</pre>

*----------------------------------------------------------------------*/
void ls_finalize()
{
#ifdef DEBUG 
  dstrc_enter("ls_finalize");
#endif
/*----------------------------------------------------------------------*/

  /* delete the arrays */
  amdel(&lset00_a);
  amdel(&lset01_a);
  amdel(&xyze_a);
  amdel(&node_to_node_a);
  amdel(&edge_to_edge_a);
  amdel(&edge_to_node_tri_a);
  amdel(&node_to_edge_tri_a);
  amdel(&edge_to_triangle_a);

  /* close the files */
  fclose(file01);
  fclose(file02);
  fclose(file03);
  fclose(file04);
  fclose(file04_);
  fclose(file05);
  fclose(file06);
  fclose(file07);
  fclose(file08);
  fclose(file09);
  fclose(file10);
  fclose(file11);
  fclose(file12);
    
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_finalize */



/*!----------------------------------------------------------------------
\brief print element data onto screen

<pre>                                                            irhan 05/04
print element data
</pre>

*----------------------------------------------------------------------*/
void ls_printelinfo(
  ELEMENT *ele
  )
{
  INT     i;

#ifdef DEBUG 
  dstrc_enter("ls_printelinfo");
#endif
/*----------------------------------------------------------------------*/

  /* print element Id */    
  printf("\nElement=%6d\n\n",ele->Id+1);

  /* print connectivity */  
  printf("con=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      printf("%3d ",ele->lm[i]+1);
  }
  printf("\n]\n");

  /* print coordinates */
  printf("\nxy=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      printf("%10.5E %10.5E\n",ele->node[i]->x[0],ele->node[i]->x[1]);
  }
  printf("]\n");  

  /* access to the interface number */
  for (i=0; i<ele->e.ls2->intcnt; i++)
  {
    printf("\nedge=\n[\n%6d %6d\n]\n\n",ele->e.ls2->intdata[i].edge[0],ele->e.ls2->intdata[i].edge[1]);
  }
  printf("is_elcut=%6d\n\n",ele->e.ls2->is_elcut);
  /* print level set profile */  
  printf("lset01=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      printf("%10.5E ",lset01[i]);
  }
  printf("\n]");
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_printelinfo */



/*!----------------------------------------------------------------------
\brief print element data into file

<pre>                                                            irhan 15/04
print element data into file
</pre>

*----------------------------------------------------------------------*/
void ls_printelinfo_to_file(
  ELEMENT *ele
  )
{
  INT     i;

#ifdef DEBUG 
  dstrc_enter("ls_printelinfo_to_file");
#endif
/*----------------------------------------------------------------------*/

  /* print element Id */    
  fprintf(file12,"\nElement=%6d\n\n",ele->Id+1);

  /* print connectivity */  
  fprintf(file12,"con=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      fprintf(file12,"%3d ",ele->lm[i]+1);
  }
  fprintf(file12,"\n]\n");

  /* print coordinates */
  fprintf(file12,"\nxy=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      fprintf(file12,"%10.5E %10.5E\n",ele->node[i]->x[0],ele->node[i]->x[1]);
  }
  fprintf(file12,"]\n");  

  /* access to the interface number */
  for (i=0; i<ele->e.ls2->intcnt; i++)
  {
    fprintf(file12,"\nedge=\n[\n%6d %6d\n]\n\n",ele->e.ls2->intdata[i].edge[0],ele->e.ls2->intdata[i].edge[1]);
  }
  fprintf(file12,"is_elcut=%6d\n\n",ele->e.ls2->is_elcut);
  /* print level set profile */  
  fprintf(file12,"lset01=\n[\n");
  for (i=0; i<ele->numnp; i++)
  {
      fprintf(file12,"%10.5E ",lset01[i]);
  }
  fprintf(file12,"\n]");
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_printelinfo */



/*!----------------------------------------------------------------------
\brief print node data to screen

<pre>                                                            irhan 05/04
print node data to screen
</pre>

*----------------------------------------------------------------------*/
void ls_printnodeinfo(
  NODE *node
  )
{

#ifdef DEBUG 
  dstrc_enter("ls_printnodeinfo");
#endif
/*----------------------------------------------------------------------*/
  
  printf("\n\nNode=%6d",node->Id+1);
  printf("\nx=\n[\n%10.5E %10.5E\n]",node->x[0],node->x[1]);
  printf("\nval = %10.5E",node->sol_increment.a.da[1][0]);
  printf("\nis_node_active =%2d",node->gnode->is_node_active);
  printf("\nis_node_inside =%2d",node->gnode->is_node_inside);
  printf("\nNumInterface=%2d\n\n",lsupdt.numinterface);  
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_printnodeinfo */



/*!----------------------------------------------------------------------
\brief print node data to file

<pre>                                                            irhan 05/04
print node data to file
</pre>

*----------------------------------------------------------------------*/
void ls_printnodeinfo_to_file(
  NODE *node
  )
{

#ifdef DEBUG 
  dstrc_enter("ls_printnodeinfo_to_file");
#endif
/*----------------------------------------------------------------------*/
  
  fprintf(file12,"\n\nNode=%6d",node->Id+1);
  fprintf(file12,"\nx=\n[\n%10.5E %10.5E\n]",node->x[0],node->x[1]);
  fprintf(file12,"\nval = %10.5E",node->sol_increment.a.da[1][0]);
  fprintf(file12,"\nis_node_active =%2d",node->gnode->is_node_active);
  fprintf(file12,"\nis_node_inside =%2d",node->gnode->is_node_inside);
  fprintf(file12,"\nNumInterface=%2d\n\n",lsupdt.numinterface);  
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_printnodeinfo_to_file */



/*!----------------------------------------------------------------------
\brief initialize the parameters used in level set update

<pre>                                                            irhan 05/04
initialize the parameters used in level set update
</pre>

*----------------------------------------------------------------------*/
void ls_setdata()
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("ls_setdata");
#endif
/*----------------------------------------------------------------------*/
  
  lset00              = amdef("lset00"          , &lset00_a          , MAXNOD_LS2, ONE, "DV");
  lset01              = amdef("lset01"          , &lset01_a          , MAXNOD_LS2, ONE, "DV");
  xyze                = amdef("xyze"            , &xyze_a            , TWO, MAXNOD_LS2, "DA");
  node_to_node        = amdef("node_to_node"    , &node_to_node_a    , 2, 3, "IA");
  edge_to_edge        = amdef("edge_to_edge"    , &edge_to_edge_a    , 2, 3, "IA");
  edge_to_node_tri    = amdef("edge_to_node_tri", &edge_to_node_tri_a, 2, 3, "IA");
  node_to_edge_tri    = amdef("node_to_edge_tri", &node_to_edge_tri_a, 2, 3, "IA");
  node_to_node_tri    = amdef("node_to_node_tri", &node_to_node_tri_a, 2, 3, "IA");  
  edge_to_triangle    = amdef("edge_to_tringle" , &edge_to_triangle_a, 4, 1, "IV");

/*----------------------------------------------------------------------*/
/*-------------------------------- connectivities rectangle => triangle */
/*----------------------------------------------------------------------*/
/* prepare data connectivities useful in coding */
/**
 * FIGURE =>
 *
 *       4          |3|            3
 *        -------------------------
 *        |*                      |
 *        |   *              _    |
 *        |      *          |2|   |
 *        |         *  |-1|       |
 *     |4||     _      *          ||2|
 *        |    |1|        *       |
 *        |                  *    |
 *        |                     * |
 *        -------------------------
 *       1          |1|            2
 */
/*
 * DESCRIPTION =>
 * node_to_node[I-1][J-1] returns the local "node" number of the rectangle
 * corresponding to "Jth" local node of "Ith" triangle 
 */
  node_to_node[0][0] = 4;
  node_to_node[0][1] = 1;
  node_to_node[0][2] = 2;
  node_to_node[1][0] = 2;
  node_to_node[1][1] = 3;
  node_to_node[1][2] = 4;
/*
 * DESCRIPTION =>
 * edge_to_edge[I-1][J-1] returns the local "edge" number of the rectangle
 * corresponding to "Jth" local edge of "Ith" triangle 
 */
  edge_to_edge[0][0] = 4;
  edge_to_edge[0][1] = 1;
  edge_to_edge[0][2] = -1;
  edge_to_edge[1][0] = 2;
  edge_to_edge[1][1] = 3;
  edge_to_edge[1][2] = -1;
/*
 * DESCRIPTION =>
 * edge_to_triangle[I-1] returns the number of traingle to which
 * "Ith" local "edge" of the rectangle corresponds 
 */
  edge_to_triangle[0] = 1;
  edge_to_triangle[1] = 2;
  edge_to_triangle[2] = 2;
  edge_to_triangle[3] = 1;
/*----------------------------------------------------------------------*/
/*------------------------------------------ connectivities => triangle */
/*----------------------------------------------------------------------*/
/*
 * FIGURE =>
 *
 *       3
 *        
 *        |*         
 *        |   *      
 *        |      *   
 *        |         *  |2|
 *     |3||            *  
 *        |               * 
 *        |                  *
 *        |                     *
 *        -------------------------
 *       1          |1|            2
 */
/*
 * DESCRIPTION =>
 * edge_to_node_tri[I-1][J-1] returns local node number
 * corresponding to "Ith" node of "Jth" local edge of triangle 
 */
  edge_to_node_tri[0][0] = 1;
  edge_to_node_tri[1][0] = 2;
  edge_to_node_tri[0][1] = 2;
  edge_to_node_tri[1][1] = 3;
  edge_to_node_tri[0][2] = 3;
  edge_to_node_tri[1][2] = 1;
/*
 * DESCRIPTION =>
 * node_to_node_tri[I-1][J-1] returns local node number
 * corresponding to "Ith" node following the "Jth" local
 * node of the triangle
 */
  node_to_node_tri[0][0] = 2;
  node_to_node_tri[1][0] = 3;
  node_to_node_tri[0][1] = 3;
  node_to_node_tri[1][1] = 1;
  node_to_node_tri[0][2] = 1;
  node_to_node_tri[1][2] = 2;
/*
 * DESCRIPTION =>
 * node_to_edge_tri[I-1][J-1] returns local edge number 
 * corresponding to "Ith" edge of "Jth" local node of triangle 
 */
  node_to_edge_tri[0][0] = 3;
  node_to_edge_tri[1][0] = 1;
  node_to_edge_tri[0][1] = 1;
  node_to_edge_tri[1][1] = 2;
  node_to_edge_tri[0][2] = 2;
  node_to_edge_tri[1][2] = 3;  
  
/* set the flag to check successfullness of level set update */
  completionflag = 0;

/* open files to write some useful information */
  file01 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_levelsetzero","w");
  file02 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_levelsetprofile","w");
  file03 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_lael","w");  
  file04 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nael","w");
  file04_ = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_nael","w");
  file05 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nint","w");
  file06 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nand","w");
  file07 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_land","w");
  file08 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nseg","w");
  file09 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_firstinlist","w");
  file10 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_layer","w");
  file11 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nlayer","w");
  /* this file is used for debugging purposes */
  file12 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_debug","w");  
  
  /*initialize the lsupdt */
  lsupdt.nael = (INT*)CCACALLOC(50,sizeof(INT));
  lsupdt.typeinterface = (INT*)CCACALLOC(50,sizeof(INT));

  /* set actfield */
  actfield = &(field[genprob.numls]);

  /* perform check for the field type */
  if (actfield->fieldtyp!=levelset) dserror("field != levelset");  

  /* initialize the flag for node activation */
  for (i=0; i<actfield->dis[0].ngnode; i++)
    actfield->dis[0].gnode[i].is_node_active = 0;

  /* set lsdyn */
  lsdyn = alldyn[genprob.numls].lsdyn;

  /* initialise level set field */
  ls_init(actfield,lsdyn,3);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_setdata */



/*!----------------------------------------------------------------------
\brief print all the elements into a file

<pre>                                                            irhan 05/04
print all the elements into a file
</pre>

*----------------------------------------------------------------------*/
void ls_printelements()
{
  INT          i;
  ELEMENT     *actele;
  
#ifdef DEBUG
  dstrc_enter("ls_printelements");
#endif
/*----------------------------------------------------------------------*/

  /*loop over the nodes of the discretization associated with field */ 
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    /* access to element */
    actele = &(actfield->dis[0].element[i]);
    /* print the element information */
    ls_printelement(actele);
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_printelements */



/*!----------------------------------------------------------------------
\brief print an element into a file

<pre>                                                            irhan 05/04
print an element into a file
</pre>

*----------------------------------------------------------------------*/
void ls_printelement(
  ELEMENT* ele
  )
{
  INT       i;
  NODE     *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_printelement");
#endif
/*----------------------------------------------------------------------*/

  /*loop over the nodes of ele */
  for (i=0; i<ele->numnp; i++)
  {
    actnode = ele->node[i];
    /* print out the coordinates of the node */
    fprintf(file01,"%10.5E %10.5E\n",actnode->x[0],actnode->x[1]);
  }
  fprintf(file01,"\n");
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_printelement */



/*!----------------------------------------------------------------------
\brief localize

<pre>                                                            irhan 05/04
set the elements and nodes located in one of the layers surrounding the
interface. The total number of layers for localization is set by
'lsdyn->lsdata->numlayer'. 
</pre>

*----------------------------------------------------------------------*/
void ls_localize()
{
  INT          ii,i,j;
  INT         *flagactive,*flaginside;
  NODE        *actnode;
  ELEMENT     *actele;
  
#ifdef DEBUG 
  dstrc_enter("ls_localize");
#endif
/*----------------------------------------------------------------------*/
  
  if (lsdyn->lsdata->numlayer>2)
  {
    for (ii=2; ii<lsdyn->lsdata->numlayer; ii++)
    {
      /* add elements */
      for (i=0; i<actfield->dis[0].numele; i++)
      {
        actele = &(actfield->dis[0].element[i]);
        /* check */ 
        if (actele->e.ls2->nlayer!=0) continue;
        /* loop */
        for (j=0; j<actele->numnp; j++)
        {
          actnode =  actele->node[j];
          flagactive = &(actnode->gnode->is_node_active);
          if (*flagactive!=ii) continue;
          /* add element to second layer */
          actele->e.ls2->nlayer = ii+1;      
        }
      }
      /* add nodes */
      for (i=0; i<actfield->dis[0].numele; i++)
      {
        actele = &(actfield->dis[0].element[i]);
        if (actele->e.ls2->nlayer==ii+1)
        {
          for (j=0; j<actele->numnp; j++)
          {
            actnode =  actele->node[j];
            flagactive = &(actnode->gnode->is_node_active);
            if (*flagactive==ii)
            {
              flaginside = &(actnode->gnode->is_node_inside);
              break;
            }
          }
          for (j=0; j<actele->numnp; j++)
          {
            actnode =  actele->node[j];
            flagactive = &(actnode->gnode->is_node_active);
            if (*flagactive==0)
            {
              /* activate the node */
              *flagactive = ii+1;
              /* set whether node is inside or outside */
              actnode->gnode->is_node_inside = *flaginside;          
              /* increment the nand */
              lsupdt.nand++;
            }
          }
        }
      }
    }
  }

  /* print to screen */
  printf("\n**WARNING** LEVELSET LAYER ACTIVATION FOR LOCALIZATION COMPLETED!");
  /* print to file */
  fprintf(file12,"\n**WARNING** LEVELSET LAYER ACTIVATION FOR LOCALIZATION COMPLETED!");
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_localize */



/*!----------------------------------------------------------------------
\brief write mesh dependent data into files

<pre>                                                            irhan 05/04
write mesh dependent data into files
</pre>

*----------------------------------------------------------------------*/
void ls_to_matlab()
{
  INT          i,j,k;
  ELEMENT     *actele;
  NODE        *actnode;
  FILE        *f1;
  
#ifdef DEBUG 
  dstrc_enter("ls_printelement");
#endif
/*----------------------------------------------------------------------*/

  /* I need the connectivity */
  f1 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_connectivity","w");

  /* write connectivity to file */
  for (i=0; i<genprob.numfld; i++)
  {
    for (j=0; j<field[i].dis[0].numele; j++)
    {
      actele = &(field[i].dis[0].element[j]);
      fprintf(f1,"%5d ",i+1);      
      for (k=0; k<actele->numnp; k++)
      {
        fprintf(f1,"%5d ",actele->lm[k]+1);
      }
      for (k=actele->numnp; k<10; k++)
      {
        fprintf(f1,"%5d ",0);
      }
      fprintf(f1,"\n");
    }
  }
  fclose(f1);

  /* I need the node coordinates */
  f1 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nodecoordinates","w");

  /* write node coordinates to file */
  for (i=0; i<genprob.numfld; i++)
  {
    for (j=0; j<field[i].dis[0].numnp; j++)
    {
      actnode = &(field[i].dis[0].node[j]);
      for (k=0; k<2; k++)
      {
        fprintf(f1,"%10.5E ",actnode->x[k]);
      }
      fprintf(f1,"%5d ",actnode->Id+1);      
      fprintf(f1,"\n");
    }
  }
  fclose(f1);

  /* I need some general information */
  if (actfield->fieldtyp==levelset)
  {
    f1 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_general","w");
    fprintf(f1,"%5d\n",actfield->dis[0].numnp);
    fprintf(f1,"%5d\n",actfield->dis[0].numele);
    fprintf(f1,"%5d\n",lsdyn->nstep+1);
    fclose(f1);
  }
  else
  {
    dserror("\n**ERROR** in ls_to_matlab!\nno levelset field assigned!\n");    
  }

  /* print number of nodes per element */
  f1 = fopen("src/ls/to_matlab/ls_to_matlab/ls_to_matlab_nodeperelement","w");
  fprintf(f1,"%5d",field[1].dis[0].element[0].numnp);
  fclose(f1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_to_matlab */



/*!----------------------------------------------------------------------
\brief construct polygon connectivity

<pre>                                                            irhan 05/04
The polygon connectivity data is used in X-FEM module.
</pre>

*----------------------------------------------------------------------*/
void ls_polygon_con(
  INT *polycon,
  INT snode,
  INT numtri
  )
{
  INT     nd1,nd2;
  
#ifdef DEBUG 
  dstrc_enter("ls_polygon_con");
#endif
/*----------------------------------------------------------------------*/

  if (numtri==-1)
  {
    polycon[0] = snode;
    polycon[1] = node_to_node_tri[0][snode-1];
    polycon[2] = node_to_node_tri[1][snode-1];
  }
  else
  {
    /* set pivot point */
    polycon[0] = node_to_node[numtri-1][snode-1];
    /* access to local node numbers following snode */
    nd1 = node_to_node_tri[0][snode-1];
    nd2 = node_to_node_tri[1][snode-1];  
    polycon[1] = node_to_node[numtri-1][nd1-1];
    polycon[2] = node_to_node[numtri-1][nd2-1];
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_polygon_con */



/*!----------------------------------------------------------------------
\brief polygonize

<pre>                                                            irhan 05/04
In this subroutine polygonization of the element domain is performed. This
polygon structure is used in refined integration which is necessary for elements
cut by the interface. There are seven subpolygons (triangles) for each triangular
domain, if it is cut by the interface. The coordinates and weights of the Gauss
points for modified integration are also computed by this subroutine call. 
</pre>

*----------------------------------------------------------------------*/
void ls_polygonize()
{
  INT              i;
  ELEMENT         *actele;
  XFEMPOLYFLAG     xfempolyflag;
  
#ifdef DEBUG 
  dstrc_enter("ls_polygonize");
#endif
/*----------------------------------------------------------------------*/

  /* loop over the elements */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==1)
    {
      /* initialize polygon */
      xfempolyflag = xfem_poly_initialize;
      xfem_polygon(xfempolyflag,actele);
      /* construct polygon */    
      xfempolyflag = xfem_poly_construct;
      xfem_polygon(xfempolyflag,actele);
      /* set material index array */    
      xfempolyflag = xfem_poly_material;
      xfem_polygon(xfempolyflag,actele);
      /* compute Gauss points and weights for each subpolygon */        
      xfempolyflag = xfem_poly_computeGP;
      xfem_polygon(xfempolyflag,actele);
      /* write polygon information to files */
      xfempolyflag = xfem_poly_write;
      xfem_polygon(xfempolyflag,actele);
    }
  }
  /* print to screen */
  printf("\n**WARNING** LEVELSET POLYGONIZATION COMPLETED!");    
  /* print to file */
  fprintf(file12,"\n**WARNING** LEVELSET POLYGONIZATION COMPLETED!\n\n\n\n");    
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_polygonize */



/*!----------------------------------------------------------------------
\brief initialize material data

<pre>                                                            irhan 05/04
In this subroutine each element is assigned to a material number. The elements
inside the closed interface are assigned to material number 2, whereas those
outside are assigned to 1. The elements located on the interface are assigned
to material number 3. Third material is used only when the enrichment is off,
i.e., it is used for the formulation where material properties are smeared
around the interface.
</pre>

*----------------------------------------------------------------------*/
void ls_init_material()
{
  INT          i,j;
  ELEMENT     *actele;
  
#ifdef DEBUG 
  dstrc_enter("ls_init_material");
#endif
/*----------------------------------------------------------------------*/

  /* loop over the elements */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==1)
    {
      actele->e.ls2->my_fluid->mat = 3;
    }
    else
    {
      /* access to the nodal values of level set profile */
      ls2_calset1(actele,1,lset01);
      if (lset01[0]<0.0)
      {
        /* check */
        for (j=0; j<actele->numnp; j++)
        {
          if (lset01[j]>0.0)
          {
            dserror("\n**ERROR** in ls_init_material()");
          }
        }
        actele->e.ls2->my_fluid->mat = 2;
      }
      else
      {
        /* check */
        for (j=0; j<actele->numnp; j++)
        {
          if (lset01[j]<0.0)
          {
            dserror("\n**ERROR** in ls_init_material()");
          }
        }
        actele->e.ls2->my_fluid->mat = 1;        
      }
    }
  }
  /* print to screen */
  printf("\n**WARNING** LEVELSET MATERIAL INITIALIZATION COMPLETED!");    
  /* print to file */
  fprintf(file12,"\n**WARNING** LEVELSET MATERIAL INITIALIZATION COMPLETED!");    
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_init_material */



/*!----------------------------------------------------------------------
\brief update material data

<pre>                                                            irhan 05/04
In this subroutine material data is updated. This update is necessary, because
interface change position time. This update process can be confined to a region
around the interface.
</pre>

*----------------------------------------------------------------------*/
void ls_updt_material()
{
  INT          i,j;
  ELEMENT     *actele;
  
#ifdef DEBUG 
  dstrc_enter("ls_updt_material");
#endif
/*----------------------------------------------------------------------*/

  /* loop over the elements */
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->e.ls2->is_elcut==0 && actele->e.ls2->is_elsearch==0)
    {
      /* no update */
      continue;
    }
    else if (actele->e.ls2->is_elcut==1)
    {
      actele->e.ls2->my_fluid->mat = 3;
    }
    else if (actele->e.ls2->intcnt==1)
    {
      /* access to the nodal values of level set profile */
      ls2_calset1(actele,1,lset01);
      if (lset01[0]<0.0)
      {
        /* check */
        for (j=0; j<actele->numnp; j++)
        {
          if (lset01[j]>0.0)
          {
            dserror("\n**ERROR** in ls_updt_material()");
          }
        }
        actele->e.ls2->my_fluid->mat = 2;
      }
      else
      {
        /* check */
        for (j=0; j<actele->numnp; j++)
        {
          if (lset01[j]<0.0)
          {
            dserror("\n**ERROR** in ls_updt_material()");
          }
        }
        actele->e.ls2->my_fluid->mat = 1;
      }
    }
  }
  /* print to screen */
  printf("\n**WARNING** LEVELSET MATERIAL UPDATE COMPLETED!");    
  /* print to file */
  fprintf(file12,"\n**WARNING** LEVELSET MATERIAL UPDATE COMPLETED!");    
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;					  
} /* end of ls_updt_material */



/*!----------------------------------------------------------------------
\brief check level set profile

<pre>                                                            irhan 05/04
In this subroutine nodal values of the levelset are checked. The nodes lying
on the interface are detected. Then they are assigned to THRESHOLD value
different from zero, therefore they are not on the interface anymore. This
perturbation is necessary because while developing the code it is assumed that
each line segment constructing the interface falls inside the element domain.
</pre>

*----------------------------------------------------------------------*/
void ls_check_profile()
{
  INT         i;
  DOUBLE      val;
  ELEMENT    *actele;
  NODE       *actnode;
  
#ifdef DEBUG
  dstrc_enter("ls_check_profile");
#endif
/*----------------------------------------------------------------------*/
  
  /* loop */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    val = actnode->sol_increment.a.da[1][0];
    /* check */
    if (fabs(val)<EPS5)
    {
      actnode->sol_increment.a.da[1][0] = NODETHRESHOLD;
      fprintf(file12,"\n\n**WARNING** DISTANCE FUNCTION FOR NODE %4d MODIFIED!",
              actnode->Id+1);
      ls_printnodeinfo_to_file(actnode);
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_check_profile */
/*! @} (documentation module close)*/
#endif
