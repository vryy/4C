/*!----------------------------------------------------------------------
\file
\brief chimera_autohole.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_CHIMERA


/*!
\addtogroup CHIMERA
*//*! @{ (documentation module open)*/



#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "chimera.h"





/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD                *field;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR               *solv;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION            *partition;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR                   par;
extern struct _CHIMERA_DATA         *chm_data;



/*!----------------------------------------------------------------------
\brief perform hole cutting automatically

<pre>                                                            irhan 09/04
In this subroutine the inner boundary of the background discretization
is automatically generated
</pre>

*----------------------------------------------------------------------*/
void chimera_automatic_hole_cutting(
  INT     n_back
  )
{
  INT              rr,mm,kk,ll,nn;    /* Laufvariablen */
  INT             *Lochknoten;
  DOUBLE           x,y;               /* tempor"ar */
  INT              n_Loecher;
  INT              chimera_parenttyp;
  INT              Randschichtdicke;
  struct _NODE    *Knoten;
  struct _ELEMENT *elm,  *parent;
  INT              indikator,Indikator_randnaehe,counter,dirich;
  
  INT              numeq;
  INT              numeq_total;
  
  FIELD           *actfield;      /* the active field */
  PARTITION       *actpart;       /* my partition of the active field */
  SOLVAR          *actsolv;       /* the active SOLVAR */
  INTRA           *actintra;      /* the field's intra-communicator */

#ifdef DEBUG
  dstrc_enter("chimera_automatic_hole_cutting");
#endif
/*----------------------------------------------------------------------*/

  /* Hauptprogramm */
  /* Festlegung der Randschichtdicke (determining thickness of boundary layer) */
  Randschichtdicke=chm_data[0].overlapping_length;
  /* Speicherplatz fuer L"ocherliste etc (allocate some memory) */
  n_Loecher = 0;
  Lochknoten = (INT*)malloc(field[0].dis[n_back].numnp*sizeof(INT));

  /* Initialisierung */
  for (rr=0; rr<field[0].dis[n_back].numnp; rr++)
  {
    Knoten = &field[0].dis[n_back].node[rr];
    Knoten[0].gnode[0].Knotentyp = frei;
  }

  /*
    preparations to be able to decide whether a calculated boundary
    is an inner boundary (two boundarys will be found, only one is
    correct)
  */
  /*
     Festlegung des Startbereiches zur Ermittlung des inneren Randes
     auf der Elementebene
  */
  for (rr=0;rr<field[0].dis[1-n_back].numnp;rr++) /* Alle Punkte aus dem Vordergrund */
  {
    Knoten = &field[0].dis[1-n_back].node[rr];
    indikator = 0;
    x = Knoten[0].x[0];
    y = Knoten[0].x[1];
    if (Knoten[0].gnode[0].chi_bndtype != empty)
    {
      parent = chimera_search(x,y,&field[0].dis[n_back],n_back,0.0001,1);
      parent[0].chimera_parenttyp = Kandidat_Innerer_Rand;
    }
  }

  /*
     every node of the background below the patch is marked as "Loch"
     and added to the list of holes
  */
  /*
    Allen Knoten der Hintergrunddiskretisierung die unterhalb der
    Oberdiskretisierung liegen wird die Eigenschaft Loch zugewiesen
    und sie werden in die L"ocherliste mit aufgenommen
  */
  for (rr=0; rr<field[0].dis[n_back].numnp; rr++) /* Alle Knoten */
  {
    x = field[0].dis[n_back].node[rr].x[0];
    y = field[0].dis[n_back].node[rr].x[1];
    parent = chimera_search(x,y,&field[0].dis[1-n_back],1-n_back,0.0001,1);
    if (parent != NULL)
    {
      field[0].dis[n_back].node[rr].gnode[0].Knotentyp = Loch_Standard;
      Lochknoten[n_Loecher] = rr;
      n_Loecher = n_Loecher+1;
    }
  }

  /*
     every node with the property "Loch_Standard" and an empty "frei"
     neighbour is marked as a boundary-layer node "Randschicht=1", if
     the corresponding parent element is "near" the chimera-boundary.
     In the end, the inner part of the hole is filled.

     Hat ein Knoten mit der Eigenschaft Loch_Standard  einen Nachbarn,
     der kein Loch ist, so wird ihm die zus"atzliche Eigenschaft
     Randschicht 1 zugewiesen. Dies geschieht nur, falls das zugeh"orige
     Element auf einem Chimera-Rand lag. Am Ende wird das verbleibende
     Loch aufgef"ullt.
  */
  for (mm=0; mm<n_Loecher; mm++) /* Schleife ueber alle Lochknoten */
  {
    indikator = 0;
    Indikator_randnaehe = 0;
    Knoten = &field[0].dis[n_back].node[Lochknoten[mm]];
    /* Schleife "uber alle umgebenden Elemente (loop over all surrounding elements) */
    for (kk=0; kk<Knoten[0].numele; kk++)
    {
      elm = &Knoten[0].element[kk][0];
      /* Aufweichen der Forderung auf Nachbarn */
      /* parent element has to be "near", e.g. neighboring a chimera_boundary */
      for (ll=0; ll<elm[0].numnp; ll++)
      {
        for (nn=0; nn<elm[0].node[ll][0].numele; nn++)
        {
          if (elm[0].node[ll][0].element[nn][0].chimera_parenttyp == Kandidat_Innerer_Rand)
          {
            Indikator_randnaehe = 1;
          }
        }
        if (elm[0].node[ll][0].gnode[0].Knotentyp == frei)
        {
          indikator = 1;
        }
      }
    }

    /* Zuweisen der Randschichteigenschaft */
    /* mark as boundary layer */
    if (indikator==1 && Indikator_randnaehe==1)
    {
      Knoten[0].gnode[0].Randmarker = 1;
    }
  }

  /* Schleife "uber Dicke der Randschicht*/
  /* (loop over the thickness of the overlapping region) */
  for (rr=0; rr<Randschichtdicke-1; rr++)
  {
    /*
      Erkl"are Knoten der Randschicht zu Interpolationsknoten oder setze ihre
      Eigenschaft von Loch auf frei. (Bis Uberlappungsdicke erreicht)
    */
    /*
      the nodes of the boundary-layer are set as interpolation nodes if
      the size of the overlapping region is reached. Otherwise they are
      set as "frei" and the process is repeated.
    */
    for (mm=0; mm<n_Loecher; mm++) /* Schleife ueber alle Lochknoten */
    {
      Knoten = &field[0].dis[n_back].node[Lochknoten[mm]];
      if (Knoten[0].gnode[0].Randmarker==1 && Knoten[0].gnode[0].Knotentyp==Loch_Standard)
      {
        if (rr>0)
        {
          Knoten[0].gnode[0].Knotentyp = frei_unter_patch;
        }
        else
        {
          Knoten[0].gnode[0].Knotentyp = frei;
        }
      }
    }
    for (mm=0; mm<n_Loecher; mm++) /* Schleife ueber alle Lochknoten */
    {
      Knoten = &field[0].dis[n_back].node[Lochknoten[mm]];
      if (Knoten[0].gnode[0].Knotentyp == Loch_Standard)
      {
        indikator = 0;
        Knoten = &field[0].dis[n_back].node[Lochknoten[mm]];
        /* Schleife "uber alle umgebenden Elemente */
        for (kk=0; kk<Knoten[0].numele; kk++)
        {
          elm = &Knoten[0].element[kk][0];
          /* Aufweichen der Forderung auf Nachbarn */
          for (ll=0; ll<elm[0].numnp; ll++)
          {
            if ((elm[0].node[ll][0].gnode[0].Knotentyp==frei ||
                 elm[0].node[ll][0].gnode[0].Knotentyp==frei_unter_patch ) &&
                elm[0].node[ll][0].gnode[0].Randmarker==1)
            {
              indikator = 1;
            }
          }
        }
        /* Zuweisen der Randschichteigenschaft */
        if (indikator == 1)
        {
          Knoten[0].gnode[0].Randmarker = 1;
        }
      }
    }
  }
  
  for (mm=0; mm<n_Loecher; mm++) /* Schleife ueber alle Lochknoten */
  {
    Knoten = &field[0].dis[n_back].node[Lochknoten[mm]];
    if (Knoten[0].gnode[0].Randmarker==1 && Knoten[0].gnode[0].Knotentyp==Loch_Standard)
    {
      if (Knoten[0].gnode[0].dirich==NULL)
      {
        Knoten[0].gnode[0].dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",&(Knoten[0].gnode[0].dirich->dirich_onoff),3,1,"IV");
        amdef("val",&(Knoten[0].gnode[0].dirich->dirich_val),3,1,"DV");
        amdef("curve",&(Knoten[0].gnode[0].dirich->curve),3,1,"IV");
        amzero(&(Knoten[0].gnode[0].dirich->dirich_onoff));
        amzero(&(Knoten[0].gnode[0].dirich->dirich_val));
        amzero(&(Knoten[0].gnode[0].dirich->curve));
        /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
        for (kk=0; kk<Knoten->numdf-1; kk++)
        {
          Knoten[0].gnode[0].dirich->dirich_onoff.a.iv[kk] = 1;
          Knoten[0].gnode[0].dirich->dirich_val.a.dv[kk] = 0;
        }
        Knoten[0].gnode[0].dirich->dirich_onoff.a.iv[Knoten->numdf-1] = 0;
        Knoten[0].gnode[0].chi_bndtype = Chimera_Dirichlet;
        Knoten[0].gnode[0].Knotentyp = Lochrand;
        if (par.myrank==0)
        {
          printf("Innerer Randknoten (auto):(%5f,%5f)\n",Knoten[0].x[0],Knoten[0].x[1]);
        }
      }
      else
      {
        if (par.myrank==0)
        {
          Knoten[0].gnode[0].chi_bndtype = Chimera_Dirichlet;
          /* Test */
          printf("Aeusserer Randknoten: (%5f,%5f)\n",Knoten[0].x[0],Knoten[0].x[1]);
        }
      }
    }
    else if (Knoten[0].gnode[0].Randmarker==0 && Knoten[0].gnode[0].Knotentyp==Loch_Standard)
    {
      /* ... */
    }
  }

  /* Auff"ullen des inneren Loches */
  /* fill the hole */
  indikator = 0;
  while (indikator == 0)
  {
    indikator = 1;
    for (mm=0; mm<field[0].dis[n_back].numnp; mm++) /* Schleife ueber alle Knoten */
    {
      Knoten = &field[0].dis[n_back].node[mm];
      
      if (Knoten[0].gnode[0].Knotentyp==Loch_Standard && Knoten[0].gnode[0].Randmarker==0)
      {
        /* Schleife "uber alle umgebenden Elemente */
        for (kk=0; kk<Knoten[0].numele; kk++)
        {
          elm = &Knoten[0].element[kk][0];
          for (ll=0; ll<elm[0].numnp; ll++)
          {
            if (elm[0].node[ll][0].gnode[0].Knotentyp == frei)
            {
              elm[0].node[ll][0].gnode[0].Knotentyp = Loch_Standard;
              indikator = 0;
            }
          }
        }
      }
    }
  }

  /* Ab hier Vorbereitung f"ur Rechnung */
  /* renumbering of degrees of freedom */
  counter = 0;
  for (rr=0; rr<field[0].dis[n_back].numnp; rr++) /* loop over nodes */
  {
    Knoten = &(field[0].dis[n_back].node[rr]);
    /* the node does not have any conditions */
    if (Knoten->gnode->couple==NULL && Knoten->gnode->dirich==NULL)
    {
      for (ll=0; ll<Knoten->numdf; ll++)
      {
        Knoten->dof[ll] = counter;
        counter++;
      }
    }

    /* the node does have conditions */
    /*
      Speicher f"ur Dirichletrandbedingungen --- Interpolationsknoten m"ussen
      als Randknoten vorbereitet werden
    */
    /*
      allocate memory for the Dirichlet-conditions --- prepare
      interpolation nodes as dirichlet-nodes
    */
    else
    {
      for (ll=0; ll<Knoten->numdf; ll++)
      {
        dirich = 0;
        /* dof has dirichlet condition */
        if (Knoten->gnode->dirich!=NULL && Knoten->gnode->dirich->dirich_onoff.a.iv[ll]!=0)
          dirich = 1;
        else
        {
          Knoten->dof[ll] = counter;
          counter++;
        }
        if (dirich == 1)
        {
          Knoten->dof[ll] = -1;
        }
      }
    }
  }

  /* Anzahl Gleichungen */
  field[0].dis[n_back].numeq = counter;

  /*
    Now all free dofs are numbered, so now number the dirichlet conditioned
    dofs from here on
  */
  for (rr=0; rr<field[0].dis[n_back].numnp; rr++)
  {
    Knoten = &(field[0].dis[n_back].node[rr]);
    if (Knoten->gnode->couple==NULL && Knoten->gnode->dirich==NULL) continue;
    for (ll=0; ll<Knoten->numdf; ll++)
    {
      if (Knoten->dof[ll] == -1)
      {
        Knoten->dof[ll] = counter;
        counter++;
      }
    }
  }
  field[0].dis[n_back].numdf = counter;
  /* Jetzt noch der sparse-K"ase */
  /* initialise sparse-arrays */
  actfield = &(field[0]);
  actsolv = &(solv[0]);
  actpart = &(partition[0]);
#ifdef PARALLEL
  actintra = &(par.intra[0]);
#else
  /* if we are not parallel here, we have to allocate a pseudo-intracommunicator */
  actintra = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  actintra->intra_fieldtyp = actfield->fieldtyp;
  actintra->intra_rank = 0;
  actintra->intra_nprocs = 1;
#endif
  mask_msr(actfield,actpart,actsolv,actintra,actsolv->sysarray[n_back].msr,n_back);

  /* Aufrauemen */
  /* cleaning */
  free(Lochknoten);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  
  return;
} /* end of chimera_automatic_hole_cutting */
/*! @} (documentation module close)*/
#endif
