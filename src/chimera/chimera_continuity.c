/*!----------------------------------------------------------------------
\file
\brief chimera_continuity.c

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
#include "chimera.h"





/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
  *----------------------------------------------------------------------*/
extern struct _GENPROB               genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD                *field;
extern struct _CHIMERA_DATA         *chm_data;
/*----------------------------------------------------------------------+
| This structure struct _PAR par; is defined in main_ccarat.c           |
| and the type is in partition.h                                        |
*----------------------------------------------------------------------*/
extern struct _PAR                  par;





/*!----------------------------------------------------------------------
\brief continuity check

<pre>                                                            irhan 09/04
Berechnet zu gegebenem Teilgebiet den Fluss "uber den inneren Rand.

Durchgef"uhrte Schritte:
- Triangulierungstyo erkennen
- Alle Randelemente abklappern
  * "au"sere Einheitsnormale berechnen
  * Normalprojektion der Geschwindigkeiten in den Ecken berechnen
  * Flu"s "uber Kante auswerten und aufaddieren

Mittlerweile redundant.
</pre>

*----------------------------------------------------------------------*/
INT chimera_conti_check(
  INT     n_dis
  )
{
  INT          rr,mm,kk;
  INT          n_ele_coupling_node, coupling_nodes[4],frei,indi;
  DISCRET     *actdis;
  ELEMENT     *actele;
  DOUBLE       tempvec[2],temp,norm;
  DOUBLE       unit_normal[2][2];
  DOUBLE       total_flux=0;

#ifdef DEBUG
  dstrc_enter("chimera_conti_check");
#endif
/*----------------------------------------------------------------------*/

  /* Initialisierung */
  actdis = &(field[genprob.numff].dis[n_dis]);

  /* Eigentliches Programm */

  /* Rechteckselemente*/
  if (actdis[0].element[0].distyp == quad4)
  {
    /* Loop over all elements */
    for (rr=0; rr<actdis[0].numele; rr++)
    {
      n_ele_coupling_node = 0;
      frei = 0;
      actele = &(actdis[0].element[rr]);
      indi = 0;
      for (mm=0; mm<actele[0].numnp; mm++)
      {
        if(actele[0].node[mm][0].gnode[0].Knotentyp==Loch_Standard ||
           actele[0].node[mm][0].gnode[0].Knotentyp==Lochrand)
        {
          indi++;
        }
      }
      if(indi <= 3)
      {
        for (mm=0; mm<actele[0].numnp; mm++)
        {
          if (actele[0].node[mm][0].gnode[0].chi_bndtype != empty)
          {
            /* in this case, we have a coupling node */
            coupling_nodes[n_ele_coupling_node] = mm;
            n_ele_coupling_node++;
          }
          else
          {
            frei = mm;
          }
        }
        if (n_ele_coupling_node == 2)
        {
          /* calculation of outward normal */
          temp = 0;
          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk];
            temp += tempvec[kk]*tempvec[kk];
          }
          temp = sqrt(temp);
          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = tempvec[kk]/temp;
          }
          unit_normal[0][0] = -tempvec[1];
          unit_normal[0][1] = tempvec[0];

          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = actele[0].node[frei][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk];
          }
          temp = 0;
          for (kk=0; kk<2; kk++)
          {
            temp += tempvec[kk]*unit_normal[0][kk];
          }
          if (temp > 0)
          {
            for (kk=0; kk<2; kk++)
            {
              unit_normal[0][kk] = -unit_normal[0][kk];
            }
          }
          /* Calculation of flux-contribution (assuming linear shape-fuctions) */
          for (kk=0; kk<2; kk++)
          {
            temp = actele[0].node[coupling_nodes[0]][0].sol_increment.a.da[1][kk]*
                   unit_normal[0][kk]+actele[0].node[coupling_nodes[1]][0].sol_increment.a.da[1][kk]*
                   unit_normal[0][kk];
          }
          norm = 0;
          for (kk=0; kk<2; kk++)
          {
            norm += (actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk])*
                    (actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk]);
          }
          norm = sqrt(norm);
          total_flux += 0.5*temp*norm;
        }
        if (n_ele_coupling_node == 3)
        {
          dserror("Flussberechnung fuer quad4 mit 3 Randpunkten noch nicht implementiert!\n");
          /* Calculation of outward normal (twice) */
          /* Calculation of flux-contributions (assuming linear shape-fuctions) */
        }
      }
    }
  }

  /* Dreieckselemente */
  if (field[genprob.numff].dis[n_dis].element[0].distyp == tri3)
  {
    for (rr=0; rr<actdis[0].numele; rr++)
    {
      n_ele_coupling_node = 0;
      frei = 0;
      actele = &(actdis[0].element[rr]);
      indi = 0;
      for (mm=0; mm<actele[0].numnp; mm++)
      {
        if(actele[0].node[mm][0].gnode[0].Knotentyp==Loch_Standard ||
           actele[0].node[mm][0].gnode[0].Knotentyp==Lochrand)
        {
          indi++;
        }
      }
      if(indi <= 2)
      {
        for (mm=0; mm<actele[0].numnp; mm++)
        {
          if (actele[0].node[mm][0].gnode[0].chi_bndtype != empty)
          {
            /* in this case, we have a coupling node */
            coupling_nodes[n_ele_coupling_node] = mm;
            n_ele_coupling_node++;
          }
          else
          {
            frei = mm;
          }
        }
        if (n_ele_coupling_node == 2)
        {
          /* Calculation of outward normal */
          temp = 0;
          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk];
            temp += tempvec[kk]*tempvec[kk];
          }
          temp = sqrt(temp);
          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = tempvec[kk]/temp;
          }
          unit_normal[0][0] = -tempvec[1];
          unit_normal[0][1] = tempvec[0];

          for (kk=0; kk<2; kk++)
          {
            tempvec[kk] = actele[0].node[frei][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk];
          }
          temp = 0;
          for (kk=0; kk<2; kk++)
          {
            temp += tempvec[kk]*unit_normal[0][kk];
          }
          if (temp > 0)
          {
            for (kk=0; kk<2; kk++)
            {
              unit_normal[0][kk] = -unit_normal[0][kk];
            }
          }

          /* Calculation of flux-contribution (assuming linear shape-fuctions) */
          for (kk=0; kk<2; kk++)
          {
            temp = actele[0].node[coupling_nodes[0]][0].sol_increment.a.da[1][kk]*
                   unit_normal[0][kk]+actele[0].node[coupling_nodes[1]][0].sol_increment.a.da[1][kk]*
                   unit_normal[0][kk];
          }
          norm = 0;
          for (kk=0; kk<2; kk++)
          {
            norm += (actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk])*
                    (actele[0].node[coupling_nodes[1]][0].x[kk]-actele[0].node[coupling_nodes[0]][0].x[kk]);
          }
          norm = sqrt(norm);
          total_flux += 0.5*temp*norm;
        }
      }
    }
  }
  if (par.myrank == 0)
  {
    printf("  |      |  Gesamtfluss ueber Inneren Rand: %5g\n",total_flux);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return(0);
} /* end of chimera_conti_check */





/*!----------------------------------------------------------------------
\brief continuity interpolation

<pre>                                                            irhan 09/04

This function modifies the boundary values which were calculated
by using Lagrange-interpolation. The new boundary values are the
best approximation of the Lagrange-values under the constraint of
mass conservation.

Linear shape functions are assumed.

The solution of the least-square problem uses orthogonal transformation.


Modifiziert die aus Lagrange-Interpolation gewonnenen Randwerte auf dem
inneren Rand durch eine least-square-Approximation unter der Nebenbedingung
der Kontinuit"atserhaltung der numerischen L"osung.

Es werden lineare Ansatzfunktionen vorrausgesetzt.

Die L"osung des linearen Ausgleichsproblems erfolgt mit Hilfe von
Orthogonaltransformationen (Givens-Rotationen)

</pre>

*----------------------------------------------------------------------*/
void chimera_continuity_interpolation(INT n_dis)
{

  INT        rr,mm,kk,nn;
  INT        n_IRP;               /* number of inner boundary-nodes */
  INT       *connect;             /* connection between inner boundary nodes an node IDs */
  INT        counter;             /* Z"ahler f"ur Interpolationsrandknoten und
                                     Nachbarinterpolationsknoten */
  INT        Nachbarn[2];         /* contains IDs of neighboring interpolation nodes */
  INT        pos;                 /* chooses column in sol-increment */

  INT        indi;                /* indicates whether element belongs to the inner part of the subdomain */

  INT        blub;
  INT        frei;                /* number of a free node to calculate outward normal */

  DOUBLE     Flusskontrolle = 0;
  DOUBLE    *Ausgangswerte;       /* values from the lagrangeinterpolation */

  /*                       --- --- ------ --------- ---------      */
  /* Komponentenanordnung |u_1|v_1|......|u_{n_IRP}|v_{n_IRP}|     */
  /*                       --- --- ------ --------- ---------      */
  /*                                                               */
  /*                      |-------------2*n_IRP--------------|     */

  DOUBLE    *angepasste_Werte;    /* calculated values */
                                  /* Komponentenanordnung wie oben */

  DOUBLE    *lS_Zb,*lS_Zb_alt;    /* left hand side of linear constraint */
  DOUBLE   **Q;                   /* transformation-matrix */
  DOUBLE    *Q_Zeile;             /* temp-vetor */
  DOUBLE     cphi,sphi;           /* angles for the Givens-rotations */

  DOUBLE     Nenner;              /* Nenner der Ausdr"ucke f"ur cphi, sphi */
                                  /* Zur Berechnung der "au"seren Einheitsnormale */
  DOUBLE     tempvec[2],temp;
  DOUBLE     norm[2];             /* length of the edge to the neighbouring interpolation node */

  DOUBLE     unit_normal[2][2];   /* unit normal; first index: contains number of neighbour (specifies edge) */

  DISCRET   *actdis;
  ELEMENT   *actele;
  DOUBLE     TOL;
  INT        next;                /* next element to be eliminated */
  INT        next_us;             /* next element used to eliminate element next */

  DOUBLE     steady;

#ifdef DEBUG
  dstrc_enter("chimera_continuity_interpolation");
#endif
/*----------------------------------------------------------------------*/

  /*****************/
  /* PREPROCESSING */
  /*****************/
  if (chm_data[0].conti_interpolation == 1)
  {
    /* Initialisierung */
    actdis = &(field[genprob.numff].dis[n_dis]);
    n_IRP = 0;

    /* Ermittlung ben"otigter Gr"o"sen */
    for (rr=0; rr<actdis[0].numnp; rr++)
    {
      /* Schleife "uber alle Punkte zur Ermittlung der Randpunktnummern */
      if (actdis[0].node[rr].gnode[0].chi_bndtype == Chimera_Dirichlet)
      {
        n_IRP++;
      }
    }

    /* Speicherorganisation */
    connect = (INT*)malloc(n_IRP*sizeof(INT));
    Ausgangswerte = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));
    angepasste_Werte = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));

    for (rr=0; rr<2*n_IRP; rr++)
    {
      angepasste_Werte[rr] = 0;
    }
    lS_Zb = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));
    lS_Zb_alt = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));
    Q = (DOUBLE**)malloc(2*n_IRP*sizeof(DOUBLE));
    for(rr=0; rr<2*n_IRP; rr++)
    {
      Q[rr] = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));
      for (mm=0; mm<2*n_IRP; mm++)
      {
        Q[rr][mm] = 0;
      }
      Q[rr][rr] = 1;
    }
    Q_Zeile = (DOUBLE*)malloc(2*n_IRP*sizeof(DOUBLE));

    /* connect nodes and calculation array */
    counter = 0;
    for (rr=0; rr<actdis[0].numnp; rr++)
    {
      if (actdis[0].node[rr].gnode[0].chi_bndtype == Chimera_Dirichlet)
      {
        connect[counter] = actdis[0].node[rr].Id_loc;
        counter++;
      }
    }

    /* read values into calculation array */
    /* Interpolierte Werte */
    for (rr=0; rr<n_IRP; rr++)
    {
      /* Alle Randknoten werden abgeklappert */
      /* loop over all boundary nodes */
      /* Zeitschritt n+1 */
      pos = 3;
      for (mm=0; mm<2; mm++)
      {
        /* loop over velocity-degrees of freedom */
        Ausgangswerte[2*rr+mm] = actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm];
      }
    }

    /* Linke Seite der Nebenbedingung */
    for (rr=0; rr<n_IRP; rr++) /* rr hier Knotenindex */
    {
      /* zwei Eintr"age f"ur Knoten rr */
      /* Ermittlung der benachbarten Randknoten --- nur Koordinaten und
         Anzahl sind interessant */
      counter = 0;
      for (mm=0; mm<actdis[0].node[connect[rr]].numele; mm++)
      {
        /* Schleife "uber alle umgebenden Elemente */
        actele = &(actdis[0].node[connect[rr]].element[mm][0]);
        /* Kontrolle, ob es sich um ein Inneres Element handelt */
        /* control whether element is situated in the interior of the domain */
        indi = 0;
        for (nn=0; nn<actele[0].numnp; nn++)
        {
          if (actele[0].node[nn][0].gnode[0].Knotentyp==Loch_Standard ||
              actele[0].node[nn][0].gnode[0].Knotentyp==Lochrand)
          {
            indi++;
          }
        }

        if (indi < actele[0].numnp)
        {
          for (kk=0; kk<actele[0].numnp; kk++)
          {
            /* Schleife "uber Knoten des Elements */
            if (actele[0].distyp == quad4)
            {
              if (actele[0].node[kk][0].gnode[0].chi_bndtype == Chimera_Dirichlet)
              {
                TOL = 0.000000001;
                if(actele[0].node[kk][0].Id_loc != connect[rr])
                {
                  if ((actele[0].node[kk][0].x[0]-actdis[0].node[connect[rr]].x[0]<TOL &&
                       actele[0].node[kk][0].x[0]-actdis[0].node[connect[rr]].x[0]>-TOL )||
                      (actele[0].node[kk][0].x[1]-actdis[0].node[connect[rr]].x[1]<TOL &&
                       actele[0].node[kk][0].x[1]-actdis[0].node[connect[rr]].x[1]>-TOL))
                  {
                    Nachbarn[counter] = actele[0].node[kk][0].Id_loc;
                    counter++;
                  }
                }
              }
            }
            else if (actele[0].distyp == tri3)
            {
              if (actele[0].node[kk][0].gnode[0].chi_bndtype == Chimera_Dirichlet)
              {
                if(actele[0].node[kk][0].Id_loc != connect[rr])
                {
                  Nachbarn[counter] = actele[0].node[kk][0].Id_loc;
                  counter++;
                }
              }
            }
            else
            {
              dserror("elementtyp not implemented in chimera_continuity");
            }
          }
        }
      }

      /* Berechnung der Normalenvektoren */
      /* calculation of normal vectors */
      for(mm=0; mm<counter; mm++)
      {
        /* F"ur die Kante zum mm-ten Nachbarn */
        /* for the edge to the neighbor mm */
        temp = 0;
        for(kk=0; kk<2; kk++)
        {
          tempvec[kk] = actdis[0].node[Nachbarn[mm]].x[kk]-actdis[0].node[connect[rr]].x[kk];
          temp += tempvec[kk]*tempvec[kk];
        }
        norm[mm] = sqrt(temp);/* Entspricht Berechnung der Kantenl"angen */
        temp = 1.0/norm[mm];
        unit_normal[mm][0] = -tempvec[1]*temp;
        unit_normal[mm][1] = tempvec[0]*temp;
        /* Jetzt mu"s die Orientierung noch angepa"st werden */
        /* Suche dazu "inneres" Element bzw. "inneren" Punkt */
        /* select orientation (outward normal) */
        /* for that purpose, an inner element is searched */
        for (kk=0; kk<actdis[0].node[connect[rr]].numele; kk++)
        {
          actele = &(actdis[0].node[connect[rr]].element[kk][0]);
          indi = 0;
          blub = 0;
          for (nn=0; nn<actele[0].numnp; nn++)
          {
            if (actele[0].node[nn][0].gnode[0].Knotentyp==Loch_Standard ||
                actele[0].node[nn][0].gnode[0].Knotentyp==Lochrand)
            {
              indi++;
            }
            else if (actele[0].node[nn][0].gnode[0].chi_bndtype != Chimera_Dirichlet)
            {
              frei = actele[0].node[nn][0].Id_loc;
            }
            if (actele[0].node[nn][0].Id_loc == Nachbarn[mm])
            {
              blub = 1;
            }
          }
          if (indi<actele[0].numnp && blub==1)
          {
            /* Das Element hat inneren Punkt frei und der Nachbar geh"ort auch dazu */
            break;
          }
        }
        /* Berechne dann Skalarprodukt aus Normalenvektor und Vektor vom
           Knoten zum inneren Punkt */
        for (nn=0; nn<2; nn++)
        {
          tempvec[nn] = actdis[0].node[frei].x[nn]-actdis[0].node[connect[rr]].x[nn];
        }
        temp = 0;
        for (nn=0; nn<2; nn++)
        {
          temp += tempvec[nn]*unit_normal[mm][nn];
        }
        if (temp>0)
        {
          for (nn=0; nn<2; nn++)
          {
            unit_normal[mm][nn] = -unit_normal[mm][nn];
          }
        }
      }

      for (mm=0; mm<2; mm++)
      {
        lS_Zb[2*rr+mm] = 0;
      }

      for (kk=0; kk<counter; kk++) /* kk hier Kantenindex */
      {
        for (mm=0; mm<2; mm++)     /* mm hier Raumrichtungsindex */
        {
          lS_Zb[2*rr+mm] = lS_Zb[2*rr+mm]+norm[kk]*unit_normal[kk][mm];
          lS_Zb_alt[2*rr+mm] = lS_Zb[2*rr+mm];
        }
      }
    }

    /* Berechnung des Gesamtflusses (zu Kontrollzwecken) */
    /* compute flux for control purpose */
    for (rr=0; rr<n_IRP; rr++)
    {
      for (mm=0; mm<2; mm++)
      {
        Flusskontrolle += lS_Zb[2*rr+mm]*Ausgangswerte[2*rr+mm];
      }
    }
    Flusskontrolle *= 0.5;

    if (par.myrank == 0)
    {
      printf("  |      |  Flusskontrolle fuer Lagrangewerte: %5g\n",Flusskontrolle);
    }

    /*********************************/
    /* L"OSEN DES AUSGLEICHSPROBLEMS */
    /*********************************/

    /* Aufstellen der orthogonalen Transformationsmatrix Q mit Hilfe von
       Givens-Rotationen */
    /* construction of transformation matrices. Use Givens-rotations! */
    next = 0;
    next_us = 0;
    TOL = 0.0000001;
    while (next < 2*n_IRP-1)
    {
      /* Solange muss eleminiert werden */
      /* Finde n"achsten zu eliminierenden Eintrag in der linken Seite */
      /* search for the next entry which has to be eliminated */
      while (-TOL<lS_Zb[next] && lS_Zb[next]<TOL)
      {
        next++;
      }
      /* next ist jetzt die erste Zeile mit einem Eintrag != 0 */
      /* next is a row with an entry not equal 0 */
      /* Finde darauf folgenden Eintrag!=0, der zum eliminieren verwendet
         werden kann */
      /* locate the next entry!=0. Use it to eliminate entry "next" */

      next_us = next + 1;
      while (-TOL<lS_Zb[next_us] && lS_Zb[next_us]<TOL)
      {
        next_us++;
      }
      /* Es muss nur weitergemacht werden falls noch ein weiterer Eintrag
         existiert */
      /* continue as long as there exist entries!=0 */
      if (next_us < 2*n_IRP)
      {
        /* Berechne cphi und sphi */
        /* compute cphi and sphi (this determines the rotation) */
        Nenner = lS_Zb[next]*lS_Zb[next]+lS_Zb[next_us]*lS_Zb[next_us];
        Nenner = sqrt(Nenner);
        Nenner = 1.0/Nenner;
        cphi = Nenner*lS_Zb[next_us];
        sphi = Nenner*lS_Zb[next];
        /* Bringe die Q-Zeile in Sicherheit */
        /* save one row of Q */
        for (rr=0; rr<2*n_IRP; rr++)
        {
          Q_Zeile[rr] = Q[next][rr];
        }
        /* Update der Transformationsmatrix Q */
        /* update Q */
        /* Es m"ussen nur zwei Zeilen ge"ander werden */
        /* to rows have to be modified*/
        for (rr=0; rr<2*n_IRP; rr++)
        {
          Q[next][rr] = Q_Zeile[rr]*cphi-Q[next_us][rr]*sphi;
        }
        for (rr=0; rr<2*n_IRP; rr++)
        {
          Q[next_us][rr] = Q_Zeile[rr]*sphi+Q[next_us][rr]*cphi;
        }
        /* Update der linken Seite */
        /* update the lefthandside */
        lS_Zb[next_us] = sphi*lS_Zb[next]+cphi*lS_Zb[next_us];
        lS_Zb[next] = 0;
      }
      else
      {
        break;
      }
    }

    /* Berechne das erste Element von Q*Ausgangswerte */
    /* compute first entry of  Q*Ausgangswerte */
    temp = 0;
    for (rr=0; rr<2*n_IRP; rr++)
    {
      temp += Q[next][rr]*Ausgangswerte[rr];
    }
    /* Berechne die angepassten Werte */
    /* compute the corrected values */
    for (rr=0; rr<2*n_IRP; rr++)
    {
      angepasste_Werte[rr] = Ausgangswerte[rr] - temp*Q[next][rr];
    }

    /***********************************************************/
    /* Zuweisen im Fall kontinuit"atserhaltender Interpolation */
    /***********************************************************/
    /* In the case of continuity-preserving interpolation write the corrected
       values to the nodes */

    /* Berechnung des Gesamtflusses (zu Kontrollzwecken) */
    /* compute flux */
    Flusskontrolle = 0;
    for (rr=0; rr<n_IRP; rr++)
    {
      for (mm=0; mm<2; mm++)
      {
        Flusskontrolle += lS_Zb_alt[2*rr+mm]*angepasste_Werte[2*rr+mm];
      }
    }
    Flusskontrolle *= 0.5;
    if (par.myrank == 0)
    {
      printf("  |      |  Flusskontrolle fuer angepasste Werte: %5g\n",Flusskontrolle);
      printf("  |      |  (Werden im naechsten Zeitschritt vorgegeben werden.)\n");
    }
    /* Auslesen aus Rechenarray in Werte auf den Knoten */
    /* calculated values are written to the nodes */
    for (rr=0; rr<n_IRP; rr++)
    {
      /* Alle Randknoten werden abgeklappert */
      pos = 3;            /* Zeitschritt n+1 */
      for (mm=0; mm<2; mm++)
      {
        /* Schleife "uber Geschwindigkeitsfreiheitsgrade */
        actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm] = angepasste_Werte[2*rr+mm];
      }
    }

    /*****************************************************************************/
    /* Manuelle Relaxation
       for (rr=0;rr<n_IRP;rr++)
       {
       for (mm=0;mm<2;mm++)
       {
       actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]=0.4*actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm];
       actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]=actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]+0.6*actdis[0].node[connect[rr]].sol_increment.a.da[pos+1][mm];
       }
       }*/
    /*****************************************************************************/

    /* Modifizieren des Abbruchkriteriums der Gebietsiteration */
    /* Ueberschreibt die Information aus chimera_boundary_update */
    chm_data[0].abs_err[n_dis] = 0;
    chm_data[0].rel_err[n_dis] = 0;

    /* Konvergenz"uberpr"ufung */
    for (rr=0; rr<n_IRP; rr++)
    {/* Alle Randknoten werden abgeklappert */
      pos = 3;           /* Zeitschritt n+1 */
      /* absouluter Zuwachs */
      /* absolute increment */
      steady = 0;
      temp = 0;
      for (mm=0; mm<2; mm++)
      {
        /* Schleife "uber Geschwindigkeitsfreiheitsgrade */
        /* loop over velocity-degrees of freedom */
        temp = (actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]-
                actdis[0].node[connect[rr]].sol_increment.a.da[pos+1][mm]);
        temp = temp*temp;
        steady = steady+temp;
      }
      steady = sqrt(steady);
      if (steady > chm_data[0].abs_err[n_dis])
      {
        chm_data[0].abs_err[n_dis] = steady;
      }
      /* "relativer" Zuwachs */
      /* relative increment */
      steady = 0;
      temp = 0;
      Nenner = 0;
      for (mm=0; mm<2; mm++)
      {
        /* Schleife "uber Geschwindigkeitsfreiheitsgrade */
        /* loop over velocity-degrees of freedom */
        temp = (actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]-
                actdis[0].node[connect[rr]].sol_increment.a.da[pos+1][mm]);
        temp = temp*temp;
        Nenner += actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm]*
                  actdis[0].node[connect[rr]].sol_increment.a.da[pos][mm];
        steady += temp;
      }
      if(Nenner != 0)
      {
        steady /= Nenner;
      }
      steady = sqrt(steady);
      if (steady > chm_data[0].rel_err[n_dis])
      {
        chm_data[0].rel_err[n_dis] = steady;
      }
    }

    /* Freigabe des Speichers */
    /* free memory */
    free(connect);
    free(Ausgangswerte);
    free(angepasste_Werte);
    free(lS_Zb);
    for(rr=0; rr<2*n_IRP; rr++)
    {
      free(Q[rr]);
    }
    free(Q);
    free(lS_Zb_alt);
    free(Q_Zeile);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_continuity_interpolation */
#endif
