#ifdef D_CHIMERA
/*
Peter Gamnitzer
"Anderungsdatum: 2004-01-30
"Anderungsdatum: 2004-02-16

Enth"alt Suchalgorithmen zur Bestimmung von parent-Elementen
0) brute force
1) quadtree
------------------------------------------------------
author: Peter Gamnitzer
last change on 2004-02-16

contains 7 functions

abstract: Provides search algorithms for parent-sarch:
          - brute force
          - quadtree

1] chimera_search
   --------------
   is  called from main program. Calls brute force/quadtree-search


2] chimera_search_brute_force
   --------------------------
   brute force search. Direct usage or indirect usage as a soubroutine
   to chimera_search_quadtree_search possible.

3] chimera_search_quadtree_init
   ----------------------------
   initializes quadtree-search (preprocessing part)

4] quadtree_recursion
   ------------------
   soubroutine used by chimera_search_quadtree_init to recursively
   generate the quadtree structure

5] chimera_search_quadtree_search
   ------------------------------
   requires quadtree-structure. Determines parent in quadtree-structure

6] mem_quadtree_free
   -----------------
   subroutine used by chimera_search_quadtree_free

7] chimera_search_quadtree_free
   ----------------------------
   clean up quadtree-structure (postprocessing part)


*/
#include "../headers/standardtypes.h"
#include "chimera.h"


extern struct _ELEMENT              *parent;
extern struct quadtree_node_struct  *quad_root;
extern struct _CHIMERA_DATA         *chm_data;

void quadtree_recursion (DISCRET *disrec,int nmax_pro_Blatt, struct quadtree_node_struct* Knoten);
void mem_quadtree_free(
    quadtree_node_type *Knoten
    );

/*--------------------------------------------------*/
/*Liefert parent aus quadtree oder brute-force-Suche*/
/*--------------------------------------------------*/
ELEMENT* chimera_search(
                        double      x,
                        double      y,
                        DISCRET    *dis,
			int         ndis,
                        double      TOL,
                        int         Typ_des_suchalgorithmus
                        )
{
/*Variablendeklarationen*/
    ELEMENT *parent;

#ifdef DEBUG
    dstrc_enter("chimera_search");
#endif
/*Aufruf der brute-force-Suche*/
    if (Typ_des_suchalgorithmus==0)
    {
	parent=chimera_search_brute_force(x,y,dis,TOL);
    }
/*Aufruf der quadtree-Suche*/
    if (Typ_des_suchalgorithmus==1)
    {
	parent=chimera_search_quadtree_search(x,y,dis,ndis,TOL);
    }


#ifdef DEBUG
    dstrc_exit();
#endif
    return(parent);
}

/*--------------------------------------------------*/
/* Durchsucht eine gegebene Menge von Elementen     */
/*--------------------------------------------------*/
ELEMENT* chimera_search_brute_force(
                        double      x,
                        double      y,
                        DISCRET    *dis,
                        double      TOL
                        )
{
/*Variablendeklarationen*/
    ELEMENT *sbparent;
    double   xi,zeta,temp;
    DOUBLE   N0[2],N1[2],N2[2],N3[2];
    int      i,j,rr;
    int      indi;

#ifdef DEBUG
    dstrc_enter("chimera_search_brute_force");
#endif
/****************************************************** Dreieckselemente*/
    if (dis[0].element[0].distyp == tri3)
    {
	indi=0;
	for (j=0;j<dis->numele;j++)
	{
/* loop over all background-elements*/
	    sbparent=&(dis[0].element[j]);
	    for (rr=0;rr<2;rr++){
		N0[rr]=(sbparent->node[0])->x[rr];
		N1[rr]=(sbparent->node[1])->x[rr];
		N2[rr]=(sbparent->node[2])->x[rr];
	    }
	    temp=1/((N1[0]-N0[0])*(N1[1]-N2[1])-(N1[0]-N2[0])*(N1[1]-N0[1]));
	    xi = (x*(N0[1]-N2[1])+y*(N2[0]-N0[0])+N0[0]*N2[1]-N2[0]*N0[1])*temp;
	    temp=1/((N2[0]-N0[0])*(N2[1]-N1[1])-(N2[0]-N1[0])*(N2[1]-N0[1]));
	    zeta=(x*(N0[1]-N1[1])+y*(N1[0]-N0[0])+N0[0]*N1[1]-N1[0]*N0[1])*temp;

	    if (0<=xi+TOL && xi<=1+TOL)
	    {
		if (0<=zeta+TOL && zeta<=1+TOL)
		{
		    temp =xi+zeta;
		    if (0<=temp+TOL && temp<=1+TOL)
		    {
			indi=1000;
			break;
		    }
		}
	    }
	}
	i=j;
    }
/*************************************************** Rechteckselemente*/
    if (dis[0].element[0].distyp == quad4)
    {
	indi=0;
	for (j=0;j<dis->numele;j++){/* loop over all background-elements*/
	    sbparent=&(dis->element[j]);
	    for (rr=0;rr<2;rr++){
		N0[rr]=(sbparent->node[0])->x[rr];
		N1[rr]=(sbparent->node[1])->x[rr];
		N3[rr]=(sbparent->node[3])->x[rr];
	    }
	    temp=(N1[0]-N0[0])*(N1[0]-N0[0])+(N1[1]-N0[1])*(N1[1]-N0[1]);
	    xi=((x-N0[0])*(N1[0]-N0[0])+(y-N0[1])*(N1[1]-N0[1]))/temp;
	    temp=(N3[0]-N0[0])*(N3[0]-N0[0])+(N3[1]-N0[1])*(N3[1]-N0[1]);
	    zeta=((x-N0[0])*(N3[0]-N0[0])+(y-N0[1])*(N3[1]-N0[1]))/temp;
	    if (0<=xi+TOL && xi<=1+TOL)
	    {
		if (0<=zeta+TOL && zeta<=1+TOL)
		{
		    indi=1000;
		    break;
		}
	    }
	}
	i=j;
    }
/*Plausibilit"ats"uberpr"ufung*/
    if (indi!=1000 && i==dis->numele)
    {
/*	printf("-------------------------------------error\n");*/
	/*dserror("CHIMERA_SEARCH no parent found --- probably to TOL to small!\n");*/
	sbparent=NULL;
    }
#ifdef DEBUG
    dstrc_exit();
#endif
    return(&sbparent[0]);
}

/*--------------------------------------------------*/
/* Baut quadtree fuer Diskretisierung dis auf. Ein
Zeiger auf die Wurzel wird in chimera_dyn gespeichert
(Eintrag .quadtree_root[ndis])                      */
/*--------------------------------------------------*/
void chimera_search_quadtree_init(
                        DISCRET    *dis,
			int         ndis,
                        int         nmax_pro_Blatt
                        )
{
/*--------------------------------------------------*/
/*                Variablendeklarationen            */
    int rr,mm;
    double **mEcken;
    quadtree_node_type *quad_root;

#ifdef DEBUG
    dstrc_enter("chimera_search_quadtree_init");
#endif

/*--------------------------------------------------*/
/*               Speicherorganisation               */
    mEcken=(double**)malloc(4*sizeof(double*));
    for (rr=0;rr<4;rr++){
	mEcken[rr]=(double*)malloc(2*sizeof(double));
	for(mm=0;mm<2;mm++){
/* Vorinitialisierung mit einer Ecke eines beliebigen
    Elements*/
	    mEcken[rr][mm]=dis[0].element[0].node[0][0].x[mm];
	};
    };
    quad_root=(quadtree_node_type*)malloc(sizeof(quadtree_node_type));
/*--------------------------------------------------*/
/*       Berechnung der Gr"o"se der Ausgangsbox     */
    for (rr=0;rr<dis[0].numnp;rr++)
    {
	if(dis[0].node[rr].x[0]<=mEcken[0][0])/*sucht kleinste x-Koordinate*/
	{
	    mEcken[0][0]=dis[0].node[rr].x[0];
	    mEcken[3][0]=mEcken[0][0];
	}
	if(dis[0].node[rr].x[1]<=mEcken[0][1])/*sucht kleinste auftretende y-Koordinate*/
	{
	    mEcken[0][1]=dis[0].node[rr].x[1];
	    mEcken[1][1]=mEcken[0][1];
	}
	if(dis[0].node[rr].x[0]>=mEcken[2][0])/*sucht gr"o"ste x-Koordinate*/
	{
	    mEcken[2][0]=dis[0].node[rr].x[0];
	    mEcken[1][0]=mEcken[2][0];
	}
	if(dis[0].node[rr].x[1]>=mEcken[2][1])/*sucht gr"o"ste auftretende y-Koordinate*/
	{
	    mEcken[2][1]=dis[0].node[rr].x[1];
	    mEcken[3][1]=mEcken[2][1];
	}
    };
/*--------------------------------------------------*/
/*         Festlegen der Wurzel des Baumes          */
    quad_root->n_enthaltene_Elemente = dis[0].numele;
    quad_root->node_Ecken=malloc(4*sizeof(double*));
    for (rr=0;rr<4;rr++){
      quad_root->node_Ecken[rr]=malloc(2*sizeof(double));
      quad_root->node_Ecken[rr][0]=mEcken[rr][0];
      quad_root->node_Ecken[rr][1]=mEcken[rr][1];
    };
    quad_root->enthaltene_Elemente=(int*)malloc(quad_root->n_enthaltene_Elemente*sizeof(int));
    for (rr=0;rr<dis[0].numele;rr++){
	quad_root->enthaltene_Elemente[rr]=dis[0].element[rr].Id_loc;/*Id_loc*/
    }
/*--------------------------------------------------*/
/* Erzeugung der quadtree-Struktur durch rekursives
   anh"angen neuer Knoten                           */
    quadtree_recursion (dis,nmax_pro_Blatt, quad_root);
/*--------------------------------------------------*/
/*       Abspeichern der Wurzel im dynamic-struct   */
    chm_data[0].quad_root[ndis]=quad_root;
#ifdef DEBUG
    dstrc_exit();
#endif
    return;
}
/*--------------------------------------------------*/
/*         quadtree-Rekursion f"ur Erzeugung        */
/*BESCHREIBUNG:

"Ubergeben wird im Prinzip ein Knoten der quadtree-Zerlegung. Es werden
folgende Schritte durchgef"uhrt:
 1) "Uberpr"ufung ob weitere Zerlegung notwendig
 2) JA:   Vier neue Knoten werden erzeugt und die Funktion auf diese rekursiv
          angewendet
 3) NEIN: Alle weiterzeiger werden auf NULL gesetzt. Blatt erreicht!

Problem: Unter wiedrigen Umst"anden kann es passieren, da"s einzelne Elemente
-------  nicht in den Baum aufgenommen werden (dies ist der Fall, falls das
	 Element keinen Knoten mit dem Rechteck gemeinsam hat)

	 Weiter: Wird nmax_pro_Blatt zu klein gew"ahlt, so mu"s ein Fehler
	         auftreten sobald mehr Elemente an einem Knoten
		 zusammenstossen wie nmax_pro_Blatt gro"s ist.

	 Hier mu"s im Suchalgorithmus eine brute-force-Suche durchgef"uhrt
	 werden.*/

/*--------------------------------------------------
quadtree recursion
------------------

description:
a quadtree node is passed to this subroutine. In a first step, the number of elements in the
node  is calculated to decide whether to split further or not.

If the node has to be split further, four new children are created and the subroutine is
called for them.

Known problems:
1) nmax_pro_Blatt has to be chosen big enough (at least the maximum number of elements belonging
   to one node)
2) quadtree search might fail --- in that case, a brute-force-search has to be performed

*/
void quadtree_recursion(DISCRET *disrec,int nmax_pro_Blatt,struct quadtree_node_struct *Knoten)
{
/* Variablendeklaration */
    int rr;                                     /* Laufvariablen          */
    quadtree_node_type *nklu,*nkru,*nkro,*nklo; /* Zeiger auf neue Knoten */
    double d[3][2],a[3],b[3],qTOL;              /* Abk"urzungen           */


    qTOL=chm_data[0].quadtree_TOL;
/*    printf("quadtree------------------n_enthaltene_Elemente %d\n",Knoten->n_enthaltene_Elemente);
    printf("(%5f,%5f)          (%5f,%5f)\n",Knoten->node_Ecken[3][0],Knoten->node_Ecken[3][1],Knoten->node_Ecken[2][0],Knoten->node_Ecken[2][1]);
    printf("(%5f,%5f)          (%5f,%5f)\n",Knoten->node_Ecken[0][0],Knoten->node_Ecken[0][1],Knoten->node_Ecken[1][0],Knoten->node_Ecken[1][1]);
    printf("quadtree------------------\n");*/
/* Fall 1: Knoten ist Blatt*/
    if (Knoten->n_enthaltene_Elemente<=nmax_pro_Blatt)
    {
	Knoten->Nachfolger_lu=NULL;
	Knoten->Nachfolger_ru=NULL;
	Knoten->Nachfolger_ro=NULL;
	Knoten->Nachfolger_lo=NULL;
	return;
    }
/* Fall 2: Knoten muss geteilt werden*/
    if (Knoten->n_enthaltene_Elemente>nmax_pro_Blatt)
    {
        /* Speicher f"ur neue Knoten*/
	nklu=(quadtree_node_type*)malloc(sizeof(quadtree_node_type));
	nkru=(quadtree_node_type*)malloc(sizeof(quadtree_node_type));
 	nkro=(quadtree_node_type*)malloc(sizeof(quadtree_node_type));
	nklo=(quadtree_node_type*)malloc(sizeof(quadtree_node_type));
	/*------------------------------------*/
	nklu->node_Ecken=(double**)malloc(4*sizeof(double*));
	nkru->node_Ecken=(double**)malloc(4*sizeof(double*));
	nkro->node_Ecken=(double**)malloc(4*sizeof(double*));
	nklo->node_Ecken=(double**)malloc(4*sizeof(double*));

	for (rr=0;rr<4;rr++){
	    nklu->node_Ecken[rr]=(double*)malloc(2*sizeof(double));
	    nkru->node_Ecken[rr]=(double*)malloc(2*sizeof(double));
	    nkro->node_Ecken[rr]=(double*)malloc(2*sizeof(double));
	    nklo->node_Ecken[rr]=(double*)malloc(2*sizeof(double));
	};
	/*------------------------------------ (erstmal zu gro"s)*/
	nklu->enthaltene_Elemente=(int*)malloc(Knoten->n_enthaltene_Elemente*sizeof(int));
	nkru->enthaltene_Elemente=(int*)malloc(Knoten->n_enthaltene_Elemente*sizeof(int));
	nkro->enthaltene_Elemente=(int*)malloc(Knoten->n_enthaltene_Elemente*sizeof(int));
	nklo->enthaltene_Elemente=(int*)malloc(Knoten->n_enthaltene_Elemente*sizeof(int));

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* Neue Nachfolger werden angeh"angt*/
	Knoten->Nachfolger_lu=nklu;
	Knoten->Nachfolger_ru=nkru;
	Knoten->Nachfolger_lo=nklo;
	Knoten->Nachfolger_ro=nkro;

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* Ecken werden zugeordnet*/
	Knoten->Nachfolger_lu->node_Ecken[0][0]=Knoten->node_Ecken[0][0];
	Knoten->Nachfolger_lu->node_Ecken[0][1]=Knoten->node_Ecken[0][1];

	Knoten->Nachfolger_lu->node_Ecken[1][0]=0.5*(Knoten->node_Ecken[0][0]+Knoten->node_Ecken[1][0]);
	Knoten->Nachfolger_lu->node_Ecken[1][1]=Knoten->node_Ecken[0][1];

	Knoten->Nachfolger_lu->node_Ecken[2][0]=Knoten->Nachfolger_lu->node_Ecken[1][0];
	Knoten->Nachfolger_lu->node_Ecken[2][1]=0.5*(Knoten->node_Ecken[0][1]+Knoten->node_Ecken[3][1]);

	Knoten->Nachfolger_lu->node_Ecken[3][0]=Knoten->node_Ecken[0][0];
	Knoten->Nachfolger_lu->node_Ecken[3][1]=Knoten->Nachfolger_lu->node_Ecken[2][1];
	/*------------------------------------*/
	Knoten->Nachfolger_ru->node_Ecken[0][0]=Knoten->Nachfolger_lu->node_Ecken[1][0];
	Knoten->Nachfolger_ru->node_Ecken[0][1]=Knoten->Nachfolger_lu->node_Ecken[1][1];

	Knoten->Nachfolger_ru->node_Ecken[1][0]=Knoten->node_Ecken[1][0];
	Knoten->Nachfolger_ru->node_Ecken[1][1]=Knoten->node_Ecken[1][1];

	Knoten->Nachfolger_ru->node_Ecken[2][0]=Knoten->node_Ecken[1][0];
	Knoten->Nachfolger_ru->node_Ecken[2][1]=Knoten->Nachfolger_lu->node_Ecken[2][1];

	Knoten->Nachfolger_ru->node_Ecken[3][0]=Knoten->Nachfolger_lu->node_Ecken[1][0];
	Knoten->Nachfolger_ru->node_Ecken[3][1]=Knoten->Nachfolger_lu->node_Ecken[2][1];
	/*------------------------------------*/
	Knoten->Nachfolger_ro->node_Ecken[0][0]=Knoten->Nachfolger_ru->node_Ecken[3][0];
	Knoten->Nachfolger_ro->node_Ecken[0][1]=Knoten->Nachfolger_ru->node_Ecken[3][1];

	Knoten->Nachfolger_ro->node_Ecken[1][0]=Knoten->node_Ecken[1][0];
	Knoten->Nachfolger_ro->node_Ecken[1][1]=Knoten->Nachfolger_lu->node_Ecken[2][1];

	Knoten->Nachfolger_ro->node_Ecken[2][0]=Knoten->node_Ecken[2][0];
	Knoten->Nachfolger_ro->node_Ecken[2][1]=Knoten->node_Ecken[2][1];

	Knoten->Nachfolger_ro->node_Ecken[3][0]=Knoten->Nachfolger_ro->node_Ecken[0][0];
	Knoten->Nachfolger_ro->node_Ecken[3][1]=Knoten->node_Ecken[2][1];
	/*------------------------------------*/
	Knoten->Nachfolger_lo->node_Ecken[0][0]=Knoten->Nachfolger_lu->node_Ecken[3][0];
	Knoten->Nachfolger_lo->node_Ecken[0][1]=Knoten->Nachfolger_lu->node_Ecken[3][1];

	Knoten->Nachfolger_lo->node_Ecken[1][0]=Knoten->Nachfolger_lu->node_Ecken[2][0];
	Knoten->Nachfolger_lo->node_Ecken[1][1]=Knoten->Nachfolger_lu->node_Ecken[2][1];

	Knoten->Nachfolger_lo->node_Ecken[2][0]=Knoten->Nachfolger_ro->node_Ecken[3][0];
	Knoten->Nachfolger_lo->node_Ecken[2][1]=Knoten->Nachfolger_ro->node_Ecken[3][1];

	Knoten->Nachfolger_lo->node_Ecken[3][0]=Knoten->node_Ecken[3][0];
	Knoten->Nachfolger_lo->node_Ecken[3][1]=Knoten->node_Ecken[3][1];

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
        /* Es werden den neuen Knoten Elemente zugeordnet nach folgendem Schema*/
	/* Enth"alt Rechteck einen Koordinatenpunkt des Elements, dann zuordnen*/

	/* the elements are attached to the quadtree-node as follows:
           if the box contains one node of the element, the element belongs to the box*/
	Knoten->Nachfolger_lu->n_enthaltene_Elemente=0;
	Knoten->Nachfolger_ru->n_enthaltene_Elemente=0;
	Knoten->Nachfolger_lo->n_enthaltene_Elemente=0;
	Knoten->Nachfolger_ro->n_enthaltene_Elemente=0;

	for (rr=0;rr<(Knoten->n_enthaltene_Elemente);rr++)
	{
	    a[0]=Knoten->Nachfolger_lu->node_Ecken[0][0];
	    a[1]=Knoten->Nachfolger_ru->node_Ecken[0][0];
	    a[2]=Knoten->Nachfolger_ru->node_Ecken[1][0];
	    b[0]=Knoten->Nachfolger_lu->node_Ecken[0][1];
	    b[1]=Knoten->Nachfolger_lu->node_Ecken[3][1];
	    b[2]=Knoten->Nachfolger_lo->node_Ecken[3][1];

	    if (disrec[0].element[0].distyp == tri3)
	    {
		d[0][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[0][0].x[0];
		d[0][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[0][0].x[1];
		d[1][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[1][0].x[0];
		d[1][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[1][0].x[1];
		d[2][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[2][0].x[0];
		d[2][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[2][0].x[1];
		/* d sind nun die Koordinaten der Eckpunkte des Elements*/
		/* d are the coordinates of the nodes of the element*/
		if (((a[0]-qTOL<=d[0][0] && d[0][0]<=a[1]+qTOL) && (b[0]-qTOL<=d[0][1] && d[0][1]<=b[1]+qTOL)) || ((a[0]-qTOL<=d[1][0] && d[1][0]<=a[1]+qTOL) && (b[0]-qTOL<=d[1][1] && d[1][1]<=b[1]+qTOL)) || ((a[0]-qTOL<=d[2][0] && d[2][0]<=a[1]+qTOL) && (b[0]-qTOL<=d[2][1] && d[2][1]<=b[1]+qTOL))){/* linke Seite unten */
		    Knoten->Nachfolger_lu->enthaltene_Elemente[Knoten->Nachfolger_lu->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_lu->n_enthaltene_Elemente++;
		}

		/*------------------------------------*/

		if (((a[0]-qTOL<=d[0][0] && d[0][0]<=a[1]+qTOL) && (b[1]-qTOL<=d[0][1] && d[0][1]<=b[2]+qTOL)) || ((a[0]-qTOL<=d[1][0] && d[1][0]<=a[1]+qTOL) && (b[1]-qTOL<=d[1][1] && d[1][1]<=b[2]+qTOL)) || ((a[0]-qTOL<=d[2][0] && d[2][0]<=a[1]+qTOL) && (b[1]-qTOL<=d[2][1] && d[2][1]<=b[2]+qTOL))){/* linke Seite oben */
		    Knoten->Nachfolger_lo->enthaltene_Elemente[Knoten->Nachfolger_lo->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_lo->n_enthaltene_Elemente++;
		}
		/*------------------------------------*/

		if (((a[1]-qTOL<=d[0][0] && d[0][0]<=a[2]+qTOL) && (b[0]-qTOL<=d[0][1] && d[0][1]<=b[1]+qTOL)) || ((a[1]-qTOL<=d[1][0] && d[1][0]<=a[2]+qTOL) && (b[0]-qTOL<=d[1][1] && d[1][1]<=b[1]+qTOL)) || ((a[1]-qTOL<=d[2][0] && d[2][0]<=a[2]+qTOL) && (b[0]-qTOL<=d[2][1] && d[2][1]<=b[1]+qTOL))){/* rechts unten */
		    Knoten->Nachfolger_ru->enthaltene_Elemente[Knoten->Nachfolger_ru->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_ru->n_enthaltene_Elemente++;
		}
		/*------------------------------------*/

		if (((a[1]-qTOL<=d[0][0] && d[0][0]<=a[2]+qTOL) && (b[1]-qTOL<=d[0][1] && d[0][1]<=b[2]+qTOL)) || ((a[1]-qTOL<=d[1][0] && d[1][0]<=a[2]+qTOL) && (b[1]-qTOL<=d[1][1] && d[1][1]<=b[2]+qTOL)) || ((a[1]-qTOL<=d[2][0] && d[2][0]<=a[2]+qTOL) && (b[1]-qTOL<=d[2][1] && d[2][1]<=b[2]+qTOL))){/* rechts oben */
		    Knoten->Nachfolger_ro->enthaltene_Elemente[Knoten->Nachfolger_ro->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_ro->n_enthaltene_Elemente++;
		}
	    }
	    if (disrec[0].element[0].distyp == quad4)
	    {
		d[0][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[0][0].x[0];
		d[0][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[0][0].x[1];
		d[1][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[1][0].x[0];
		d[1][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[1][0].x[1];
		d[2][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[2][0].x[0];
		d[2][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[2][0].x[1];
		d[3][0]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[3][0].x[0];
		d[3][1]=disrec[0].element[Knoten->enthaltene_Elemente[rr]].node[3][0].x[1];
		/* d sind nun die Koordinaten der Eckpunkte des Elements*/
		/* d are the coordinates of the nodes of the element*/

		if (((a[0]<=d[0][0]+qTOL && d[0][0]<=a[1]+qTOL) && (b[0]<=d[0][1]+qTOL && d[0][1]<=b[1]+qTOL)) || ((a[0]<=d[1][0]+qTOL && d[1][0]<=a[1]+qTOL) && (b[0]<=d[1][1]+qTOL && d[1][1]<=b[1]+qTOL)) || ((a[0]<=d[2][0]+qTOL && d[2][0]<=a[1]+qTOL) && (b[0]<=d[2][1]+qTOL && d[2][1]<=b[1]+qTOL)) || ((a[0]<=d[3][0]+qTOL && d[3][0]<=a[1]+qTOL) && (b[0]<=d[3][1]+qTOL && d[3][1]<=b[1]+qTOL))){/* linke Seite unten */
		    Knoten->Nachfolger_lu->enthaltene_Elemente[Knoten->Nachfolger_lu->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_lu->n_enthaltene_Elemente++;
		}
		/*------------------------------------*/

		if (((a[0]<=d[0][0]+qTOL && d[0][0]<=a[1]+qTOL) && (b[1]<=d[0][1]+qTOL && d[0][1]<=b[2]+qTOL)) || ((a[0]<=d[1][0]+qTOL && d[1][0]<=a[1]+qTOL) && (b[1]<=d[1][1]+qTOL && d[1][1]<=b[2]+qTOL)) || ((a[0]<=d[2][0]+qTOL && d[2][0]<=a[1]+qTOL) && (b[1]<=d[2][1]+qTOL && d[2][1]<=b[2]+qTOL))||((a[0]<=d[3][0]+qTOL && d[3][0]<=a[1]+qTOL) && (b[1]<=d[3][1]+qTOL && d[3][1]<=b[2]+qTOL))){/* linke Seite oben */
		    Knoten->Nachfolger_lo->enthaltene_Elemente[Knoten->Nachfolger_lo->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_lo->n_enthaltene_Elemente++;
		}
		/*------------------------------------*/

		if (((a[1]<=d[0][0]+qTOL && d[0][0]<=a[2]+qTOL) && (b[0]<=d[0][1]+qTOL && d[0][1]<=b[1]+qTOL)) || ((a[1]<=d[1][0]+qTOL && d[1][0]<=a[2]+qTOL) && (b[0]<=d[1][1]+qTOL && d[1][1]<=b[1]+qTOL)) || ((a[1]<=d[2][0]+qTOL && d[2][0]<=a[2]+qTOL) && (b[0]<=d[2][1]+qTOL && d[2][1]<=b[1]+qTOL))|| ((a[1]<=d[3][0]+qTOL && d[3][0]<=a[2]+qTOL) && (b[0]<=d[3][1]+qTOL && d[3][1]<=b[1]+qTOL))){/* rechts unten */
		    Knoten->Nachfolger_ru->enthaltene_Elemente[Knoten->Nachfolger_ru->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_ru->n_enthaltene_Elemente++;
		}
		/*------------------------------------*/

		if (((a[1]<=d[0][0]+qTOL && d[0][0]<=a[2]+qTOL) && (b[1]<=d[0][1]+qTOL && d[0][1]<=b[2]+qTOL)) || ((a[1]<=d[1][0]+qTOL && d[1][0]<=a[2]+qTOL) && (b[1]<=d[1][1]+qTOL && d[1][1]<=b[2]+qTOL)) || ((a[1]<=d[2][0]+qTOL && d[2][0]<=a[2]+qTOL) && (b[1]<=d[2][1]+qTOL && d[2][1]<=b[2]+qTOL))|| ((a[1]<=d[3][0]+qTOL && d[3][0]<=a[2]+qTOL) && (b[1]<=d[3][1]+qTOL && d[3][1]<=b[2]+qTOL))){/* rechts oben */
		    Knoten->Nachfolger_ro->enthaltene_Elemente[Knoten->Nachfolger_ro->n_enthaltene_Elemente]=Knoten->enthaltene_Elemente[rr];
		    Knoten->Nachfolger_ro->n_enthaltene_Elemente++;
		}
	    }
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* Rekursion*/
	quadtree_recursion (disrec,nmax_pro_Blatt,Knoten->Nachfolger_lu);
	quadtree_recursion (disrec,nmax_pro_Blatt,Knoten->Nachfolger_ru);
	quadtree_recursion (disrec,nmax_pro_Blatt,Knoten->Nachfolger_ro);
	quadtree_recursion (disrec,nmax_pro_Blatt,Knoten->Nachfolger_lo);
    }
}

/*--------------------------------------------------*/
/*F"uhrt Suche zu gegebenem quadtree durch          */
/*--------------------------------------------------*/
ELEMENT* chimera_search_quadtree_search(
                        double      x,
                        double      y,
                        DISCRET    *dis,
			int         ndis,
                        double      TOL
                        )
{
    ELEMENT            *parent;
    quadtree_node_type *Knoten;         /* Zeiger auf neue Knoten */
    double              a[3],b[3];                   /* Abk"urzungen           */
    int                 indi;
    double              xi,zeta,temp;
    DOUBLE              N0[2],N1[2],N2[2],N3[2];
    int                 i,j,rr;

#ifdef DEBUG
    dstrc_enter("chimera_search_quadtree_search");
#endif
    parent=NULL;
/*--------------------------------------------------*/
/*     Bestimme das Blatt in welchem (x,y) liegt    */
    Knoten=chm_data[0].quad_root[ndis];
    while (Knoten->Nachfolger_lu != NULL){
	a[0]=Knoten->Nachfolger_lu->node_Ecken[0][0];
	a[1]=Knoten->Nachfolger_lu->node_Ecken[1][0];
	a[2]=Knoten->Nachfolger_ru->node_Ecken[2][0];
	b[0]=Knoten->Nachfolger_lu->node_Ecken[0][1];
	b[1]=Knoten->Nachfolger_lu->node_Ecken[3][1];
	b[2]=Knoten->Nachfolger_lo->node_Ecken[3][1];
	indi=0;
	if((a[0]-TOL<= x) && (x<a[1])){
	    if ((b[0]-TOL<=y) && (y<b[1])){
		Knoten=Knoten->Nachfolger_lu;
		indi=1;
	    }
	    if ((b[1]<=y) && (y<=b[2]+TOL)){
		Knoten=Knoten->Nachfolger_lo;
		indi=1;
	    }
	}
	if((a[1]<= x) && (x<=a[2]+TOL)){
	    if ((b[0]-TOL<=y) && (y<b[1])){
		Knoten=Knoten->Nachfolger_ru;
		indi=1;
	    }
	    if ((b[1]<=y) && (y<=b[2]+TOL)){
		Knoten=Knoten->Nachfolger_ro;
		indi=1;
	    }
	}
	if (indi == 0)
	{
	    return(NULL);
/*    	    printf("CHIMERA_SEARCH quadtree-error --- Kein Blatt zu finden!\n");*/
	}
    };
    /* Jetzt ist Knoten der Zeiger auf das entsprechende Blatt */
/*--------------------------------------------------*/
/*            brute-force-Suche im Blatt            */
/****************************************************** Dreieckselemente*/
    if (dis[0].element[0].distyp == tri3)
    {
	indi=0;
	for (j=0;j<Knoten[0].n_enthaltene_Elemente;j++)
	{
/* loop over all background-elements*/
	    parent=&(dis[0].element[Knoten[0].enthaltene_Elemente[j]]);
	    for (rr=0;rr<2;rr++){
		N0[rr]=(parent->node[0])->x[rr];
		N1[rr]=(parent->node[1])->x[rr];
		N2[rr]=(parent->node[2])->x[rr];
	    }
	    temp=1/((N1[0]-N0[0])*(N1[1]-N2[1])-(N1[0]-N2[0])*(N1[1]-N0[1]));
	    xi = (x*(N0[1]-N2[1])+y*(N2[0]-N0[0])+N0[0]*N2[1]-N2[0]*N0[1])*temp;
	    temp=1/((N2[0]-N0[0])*(N2[1]-N1[1])-(N2[0]-N1[0])*(N2[1]-N0[1]));
	    zeta=(x*(N0[1]-N1[1])+y*(N1[0]-N0[0])+N0[0]*N1[1]-N1[0]*N0[1])*temp;

	    if (0<=xi+TOL && xi<=1+TOL)
	    {
		if (0<=zeta+TOL && zeta<=1+TOL)
		{
		    temp =xi+zeta;
		    if (0<=temp+TOL && temp<=1+TOL)
		    {
			indi=1000;
			break;
		    }
		}
	    }
	}
	i=j;
    }
/*************************************************** Rechteckselemente*/
    if (dis[0].element[0].distyp == quad4)
    {
	indi=0;
	for (j=0;j<Knoten[0].n_enthaltene_Elemente;j++){/* loop over all background-elements*/
	    parent=&(dis[0].element[Knoten[0].enthaltene_Elemente[j]]);
	    for (rr=0;rr<2;rr++){
		N0[rr]=(parent->node[0])->x[rr];
		N1[rr]=(parent->node[1])->x[rr];
		N3[rr]=(parent->node[3])->x[rr];
	    }
	    temp=(N1[0]-N0[0])*(N1[0]-N0[0])+(N1[1]-N0[1])*(N1[1]-N0[1]);
	    xi=((x-N0[0])*(N1[0]-N0[0])+(y-N0[1])*(N1[1]-N0[1]))/temp;
	    temp=(N3[0]-N0[0])*(N3[0]-N0[0])+(N3[1]-N0[1])*(N3[1]-N0[1]);
	    zeta=((x-N0[0])*(N3[0]-N0[0])+(y-N0[1])*(N3[1]-N0[1]))/temp;
	    if (0<=xi+TOL && xi<=1+TOL)
	    {
		if (0<=zeta+TOL && zeta<=1+TOL)
		{
		    indi=1000;
		    break;
		}
	    }
	}
	i=j;
    }
/*Plausibilit"ats"uberpr"ufung*/
    if (indi!=1000 && i==Knoten[0].n_enthaltene_Elemente)
    {
/*	printf("-------------------------------------error\n");
	printf("no parent found --- quadtree-miss\n");
	printf("Punkt %5f   %5f \n",x,y);*/
	parent=chimera_search_brute_force(x,y,dis,TOL);
	if (parent!=NULL)
	{
/*	    printf("Brute-force lieferte: parent-Element %d\n",parent[0].Id );
	    printf("Die Koordinaten waren (%5f,%5f)   (%5f,%5f)   (%5f,%5f)\n",parent[0].node[0][0].x[0],parent[0].node[0][0].x[1],parent[0].node[0][1].x[0],parent[0].node[0][1].x[1],parent[0].node[0][2].x[0],parent[0].node[0][2].x[1]);*/
	}
	else
	{
/*	    printf("Brute-force lieferte kein parent-Element\n");   */
	}
    }

#ifdef DEBUG
    dstrc_exit();
#endif
    return(parent);
}
/*--------------------------------------------------*/
/* Gibt den Speicher aus quadtree ndis wieder frei  */
/*--------------------------------------------------*/
void mem_quadtree_free(
                         quadtree_node_type *Knoten
                        )
{
#ifdef DEBUG
    dstrc_enter("mem_quadtree_free");
#endif
/*--------------------------------------------------*/
/*             Variablendeklarationen               */
/*--------------------------------------------------*/
    int rr; /* Laufvariable */
/*--------------------------------------------------*/
/* Rekursive Aufrufe und Freigabe der Strukturen    */
/*--------------------------------------------------*/
    if (Knoten->Nachfolger_lu != NULL){/* noch Nachfolger da  (kein Blatt)*/
	mem_quadtree_free(Knoten->Nachfolger_lu);
	free(Knoten->Nachfolger_lu);
	Knoten->Nachfolger_lu=NULL;
    };
    if (Knoten->Nachfolger_ru != NULL){/* noch Nachfolger da  (kein Blatt)*/
	mem_quadtree_free(Knoten->Nachfolger_ru);
	free(Knoten->Nachfolger_ru);
	Knoten->Nachfolger_ru=NULL;
    };
    if (Knoten->Nachfolger_ro != NULL){/* noch Nachfolger da  (kein Blatt)*/
	mem_quadtree_free(Knoten->Nachfolger_ro);
	free(Knoten->Nachfolger_ro);
	Knoten->Nachfolger_ro=NULL;
    };
    if (Knoten->Nachfolger_lo != NULL){/* noch Nachfolger da  (kein Blatt)*/
	mem_quadtree_free(Knoten->Nachfolger_lo);
	free(Knoten->Nachfolger_lo);
	Knoten->Nachfolger_lo=NULL;
    };
/*--------------------------------------------------*/
/*   Freigeben des Speichers am aktuellen Knoten    */
/*--------------------------------------------------*/
    for (rr=0;rr<4;rr++){
	free((Knoten->node_Ecken)[rr]);
    }
    free(Knoten->node_Ecken);
    free(Knoten->enthaltene_Elemente);
#ifdef DEBUG
  dstrc_exit();
#endif
    return;
}
void chimera_search_quadtree_free(
                        INT         ndis
    )
{
#ifdef DEBUG
    dstrc_enter("chimera_search_quadtree_free");
#endif
    mem_quadtree_free(chm_data[0].quad_root[ndis]);
    free(chm_data[0].quad_root[ndis]);
#ifdef DEBUG
    dstrc_exit();
#endif
    return;
}
#endif
