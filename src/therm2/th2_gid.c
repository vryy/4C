/*======================================================================*/
/*!
\file
\brief Print THERM2 element data to Ccarat output file

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/

/*----------------------------------------------------------------------*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"
#include "therm2.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/06
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 10/06
*/
extern FIELD *field;


/*======================================================================*/
/*!
\brief Set required Gauss point sets for Gid output file

\param   *actele        ELEMENT     (i)   pointer to current element
\param   *actgid        GIDSET      (o)   Gid data (incl. Gauss point 
                                          sets

\return void

\author bborn
\date 08/06
*/
void th2_gid_init(ELEMENT *actele,
                  GIDSET *actgid)
{
  
  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_gid_init");
#endif

  /*--------------------------------------------------------------------*/
  switch (actele->distyp)
  {
    /* linear quadrilaterals */
    case quad4:
      if (actele->e.th2->nGP[0] == 1)
      {
        actgid->is_therm2_q4_11 = 1;
        actgid->therm2_q4_11_name = "therm2_q4_11";
      }
      if (actele->e.th2->nGP[0] == 2)
      {
        actgid->is_therm2_q4_22 = 1;
        actgid->therm2_q4_22_name = "therm2_q4_22";
      }
    /* quadratic quadrilaterals */
    case quad8:
      if (actele->e.th2->nGP[0] == 3)
      {
        actgid->is_therm2_q8_33 = 1;
        actgid->therm2_q8_33_name = "therm2_q8_33";
      }
      break;
    /* quadratic quadrilaterals */
    case quad9:
      if (actele->e.th2->nGP[0] == 3)
      {
        actgid->is_therm2_q9_33 = 1;
        actgid->therm2_q9_33_name = "therm2_q9_33";
      }
      break;
    /* linear triangles */
    case tri3:
      if (actele->e.th2->nGP[0] == 1)
      {
        actgid->is_therm2_t3_1 = 1;
        actgid->therm2_t3_1_name = "therm2_t3_1";
      }
      break;
    /* quadratic triangles */
    case tri6:
      if (actele->e.th2->nGP[0]==3)
      {
        actgid->is_therm2_t6_3 = 1;
        actgid->therm2_t6_3_name = "therm2_t6_3";
      }
      break;
    /* catch errors */
    default:
      dserror("Discretisation type is impossible!");
      break;
  }  /* end switch */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_gid_init */



/*======================================================================*/
/*!
\brief Print header of Gauss point set to Gid output file

\param   *actfield      FIELD       (i)   current field
\param    ldis          INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th2_gid_msh(FIELD *actfield,
                 INT ldis,
                 GIDSET *actgid,
                 FILE *out)
{
  INT nelenod;
  ELEMENT *actele;
  INT j, k;

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_gid_msh");
#endif
  
  /*----------------------------------------------------------------*/
  /* quadrilateral element with 4 nodes and 2x2 Gauss points */
  if (actgid->is_therm2_q4_22)
  {
    nelenod = 4;
    /* print header */
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# MESH %s FOR FIELD %s, DIS %1i: THERM2 QUAD4 2x2 GP\n",
            actgid->therm2_q4_22_name, actgid->fieldname, ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 2 ELEMTYPE Quadrilateral NNODE %d\n",
            actgid->therm2_q4_22_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 4) )
      {
        fprintf(out, " %6d ", actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out, "%6d ", actele->node[k]->Id+1);
        }
        fprintf(out, "\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*----------------------------------------------------------------*/
  /* quadrilateral element with 8 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q8_33)
  {
    nelenod = 8;
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM2 QUAD8 3x3 GP\n",
            actgid->therm2_q8_33_name, actgid->fieldname, ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d\n",
            actgid->therm2_q8_33_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp >= 8) )
      {
        fprintf(out, " %6d ", actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out, "%6d ", actele->node[k]->Id+1);
        }
        fprintf(out,"\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*----------------------------------------------------------------*/
  /* quadrilateral element with 8 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q9_33)
  {
    nelenod = 9;
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM2 QUAD9 3x3 GP\n",
            actgid->therm2_q8_33_name, actgid->fieldname, ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d\n",
            actgid->therm2_q8_33_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp >= 8) )
      {
        fprintf(out, " %6d ", actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out, "%6d ", actele->node[k]->Id+1);
        }
        fprintf(out,"\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*----------------------------------------------------------------*/
  /* triangular element with 3 nodes and 1 Gauss point */
  if (actgid->is_therm2_t3_1)
  {
    nelenod = 3;
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM2 TRI3 1 GP\n",
            actgid->therm2_t3_1_name, actgid->fieldname, ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 2 ELEMTYPE Triangle NNODE %d\n",
            actgid->therm2_t3_1_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 3) )
      {
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        }
        fprintf(out,"\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*----------------------------------------------------------------*/
  /* triangular element with 6 nodes and 3 Gauss points */
  if (actgid->is_therm2_t6_3)
  {
    nelenod = 6;
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM2 TRI6 3 GP\n",
            actgid->therm2_t6_3_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 2 ELEMTYPE Triangle NNODE %d\n",
            actgid->therm2_t6_3_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 6) )
      {
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        }
        fprintf(out,"\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_gid_msh */


/*======================================================================*/
/*!
\brief Print header of Gauss point set to Gid output file

\param    jdis          INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data (incl. GP sets)
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th2_gid_gpset(INT jdis,
                   GIDSET *actgid,
                   FILE *out)
{
  
  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_gid_gpset");
#endif

   /*-------------------------------------------------------------------*/
   /* quadrilateral element with 4 nodes and 2x2 Gauss points */
   if (actgid->is_therm2_q4_22)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM2 QUAD4 2x2 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Quadrilateral \"%s_dis_%1i\"\n",
             actgid->therm2_q4_22_name, jdis,
             actgid->therm2_q4_22_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 4\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* quadrilateral element with 8/9 nodes and 3x3 Gauss points */
   if (actgid->is_therm2_q8_33)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM2 QUAD8 3x3 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Quadrilateral \"%s_dis_%1i\"\n",
             actgid->therm2_q8_33_name, jdis,
             actgid->therm2_q8_33_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 9\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* quadrilateral element with 8/9 nodes and 3x3 Gauss points */
   if (actgid->is_therm2_q9_33)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM2 QUAD9 3x3 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Quadrilateral \"%s_dis_%1i\"\n",
             actgid->therm2_q9_33_name, jdis,
             actgid->therm2_q9_33_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 9\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* triangular element with 3 nodes and 1 Gauss point */
   if (actgid->is_therm2_t3_1)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM2 TRI3 1 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Triangle \"%s_dis_%1i\"\n",
             actgid->therm2_t3_1_name, jdis,
             actgid->therm2_t3_1_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 1\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* triangular element with 6 nodes and 3 Gauss points */
   if (actgid->is_therm2_t6_3)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM2 TRI6 3 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Triangle \"%s_dis_%1i\"\n",
             actgid->therm2_t6_3_name, jdis, 
             actgid->therm2_t6_3_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 3\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th2_gid_gpset */



/*======================================================================*/
/*!
\brief Print header for domain section

\param    disnum        INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data (incl. Gauss point 
                                          sets
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th2_gid_dom(FIELD *actfield,
                 INT disnum,
                 GIDSET *actgid,
                 FILE *out)
{
  INT ngauss;
  ELEMENT *actele;
  INT i, j;

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_gid_domhdr");
#endif


  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 4 nodes and 2x2 Gauss points */
  if (actgid->is_therm2_q4_22)
  {
    ngauss = 4;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm2_q4_22_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm2_q4_22_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 4) )
      {
        fprintf(out, "    %6d  %18.5E\n",
                actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)  /* loop number of GPs */
        {
          fprintf(out, "            %18.5E\n",
                  (DOUBLE)actele->proc);
        }
      }      
    }
    fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 8/9 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q8_33)
  {
    ngauss = 9;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm2_q8_33_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm2_q8_33_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp >= 8) )
      {
        fprintf(out, "    %6d  %18.5E\n",
                actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)  /* loop number of GPs */
        {
          fprintf(out, "            %18.5E\n",
                  (DOUBLE)actele->proc);
        }
      }      
    }
    fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 8/9 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q9_33)
  {
    ngauss = 9;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm2_q9_33_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm2_q9_33_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp >= 8) )
      {
        fprintf(out, "    %6d  %18.5E\n",
                actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)  /* loop number of GPs */
        {
          fprintf(out, "            %18.5E\n",
                  (DOUBLE)actele->proc);
        }
      }      
    }
    fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* triangular element with 3 nodes and 1 Gauss point */
  if (actgid->is_therm2_t3_1)
  {
    ngauss = 1;
    /* print headers */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm2_t3_1_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm2_t3_1_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 3) )
      {
        fprintf(out, "    %6d  %18.5E\n", 
                actele->Id+1, (DOUBLE)actele->proc);
      }
    }
    fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* triangular element with 6 nodes and 3 Gauss point */
  if (actgid->is_therm2_t6_3)
  {
    ngauss = 3;
    /* print headers */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm2_t3_1_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm2_t3_1_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm2) && (actele->numnp == 3) )
      {
        fprintf(out, "    %6d  %18.5E\n", 
                actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)  /* loop number of GPs */
        {
          fprintf(out, "            %18.5E\n",
                  (DOUBLE)actele->proc);
        }
      }
    }
    fprintf(out, "END VALUES\n");
  }

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th2_gid_dom */



/*======================================================================*/
/*!
\brief Print header of heat flux section in Gid output file

\param    resstring     char        (i)   result string
\param   *actfield      FIELD       (i)   current field
\param    disnum        INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th2_gid_hflux(char resstring[],
                   FIELD *actfield,
                   INT disnum,
                   INT step,
                   GIDSET *actgid,
                   FILE *out)
{
  const INT place = 0;

  INT ngauss;
  char *resulttype;
  char *resultplace;
  char *gpset = NULL;
  char *rangetable;
  INT  ncomponent;
  char *componentnames[3];

  ELEMENT *actele;
  DOUBLE **heatflux;
  INT i, j;

  INT gaussperm3[3] = {0,1,2};  /* verify! */
  INT gaussperm4[4] = {3,1,0,2};
  INT gaussperm9[9] = {8,2,0,6,5,1,3,7,4};

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th2_gid_hflux");
#endif

  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 4 nodes and 2x2 Gauss points */
  if (actgid->is_therm2_q4_22)
  {
      ngauss            = 4;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm2_q4_22_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
      /* print title */
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
              resstring, actgid->fieldname, disnum);
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
              resstring,
              step,
              resulttype,
              resultplace,
              gpset,
              disnum);
      fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\"\n",
              componentnames[0],
              componentnames[1],
              componentnames[2]);
      /* print values */
      fprintf(out, "VALUES\n");
      for (i=0; i<actfield->dis[disnum].numele; i++)
      {
        /* pointer to current element */
        actele = &(actfield->dis[disnum].element[i]);
        /* print heat fluxes */
        if ( (actele->eltyp == el_therm2) && (actele->numnp == 4) )
        {
          heatflux = actele->e.th2->hflux_gp.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out, " %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][gaussperm4[0]],
                  heatflux[1][gaussperm4[0]],
                  heatflux[2][gaussperm4[0]]);
          /* print values of second to ... Gauss point */
          for (j=1; j<ngauss; j++)
          {
            fprintf(out, "        %18.5E %18.5E %18.5E \n",
                    heatflux[0][gaussperm4[j]],
                    heatflux[1][gaussperm4[j]],
                    heatflux[2][gaussperm4[j]]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 8 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q8_33)
  {
      ngauss            = 9;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm2_q8_33_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
      /*----------------------------------------------------------------*/
      /* print title */
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
              resstring, actgid->fieldname, disnum);
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
              resstring,
              step,
              resulttype,
              resultplace,
              gpset, disnum);
      fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\"\n",
              componentnames[0],
              componentnames[1],
              componentnames[2]);
      /* print values */
      fprintf(out, "VALUES\n");
      for (i=0; i<actfield->dis[disnum].numele; i++)
      {
        /* pointer to current element */
        actele = &(actfield->dis[disnum].element[i]);
        /* print heat fluxes */
        if ( (actele->eltyp == el_therm2) && (actele->numnp == 8) )
        {
          heatflux = actele->e.th2->hflux_gp.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][gaussperm9[0]],
                  heatflux[1][gaussperm9[0]],
                  heatflux[2][gaussperm9[0]]);
          /* print values of second to ... Gauss point */
          for (j=1; j<ngauss; j++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[0][gaussperm9[j]],
                    heatflux[1][gaussperm9[j]],
                    heatflux[2][gaussperm9[j]]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* quadrilateral element with 9 nodes and 3x3 Gauss points */
  if (actgid->is_therm2_q9_33)
  {
      ngauss            = 9;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm2_q9_33_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
      /*----------------------------------------------------------------*/
      /* print title */
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
              resstring, actgid->fieldname, disnum);
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
              resstring,
              step,
              resulttype,
              resultplace,
              gpset, disnum);
      fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\"\n",
              componentnames[0],
              componentnames[1],
              componentnames[2]);
      /* print values */
      fprintf(out, "VALUES\n");
      for (i=0; i<actfield->dis[disnum].numele; i++)
      {
        /* pointer to current element */
        actele = &(actfield->dis[disnum].element[i]);
        /* print heat fluxes */
        if ( (actele->eltyp == el_therm2) && (actele->numnp == 9) )
        {
          heatflux = actele->e.th2->hflux_gp.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][gaussperm9[0]],
                  heatflux[1][gaussperm9[0]],
                  heatflux[2][gaussperm9[0]]);
          /* print values of second to ... Gauss point */
          for (j=1; j<ngauss; j++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[0][gaussperm9[j]],
                    heatflux[1][gaussperm9[j]],
                    heatflux[2][gaussperm9[j]]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* triangular element with 3 nodes and 1 Gauss point */
  if (actgid->is_therm2_t3_1)
  {
      ngauss            = 1;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm2_t3_1_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
      /*----------------------------------------------------------------*/
      /* print title */
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
              resstring, actgid->fieldname, disnum);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
              resstring,
              step,
              resulttype,
              resultplace,
              gpset, disnum);
      fprintf(out,"COMPONENTNAMES \"%s\",\"%s\",\"%s\"\n",
              componentnames[0],
              componentnames[1],
              componentnames[2]);
      /* print values */
      fprintf(out, "VALUES\n");
      for (i=0; i<actfield->dis[disnum].numele; i++)
      {
        /* pointer to current element */
        actele = &(actfield->dis[disnum].element[i]);
        /* print heat fluxes */
        if ( (actele->eltyp == el_therm2) && (actele->numnp == 3) )
        {
          heatflux = actele->e.th2->hflux_gp.a.d3[place];
          /* print values at the Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][0],
                  heatflux[1][0],
                  heatflux[2][0]);
        }
      }
      fprintf(out, "END VALUES\n");
  }  /* end of if (actgid->is_therm2_t3_1) */
  /*--------------------------------------------------------------------*/
  /* triangular element with 6 nodes and 3 Gauss points */
  if (actgid->is_therm2_t6_3)
  {
      ngauss            = 3;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm2_t6_3_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
      /*----------------------------------------------------------------*/
      /* print title */
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
              resstring, actgid->fieldname, disnum);
      fprintf(out, "#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
              resstring,
              step,
              resulttype,
              resultplace,
              gpset,disnum);
      fprintf(out,"COMPONENTNAMES \"%s\",\"%s\",\"%s\"\n",
              componentnames[0],
              componentnames[1],
              componentnames[2]);
      /* print values */
      fprintf(out, "VALUES\n");
      for (i=0; i<actfield->dis[disnum].numele; i++)
      {
        /* pointer to current element */
        actele = &(actfield->dis[disnum].element[i]);
        /* print heat fluxes */
        if ( (actele->eltyp == el_therm2) && (actele->numnp == 6) )
        {
          heatflux = actele->e.th2->hflux_gp.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][gaussperm3[0]],
                  heatflux[1][gaussperm3[0]],
                  heatflux[2][gaussperm3[0]]);
          /* print values of second to ... Gauss point */
          for (j=1; j<ngauss; j++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[0][gaussperm3[j]],
                    heatflux[1][gaussperm3[j]],
                    heatflux[2][gaussperm3[j]]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }  /* end of if (actgid->is_therm2_q8_33) */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th2_gid_hflux */



#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
