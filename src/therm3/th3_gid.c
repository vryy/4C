/*======================================================================*/
/*!
\file
\brief Print THERM3 element data to Ccarat output file

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/

#ifndef CCADISCRET
/*----------------------------------------------------------------------*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"
#include "therm3.h"

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
void th3_gid_init(ELEMENT *actele,
                  GIDSET *actgid)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_gid_init");
#endif

  /*--------------------------------------------------------------------*/
  switch (actele->distyp)
  {
    /* linear hexahedra */
    case hex8:
      if (actele->e.th3->gpnum[0] == 2)
      {
        actgid->is_therm3_h8_222 = 1;
        actgid->therm3_h8_222_name = "therm3_h8_222";
      }
    /* quadratic hexahedra */
    case hex20:
      if (actele->e.th3->gpnum[0] == 3)
      {
        actgid->is_therm3_h20_333 = 1;
        actgid->therm3_h20_333_name = "therm3_h20_333";
      }
      break;
    /* quadratic hexahedra */
    case hex27:
      if (actele->e.th3->gpnum[0] == 3)
      {
        actgid->is_therm3_h27_333 = 1;
        actgid->therm3_h27_333_name = "therm3_h27_333";
      }
      break;
    /* linear tetrahedra */
    case tet4:
      if (actele->e.th3->gpnum[0] == 1)
      {
        actgid->is_therm3_t4_1 = 1;
        actgid->therm3_t4_1_name = "therm3_t4_1";
      }
      break;
    /* quadratic tetrahedra */
    case tet10:
      if (actele->e.th3->gpnum[0] == 4)
      {
        actgid->is_therm3_t10_4 = 1;
        actgid->therm3_t10_4_name = "therm3_t10_4";
      }
      break;
    /* catch errors */
    default:
      dserror("Discretisation type is impossible!");
      break;
  }  /* end switch */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_gid_init */

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
void th3_gid_msh(FIELD *actfield,
                 INT ldis,
                 GIDSET *actgid,
                 FILE *out)
{
  INT nelenod;
  ELEMENT *actele;
  INT j, k;  /* loop index */

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("th3_gid_msh");
#endif

  /*----------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_therm3_h8_222)
  {
    nelenod = 8;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM3 2x2x2 GP\n",
            actgid->therm3_h8_222_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE %d\n",
            actgid->therm3_h8_222_name, ldis, nelenod);
    /* print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm3) && (actele->numnp == 8) )
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
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h20_333)
  {
    nelenod = 20;
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM3 3x3x3 GP\n",
            actgid->therm3_h20_333_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE %d\n",
            actgid->therm3_h20_333_name, ldis, nelenod);
    fprintf(out,"ELEMENTS\n");
    /* print elements */
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm3) && (actele->numnp == 20) )
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
  /* hexahedron element with 27 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h20_333)
  {
    nelenod = 27;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM3 3x3x3 GP\n",
            actgid->therm3_h27_333_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE %d\n",
            actgid->therm3_h27_333_name, ldis, nelenod);
    /* print elements */
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm3) && (actele->numnp == 27) )
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
  /* tetrahedron element with 4 nodes and 1 Gauss point */
  if (actgid->is_therm3_t4_1)
  {
    nelenod = 4;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM3 1 GP\n",
            actgid->therm3_t4_1_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Tetrahedra NNODE %d\n",
            actgid->therm3_t4_1_name, ldis, nelenod);
    /* print elements */
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm3) && (actele->numnp == 4) )
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
  /* tetrahedron element with 10 nodes and 4 Gauss points */
  if (actgid->is_therm3_t10_4)
  {
    nelenod = 10;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: THERM3 4 GP\n",
            actgid->therm3_t10_4_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Tetrahedra NNODE %d\n",
            actgid->therm3_t10_4_name, ldis, nelenod);
    /* print elements */
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_therm3) && (actele->numnp == 10) )
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
}  /* end of th3_gid_msh */


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
void th3_gid_gpset(INT jdis,
                   GIDSET *actgid,
                   FILE *out)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_gid_gpset");
#endif

   /*-------------------------------------------------------------------*/
   /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
   if (actgid->is_therm3_h8_222)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM3 2x2x2 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Hexahedra \"%s_dis_%1i\"\n",
             actgid->therm3_h8_222_name, jdis,
             actgid->therm3_h8_222_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 8\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
   if (actgid->is_therm3_h20_333)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM3 3x3x3 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Hexahedra \"%s_dis_%1i\"\n",
             actgid->therm3_h20_333_name, jdis,
             actgid->therm3_h20_333_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 27\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* hexahedron element with 27 nodes and 3x3x3 Gauss points */
   if (actgid->is_therm3_h27_333)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM3 3x3x3 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Hexahedra \"%s_dis_%1i\"\n",
             actgid->therm3_h27_333_name, jdis,
             actgid->therm3_h27_333_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 27\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* tetrahedron element with 4 nodes and 1 Gauss point */
   if (actgid->is_therm3_t4_1)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM3 1 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Tetrahedra \"%s_dis_%1i\"\n",
             actgid->therm3_t4_1_name, jdis,
             actgid->therm3_t4_1_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 1\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------------------------*/
   /* tetrahedron element with 10 nodes and 4 Gauss points */
   if (actgid->is_therm3_t10_4)
   {
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i THERM3 1 GP\n",
             actgid->fieldname, jdis);
     fprintf(out, "#-------------------------------------------------------------------------------\n");
     fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Tetrahedra \"%s_dis_%1i\"\n",
             actgid->therm3_t10_4_name, jdis,
             actgid->therm3_t10_4_name, jdis);
     fprintf(out, "NUMBER OF GAUSS POINTS: 4\n");
     fprintf(out, "NATURAL COORDINATES: Internal\n");
     fprintf(out, "END GAUSSPOINTS\n");
   }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th3_gid_gpset */



/*======================================================================*/
/*!
\brief Print header for domain section

\param    disnum        INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th3_gid_dom(FIELD *actfield,
                 INT disnum,
                 GIDSET *actgid,
                 FILE *out)
{
  INT ngauss;
  INT i, j;  /* loop index */
  ELEMENT *actele;  /* current element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_gid_dom");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_therm3_h8_222)
  {
    ngauss = 8;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm3_h8_222_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm3_h8_222_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm3) && (actele->distyp == hex8) )
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
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h20_333)
  {
    ngauss = 27;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm3_h20_333_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm3_h20_333_name, disnum);
    /* print values */
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm3) && (actele->distyp == hex20) )
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
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h27_333)
  {
    ngauss = 27;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm3_h27_333_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm3_h27_333_name, disnum);
    /* print values */
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm3) && (actele->distyp == hex27) )
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
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron element with 4 nodes and 1 Gauss point */
  if (actgid->is_therm3_t4_1)
  {
    ngauss = 1;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm3_t4_1_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm3_t4_1_name, disnum);
    /* print values */
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm3) && (actele->distyp == tet4) )
      {
        fprintf(out, "    %6d  %18.5E\n",
                actele->Id+1, (DOUBLE)actele->proc);
      }
    }
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron element with 10 nodes and 4 Gauss point */
  if (actgid->is_therm3_t10_4)
  {
    ngauss = 4;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->therm3_t10_4_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->therm3_t10_4_name, disnum);
    /* print values */
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_therm3) && (actele->distyp == tet10) )
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
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th3_gid_dom */




/*======================================================================*/
/*!
\brief Print header of heat flux section in Gid output file

\param    resstring     char        (i)   result string
\param   *actfield      FIELD       (i)   current field
\param    disnum        INT         (i)   discretisation index
\param    step          INT         (i)   time step index
\param   *actgid        GIDSET      (i)   Gid data
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void th3_gid_hflux(char resstring[],
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
  INT i;

  INT gperm_h_222[8] = {0,4,6,2,1,5,7,3};
  INT gperm_h_333[27] = {0,10,24,6,2,20,26,8,
                         9,21,15,3,1,19,25,7,11,23,17,5,
                         12,10,22,16,4,14,13};
  INT gperm_t_4[4] = {0,1,2,3};
  /* INT gperm_t_10[10] = {0,1,2,3,4,5,6,7,8,9}; */  /* check this */

  DOUBLE **heatflux;

  INT jgp, igp;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_gid_hfluxhdr");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_therm3_h8_222)
  {
      ngauss            = 8;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm3_h8_222_name;
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
        if ( (actele->eltyp == el_therm3) && (actele->distyp == hex8) )
        {
          heatflux = actele->e.th3->hflux_gp_xyz.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          igp = gperm_h_222[0];
          fprintf(out, " %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[igp][0],  /* heatflux x-component */
                  heatflux[igp][1],  /* heatflux y-component */
                  heatflux[igp][2]);  /* heatflux z-component */
          /* print values of second to ... Gauss point */
          for (jgp=1; jgp<ngauss; jgp++)
          {
            igp = gperm_h_222[jgp];
            fprintf(out, "        %18.5E %18.5E %18.5E \n",
                    heatflux[igp][0],  /* heatflux x-component */
                    heatflux[igp][1],  /* heatflux y-component */
                    heatflux[igp][2]);  /* heatflux z-component */
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h20_333)
  {
      ngauss            = 27;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm3_h20_333_name;
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
        if ( (actele->eltyp == el_therm3) && (actele->distyp == hex20) )
        {
          heatflux = actele->e.th3->hflux_gp_xyz.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[gperm_h_333[0]][0],
                  heatflux[gperm_h_333[0]][1],
                  heatflux[gperm_h_333[0]][2]);
          /* print values of second to ... Gauss point */
          for (jgp=1; jgp<ngauss; jgp++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[gperm_h_333[jgp]][0],
                    heatflux[gperm_h_333[jgp]][1],
                    heatflux[gperm_h_333[jgp]][2]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 27 nodes and 3x3x3 Gauss points */
  if (actgid->is_therm3_h27_333)
  {
      ngauss            = 27;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm3_h27_333_name;
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
        if ( (actele->eltyp == el_therm3) && (actele->distyp == hex27) )
        {
          heatflux = actele->e.th3->hflux_gp_xyz.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[gperm_h_333[0]][0],
                  heatflux[gperm_h_333[0]][1],
                  heatflux[gperm_h_333[0]][2]);
          /* print values at 2nd to last Gauss point */
          for (jgp=1; jgp<ngauss; jgp++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[gperm_h_333[jgp]][0],
                    heatflux[gperm_h_333[jgp]][1],
                    heatflux[gperm_h_333[jgp]][2]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron element with 4 nodes and 1 Gauss point */
  if (actgid->is_therm3_t4_1)
  {
      ngauss            = 1;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm3_t4_1_name;
      rangetable        = actgid->standardrangetable;
      ncomponent        = 3;
      componentnames[0] = "heatflux-x";
      componentnames[1] = "heatflux-y";
      componentnames[2] = "heatflux-z";
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
        if ( (actele->eltyp == el_therm3) && (actele->distyp == tet4) )
        {
          heatflux = actele->e.th3->hflux_gp_xyz.a.d3[place];
          /* print values at the Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[0][0],
                  heatflux[0][1],
                  heatflux[0][2]);
        }
      }
      fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron element with 10 nodes and 4 Gauss points */
  if (actgid->is_therm3_t10_4)
  {
      ngauss            = 4;
      resulttype        = "VECTOR";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->therm3_t10_4_name;
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
        if ( (actele->eltyp == el_therm3) && (actele->distyp == tet10) )
        {
          heatflux = actele->e.th3->hflux_gp_xyz.a.d3[place];
          /* print values at first Gauss point _with_ element index */
          fprintf(out," %6d %18.5E %18.5E %18.5E \n",
                  actele->Id+1,
                  heatflux[gperm_t_4[0]][0],
                  heatflux[gperm_t_4[0]][1],
                  heatflux[gperm_t_4[0]][2]);
          /* print values of second to ... Gauss point */
          for (jgp=1; jgp<ngauss; jgp++)
          {
            fprintf(out,"        %18.5E %18.5E %18.5E \n",
                    heatflux[gperm_t_4[jgp]][0],
                    heatflux[gperm_t_4[jgp]][1],
                    heatflux[gperm_t_4[jgp]][2]);
          }
        }
      }
      fprintf(out, "END VALUES\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void th3_gid_hflux */



#endif /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
#endif
