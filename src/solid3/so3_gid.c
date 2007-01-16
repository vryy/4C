/*======================================================================*/
/*!
\file
\brief Print SOLID3 element data to Ccarat Gid output file

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 01/07
*/

/*----------------------------------------------------------------------*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"
#include "solid3.h"

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
void so3_gid_init(ELEMENT *actele,
                  GIDSET *actgid)
{
  
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_init");
#endif

  /*--------------------------------------------------------------------*/
  switch (actele->distyp)
  {
    /* linear hexahedra */
    case hex8:
      if (actele->e.so3->gptot == 8)
      {
        actgid->is_solid3_h8_222 = 1;
        actgid->solid3_h8_222_name = "solid3_h8_222";
      }
    /* quadratic hexahedra */
    case hex20:
      if (actele->e.so3->gptot == 27)
      {
        actgid->is_solid3_h20_333 = 1;
        actgid->solid3_h20_333_name = "solid3_h20_333";
      }
      break;
    /* quadratic hexahedra */
    case hex27:
      if (actele->e.so3->gptot == 27)
      {
        actgid->is_solid3_h27_333 = 1;
        actgid->solid3_h27_333_name = "solid3_h27_333";
      }
      break;
    /* linear tetrahedra */
    case tet4:
      if (actele->e.so3->gptot == 1)
      {
        actgid->is_solid3_t4_1 = 1;
        actgid->solid3_t4_1_name = "solid3_t4_1";
      }
      break;
    /* quadratic tetrahedra */
    case tet10:
      if (actele->e.so3->gptot == 4)
      {
        actgid->is_solid3_t10_4 = 1;
        actgid->solid3_t10_4_name = "solid3_t10_4";
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
}  /* end of so3_gid_init */

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
void so3_gid_msh(FIELD *actfield,
                 INT ldis,
                 GIDSET *actgid,
                 FILE *out)
{
  ELEMENT *actele;
  INT j, k;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_msh");
#endif

  /*----------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_solid3_h8_222)
  {
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: SOLID3 HEX8 2x2x2 GP\n",
            actgid->solid3_h8_222_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",
            actgid->solid3_h8_222_name, ldis);
    /* print elements */
    fprintf(out, "ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_solid3) && (actele->numnp == 8) )
      {
        fprintf(out, " %6d ", actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out, "%6d ", actele->node[k]->Id+1);
        fprintf(out, "\n");
      }
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_solid3_h20_333)
  {
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# MESH %s FOR FIELD %s, DIS %1i: SOLID3 HEX20 3x3x3 GP\n",
            actgid->solid3_h20_333_name, actgid->fieldname, ldis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE %d\n",
            actgid->solid3_h20_333_name, ldis, actele->numnp);
    /* print elements */
    fprintf(out, "ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if ( (actele->eltyp == el_solid3) && (actele->numnp == 20) )
      {
        fprintf(out, " %6d ", actele->Id+1);
        for (k=0; k<actele->numnp; k++)
        {
          fprintf(out, "%6d ", actele->node[k]->Id+1);
        }
        fprintf(out, "\n");
      }
    }
    fprintf(out, "END ELEMENTS\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_gid_msh */


/*======================================================================*/
/*!
\brief Print header of Gauss point set to Gid output file

\param   *actele        ELEMENT     (i)   pointer to current element
\param    jdis          INT         (i)   discretisation index
\param   *actgid        GIDSET      (i)   Gid data (incl. GP sets)
\param   *out           FILE        (o)   Gid output file

\return void

\author bborn
\date 10/06
*/
void so3_gid_gpset(INT jdis,
                   GIDSET *actgid,
                   FILE *out)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_gpset");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_solid3_h8_222)
  {
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# GAUSSPOINTSET FOR FIELD %s DIS %1i SOLID3 2x2x2 GP\n",
            actgid->fieldname, jdis);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Hexahedra \"%s_dis_%1i\"\n",
            actgid->solid3_h8_222_name, jdis, 
            actgid->solid3_h8_222_name, jdis);
    fprintf(out, "NUMBER OF GAUSS POINTS: 8\n");
    fprintf(out, "NATURAL COORDINATES: Internal\n");
    fprintf(out, "END GAUSSPOINTS\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_solid3_h20_333)
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# GAUSSPOINTSET FOR FIELD %s DIS %1i SOLID3 3x3x3 GP\n",
            actgid->fieldname, jdis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"GAUSSPOINTS \"%s_dis_%1i\" ELEMTYPE Hexahedra \"%s_dis_%1i\"\n",
            actgid->solid3_h20_333_name, jdis,
            actgid->solid3_h20_333_name, jdis);
    fprintf(out, "NUMBER OF GAUSS POINTS: 27\n");
    fprintf(out, "NATURAL COORDINATES: Internal\n");
    fprintf(out, "END GAUSSPOINTS\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void so3_gid_gpset */



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
void so3_gid_dom(FIELD *actfield,
                 INT disnum,
                 GIDSET *actgid,
                 FILE *out)
{
  INT ngauss;
  INT i, j;  /* loop index */
  ELEMENT *actele;  /* current element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_dom");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_solid3_h8_222)
  {
    ngauss = 8;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->solid3_h8_222_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->solid3_h8_222_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_solid3) && (actele->numnp == 8) )
      {
        fprintf(out, "    %6d  %18.5E\n", actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)
        {
          fprintf(out, "            %18.5E\n", (DOUBLE)actele->proc);
        }
      }
    }
    fprintf(out, "END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_solid3_h20_333)
  {
    ngauss = 27;
    /* print header */
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "# RESULT Domains on DIS %1i, MESH %s\n",
            disnum, actgid->solid3_h20_333_name);
    fprintf(out, "#-------------------------------------------------------------------------------\n");
    fprintf(out, "RESULT \"Domains\" \"ccarat\" 0 SCALAR ONGAUSSPOINTS \"%s_dis_%1i\"\n",
            actgid->solid3_h20_333_name, disnum);
    /* print values */
    fprintf(out, "VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if ( (actele->eltyp == el_solid3) && (actele->numnp == 20) )
      {
        fprintf(out, "    %6d  %18.5E\n", actele->Id+1, (DOUBLE)actele->proc);
        for (j=1; j<ngauss; j++)
        {
          fprintf(out, "            %18.5E\n", (DOUBLE)actele->proc);
        }
      }
    }
    fprintf(out,"END VALUES\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void so3_gid_dom */




/*======================================================================*/
/*!
\brief Print stress in Gid output file

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
void so3_gid_stress(CHAR resstring[],
                    FIELD *actfield,
                    INT disnum,
                    INT step,
                    GIDSET *actgid,
                    FILE *out)
{
  const INT place = 0;

  INT ngauss;  /* number of Gauss points in element domain */
  INT nelenod;  /* number of element nodes */
  CHAR *resultname;
  CHAR *resulttype;
  CHAR *resultplace;
  CHAR *gpset = NULL;
  CHAR *rangetable;
  INT ncomponent;  /* number of output components */
  CHAR *componentnames[6];
  ELEMENT *actele;

  /* Gauss point permutations */
  INT gperm_h_222[8] = {0,4,6,2,1,5,7,3};
  INT gperm_h_333[27] = {0,10,24,6,2,20,26,8,
                         9,21,15,3,1,19,25,7,11,23,17,5,
                         12,10,22,16,4,14,13};
  INT gperm_t_4[4] = {0,1,2,3};
  INT gperm_t_10[10] = {0,1,2,3,4,5,6,7,8,9};  /* check this */
  INT *gperm;  /* pointer to current permutation */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* SOLID3s have 6 stress - use 3D matrix */
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_solid3_h8_222)
  {
    /* check only first element and assume, that the others are the same */
    actele = &(actfield->dis[disnum].element[0]);
    /* assure: we deal with SOLID3s */
    if (actele->eltyp != el_solid3)
    {
      dserror("All elements in discretisation must be of SOLID3 kind, "
              "otherwise Gid output fails.\n");
    }
    /* distinguish stress type */
    switch(actele->e.so3->stresstype)
    {
      /* equivalent stresses at nodes */
      case so3_stress_ndeqv:
        /* ~~~ */
        dserror("Equivalent stresses at nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "SCALAR";
        resultplace       = "ONNODES";  /* extrapolated to nodal points! */
        gpset             = actgid->solid3_h8_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 1;
        componentnames[0] = "equivStress";
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out, "COMPONENTNAMES \"%s\"\n", componentnames[0]);
        fprintf(out, "VALUES\n");
        /* extrapolated to nodal points! */
        /* so3_out_gid_sol_str(out,actfield,disnum,place,0); */
        fprintf(out, "END VALUES\n");
        break;
      /* xyz-oriented stresses at nodes */
      case so3_stress_ndxyz:
        /* ~~~ */
        dserror("XYZ-oriented stresses at nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->solid3_h8_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-XX";
        componentnames[1] = "Stress-YY";
        componentnames[2] = "Stress-ZZ";
        componentnames[3] = "Stress-XY";
        componentnames[4] = "Stress-YZ";
        componentnames[5] = "Stress-ZX";
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        fprintf(out, "VALUES\n");
        /* extrapolated to nodal points! */
        /* so3_out_gid_sol_str(out,actfield,disnum,place,0); */
        fprintf(out,"END VALUES\n");
        break;
      /* parameter-space rst-oriented stresses at nodes */
      case so3_stress_ndrst:
        /* ~~~ */
        dserror("rst-oriented stresses at nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /* extrapolated to nodal points! */
        gpset             = actgid->solid3_h8_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-rr";
        componentnames[1] = "Stress-ss";
        componentnames[2] = "Stress-tt";
        componentnames[3] = "Stress-rs";
        componentnames[4] = "Stress-st";
        componentnames[5] = "Stress-tr";
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        fprintf(out, "VALUES\n");
        /*extrapolated to nodal points!*/
        /*so3_out_gid_sol_str(out,actfield,disnum,place,0); */
        fprintf(out, "END VALUES\n");
        break;
      /* xyz-oriented stresses at Gauss points */
      case so3_stress_gpxyz:
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONGAUSSPOINTS";
        gpset             = actgid->solid3_h8_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-XX";
        componentnames[1] = "Stress-YY";
        componentnames[2] = "Stress-ZZ";
        componentnames[3] = "Stress-XY";
        componentnames[4] = "Stress-YZ";
        componentnames[5] = "Stress-ZX";
        nelenod           = 8;
        gperm             = &(gperm_h_222[0]);
        ngauss            = 8;
        so3_gid_stress_gp(actfield, disnum, actele, actgid, 
                          so3_stress_gpxyz, resultname, step,
                          resulttype, resultplace, gpset,
                          ncomponent, componentnames,
                          nelenod, gperm, ngauss, place, out);
        break;
      /* catch the remaining stress types */
      default:
        fprintf(out,"no stresses available\n");
        break;
    }  /* end switch */
  }  /* end if */
  /*--------------------------------------------------------------------*/
  /* hexahderon element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_solid3_h20_333)
  {
    /*check only first element and assume, that the others are the same*/
    actele = &(actfield->dis[disnum].element[0]);
    switch(actele->e.so3->stresstype)
    {
      /* equivalent stress */
      case so3_stress_ndeqv:
        /* ~~~ */
        dserror("Equivalent stresses at element nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "SCALAR";
        resultplace       = "ONNODES";  /* extrapolated to nodal points! */
        gpset             = actgid->solid3_h20_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 1;
        componentnames[0] = "equivStress";
        /* print header */
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out, "COMPONENTNAMES \"%s\"\n", componentnames[0]);
        /* print values */
        fprintf(out, "VALUES\n");
         /*extrapolated to nodal points!*/
        /* so3_out_gid_sol_str(out,actfield,disnum,place,0);*/
        fprintf(out, "END VALUES\n");
        break;
      /* xyz-oriented stresses at nodes */
      case so3_stress_ndxyz:
        /* ~~~ */
        dserror("xyz-oriented stresses at global nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->solid3_h20_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-xx";
        componentnames[1] = "Stress-yy";
        componentnames[2] = "Stress-zz";
        componentnames[3] = "Stress-xy";
        componentnames[4] = "Stress-xz";
        componentnames[5] = "Stress-yz";
        /* print header */
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out, "#-------------------------------------------------------------------------------\n");
        fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out, "COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        /* print values */
        fprintf(out, "VALUES\n");
         /*extrapolated to nodal points!*/
        /*so3_out_gid_sol_str(out,actfield,disnum,place,0);*/
        fprintf(out, "END VALUES\n");
        break;
      /* parameter-space rst-oriented stresses at nodes */
      case so3_stress_ndrst:
        /* ~~~ */
        dserror("rst-oriented stresses at global nodes are not available!\n");
        /* ~~~ */
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /* extrapolated to nodal points! */
        gpset             = actgid->solid3_h20_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-rr";
        componentnames[1] = "Stress-ss";
        componentnames[2] = "Stress-tt";
        componentnames[3] = "Stress-rs";
        componentnames[4] = "Stress-st";
        componentnames[5] = "Stress-tr";
        /* print header */
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out,"COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        /* print values */
        fprintf(out,"VALUES\n");
        /*extrapolated to nodal points!*/
        /*so3_out_gid_sol_str(out,actfield,disnum,place,0);*/
        fprintf(out,"END VALUES\n");
        break;
      /* xyz-oriented stresses at Gauss points */
      case so3_stress_gpxyz:
        resultname        = "solid3_stress";
        resulttype        = "MATRIX";
        resultplace       = "ONGAUSSPOINTS";
        gpset             = actgid->solid3_h20_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-XX";
        componentnames[1] = "Stress-YY";
        componentnames[2] = "Stress-ZZ";
        componentnames[3] = "Stress-XY";
        componentnames[4] = "Stress-YZ";
        componentnames[5] = "Stress-ZX";
        nelenod           = 20;
        gperm             = &(gperm_h_333[0]);
        ngauss            = 27;
        /* print */
        so3_gid_stress_gp(actfield, disnum, actele, actgid, 
                          so3_stress_gpxyz, resultname, step,
                          resulttype, resultplace, gpset,
                          ncomponent, componentnames,
                          nelenod, gperm, ngauss, place, out);
        break;
      /* catch the remains */
      default:
        fprintf(out,"no stresses available\n");
        break;
    }  /* end switch */
  }  /* end if */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void so3_gid_stress */


/*======================================================================*/
/*!
\brief Print stress in Gid output file

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
void so3_gid_stress_gp(FIELD *actfield,
                       INT disnum,
                       ELEMENT *actele,
                       GIDSET *actgid,
                       SO3_STRESSOUT stresstype,
                       CHAR *resultname,
                       INT step,
                       CHAR *resulttype,
                       CHAR *resultplace,
                       CHAR *gpset,
                       INT ncomp,
                       CHAR *componentnames[],
                       INT nelenod,
                       INT *gperm,
                       INT ngauss,
                       INT place,
                       FILE *out)
{
  INT icomp;
  INT iele;
  INT igp, jgp;
  DOUBLE **stress;
  

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_gid_stress_gpxyz");
#endif

  /*--------------------------------------------------------------------*/
  /* print header */
  fprintf(out, "#-------------------------------------------------------------------------------\n");
  fprintf(out, "# RESULT %s on FIELD %s, DIS %1i\n",
          resultname, actgid->fieldname, disnum);
  fprintf(out, "#-------------------------------------------------------------------------------\n");
  fprintf(out, "RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
          resultname, step, resulttype, resultplace, gpset,disnum);
  fprintf(out, "COMPONENTNAMES ");
  for (icomp=0; icomp<ncomp; ncomp++)
  {
    if (icomp == ncomp-1)
    {
      fprintf(out, "\"%s\"\n", componentnames[icomp]);
    }
    else
    {
      fprintf(out, "\"%s\",", componentnames[icomp]);
    }
  }
  /* print values */
  fprintf(out, "VALUES\n");
  for (iele=0; iele<actfield->dis[disnum].numele; iele++)
  {
    /* pointer to current element */
    actele = &(actfield->dis[disnum].element[iele]);
    /* print heat fluxes */
    if ( (actele->eltyp == el_solid3) && (actele->numnp == nelenod) )
    {
      if (stresstype == so3_stress_gpxyz)
      {
        stress = actele->e.so3->stress_gpxyz.a.d3[place];
      }
      else if (stresstype == so3_stress_gprst)
      {
        stress = actele->e.so3->stress_gprst.a.d3[place];
      }
      else
      {
        dserror("stresstype is not printable to Gid\n");
      }
      /* print values at first Gauss point _with_ element index */
      igp = gperm[0];  /* get corresponding Gauss point index
                        * of 0th Gid Gauss point in
                        * SOLID3 Gauss point numbering */
      fprintf(out, " %6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
              actele->Id+1,
              stress[igp][0],  /* stress xx-component */
              stress[igp][1],  /* stress yy-component */
              stress[igp][2],  /* stress zz-component */
              stress[igp][3],  /* stress xy-component */
              stress[igp][4],  /* stress yz-component */
              stress[igp][5]);  /* stress zx-component */
      /* print values at 2nd to last Gauss point */
      for (jgp=1; jgp<ngauss; jgp++)
      {
        igp = gperm[jgp];
        fprintf(out, "        %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
                stress[igp][0],  /* stress xx-compo */
                stress[igp][1],  /* stress yy-compo */
                stress[igp][2],  /* stress zz-compo */
                stress[igp][3],  /* stress xy-compo */
                stress[igp][4],  /* stress yz-compo */
                stress[igp][5]);  /* stress zx-compo */
      }
    }
  }
  fprintf(out, "END VALUES\n");

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void so3_gid_stress_gp */

#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
