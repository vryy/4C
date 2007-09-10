/*======================================================================*/
/*!
\file
\brief Print BRICK1 element data to Ccarat Gid output file

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
#ifndef CCADISCRET
#ifdef D_BRICK1

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "../output/gid.h"
#include "brick1.h"

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
void c1_gid_init(ELEMENT *actele,
                 GIDSET *actgid)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_gid_init");
#endif

  /*--------------------------------------------------------------------*/
  if (actele->numnp==8)
  {
    actgid->is_brick1_222   = 1;
    actgid->brick1_222_name = "brick1_222";
  }
  if (actele->numnp==20 || actele->numnp==27)
  {
    actgid->is_brick1_333   = 1;
    actgid->brick1_333_name = "brick1_333";
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of c1_gid_init */

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
void c1_gid_msh(FIELD *actfield,
                INT ldis,
                GIDSET *actgid,
                FILE *out)
{
  ELEMENT *actele;
  INT j, k;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_gid_msh");
#endif

  /*----------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_brick1_222)
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# MESH %s FOR FIELD %s, DIS %1i: BRICK1 2x2x2 GP\n",
            actgid->brick1_222_name,actgid->fieldname,ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",
            actgid->brick1_222_name,ldis);
    /*------------------------------------------------ print elements */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if (actele->eltyp != el_brick1 || actele->numnp !=8) continue;
      fprintf(out," %6d ",actele->Id+1);
      for (k=0; k<actele->numnp; k++)
        fprintf(out,"%6d ",actele->node[k]->Id+1);
      fprintf(out,"\n");
    }
    fprintf(out,"END ELEMENTS\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_brick1_333)
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# MESH %s FOR FIELD %s, DIS %1i: BRICK1 3x3x3 GP as HEX8!\n",
            actgid->brick1_333_name,actgid->fieldname,ldis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    /*fprintf(out,"MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->brick1_333_name,ldis);*/
    fprintf(out,"MESH %s_dis_%1i DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",
            actgid->brick1_333_name,ldis);
    /*------------------------------------------------ print elements */
    /* gid shows Hex20 as Hex with eight nodes only,
       so the output is reduced */
    fprintf(out,"ELEMENTS\n");
    for (j=0; j<actfield->dis[ldis].numele; j++)
    {
      actele = &(actfield->dis[ldis].element[j]);
      if (actele->eltyp != el_brick1 || actele->numnp !=20) continue;
      fprintf(out," %6d ",actele->Id+1);
      /*for (k=0; k<actele->numnp; k++)*/
      for (k=0; k<8; k++)
        fprintf(out,"%6d ",actele->node[k]->Id+1);
      fprintf(out,"\n");
    }
    fprintf(out,"END ELEMENTS\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of c1_gid_msh */


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
void c1_gid_gpset(INT jdis,
                  GIDSET *actgid,
                  FILE *out)
{
  const char sign='"';

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_gid_gpset");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_brick1_222)
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# GAUSSPOINTSET FOR FIELD %s DIS %1i BRICK1 2x2x2 GP\n",actgid->fieldname,jdis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"GAUSSPOINTS %c%s_dis_%1i%c ELEMTYPE Hexahedra %c%s_dis_%1i%c\n",
            sign,actgid->brick1_222_name,jdis,sign,
            sign,actgid->brick1_222_name,jdis,sign);
    fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
    fprintf(out,"NATURAL COORDINATES: Internal\n");
    fprintf(out,"END GAUSSPOINTS\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_brick1_333)
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# GAUSSPOINTSET FOR FIELD %s DIS %1i BRICK1 3x3x3 GP\n",actgid->fieldname,jdis);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"GAUSSPOINTS %c%s_dis_%1i%c ELEMTYPE Hexahedra %c%s_dis_%1i%c\n",
            sign,actgid->brick1_333_name,jdis,sign,
            sign,actgid->brick1_333_name,jdis,sign);
    fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
    fprintf(out,"NATURAL COORDINATES: Internal\n");
    fprintf(out,"END GAUSSPOINTS\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void c1_gid_gpset */



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
void c1_gid_dom(FIELD *actfield,
                INT disnum,
                GIDSET *actgid,
                FILE *out)
{
  const char sign='"';
  INT ngauss;
  INT i, j;  /* loop index */
  ELEMENT *actele;  /* current element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_gid_dom");
#endif

  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_brick1_222)
  {
    ngauss = 8;
    /* print header */
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT Domains on DIS %1i, MESH %s\n"
            ,disnum,actgid->brick1_222_name);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s_dis_%1i%c\n",
            sign,sign,sign,sign,sign,actgid->brick1_222_name,disnum,sign);
    /* print values */
    fprintf(out,"VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if (actele->eltyp != el_brick1 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<ngauss; j++)
        fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
    }
    fprintf(out,"END VALUES\n");
  }
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_brick1_333)
  {
    ngauss = 27;
    /* print header */
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT Domains on DIS %1i, MESH %s\n",disnum,actgid->brick1_333_name);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s_dis_%1i%c\n",
            sign,sign,sign,sign,sign,actgid->brick1_333_name,disnum,sign);
    /* print values */
    fprintf(out,"VALUES\n");
    for (i=0; i<actfield->dis[disnum].numele; i++)
    {
      actele = &(actfield->dis[disnum].element[i]);
      if (actele->eltyp != el_brick1 || (actele->numnp != 20 || actele->numnp != 27)) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<ngauss; j++)
        fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
    }
    fprintf(out,"END VALUES\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void c1_gid_dom */




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
void c1_gid_stress(char resstring[],
                   FIELD *actfield,
                   INT disnum,
                   INT step,
                   GIDSET *actgid,
                   FILE *out)
{
  const INT place = 0;

  INT ngauss;
  INT ist;
  char *resultname;
  char *resulttype;
  char *resultplace;
  char *gpset = NULL;
  char *rangetable;
  INT  ncomponent;
  char *componentnames[6];
  ELEMENT *actele;
  INT i;

  /* INT gaussperm8[8] = {0,4,2,6,1,5,3,7}; */
  INT gperm_222[8] = {4,6,2,0,5,7,3,1};  /* verified */
  /* INT gaussperm27[27] = {0,9,18,3,12,21,6,15,24,
                         1,10,19,4,13,22,7,16,25,
                         2,11,20,5,14,23,8,17,26}; */
  INT gperm_333[27] = {18,24,6,0,20,26,8,2,
                       21,15,3,9,19,25,7,1,23,17,5,11,
                       12,22,16,4,10,14,13};

  DOUBLE **stress;

  INT jgp;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_gid_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* bricks have 6 stress - use 3D matrix */
  /*--------------------------------------------------------------------*/
  /* hexahedron element with 8 nodes and 2x2x2 Gauss points */
  if (actgid->is_brick1_222)
  {
    /*check only first element and assume, that the others are the same*/
    actele = &(actfield->dis[disnum].element[0]);
    /* assure: we deal with bricks */
    if (actele->eltyp != el_brick1)
    {
      dserror("All elements in discretisation must be of BRICK1 kind, "
              "otherwise Gid output fails.\n");
    }
    /* distinguish stress type */
    switch(actele->e.c1->stresstyp)
    {
      /* equivalent stresses at nodes */
      case c1_npeqs:
        resultname        = "brick1_forces";
        resulttype        = "SCALAR";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 1;
        componentnames[0] = "equivStress";
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out,"COMPONENTNAMES \"%s\"\n", componentnames[0]);
        fprintf(out,"VALUES\n");
        /*extrapolated to nodal points!*/
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* xyz-oriented stresses at nodes */
      case c1_npxyz:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-xx";
        componentnames[1] = "Stress-yy";
        componentnames[2] = "Stress-zz";
        componentnames[3] = "Stress-xy";
        componentnames[4] = "Stress-xz";
        componentnames[5] = "Stress-yz";
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out,"COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        fprintf(out,"VALUES\n");
        /*extrapolated to nodal points!*/
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* parameter-space rst-oriented stresses at nodes */
      case c1_nprst:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-rr";
        componentnames[1] = "Stress-ss";
        componentnames[2] = "Stress-tt";
        componentnames[3] = "Stress-rs";
        componentnames[4] = "Stress-st";
        componentnames[5] = "Stress-tr";
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname,disnum);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out,"COMPONENTNAMES \"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
                componentnames[0], componentnames[1], componentnames[2],
                componentnames[3], componentnames[4], componentnames[5]);
        fprintf(out,"VALUES\n");
        /*extrapolated to nodal points!*/
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* xyz-oriented stresses at Gauss points */
      case c1_gpxyz:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONGAUSSPOINTS";
        gpset             = actgid->brick1_222_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-xx";
        componentnames[1] = "Stress-yy";
        componentnames[2] = "Stress-zz";
        componentnames[3] = "Stress-xy";
        componentnames[4] = "Stress-xz";
        componentnames[5] = "Stress-yz";
        ngauss            = 8;
        ist               = 6;  /* first stress index in stress_GP array */
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
        for (i=0; i<actfield->dis[disnum].numele; i++)
        {
          /* pointer to current element */
          actele = &(actfield->dis[disnum].element[i]);
          /* print heat fluxes */
          if ( (actele->eltyp == el_brick1) && (actele->numnp == 8) )
          {
            stress = actele->e.c1->stress_GP.a.d3[place];
            /* print values at first Gauss point _with_ element index */
            fprintf(out, " %6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
                    actele->Id+1,
                    stress[ist+0][gperm_222[0]],  /* stress xx-component */
                    stress[ist+1][gperm_222[0]],  /* stress yy-component */
                    stress[ist+2][gperm_222[0]],  /* stress zz-component */
                    stress[ist+3][gperm_222[0]],  /* stress xy-component */
                    stress[ist+4][gperm_222[0]],  /* stress yz-component */
                    stress[ist+5][gperm_222[0]]);  /* stress zx-component */
            /* print values at 2nd to last Gauss point */
            for (jgp=1; jgp<ngauss; jgp++)
            {
              fprintf(out, "        %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
                      stress[ist+0][gperm_222[jgp]],  /* stress xx-compo */
                      stress[ist+1][gperm_222[jgp]],  /* stress yy-compo */
                      stress[ist+2][gperm_222[jgp]],  /* stress zz-compo */
                      stress[ist+3][gperm_222[jgp]],  /* stress xy-compo */
                      stress[ist+4][gperm_222[jgp]],  /* stress yz-compo */
                      stress[ist+5][gperm_222[jgp]]);  /* stress zx-compo */
            }
          }
        }
        fprintf(out, "END VALUES\n");
        break;
      /* catch the remaining stress types */
      default:
        fprintf(out,"no stresses available\n");
        break;
    }  /* end switch */
  }  /* end if */
  /*--------------------------------------------------------------------*/
  /* hexahderon element with 20 nodes and 3x3x3 Gauss points */
  if (actgid->is_brick1_333)
  {
    /*check only first element and assume, that the others are the same*/
    actele = &(actfield->dis[disnum].element[0]);
    switch(actele->e.c1->stresstyp)
    {
      /* equivalent stress */
      case c1_npeqs:
        resultname        = "brick1_forces";
        resulttype        = "SCALAR";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 1;
        componentnames[0] = "equivStress";
        /* print header */
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT %s on FIELD %s, DIS %1i\n",
                resultname, actgid->fieldname, disnum);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT \"%s\" \"ccarat\" %d %s %s \"%s_dis_%1i\"\n",
                resultname, step, resulttype, resultplace, gpset,disnum);
        fprintf(out,"COMPONENTNAMES \"%s\"\n", componentnames[0]);
        /* print values */
        fprintf(out,"VALUES\n");
         /*extrapolated to nodal points!*/
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* xyz-oriented stresses at nodes */
      case c1_npxyz:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-xx";
        componentnames[1] = "Stress-yy";
        componentnames[2] = "Stress-zz";
        componentnames[3] = "Stress-xy";
        componentnames[4] = "Stress-xz";
        componentnames[5] = "Stress-yz";
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
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* parameter-space rst-oriented stresses at nodes */
      case c1_nprst:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONNODES";  /*extrapolated to nodal points!*/
        gpset             = actgid->brick1_333_name;
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
        c1_out_gid_sol_str(out,actfield,disnum,place,0);
        fprintf(out,"END VALUES\n");
        break;
      /* xyz-oriented stresses at Gauss points */
      case c1_gpxyz:
        resultname        = "brick1_forces";
        resulttype        = "MATRIX";
        resultplace       = "ONGAUSSPOINTS";
        gpset             = actgid->brick1_333_name;
        rangetable        = actgid->standardrangetable;
        ncomponent        = 6;
        componentnames[0] = "Stress-xx";
        componentnames[1] = "Stress-yy";
        componentnames[2] = "Stress-zz";
        componentnames[3] = "Stress-xy";
        componentnames[4] = "Stress-xz";
        componentnames[5] = "Stress-yz";
        ngauss            = 27;
        ist               = 6;  /* 1st stress index in stress_GP array */
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
        for (i=0; i<actfield->dis[disnum].numele; i++)
        {
          /* pointer to current element */
          actele = &(actfield->dis[disnum].element[i]);
          /* print heat fluxes */
          if ( (actele->eltyp == el_brick1) && (actele->numnp == 20) )
          {
            stress = actele->e.c1->stress_GP.a.d3[place];
            /* print values at first Gauss point _with_ element index */
            fprintf(out, " %6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
                    actele->Id+1,
                    stress[ist+0][gperm_333[0]],  /* stress xx-component */
                    stress[ist+1][gperm_333[0]],  /* stress yy-component */
                    stress[ist+2][gperm_333[0]],  /* stress zz-component */
                    stress[ist+3][gperm_333[0]],  /* stress xy-component */
                    stress[ist+4][gperm_333[0]],  /* stress yz-component */
                    stress[ist+5][gperm_333[0]]);  /* stress zx-component */
            /* print values at 2nd to last Gauss point */
            for (jgp=1; jgp<ngauss; jgp++)
            {
              fprintf(out, "        %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E\n",
                      stress[ist+0][gperm_333[jgp]],  /* stress xx-compo */
                      stress[ist+1][gperm_333[jgp]],  /* stress yy-compo */
                      stress[ist+2][gperm_333[jgp]],  /* stress zz-compo */
                      stress[ist+3][gperm_333[jgp]],  /* stress xy-compo */
                      stress[ist+4][gperm_333[jgp]],  /* stress yz-compo */
                      stress[ist+5][gperm_333[jgp]]);  /* stress zx-compo */
            }
          }
        }
        fprintf(out, "END VALUES\n");
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
}  /* end void c1_gid_stress */



#endif /* end of #ifdef D_BRICK1 */
/*! @} (documentation module close)*/
#endif
