/*!----------------------------------------------------------------------
\file
\brief contains the routine 'saxi_cal_stress' which evaluates the element
       stresses for the axisymmmetric shell element

 *----------------------------------------------------------------------*/
#ifdef D_AXISHELL
#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*! 
\addtogroup AXISHELL 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  Computation of the stresses for the axialsymmetric shell element 

<pre>                                                              mn 05/03 
This routine computes the stresses for the axialsymmetric shell element
after computation of the displacements.


</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          SAXI_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void                                               
\sa calling:   saxiintg, saxi_B, saxi_mat;
    called by: axishell();

*----------------------------------------------------------------------*/
void saxi_cal_stress(ELEMENT   *ele, 
                     SAXI_DATA *data, 
                     MATERIAL  *mat,
		     INT        init)
{
INT                 i,j,k;                    /* some loopers           */
INT                 nixsi;                    /* num GP in xsi direction*/
INT                 lxsi;                     /* loopers over GP        */
INT                 iel;                      /* numnp to this element  */
INT                 nd;
INT                 istore = 0;/* controls storing of new stresses to wa*/

INT                 diff,max;

const INT           numdf  = 2;
const INT           numeps = 4;

DOUBLE              thick;                    /* thickness              */
DOUBLE              xsi;                      /*GP-coords               */
DOUBLE              r,dr,dz,dl,cosa,sina;     /* infos about the element*/
DOUBLE              fac;                      /* weights at GP          */

DOUBLE              disp_global[7];/* vector of global nodal displ.     */
DOUBLE              disp_local[7]; /* vector of local nodal displa.     */
DOUBLE              S[4][7];       /* stress matrix S=D*B(xsi=0.5)      */
DOUBLE              n[5];          /* Schnittkraftvektor                */
DOUBLE              T[6][6];       /* Transformation matrix             */
DOUBLE              Bs[4][7];      /* Ableitung der Matrix B nach xsi   */
DOUBLE              DBs[7];        /* Multiplication of D and Bs        */


static ARRAY    D_a;               /* material tensor                   */     
static DOUBLE **D;  
static ARRAY    work1_a;           /* work1 = D * B                     */   
static DOUBLE **work1;                
static ARRAY    work2_a;           /* work2 = B^T * work1               */   
static DOUBLE **work2;                       
static ARRAY    B_a;               /* B-operator                        */   
static DOUBLE **B;

DOUBLE ***gp_stress;   /* pointer to array of stresses at the middle GP */
DOUBLE * statcond;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("saxi_cal_stress");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
D         = amdef("D"      ,&D_a     ,6,6             ,"DA");           
work1     = amdef("work1"  ,&work1_a ,4,7,"DA");   
work2     = amdef("work2"  ,&work2_a ,7,7,"DA");                     
B         = amdef("B"      ,&B_a     ,numeps,7,"DA");                     
goto end;
}
/*------------------------------------------- integration parameters ---*/
saxiintg(data);
/*------------------------------------------- integration parameters ---*/
nixsi   = 5; /* number of Gauss-points */
iel     = ele->numnp;
nd      = numdf * iel;
/* -------------------------datas about the element under consideration */
dr = ele->node[1]->x[0] - ele->node[0]->x[0]; /* delta r */
dz = ele->node[1]->x[1] - ele->node[0]->x[1]; /* delta z */
dl = sqrt(dr*dr+dz*dz);                       /* delta l */
cosa = dr/dl;
sina = dz/dl;
/*------------------ store node displacements (=local) on local field --*/
/* call global displacements */


/* ------------------------------------ compute  transfomation matrix T */
for (i=0; i<6; i++)
{
  for (j=0; j<6; j++)
  {
    T[i][j] = 0.0;
  }
}
T[0][0] =  cosa;
T[0][1] =  sina;
T[1][0] =  sina;
T[1][1] = -cosa;
T[2][2] =  -1.0;
T[3][3] =  cosa;
T[3][4] =  sina;
T[4][3] =  sina;
T[4][4] = -cosa;
T[5][5] =  -1.0;


/* ------------------------------------------------------------- node 0 */
if (ele->node[0]->gnode->ondesigntyp==ondnode && ele->node[0]->gnode->d.dnode->cos_type == 1)
{
  /* local cos */
  disp_local[0] = ele->node[0]->sol.a.da[0][0]; 
  disp_local[1] = ele->node[0]->sol.a.da[0][1];
  disp_local[2] = ele->node[0]->sol.a.da[0][2];
  /* calculate global disps */
  for (i=0; i<3; i++)
  {
    disp_global[i] = 0.0;
    for (j=0; j<3; j++)
    {
      disp_global[i] += T[j][i] * disp_local[j];
    }
  }
}
else
{
  /* global cos */
  disp_global[0] = ele->node[0]->sol.a.da[0][0]; 
  disp_global[1] = ele->node[0]->sol.a.da[0][1];
  disp_global[2] = ele->node[0]->sol.a.da[0][2];
  /* calculate local disps */
  for (i=0; i<3; i++)
  {
    disp_local[i] = 0.0;
    for (j=0; j<3; j++)
    {
      disp_local[i] += T[i][j] * disp_global[j];
    }
  }
}
/* ------------------------------------------------------------- node 1 */
if (ele->node[1]->gnode->ondesigntyp==ondnode && ele->node[1]->gnode->d.dnode->cos_type == 1)
{
  /* local cos */
  disp_local[3] = ele->node[1]->sol.a.da[0][0];
  disp_local[4] = ele->node[1]->sol.a.da[0][1];
  disp_local[5] = ele->node[1]->sol.a.da[0][2];
  /* calculate global disps */
  for (i=3; i<6; i++)
  {
    disp_global[i] = 0.0;
    for (j=3; j<6; j++)
    {
      disp_global[i] += T[j][i] * disp_local[j];
    }
  }
}
else
{
  /* global cos */
  disp_global[3] = ele->node[1]->sol.a.da[0][0];
  disp_global[4] = ele->node[1]->sol.a.da[0][1];
  disp_global[5] = ele->node[1]->sol.a.da[0][2];
  /* calculate local disps */
  for (i=3; i<6; i++)
  {
    disp_local[i] = 0.0;
    for (j=3; j<6; j++)
    {
      disp_local[i] += T[i][j] * disp_global[j];
    }
  }
}

/* ===================================== write global dis. to the nodes */
ele->node[0]->sol.a.da[0][0] = disp_global[0];
ele->node[0]->sol.a.da[0][1] = disp_global[1];
ele->node[0]->sol.a.da[0][2] = disp_global[2];
ele->node[1]->sol.a.da[0][0] = disp_global[3];
ele->node[1]->sol.a.da[0][1] = disp_global[4];
ele->node[1]->sol.a.da[0][2] = disp_global[5];

/* ====================================== write local dis. to the nodes */
/* ------------------------------------------ enlarge sol, if necessary */
if (1 >= ele->node[0]->sol.fdim)
{
  diff = 1 - ele->node[0]->sol.fdim;
  max  = IMAX(diff,3);
  amredef(&(ele->node[0]->sol),ele->node[0]->sol.fdim+max+1,ele->node[0]->sol.sdim,"DA");
}
if (1 >= ele->node[1]->sol.fdim)
{
  diff = 1 - ele->node[1]->sol.fdim;
  max  = IMAX(diff,3);
  amredef(&(ele->node[1]->sol),ele->node[1]->sol.fdim+max+1,ele->node[1]->sol.sdim,"DA");
}

ele->node[0]->sol.a.da[1][0] = disp_local[0];
ele->node[0]->sol.a.da[1][1] = disp_local[1];
ele->node[0]->sol.a.da[1][2] = disp_local[2];
ele->node[1]->sol.a.da[1][0] = disp_local[3];
ele->node[1]->sol.a.da[1][1] = disp_local[4];
ele->node[1]->sol.a.da[1][2] = disp_local[5];




/* compute mid-side displacement */
statcond = ele->e.saxi->statcond;
disp_local[6]= 0.0;
for (i=0; i<6; i++)
{
  disp_local[6] += statcond[i]*disp_local[i];
}
disp_local[6] -= ele->e.saxi->Sk;


/*================================== loop over GP (=Auswertungsstellen) */
for (lxsi=2; lxsi<3; lxsi++) /* Auswertung nur bei xsi=0.5 */
{
   /*==================== gaussian point, weight and thickness at it ===*/
   xsi   = data->xgr[lxsi];                  /* xsi       */
   fac   = data->wgt[lxsi];                  /* weight    */
   r     = ele->node[0]->x[0] + cosa*xsi*dl; /* radius    */
   thick = ele->e.saxi->thick[0] +           /* thickness */
           xsi/1.0*(ele->e.saxi->thick[1]-ele->e.saxi->thick[0]);
   /*------------------------------------------ calculate operator B ---*/
   amzero(&B_a);
   saxi_B(B,xsi,r,dl,cosa,sina);
   /*--------------------------------- compute constitutive matrix D ---*/
   saxi_mat(mat->m.stvenant,D,thick);
   /*------------------------------------- compute stress matrix D*B ---*/
   for (k=0; k<7; k++)
   {
     for (i=0; i<4; i++)
     {
       S[i][k] = 0.0;
       for (j=0; j<4; j++)
       {
         S[i][k] += D[i][j] * B[j][k];
       }
     }
   }
   /*----------------------------------------- compute schnittkraefte n */
   /*------------------------------------------------------------------ */
   
   /*---------------- normal forces n_s, n_deta and moments m_s, m_deta */
   for (i=0; i<4; i++)
   {    
     n[i] = 0.0;
     for (j=0; j<7; j++)
     {
       n[i] += S[i][j] * disp_local[j];
     }
   }
   /*---------------------------------- Rueckrechnung der Querkraft q_s */
   /* Ermittlung der Ableitung von ms nach ds                           */
   /* Bs[i][j]...Ableitung der Matrix B nach xsi (nur 3. und 4. Zeile)  */
   
   for (i=0; i<4; i++)
   {
     for (j=0; j<7; j++)
     {
       Bs[i][j] = 0.0;
     }
   }    
   Bs[2][1] = 12.0/dl/dl;
   Bs[2][2] = 6.0/dl;
   Bs[2][4] = -Bs[2][1];
   Bs[2][5] = 6.0/dl;
   Bs[3][1] = 6.0*(-1.0+2.0*xsi)*cosa/(dl*r) - 
              6.0*cosa*xsi*(-1.0+xsi)*dr/dl/r/r;
   Bs[3][2] = (-4.0+6.0*xsi)*cosa/r - 
              cosa/r/r*dr*dl*(1.0-4.0*xsi+3.0*xsi*xsi);
   Bs[3][4] = -B[3][1];
   Bs[3][5] = (-2.0+6.0*xsi)*cosa/r - 
              cosa/r/r*dr*dl*(-2.0*xsi+3.0*xsi*xsi);
   /* since we have to derive B with respect to s we have to divide the
      derivative of B with respect to xsi additionally bz dl !!!        */
   for (i=0; i<4; i++)
   {
     for (j=0; j<7; j++)
     {
       Bs[i][j] /= dl;
     }
   }    

   /* Multiplication of D with Bs */
   for (i=0; i<7; i++)
   {
     DBs[i] = 0.0;
     for (j=2; j<4; j++) /* since D[2][0] and D[2][1] are 0 */
     {
        DBs[i] += D[2][j] * Bs[j][i];
     }
   }    
   /* Computation of the first part of q_s = DBs * disp_local */
   n[4] = 0.0;
   for (i=0; i<7; i++) 
   {
      n[4] += DBs[i] * disp_local[i];
   }
   /* final solution of q_s */
   n[4] += (n[2]-n[3])*cosa/r;
}/*=========================== end of loop over GP (=Auswertungsstellen) */

/* ---------store the local n-vector onto stress_GP for postprocessing */
gp_stress = ele->e.saxi->stress_GP.a.d3;
for (i=0; i<5; i++) 
{
   gp_stress[0][i][0] = n[i];
}

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of saxi_cal_stress */

/*----------------------------------------------------------------------*/
#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
