#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | compute Kroneker-Delta (2-stufiger Einheitstensor)      he    04/03   |
 *----------------------------------------------------------------------*/
void w1_kroneker(DOUBLE delta[3][3])
{
INT    i,j;

#ifdef DEBUG
dstrc_enter("w1_kroneker");
#endif
/*----------------------------------------------------------------------*/

for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  if (i==j) delta[i][j] = 1.0;
  if (i!=j) delta[i][j] = 0.0;
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_kroneker */

/*----------------------------------------------------------------------*
 | compute ela. stresses with Hook                        he    04/03   |
 *----------------------------------------------------------------------*/
void w1_stress_ela(DOUBLE ym,
                   DOUBLE pv,
                   DOUBLE epsilon[3][3],
                   DOUBLE sigma_el[3][3],
                   DOUBLE delta[3][3]
                   )
{
INT    i,j;
DOUBLE mu,lam,trace;
#ifdef DEBUG
dstrc_enter("w1_stress_ela");
#endif
/*-------------------- elast. variables and trace from epsilon --------*/
mu    = ym/(2.0*(1.0+ pv));
lam   = ym*pv/((1.0-2.0*pv)*(1.0+pv));
trace = epsilon[0][0]+epsilon[1][1]+epsilon[2][2];

for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  sigma_el[i][j] = 2.0*mu*epsilon[i][j] + lam*trace*delta[i][j];
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_stress_ela */


/*----------------------------------------------------------------------*
 | compute ela. stiffness tensor                           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_mat_ela(DOUBLE   ym,
                DOUBLE   pv,
                DOUBLE   delta[3][3],
                DOUBLE   c_el[3][3][3][3]
                )
{
INT    i,j,k,l;
DOUBLE mu,lam;
#ifdef DEBUG
dstrc_enter("w1_mat_ela");
#endif
/*-------------------- elast. variables and trace from epsilon --------*/
mu    = ym/(2.0*(1.0+ pv));
lam   = ym*pv/((1.0-2.0*pv)*(1.0+pv));

for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  for (k=0; k<3; k++)
  {
   for (l=0; l<3; l++)
   {
    c_el[i][j][k][l] = mu*(delta[i][k]*delta[j][l]+delta[i][l]*delta[j][k])+
                       lam*delta[i][j]*delta[k][l];
   }
  }
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_mat_ela */

/*----------------------------------------------------------------------*
 | compute actuel equivalent strains and                  he    04/03   |
 | their global derivatives of epsilon                                  |
 | Equival = 1 -  ENERGY RELEASE RATE CONCEPT - THERMODYNAMICS          |
 | Equival = 2 -  SIMO & JU (1987)                                      |
 | Equival = 3 -  JU (1989)                                             |
 | Equival = 4 -  de VREE (1995)                                        |
 *----------------------------------------------------------------------*/
void w1_equi_eps(DOUBLE   epsilon[3][3],
                 DOUBLE   sigma_el[3][3],
                 DOUBLE   delta[3][3],
                 INT      Equival,
                 DOUBLE   *eta,
                 DOUBLE   eta_der[3][3],
                 DOUBLE   pv,
                 DOUBLE   k
                 )
{
INT    i,j;
DOUBLE depsilon[3][3];
DOUBLE I_1,J_2,epseps;
DOUBLE K0,K1,K2,K3,ROOT;
#ifdef DEBUG
dstrc_enter("w1_equi_eps");
#endif
/*----------------------------------------------------------------------*/
*eta = 0.0;

/*------------------- EQUIVALENT STRAIN - BASED ON THERMODYNAMIX --------*/
if (Equival == 1)
{
 for(i=0; i<3; i++)
 {
  for(j=0; j<3; j++)
  {
   *eta += 0.5*epsilon[i][j]*sigma_el[i][j];
   eta_der[i][j] = sigma_el[i][j];
  }
 }
} /* if */

/*----------------- EQUIVALENT STRAIN - BASED ON SIMO & JU (1987) ------*/
if (Equival == 2)
{
 for(i=0; i<3; i++)
 {
  for(j=0; j<3; j++)
  {
   *eta += epsilon[i][j]*sigma_el[i][j];
  }
 }
 *eta = sqrt(*eta);

 for(i=0; i<3; i++)
 {
  for(j=0; j<3; j++)
  {
   eta_der[i][j] = sigma_el[i][j]/(*eta);
  }
 }
} /* if */

/*----------------- EQUIVALENT STRAIN - BASED ON  JU (1989) -----------*/
if (Equival == 3)
{
 for(i=0; i<3; i++)
 {
  for(j=0; j<3; j++)
  {
   *eta += 0.5*epsilon[i][j]*sigma_el[i][j];
  }
 }
 *eta = sqrt(*eta);

 for(i=0; i<3; i++)
 {
  for(j=0; j<3; j++)
  {
   eta_der[i][j] = 0.5*sigma_el[i][j]/(*eta);
  }
 }
} /* if */

/*------------ EQUIVALENT STRAIN - BASED ON  de VREE (1995) -----------*/
if (Equival == 4)
{
/*----------------------- calculate Invariant I_1 = tr(eps) -----------*/
  I_1 = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];

/*--- calculate Invariant J_2 = 0.5[eps:eps - 1/3*tr(eps)*tr(eps)] ----*/
  epseps = 0.0;
  for(i=0; i<3; i++)
  {
   for(j=0; j<3; j++)
   {
    epseps += epsilon[i][j]*epsilon[i][j];
   }
  }
  J_2 = 0.5*(epseps - I_1*I_1/3.0);

/*------------------------------------------- calculate constants ----*/
  K0 = (k-1.0)/(2.0*k*(1.0-2.0*pv));
  K1 = (k-1.0)/(1.0-2.0*pv);
  K2 = (12.0*k)/(1.0+pv)/(1.0+pv);
  K3 = (1.0)/(2.0*k);

  ROOT = sqrt(K1*K1*I_1*I_1+K2*J_2);
 *eta  = K0*I_1+K3*ROOT;

  for(i=0; i<3; i++)
  {
   for(j=0; j<3; j++)
   {
    depsilon[i][j] = epsilon[i][j];
    if(i==j) depsilon[i][j]=depsilon[i][j]-I_1/3.0;
   }
  }

  for(i=0; i<3; i++)
  {
   for(j=0; j<3; j++)
   {
    if(ROOT != 0.0)
    {
     eta_der[i][j] = K0*delta[i][j]+K3/(2.0*ROOT)
                    *(2.0*K1*K1*I_1*delta[i][j]+K2*depsilon[i][j]);
    }
    else
    {
     eta_der[i][j] = 0.0;
    }
   }
  }

} /* if */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_equi_eps */


/*----------------------------------------------------------------------*
 | create tensor from vector !!!!!!only for strains!!!!   he    04/03   |
 *----------------------------------------------------------------------*/
void w1_4to9(DOUBLE  *vector,
             DOUBLE   tensor[3][3]
             )
{
#ifdef DEBUG
dstrc_enter("w1_4to9");
#endif
/*----------------------------------------------------------------------*/

   tensor[0][0] = vector[0];
   tensor[0][1] = vector[2]/2.0;
   tensor[0][2] = 0.0;
   tensor[1][0] = vector[2]/2.0;
   tensor[1][1] = vector[1];
   tensor[1][2] = 0.0;
   tensor[2][0] = 0.0;
   tensor[2][1] = 0.0;
   tensor[2][2] = vector[3];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_4to9 */

/*----------------------------------------------------------------------*
 | reduce tensor to vector                                 he    04/03   |
 *----------------------------------------------------------------------*/
void w1_9to4(DOUBLE   tensor[3][3],
             DOUBLE  *vector
             )
{
#ifdef DEBUG
dstrc_enter("w1_9to4");
#endif
/*----------------------------------------------------------------------*/
   vector[0] = tensor[0][0];
   vector[1] = tensor[1][1];
   vector[2] = tensor[0][1];
   vector[3] = tensor[2][2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_9to4 */

/*----------------------------------------------------------------------*
 | reduce 4-stufigen tensor to 2-stufigen tensor           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_81to16(DOUBLE   tensor4[3][3][3][3],
               DOUBLE **tensor2
              )
{
#ifdef DEBUG
dstrc_enter("w1_81to16");
#endif
/*----------------------------------------------------------------------*/

tensor2[0][0] = tensor4[0][0][0][0];
tensor2[0][1] = tensor4[0][0][1][1];
tensor2[0][2] = tensor4[0][0][0][1];
tensor2[0][3] = tensor4[0][0][2][2];

tensor2[1][0] = tensor4[1][1][0][0];
tensor2[1][1] = tensor4[1][1][1][1];
tensor2[1][2] = tensor4[1][1][0][1];
tensor2[1][3] = tensor4[1][1][2][2];

tensor2[2][0] = tensor4[0][1][0][0];
tensor2[2][1] = tensor4[0][1][1][1];
tensor2[2][2] = tensor4[0][1][0][1];
tensor2[2][3] = tensor4[0][1][2][2];

tensor2[3][0] = tensor4[2][2][0][0];
tensor2[3][1] = tensor4[2][2][1][1];
tensor2[3][2] = tensor4[2][2][0][1];
tensor2[3][3] = tensor4[2][2][2][2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_81to16 */

/*----------------------------------------------------------------------*
 | reduce 4-stufigen tensor to 2-stufigen tensor           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_81to16_1(DOUBLE   tensor4[3][3][3][3],
                 DOUBLE   tensor2[4][4]
              )
{
#ifdef DEBUG
dstrc_enter("w1_81to16_1");
#endif
/*----------------------------------------------------------------------*/

tensor2[0][0] = tensor4[0][0][0][0];
tensor2[0][1] = tensor4[0][0][1][1];
tensor2[0][2] = tensor4[0][0][0][1];
tensor2[0][3] = tensor4[0][0][2][2];

tensor2[1][0] = tensor4[1][1][0][0];
tensor2[1][1] = tensor4[1][1][1][1];
tensor2[1][2] = tensor4[1][1][0][1];
tensor2[1][3] = tensor4[1][1][2][2];

tensor2[2][0] = tensor4[0][1][0][0];
tensor2[2][1] = tensor4[0][1][1][1];
tensor2[2][2] = tensor4[0][1][0][1];
tensor2[2][3] = tensor4[0][1][2][2];

tensor2[3][0] = tensor4[2][2][0][0];
tensor2[3][1] = tensor4[2][2][1][1];
tensor2[3][2] = tensor4[2][2][0][1];
tensor2[3][3] = tensor4[2][2][2][2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_81to16_1 */

/*----------------------------------------------------------------------*
 | BESTIMMUNG DER SCHAEDIGUNG UND IHRER PARTIELLEN        he    04/03   |
 | ABLEITUNG NACH KAPPA JE NACH DAMTYPE                                 |
 | DAMTYPE=1 - LINEAR-ENTFESTIGENDES DAMAGING                           |
 | DAMTYPE=2 - EXPONENTIELLES DAMAGING                                  |
 *----------------------------------------------------------------------*/
void w1_dam_typ(DOUBLE *damage,
                DOUBLE *dam_deriv,
                DOUBLE kappa,
                DOUBLE Kappa_0,
                DOUBLE Kappa_m,
                DOUBLE Alpha,
                DOUBLE Beta,
                INT    Damtyp
                )
{
#ifdef DEBUG
dstrc_enter("w1_dam_typ");
#endif
/*---------------- DAMTYPE=1 - LINEAR-ENTFESTIGENDES DAMAGING -----------*/
if (Damtyp == 1)
{
 if (kappa < Kappa_m)
 {
  *damage    = Kappa_m/kappa*(kappa-Kappa_0)/(Kappa_m-Kappa_0);
  *dam_deriv = Kappa_m*Kappa_0/(kappa*kappa*(Kappa_m-Kappa_0));
 }
 if (kappa >= Kappa_m)
 {
  *damage    = 0.999;
  *dam_deriv = 0.001;
 }
} /* end if Damtyp */

/*---------------- DAMTYPE=2 -  EXPONENTIELLES DAMAGING -----------------*/
if (Damtyp == 2)
{
  *damage    = 1.0 - Kappa_0/kappa*(1.0-Alpha+Alpha*exp(Beta*(Kappa_0-kappa)));
  *dam_deriv = Kappa_0*(1.0-Alpha)/(kappa*kappa)+Kappa_0*Alpha/kappa*
               exp(Beta*(Kappa_0-kappa))*(1.0/kappa+Beta);
} /* end if Damtyp */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_dam_typ */

/*----------------------------------------------------------------------*
 | BESTIMMUNG DES SEKANTENTENSORS                         he    04/03   |
 *----------------------------------------------------------------------*/
void w1_sec(DOUBLE damage,
            DOUBLE c_el[3][3][3][3],
            DOUBLE c_sec[3][3][3][3]
            )
{
INT    i,j,k,l;
#ifdef DEBUG
dstrc_enter("w1_sec");
#endif
/*----------------------------------------------------------------------*/

for(i=0; i<3; i++)
{
 for(j=0; j<3; j++)
 {
  for(k=0; k<3; k++)
  {
   for(l=0; l<3; l++)
   {
    c_sec[i][j][k][l] = (1.0-damage)*c_el[i][j][k][l];
   }
  }
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_sec */


/*----------------------------------------------------------------------*
 | BESTIMMUNG DES SPANNUNGSTENSORS                        he    04/03   |
 *----------------------------------------------------------------------*/
void w1_stress(DOUBLE c_sec[3][3][3][3],
               DOUBLE epsilon[3][3],
               DOUBLE sigma[3][3]
              )
{
INT    i,j,k,l;
#ifdef DEBUG
dstrc_enter("w1_stress");
#endif
/*----------------------------------------------------------------------*/

for(i=0; i<3; i++)
{
 for(j=0; j<3; j++)
 {
  sigma[i][j] = 0.0;
  for(k=0; k<3; k++)
  {
   for(l=0; l<3; l++)
   {
    sigma[i][j] = sigma[i][j] + c_sec[i][j][k][l]*epsilon[k][l];
   }
  }
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_stress */

/*----------------------------------------------------------------------*
 | BESTIMMUNG DER SPANNUNGSKOMPONENTEN sig                he    04/03   |
 *----------------------------------------------------------------------*/
void w1_cond(DOUBLE sig[4],
             DOUBLE **d
              )
{
INT    i,j;
DOUBLE sig_con[4];
DOUBLE d_con[4][4];
#ifdef DEBUG
dstrc_enter("w1_cond");
#endif
/*----------------------------------------------------------------------*/

if (d[3][3] != 0.0)
{
 for(i=0; i<3; i++)
 {
  sig_con[i]= sig[i]-sig[3]*d[i][3]/d[3][3];
 }

 for(i=0; i<3; i++)
 {
  sig[i] = sig_con[i];
 }

 for(i=0; i<3; i++)
 {
  for (j=0; j<3; j++)
  {
   d_con[i][j]=d[i][j]-d[i][3]*d[3][j]/d[3][3];
  }
 }

 for(i=0; i<3; i++)
 {
  for (j=0; j<3; j++)
  {
   d[i][j]=d_con[i][j];
  }
 }

 for(i=0; i<4; i++)
 {
  d[i][3] = 0.0;
  d[3][i] = 0.0;
 }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_cond */
#endif
