/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_tvma: which integrates the material law and stresses in thickness 
            direction 

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief integrate material law and stresses in thickness direction of shell                                      

<pre>                     m.gee 6/01              modified by    sh 11/02
This routine integrates the material law and stresses in thickness 
direction of multilayer shell element.
</pre>
\param  double **D       (o)  konstitutive matrix combining kinematic and static variables 
                              of shell element (-> dimension [12][12])
\param  double **C       (i)  konstitutive matrix from material law (-> dimension [6][6])
\param  double  *stress  (i)  PK_II stresses from konsitutive law (-> dimenesion [6])
\param  double  *stress_r(o)  stress resultants ("Schnittgroessen") (-> dimension [12])
\param  double   e3      (i)  zeta of aktuel material layer
\param  double   fact    (i)  factor including the shell shifter
\param  double   h       (i)  total thickness of this element
\param  double  *klayhgt (i)  hight of kin layer in % of total thickness of shell 
\param  double  *mlayhgt (i)  hight of mat layer in % of adjacent kin layer 
\param  int      num_klay(i)  number of kin layers to this element  
\param  int      num_mlay(i)  number of mat layers to this kin layer 
\param  int      klay    (i)  actual kin layer 
\param  int      mlay    (i)  actual mat layer of this kin layer 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_tvma(double   **D,
             double   **C,
             double    *stress,
             double    *stress_r,
             double     e3,
             double     fact,
             double     h,           /* total thickness of this element */
             double    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             double    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             int        num_klay,    /* number of kin layers to this element */  
             int        num_mlay,    /* number of mat layers to this kin layer */
             int        klay,        /* actual kin layer */
             int        mlay,        /* actual mat layer of this kin layer */
             double     condfac)
{
int i,i6,j,j6;
double deltah, h_mlay, h_kl;
double zeta_kl,zeta;
double stress_fact, C_fact;
#ifdef DEBUG 
dstrc_enter("s9_tvma");
#endif
/*----------------------------------------------------------------------*/
/*- calculate zeta_kl of kinematic layer due to local zeta_ml of material layer -> old s9tmtr --*/
h_kl   = (klayhgt[klay]/100.)*h;              /* absolute hight of the actual kinematic layer */
deltah = (mlayhgt[mlay]/100.)*h_kl;           /* absolute hight of the actual material layer */
h_mlay   = 0.0;
for (i=0; i<=mlay; i++)                        /* sum of the absolute hights of the material layers  */
{                                             /* within the actual kinematic layer up to the actual */
   h_mlay += (mlayhgt[i]/100.)*h_kl;          /* material layer */
}
zeta_kl = -1. + (-deltah*(1.-e3)+2.*h_mlay)/h_kl;  /* equals the theta3 value of the act. kinematic layer, see equation (5.45) in Dis. Braun */
fact = fact*deltah/h_kl ;
/*---- get continuity coefficient for present kinematic layer -------------*/
zeta = zeta_kl;
zeta = s9con(zeta,num_klay,klay,klay,condfac);


for (i=0; i<6; i++)
{
   i6=i+6;
   stress_fact   = stress[i]*fact;
   stress_r[i]  += stress_fact;
   stress_r[i6] +=  stress_fact * zeta;
   for (j=0; j<6; j++)
   {
      j6 = j+6;
      C_fact     = C[i][j]*fact;
      D[i][j]   += C_fact;
      D[i6][j]  += C_fact*zeta;
      D[j][i6]  += C_fact*zeta;
      D[i6][j6] += C_fact*zeta*zeta;
   }
}
/*-------------------------------------------------------- symmetrize D */
/*for (i=0; i<12; i++)
{
   for (j=i+1; j<12; j++)
   {
      D[i][j]=D[j][i];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tvma */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
 
