/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_mat_plast_mises', 'b3_blowup_stress',
'b3_compr_stress', 'b3_blowup_strain', 'b3_compr_strain', 'b3_compr_Cep',
'b3_kronecker', 'b3_cal_volCel', 'b3_cal_devCel', 'b3_cal_Cel', 'b3_condense'

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the constitutive matrix - forces for von Mises material

<pre>                                                              fh 04/03
This routine calculates the constitutive matrix and the actual stresses
within a plastic calculation for von Mises material

</pre>
\param ym              DOUBLE      (i)  Youngs Modulus
\param pv              DOUBLE      (i)  Possions ratio
\param ALFAT           DOUBLE      (i)  alpha T
\param sigy            DOUBLE      (i)  1D yield stress
\param hard            DOUBLE      (i)  1D hardening modulus
\param gf              DOUBLE      (i)  Crack energy (Gf)
\param *ele            ELEMENT    (i/o) actual element
\param *strain         DOUBLE      (i)  strain vector of actual IP
\param ip              INT         (i)  actual integration point number
\param *stress         DOUBLE      (o)  stress vector of actual IP
\param **d             DOUBLE      (o)  constitutive matrix
\param istore          INT         (i)  flag to control storing of new val to WA
\param newval          INT         (i)  flog to control evaluation of new stresses
\param init            INT         (i)  initialization (1) or calculation (2,3)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_blowup_stress() , b3_compr_stress() , b3_blowup_strain() ,
               b3_compr_strain() , b3_compr_Cep() , b3_kronecker() , 
	       b3_cal_volCel() , b3_cal_devCel() , b3_cal_Cel() , b3_condense()
    called by: beam3() , b3_call_mat()

*----------------------------------------------------------------------*/
void b3_mat_plast_mises(DOUBLE    ym,
                        DOUBLE    pv,
                        DOUBLE    ALFAT,
                        DOUBLE    sigy,
                        DOUBLE    hard,
                        DOUBLE    gf,
                        ELEMENT  *ele,
                        DOUBLE   *strain,
			INT       ip,
                        DOUBLE   *stress,       
                        DOUBLE  **d,
                        INT       istore,
                        INT       newval,
                        INT       init)
{
INT i,j,k,k1,k2,l;
INT iel;

DOUBLE yip;
DOUBLE hards;                  /* hardening module for 3-dim stress state */
INT    iupd=0;
DOUBLE bulk,shear;             /* bulk modules and shear modulus */
DOUBLE alpha_pl_old;           /* internal variable for isotropic hardening */
DOUBLE alpha_pl;               /* internal variable for isotropic hardening */
DOUBLE eps_vol;                /* volumetric strain (scalar) */
DOUBLE sig_vol;                /* volumetric stress (scalar) */
DOUBLE e_dev[3][3];            /* deviatoric strain tensor */
DOUBLE sigma_el_vol[3][3];     /* volumetric part of trial stress tensor */
DOUBLE sigma_tr_dev[3][3];     /* deviatoric part of trial stress tensor */
DOUBLE sigma_tr[3][3];         /* trial stress tensor */
DOUBLE sigma_tr_dev2;          /* sigma_tr_dev * sigma_tr_dev */
DOUBLE phi_tr;                 /* flow condition */ 
DOUBLE n_tr[3][3];             /* dphi / dsigma */
DOUBLE delta_C[3][3][3][3];    /* Cep-Cel */
DOUBLE I_dev[3][3][3][3];      /* IIdev */
DOUBLE gamma_pl;               /* plastic multiplier */
DOUBLE c_pl_1, c_pl_2;
DOUBLE dikjl,diljk,dijkl;
DOUBLE eps_it[6];              /* used for iteration to sigyy=sigzz=tauyz=0 */
DOUBLE deleps_it[3];           /* used for iteration to sigyy=sigzz=tauyz=0 */
DOUBLE sig_dum_it[3];          /* used for iteration to sigyy=sigzz=tauyz=0 */
DOUBLE det[1];                 /* used for iteration to sigyy=sigzz=tauyz=0 */

DOUBLE tol_pl = 1.0E-10;

const DOUBLE    q23 =.66666666666666666667;
const DOUBLE    q13 =.33333333333333333333;
const INT    numeps = 3;         /* epsilon xx,gamma xy,gamma xz    */
const INT    max = MAXDOFPERNODE*MAXNOD_BEAM3;
const INT    numdf = 6;
INT          nedof;

static ARRAY      sig_a; 
static DOUBLE    *sig; /* new stress vector at GP */
static ARRAY      sig_it_a;
static DOUBLE    *sig_it; /* used for iteration to sigyy=sigzz=tauyz=0 */
static ARRAY      sigma_a;
static DOUBLE   **sigma; /* new stress tensor at GP */
static ARRAY      sigma_old_a;
static DOUBLE   **sigma_old; /* old stress tensor at GP */    
static ARRAY      eps_pl_a;
static DOUBLE    *eps_pl; /* new plastic strain vector at GP */
static ARRAY      epsilon_a; 
static DOUBLE   **epsilon; /* new strain tensor at GP */ 
static ARRAY      epsilon_pl_a; 
static DOUBLE   **epsilon_pl; /* new plastic strain tensor at GP */ 
static ARRAY      epsilon_pl_old_a;
static DOUBLE   **epsilon_pl_old; /* old plastic strain tensor at GP */
static ARRAY      delta_a;
static DOUBLE   **delta; /* Kronecker-delta */
static ARRAY      C_BB_a;
static DOUBLE   **C_BB; /* used for iteration to sigyy=sigzz=tauyz=0 */
static ARRAY      C_BA_a;
static DOUBLE   **C_BA; /* used for iteration to sigyy=sigzz=tauyz=0 */
static ARRAY4D    Cep_a;
static DOUBLE ****Cep; /* elastoplastic tangent */
static ARRAY4D    Cvol_a;
static DOUBLE ****Cvol; /* volumetric part of elasticity tensor */
static ARRAY4D    Cdev_a;
static DOUBLE ****Cdev; /* deviatoric part of elasticity tensor */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("b3_mat_plast_mises");
#endif

/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/
if (init==1)
{
  sig            = amdef("sig"            ,&sig_a      ,  numdf,1 ,"DV");
  sig_it         = amdef("sig_it"         ,&sig_it_a   ,  numdf,1 ,"DV");
  eps_pl         = amdef("eps_pl"         ,&eps_pl_a   ,  numdf,1 ,"DV");
  sigma          = amdef("sigma"          ,&sigma_a         , 3,3 ,"DA");
  sigma_old      = amdef("sigma_old"      ,&sigma_old_a     , 3,3 ,"DA");
  epsilon        = amdef("epsilon"        ,&epsilon_a       , 3,3 ,"DA");
  epsilon_pl     = amdef("epsilon_pl"     ,&epsilon_pl_a    , 3,3 ,"DA");
  epsilon_pl_old = amdef("epsilon_pl_old" ,&epsilon_pl_old_a, 3,3 ,"DA");
  delta          = amdef("delta"          ,&delta_a,          3,3 ,"DA");
  C_BB           = amdef("C_BB"           ,&C_BB_a,           3,3 ,"DA");
  C_BA           = amdef("C_BA"           ,&C_BA_a,           3,3 ,"DA");
  Cep            = am4def("Cep"            ,&Cep_a  ,      3,3,3,3 ,"D4");
  Cvol           = am4def("Cvol"           ,&Cvol_a ,      3,3,3,3 ,"D4");
  Cdev           = am4def("Cdev"           ,&Cdev_a ,      3,3,3,3 ,"D4");  
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
amdel(&sig_a);
amdel(&sig_it_a);
amdel(&eps_pl_a);
amdel(&sigma_a);
amdel(&sigma_old_a);
amdel(&epsilon_a);
amdel(&epsilon_pl_a);
amdel(&epsilon_pl_old_a);
amdel(&delta_a);
amdel(&C_BB_a);
amdel(&C_BA_a);
am4del(&Cep_a);
am4del(&Cvol_a);
am4del(&Cdev_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
amzero(&sig_a);
amzero(&sig_it_a);
amzero(&eps_pl_a);
amzero(&sigma_a);
amzero(&sigma_old_a);
amzero(&epsilon_a);
amzero(&epsilon_pl_a);
amzero(&epsilon_pl_old_a);
amzero(&delta_a);
amzero(&C_BB_a);
amzero(&C_BA_a);
am4zero(&Cep_a);
am4zero(&Cvol_a);
am4zero(&Cdev_a);

bulk  = ym / (3.-6.*pv);
shear = ym / (2.+2.*pv);
iel = ele->numnp;
nedof = iel * numdf;
hards = ym * hard / (ym - hard);

/*----------------------------------------------------------------------*/
/* element working array:	0     = alpha_pl                        */
/* 				1     = yip		                */
/* 				2-7   = sigma		                */
/* 				8-13  = eps_pl		                */
/* 				14-19 = eps_it                          */
/* 				20-22 = sig_it                          */
/* 				23-40 = d[3][0]-d[5][5]_it              */
/*----------------------------------------------------------------------*/


/*----------------------------- get old values -> sig, eps , alpha------*/
alpha_pl_old = ele->e.b3->elewa.a.da[ip][0];
yip          = ele->e.b3->elewa.a.da[ip][1];
for (i=0; i<6; i++) 
{
   sig[i]       = ele->e.b3->elewa.a.da[ip][i+2];
   eps_pl[i]    = ele->e.b3->elewa.a.da[ip][i+8];
}


/*---------- get stresses from WA for calculation of M,V,N -------------*/
if (newval==1) 
{
   for (i=0; i<6; i++) stress[i]=sig[i];
   goto end;
}


/* get old values for condensation of stresses sig_yy=sig_zz=sig_yz=0   */
/*----------------------------------------------------------------------*/
/* 	eps_it = [eps_xx, gamma_xy, gamma_xz, eps_yy, eps_zz, gamma_yz]T*/
/* 	sig_it = [sig_yy, sig_zz, tau_yz]T of last iteration step       */
/*   deleps_it = actual [eps_xx, gamma_xy, gamma_xz]T - eps_it          */ 
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
{ 
   sig_it[i]    = ele->e.b3->elewa.a.da[ip][i+20];
}

for (i=0; i<6; i++)
{ 
   eps_it[i]    = ele->e.b3->elewa.a.da[ip][i+14];
   deleps_it[i] = strain[i]-eps_it[i];
}


/*-----------calculate strains eps_yy, eps_zz, gamma_yz-----------------*/
/*----------------------------------------------------------------------*/
/*           deleps_r=-1/C_BB*[C_BA*deleps_it+sig_it]                   */
/*----------------------------------------------------------------------*/

/*-----------get C_BA [20-28] and C_BB [29-37] from working array-------*/
k1=23;
k2=32;
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      C_BA[i][j]=ele->e.b3->elewa.a.da[ip][k1];
      C_BB[i][j]=ele->e.b3->elewa.a.da[ip][k2];
      k1+=1;
      k2+=1;
   }
} 
/*-----------make C_BB to 1/C_BB----------------------------------------*/
math_inv3(C_BB, det);
if ( *det==0.0 ) amzero(&C_BB_a);
/*-----------calculate C_BA*deleps_it-----------------------------------*/
math_matvecdense(sig_dum_it,C_BA,deleps_it,3,3,0,1.);
/*-----------calculate C_BA*deleps_it+sig_it----------------------------*/
for (i=0; i<3; i++) sig_dum_it[i]+=sig_it[i];
/*-----------calculate -1/C_BB*[C_BA*deleps_it+sig_it]------------------*/
math_matvecdense(deleps_it,C_BB,sig_dum_it,3,3,0,-1.0);
/*-----------calculate strain [yy,zz,yz] = eps_it+deleps_it-------------*/
for (i=0; i<3; i++) 
{
   if (deleps_it[i]!=0.) strain[i+3]=eps_it[i+3]+deleps_it[i];
}

  
b3_blowup_strain(strain,epsilon);
b3_blowup_strain(eps_pl,epsilon_pl_old);
b3_kronecker(delta);


/*----------------------------------------------------------------------|
|     yip > 0	  Stresses are available from last update               |
|         = 1.0   elastic step                                          |
|         = 2.0   plastic step                                          |
|----------------------------------------------------------------------*/
if (yip>0)
{
   if (yip==1.0)
   {    
      b3_cal_Cel(ym,pv,delta,Cep);
      b3_compr_Cep(d,Cep);
      math_matvecdense(sig,d,strain,6,6,0,1.);
      yip=-yip;
      iupd=1;
      goto stop;
   }
   else
   {
      b3_blowup_stress(sig,sigma_old);
      sig_vol=0.;							         
      for (i=0; i<3; i++)						        
      {
         for (j=0; j<3; j++) sig_vol+=sigma_old[i][j]*delta[i][j];
      }
      for (i=0; i<3; i++)
      {
         for (j=0; j<3; j++) 
	 {
	    sigma_tr_dev[i][j]=sigma_old[i][j]-q13*sig_vol*delta[i][j];
	    sigma_tr[i][j]=sigma_old[i][j];
	 }
      } 
      sigma_tr_dev2=0.;
      for (i=0; i<3; i++)
      {
         for (j=0; j<3; j++) sigma_tr_dev2+=sigma_tr_dev[i][j]*sigma_tr_dev[i][j];
      }
      phi_tr=0.;
      gamma_pl=0.;
      yip=-yip;
      iupd=1;
      goto plastic;
   }
}


/*----------calculate volumetric and deviatoric strains ----------------*/
eps_vol=0.;
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++) eps_vol+=epsilon[i][j]*delta[i][j];
}
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++) e_dev[i][j]=epsilon[i][j]-q13*eps_vol*delta[i][j];
} 

/*----------calculate vol. and dev. part of elasticity tensor ----------*/
/*b3_cal_volCel(bulk,delta,Cvol);
/*b3_cal_devCel(shear,delta,Cdev);



/*----------------------------------------------------------------------|
|           STEP 1: Trial state: predictor values                       |
|-----------------------------------------------------------------------*/

/*----------calculate trial stresses-------------------------------------*/
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {     
     sigma_el_vol[i][j]=bulk*delta[i][j]*eps_vol;
     sigma_tr_dev[i][j]=2.*shear*(e_dev[i][j]-epsilon_pl_old[i][j]);
     sigma_tr[i][j]=sigma_el_vol[i][j]+sigma_tr_dev[i][j];
   }
}

/*----------calculate trial flow condition------------------------------*/
sigma_tr_dev2=0.;
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++) sigma_tr_dev2+=sigma_tr_dev[i][j]*sigma_tr_dev[i][j];
}
phi_tr=sqrt(sigma_tr_dev2)-sqrt(2./3.)*(sigy+hards*alpha_pl_old);


/*----------------------------------------------------------------------|
|           STEP 2: elastic                                             |
|-----------------------------------------------------------------------*/
if (phi_tr<=tol_pl)
{
   yip=1;
   alpha_pl=alpha_pl_old;
   for (i=0; i<3; i++)
   {
      for (j=0; j<3; j++) 
      {
         epsilon_pl[i][j]=epsilon_pl_old[i][j];
	 sigma[i][j]=sigma_tr[i][j];
      }
   }
   b3_cal_Cel(ym,pv,delta,Cep);
}


/*----------------------------------------------------------------------|
|           STEP 3: plastic                                             |
|-----------------------------------------------------------------------*/
else
{
/*----------plastic multiplier gamma_pl---------------------------------*/
   yip=2;
   gamma_pl=phi_tr/(2.*shear+q23*hards);

plastic:
/*----------calculate stresses------------------------------------------*/   
   for (i=0; i<3; i++)
   {
      for (j=0; j<3; j++) 
      {
         n_tr[i][j]=sigma_tr_dev[i][j]/sqrt(sigma_tr_dev2);
	 sigma[i][j]=sigma_tr[i][j]-2.*shear*gamma_pl*n_tr[i][j];
      }
   }
/*----------calculate the elastoplastic tangent Cep---------------------*/
   b3_cal_Cel(ym,pv,delta,Cep);    	  	  
   c_pl_1=-2.*shear*gamma_pl/sqrt(sigma_tr_dev2);
   c_pl_2=2.*shear/(2.*shear+q23*hards)*(1-phi_tr/sqrt(sigma_tr_dev2));
   for (i=0; i<3; i++)
   {
      for (j=0; j<3; j++)
      {
         for (k=0; k<3; k++)
	 {
	    for (l=0; l<3; l++)
	    {
               dikjl=delta[i][k]*delta[j][l];
	       diljk=delta[i][l]*delta[j][k];
	       dijkl=delta[i][j]*delta[k][l];
	       I_dev[i][j][k][l]=(dikjl+diljk)/2.-q13*dijkl;
	       delta_C[i][j][k][l]=c_pl_1*2.*shear*I_dev[i][j][k][l]-c_pl_2*2.*shear*n_tr[i][j]*n_tr[k][l];
/*-----------calculate Cep----------------------------------------------*/
	       Cep[i][j][k][l]+=delta_C[i][j][k][l];
	    }
         }
      }
   }
/*-----------calculate new alpha_pl, epsilon_pl-------------------------*/
   alpha_pl=alpha_pl_old+sqrt(2./3.)*gamma_pl;
   for (i=0; i<3; i++)
   {
      for (j=0; j<3; j++) epsilon_pl[i][j]=epsilon_pl_old[i][j]+gamma_pl*n_tr[i][j];
   }
}


b3_compr_strain(eps_pl,epsilon_pl);
b3_compr_Cep(d,Cep);
b3_compr_stress(sig,sigma);


stop:
/*-----------write data for iteration sigyy=0,sigzz=0,sigyz=0 to WA ----*/
k=23;
for (i=0; i<3; i++)
{
/*-----------write sig_it to WA-----------------------------------------*/
   ele->e.b3->elewa.a.da[ip][i+20] = sig[i+3];
/*-----------write eps_it to WA-----------------------------------------*/
   ele->e.b3->elewa.a.da[ip][i+14] = strain[i];
   ele->e.b3->elewa.a.da[ip][i+17] = strain[3+i];
   for (j=0; j<3; j++) 
   {  
/*-----------write C_BA to working array [20-28]------------------------*/
      ele->e.b3->elewa.a.da[ip][k]   = d[i+3][j];
/*-----------write C_BB to working array [29-37]------------------------*/
      ele->e.b3->elewa.a.da[ip][k+9] = d[i+3][j+3];
      k+=1;
   }
}


b3_condense(d,sig);
/*-----------stresses for calculation of internal forces----------------*/
for (i=0; i<3; i++) stress[i] = sig[i];


end:
/*-----------put new values -> alpha_pl, yip, sig, eps_pl---------------*/
if (istore==1)
{   
   ele->e.b3->elewa.a.da[ip][0] = alpha_pl;
   ele->e.b3->elewa.a.da[ip][1] = yip;
   for (i=0; i<6; i++)				        	 
   {
      ele->e.b3->elewa.a.da[ip][i+2] = sig[i];
      ele->e.b3->elewa.a.da[ip][i+8] = eps_pl[i];
   }
}


/*-----------only done for new displacement step------------------------*/
if (iupd==1) ele->e.b3->elewa.a.da[ip][1] = yip;


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_mat_plast_mises */

/*!----------------------------------------------------------------------
\brief calculates stress tensor (3*3) out of stress vector (6*1)

<pre>                                                              fh 03/03
This routine calculates the stress tensor (3*3) out of stress vector (6*1)

</pre>
\param *vector   DOUBLE    (i)  stress vector
\param **matrix  DOUBLE    (o)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_blowup_stress(DOUBLE  *vector,
                      DOUBLE **matrix)
{
#ifdef DEBUG
dstrc_enter("b3_blowup_stress");
#endif

/* vector form: [sig xx, sig xy, sig xz, sig yy, sig zz, sig yz]T       */
matrix[0][0]=vector[0];
matrix[0][1]=vector[1];
matrix[0][2]=vector[2];
matrix[1][0]=vector[1];
matrix[1][1]=vector[3];
matrix[1][2]=vector[5];
matrix[2][0]=vector[2];
matrix[2][1]=vector[5];
matrix[2][2]=vector[4];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_blowup_stress */

/*!----------------------------------------------------------------------
\brief calculates stress vector (6*1) out of stress tensor (3*3)

<pre>                                                              fh 03/03
This routine calculates the stress vector (6*1) out of stress tensor (3*3)

</pre>
\param *vector   DOUBLE    (o)  stress vector
\param **matrix  DOUBLE    (i)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_stress(DOUBLE  *vector,
                     DOUBLE **matrix)
{
#ifdef DEBUG
dstrc_enter("b3_compr_stress");
#endif

/* vector form: [sig xx, sig xy, sig xz, sig yy, sig zz, sig yz]T       */
vector[0]=matrix[0][0];
vector[1]=matrix[0][1];
vector[2]=matrix[0][2];
vector[3]=matrix[1][1];
vector[4]=matrix[2][2];
vector[5]=matrix[1][2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_compr_stress */


/*!----------------------------------------------------------------------
\brief calculates strain tensor (3*3) out of strain vector (6*1)

<pre>                                                              fh 03/03
This routine calculates the strain tensor (3*3) out of strain vector (6*1)

</pre>
\param *vector   DOUBLE    (i)  strain vector
\param **matrix  DOUBLE    (o)  strain tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_blowup_strain(DOUBLE  *vector,
                      DOUBLE **matrix)
{
#ifdef DEBUG
dstrc_enter("b3_blowup_strain");
#endif

/* vector form: [eps xx, gamma xy, gamma xz, eps yy, eps zz, gamma yz]T */
matrix[0][0]=vector[0];
matrix[0][1]=vector[1]/2.;
matrix[0][2]=vector[2]/2.;
matrix[1][0]=vector[1]/2.;
matrix[1][1]=vector[3];
matrix[1][2]=vector[5]/2.;
matrix[2][0]=vector[2]/2.;
matrix[2][1]=vector[5]/2.;
matrix[2][2]=vector[4];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_blowup_strain */

/*!----------------------------------------------------------------------
\brief calculates strain vector (6*1) out of strain tensor (3*3)

<pre>                                                              fh 03/03
This routine calculates the strain vector (6*1) out of strain tensor (3*3)

</pre>
\param *vector   DOUBLE    (o)  stress vector
\param **matrix  DOUBLE    (i)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_strain(DOUBLE  *vector,
                     DOUBLE **matrix)
{
#ifdef DEBUG
dstrc_enter("b3_compr_strain");
#endif

/* vector form: [eps xx, gamma xy, gamma xz, eps yy, eps zz, gamma yz]T */
vector[0]=matrix[0][0];
vector[1]=matrix[0][1]*2.;
vector[2]=matrix[0][2]*2.;
vector[3]=matrix[1][1];
vector[4]=matrix[2][2];
vector[5]=matrix[1][2]*2.;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_compr_strain */

/*!----------------------------------------------------------------------
\brief calculates Cep matrix (6*6) out of Cep tensor (3*3*3*3)

<pre>                                                              fh 03/03
This routine calculates the elastoplastic tangent matrix (6*6) out of 
elastoplastic tangent tensor (3*3*3*3)

</pre>
\param **matrix    DOUBLE    (o)  Cep matrix
\param ****tensor  DOUBLE    (i)  Cep tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_Cep(DOUBLE   **matrix,
                  DOUBLE ****tensor)
{
#ifdef DEBUG
dstrc_enter("b3_compr_Cep");
#endif

/* stress vector : [sig xx, sig xy, sig xz, sig yy, sig zz, sig yz]T   */
/* strain vector:  [eps xx, gam xy, gam xz, eps yy, eps zz, gam yz]T   */
matrix[0][0]=tensor[0][0][0][0];
matrix[0][1]=tensor[0][0][0][1];
matrix[0][2]=tensor[0][0][0][2];
matrix[0][3]=tensor[0][0][1][1];
matrix[0][4]=tensor[0][0][2][2];
matrix[0][5]=tensor[0][0][1][2];

matrix[1][0]=tensor[0][1][0][0];
matrix[1][1]=tensor[0][1][0][1];
matrix[1][2]=tensor[0][1][0][2];
matrix[1][3]=tensor[0][1][1][1];
matrix[1][4]=tensor[0][1][2][2];
matrix[1][5]=tensor[0][1][1][2];
		  
matrix[2][0]=tensor[0][2][0][0];
matrix[2][1]=tensor[0][2][0][1];
matrix[2][2]=tensor[0][2][0][2];
matrix[2][3]=tensor[0][2][1][1];
matrix[2][4]=tensor[0][2][2][2];
matrix[2][5]=tensor[0][2][1][2];

matrix[3][0]=tensor[1][1][0][0];
matrix[3][1]=tensor[1][1][0][1];
matrix[3][2]=tensor[1][1][0][2];
matrix[3][3]=tensor[1][1][1][1];
matrix[3][4]=tensor[1][1][2][2];
matrix[3][5]=tensor[1][1][1][2];

matrix[4][0]=tensor[2][2][0][0];
matrix[4][1]=tensor[2][2][0][1];
matrix[4][2]=tensor[2][2][0][2];
matrix[4][3]=tensor[2][2][1][1];
matrix[4][4]=tensor[2][2][2][2];
matrix[4][5]=tensor[2][2][1][2];

matrix[5][0]=tensor[1][2][0][0];
matrix[5][1]=tensor[1][2][0][1];
matrix[5][2]=tensor[1][2][0][2];
matrix[5][3]=tensor[1][2][1][1];
matrix[5][4]=tensor[1][2][2][2];
matrix[5][5]=tensor[1][2][1][2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_compr_Cep */

/*!----------------------------------------------------------------------
\brief calculates Kronecker-delta for (3*3) tensor

<pre>                                                              fh 03/03
This routine calculates the Kronecker-delta for (3*3) tensor

</pre>
\param **delta    DOUBLE    (i/o)  Kronecker-delta


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_kronecker(DOUBLE **delta)
{

INT i,j; /* some loopers */

#ifdef DEBUG
dstrc_enter("b3_kronecker");
#endif

for (i=0; i<3; i++)
{
  for (j=0; j<3; j++)
  {
     if (i==j) delta[i][j]=1.;
     else delta[i][j]=0.;
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_kronecker */

/*!----------------------------------------------------------------------
\brief calculates volumetric part of elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the volumetric part of the elasticity tensor

</pre>
\param bulk       DOUBLE     (i)   bulk modulus
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_vol  DOUBLE     (o)   volumetric part of elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_volCel(DOUBLE     bulk,
                   DOUBLE   **delta,
	           DOUBLE ****c_vol)
{

INT i,j,k,l; /* some loopers */

#ifdef DEBUG
dstrc_enter("b3_cal_volCel");
#endif

for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      for (k=0; k<3; k++)
      {
         for (l=0; l<3; l++) c_vol[i][j][k][l]=bulk*delta[i][j]*delta[k][l];
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_volCel */

/*!----------------------------------------------------------------------
\brief calculates volumetric part of elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the volumetric part of the elasticity tensor

</pre>
\param bulk       DOUBLE     (i)   bulk modulus
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_dev  DOUBLE     (o)   volumetric part of elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_devCel(DOUBLE     shear,
                   DOUBLE   **delta,
	           DOUBLE ****c_dev)
{

INT i,j,k,l; /* some loopers */
DOUBLE dikjl,diljk,dijkl;
const DOUBLE q23=.666666666666666667;

#ifdef DEBUG
dstrc_enter("b3_cal_devCel");
#endif

for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      for (k=0; k<3; k++)
      {
         for (l=0; l<3; l++) 
	 { 
	    dikjl=delta[i][k]*delta[j][l];
	    diljk=delta[i][l]*delta[j][k];
	    dijkl=delta[i][j]*delta[k][l];
	    c_dev[i][j][k][l]=shear*(dikjl+diljk)-q23*shear*dijkl;
	 }
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_devCel */

/*!----------------------------------------------------------------------
\brief calculates elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the elasticity tensor

</pre>
\param ym         DOUBLE     (i)   youngs modulus
\param pv         DOUBLE     (i)   possions value
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_el   DOUBLE     (o)   elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_Cel(DOUBLE     ym,
                DOUBLE     pv,
                DOUBLE   **delta,
                DOUBLE ****c_el)
{

INT i,j,k,l; /* some loopers */
DOUBLE lambda,mu; /* lame constants */
DOUBLE dikjl,diljk,dijkl;
const DOUBLE q23=.666666666666666667;

#ifdef DEBUG
dstrc_enter("b3_cal_Cel");
#endif

lambda= ym*pv/((1.+pv)*(1.-2.*pv));
mu    = ym/(2.+2.*pv);

for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      for (k=0; k<3; k++)
      {
         for (l=0; l<3; l++) 
	 { 
	    dikjl=delta[i][k]*delta[j][l];
	    diljk=delta[i][l]*delta[j][k];
	    dijkl=delta[i][j]*delta[k][l];
	    c_el[i][j][k][l]=lambda*dijkl+mu*(dikjl+diljk);
	 }
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_cal_Cel */

/*!----------------------------------------------------------------------
\brief condensates stresses and elastoplastic tangent

<pre>                                                              fh 03/03
This routine condensates the stresses and the elastoplastic tangent

</pre>
\param **D        DOUBLE    (i/o)   elastoplastic tangent
\param *sig       DOUBLE    (i/o)   stress vector


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_condense(DOUBLE **D,
                 DOUBLE  *sig)
{

INT j,n,m; /* some loopers */
DOUBLE stfc, constress; 

#ifdef DEBUG
dstrc_enter("b3_condense");
#endif

/* sigma yy = 0 */
/* sigma zz = 0 */
/* sigma yz = 0 */

/* Adopted from CARAT */
/*----------------------------------------------------------------------*/
for (j=3; j<6; j++)
{
  stfc=D[j][j];
  if (D[j][j]!=0.)
  {	
    for (n=0; n<6; n++) D[j][n]=D[j][n]/stfc;
    /*-------------Condensation stress--------------------------------------*/
    constress=sig[j];
    for (n=0; n<6; n++) sig[n]=sig[n]-D[j][n]*constress;    
    /*-------------Condensation Kaa-----------------------------------------*/
    for (n=0; n<j; n++)
    {
      for (m=n; m<j; m++) D[n][m]=D[n][m]-D[n][j]*D[j][m];
      for (m=0; m<n; m++) D[n][m]=D[m][n];
    }	  
    /*-------------Condensation Kac and Kca---------------------------------*/
    for (n=j+1; n<6; n++)
    {
      for (m=0; m<j; m++)
      {
    	D[n][m]=D[n][m]-D[j][n]*D[m][j];
    	D[m][n]=D[n][m];
      }
    }					
    /*-------------Condensation Kcc-----------------------------------------*/
    for (n=j+1; n<6; n++)
    {
      for (m=n; m<6; m++) D[n][m]=D[n][m]-D[n][j]*D[j][m];
      for (m=j+1; m<6; m++) D[n][m]=D[m][n];
    }			      
    /*-------------Set Cond. domain equal zero------------------------------*/
    for (n=0; n<6; n++)
    {
      D[j][n]=0.;
      D[n][j]=0.;
    }
  }
  else goto out;                                   
}
goto end;
out:
for (j=0; j<6; j++)
{
  for (n=0; n<6; n++) D[j][n]=0.;
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_condense */
#endif
/*! @} (documentation module close)*/
