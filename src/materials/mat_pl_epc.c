/*!----------------------------------------------------------------------
\file
\brief contains the routine 
- mat_pl_epc_main:      which calculates the constitutive matrix for 
                        the elastoplastic concrete model described
                        in Dis. Menrath and Dis. Haufe
                        This routine is formulated in cartesian 
                        coordinate system, general 3D with the sorting 
                        [11,22,33,12,23,13]
- mat_pl_epc_radi_dp:   which performs the stress projection for a
                        standard drucker-prager cone
- mat_pl_epc_mapl_dp:   which calculates the algroithmic material tangent
                        for a standard drucker-prager cone
- mat_pl_epc_radi_ms:   performs stress projection for two coupled DP-
                        yield surfaces (including the apex -> inverted cone)
- mat_pl_epc_mapl_ms:   calculates the algorithmic material tangent for
                        two coupled DP-yield surfaces
- mat_pl_epc_radi_cap:  stress projection for spherical cap
- mat_pl_epc_mapl_cap:  algorithmic material tangent for spherical cap
- mat_pl_epc_dia:       calculates the critical internal length
- mat_pl_epc_hard:      calculates needed derivatives for the hardening law
                        in compression and tension
- mat_pl_epc_preval:    claculates needed prevalues
*----------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

/*! 
\addtogroup MAT 
*//*! @{ (documentation module open)*/ 


/*!----------------------------------------------------------------------
\brief consitutive matrix for the elastoplastic concrete model            

<pre>                                                            sh 10/03
This routine calculates the constitutive matrix and forces for the elasto
plastic concrete material model described in Dis. Menrath (see also Dis. Haufe)
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE     Ec        (i)  young's modulus of concrete             
\param  DOUBLE     vc        (i)  poisson's ratio of concrete             
\param  DOUBLE     ftm       (i)  tensile strength of concrete              
\param  DOUBLE     fcm       (i)  compressive strength of concrete              
\param  DOUBLE     gt        (i)  tensile fracture energie of concrete              
\param  DOUBLE     gc        (i)  compressive fracture energie of concrete              
\param  DOUBLE     gamma1    (i)  fitting parameter              
\param  DOUBLE     gamma2    (i)  fitting parameter              
\param  DOUBLE     gamma3    (i)  fitting parameter              
\param  DOUBLE     gamma4    (i)  fitting parameter              
\param  DOUBLE     dia       (i)  equivalent element length              
\param  DOUBLE     stress[6] (o)  vector of stresses [11,22,33,12,23,13]
\param  DOUBLE     strain[6] (i)  actual strains from displacements  [11,22,33,12,23,13]
\param  DOUBLE   **d         (o)  constitutive matrix          
\param  INT       *iupd     (i/o) controls update of new stresses to wa         
\param  INT       *yip      (i/o) flag if global predictor step an if last step was plastic/elastic               
\param  DOUBLE    *kappa_t  (i/o) uniaxial equivalent strain -> WA (tension)
\param  DOUBLE    *kappa_c  (i/o) uniaxial equivalent strain -> WA (compression)
\param  DOUBLE     sig[6]    (i)  stresses from WA  [11,22,33,12,23,13]
\param  DOUBLE     eps[6]    (i)  strains from WA  [11,22,33,12,23,13]

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_mat_plast_epc()     [s9_mat_plast_epc.c]
                             s9_mat_plast_rc()      [s9_mat_plast_rc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_main(DOUBLE   Ec,        
                     DOUBLE   vc,        
                     DOUBLE   ftm,       
                     DOUBLE   fcm,       
                     DOUBLE   gt,        
                     DOUBLE   gc, 
                     DOUBLE   gamma1,    
                     DOUBLE   gamma2,    
                     DOUBLE   gamma3,    
                     DOUBLE   gamma4,    
                     DOUBLE   dia,       
                     DOUBLE   stress[6], 
                     DOUBLE   strain[6], 
                     DOUBLE **d,         
                     INT     *iupd,      
                     INT     *yip,       
                     DOUBLE  *kappa_t,   
                     DOUBLE  *kappa_c,   
                     DOUBLE   sig[6],    
                     DOUBLE   eps[6])    
{
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE ro23;
DOUBLE dlam[2];       /*dlam[0] -> kappa_t; dlam[1] -> kappa_c*/
DOUBLE ft[4];         /*yield criterion for the 4 yield surfaces*/
DOUBLE alpha[3];
DOUBLE beta[2];
DOUBLE k_bar[3];      /*k_bar = beta * sigym*/
DOUBLE kappa_c_tr,I1_tr,dev_tr;    /*-> if projection has to be done several times*/
DOUBLE kappa_tu;
DOUBLE kappa_e;
DOUBLE kappa_cu;

DOUBLE L,R;           /*Midpoint and Radius of spheriacal cap*/
DOUBLE I1_0;          /*Schnittpunkt des Ueberganges von f2 nach f4 mit der hyd. Achse*/

DOUBLE sigma[6];        /*stress*/
DOUBLE sigma_tr[6];     /*trial stresses -> if projection has to be done several times*/
DOUBLE delsig[6],deleps[6];

DOUBLE G,K;             /*shear/bulk modulus*/
DOUBLE I1;              /*trace of tresses -> first invariant*/
DOUBLE n[6];            /*Gradient in deviatoric direction of stresses (s/|s|)*/
DOUBLE dev;             /*Norm of deviatoric stresses |s| */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_main");
#endif
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);


/*--------- check for allowed value of dia ------------------------------*/
mat_pl_epc_dia(Ec,gt,gc,ftm,fcm,&dia);

/*---- calculate some needed values that are constant within routine ---*/
K = Ec / (3. - 6.* vc);  /*bulk modulus*/
G = Ec / (2. + 2.* vc);  /*shear modlulus*/

/*-------- values for hardening/softening ----------------------------*/
/*tension*/
kappa_tu = gt / (dia*ftm); 

/*compression*/
kappa_e  = gamma4 * (fcm/Ec);
kappa_cu = (3.*gc)/(2.*dia*fcm) + kappa_e;


/*------------ alpha values ------------------------------------------*/
alpha[0] = ro23 * (gamma1*fcm - ftm)/(gamma1*fcm + ftm);
alpha[1] = ro23 * (gamma2 - 1.)/(2.*gamma2 - 1.);
alpha[2] = -1./(3.*alpha[0]);

/*------------ beta values -------------------------------------------*/
beta[0]  = (2.*gamma1*fcm)/(gamma1*fcm + ftm);
beta[1]  = gamma2/(2.*gamma2 - 1.);

/*----------------------------------------------------------------------*/
  *iupd=0;

/*----- get linear elastic isotropic material tensor -------------------*/
   mat_el_iso(Ec, vc, d); 

/*-----------------------------------------------------------------------|
|     YIP > 0  STRESSES ARE AVAILABLE FROM LAST UPDATE                   |
|         = 1  E L A S T I C                                             |
|         = 2  P L A S T I C                                             |
|     UPDATE FLAG MUST SET TO STORE CHANGE OF PARAMETER YIP              |
|     NO CHANGES HAVE BEEN MADE ON STRESS STATE                          |
|-----------------------------------------------------------------------*/
      if(*yip>0)
      {
        for (i=0; i<6; i++)
        {
          stress[i] = sig[i];
        }

        if(*yip==1)           /*elastic*/
        {
          *yip  = - *yip;
        }
        else                  /*plastic*/
        {
          dlam[0]=0.;
          dlam[1]=0.;
          
          /*--------- calculate some needed prevalues  ----------------------------*/
          mat_pl_epc_preval(stress,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                            *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

          if(*yip == 2)           /* standard DP-Mapl on 1st yield surface */
          {
              /*---- get Cep for global predictor step with dlam=0 ------*/
              mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                  ftm,fcm,alpha[0],dlam[0],Ec,vc,d,dev,n,*yip);
          }
          else if (*yip == 3)     /* Mapl for coupled yield surfaces 1 and 2*/
          {
              mat_pl_epc_mapl_ms(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                  ftm,fcm,alpha,dlam,Ec,vc,d,dev,n,*yip);
          }
          else if (*yip == 4)     /* standard DP-Mapl on 2nd yield surface */
          {
              mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                  ftm,fcm,alpha[1],dlam[1],Ec,vc,d,dev,n,*yip);
          }
          else if (*yip == 5)    /* Mapl for Cap-Region */
          {
              mat_pl_epc_mapl_cap(*kappa_c,kappa_cu,kappa_e,gamma2,gamma3,fcm,alpha[1],
                                  dlam[1],Ec,vc,R,L,stress,d);
          }

          *yip= - *yip;
        }
        *iupd=1;
        goto end;
      }
/*-----------------------------------------------------------------------|
|   1. CALCULATE INCREMENTAL STRAINS     DELEPS                          |
|   2. CALCULATE STRESS INCREMENT ASSUMING ELASTIC BEHAVIOUR             |
|   3. CALCULATE TOTAL STRESS                                            |
|   4. CHECK STRESS DEVIATOR AGAINST CURRENT YIELD SURFACE               |
|-----------------------------------------------------------------------*/
  
  for (i=0; i<6; i++) deleps[i] = strain[i] - eps[i]; 

  deleps[3] = 2. * deleps[3];
  deleps[4] = 2. * deleps[4];
  deleps[5] = 2. * deleps[5];

  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += d[i][j]*deleps[j];
  
  for (i=0; i<6; i++) sigma[i] = sig[i] + delsig[i];


/*--------- calculate some needed prevalues for trial state ----------------------*/
  mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                    *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

/*--- yield condition for all four yield surfaces --------------------------------*/
  ft[0] = dev + alpha[0]*I1 - ro23*k_bar[0];               /* 1. DP-yield surface */
  ft[1] = dev + alpha[1]*I1 - ro23*k_bar[1];               /* 2. DP-yield surface */
  ft[2] = dev + alpha[2]*I1 - ro23*k_bar[2];               /*    inverted cone    */
  ft[3] = sqrt(dev*dev + 1./9.*(I1-L)*(I1-L)) - R;         /*    spherical cap    */

/*------------- state of stress within yield surface - E L A S T I C ---*/
  if ( ft[0]<=EPS10 && ft[1]<=EPS10 && (ft[3]<EPS10 || I1 >= I1_0) )   /* case I */
  {
    *yip = 1;
    for (i=0; i<6; i++) stress[i] = sigma[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    /*-----------------------------------------------------------------------|
    |       tension- / tension-compression-region        case II             |
    |-----------------------------------------------------------------------*/
    if ( (ft[0]>EPS10 && ft[1]<=EPS10) || ft[2]<EPS10 )
    {
       *yip = 2;
       if (ft[2]>=EPS10)   /* standard DP-Projection on yield surface 1 */
       {
          /* return -> new stresses, dlam, kappa_t */
          mat_pl_epc_radi_dp(kappa_t,kappa_c,beta,k_bar,kappa_tu,kappa_cu,kappa_e,gamma3,
                             ftm,fcm,G,K,alpha[0],sigma,dlam,ft[0],n,*yip);

          /*--------- calculate some needed prevalues for updated state --------------------*/
          mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                            *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

          /* algorithmic elasto plastic tangent material matrix */
          mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                              ftm,fcm,alpha[0],dlam[0],Ec,vc,d,dev,n,*yip);
       }
       else    /* Projection in the Apex region (inverted cone) -> Multisurface with the 2 yield surfaces 1 & 3*/
       {
          /* return -> new stresses, dlam, kappa_t */
          mat_pl_epc_radi_ms(kappa_t,kappa_c,beta,k_bar,kappa_tu,kappa_cu,kappa_e,gamma3,  
                             ftm,fcm,G,K,alpha,sigma,dlam,ft[0],ft[2],n,0);    

          /*--------- calculate some needed prevalues for updated state --------------------*/
          mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                            *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

          if (dlam[0]>=EPS10) /* inverted cone is not active after return -> standard DP-Mapl*/
          {
             /* algorithmic elasto plastic tangent material matrix */
             mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                 ftm,fcm,alpha[0],dlam[0],Ec,vc,d,dev,n,*yip);
          }
          else               /* inverted cone is active after return -> Mapl for coupled yield surfaces*/
          {
             mat_pl_epc_mapl_ms(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                 ftm,fcm,alpha,dlam,Ec,vc,d,dev,n,*yip);
          }
       }
    }
    /*-----------------------------------------------------------------------| 
    |       multisurface-region (compression/tension)    case III            |
    |-----------------------------------------------------------------------*/ 
    else if ( ft[0]>EPS10 && ft[1]>EPS10 )  
    {
       /*save internal parameters and stresses if projection has to be done several times*/
       kappa_c_tr = *kappa_c;
       I1_tr      = I1;
       dev_tr     = dev;
       for (i=0; i<6; i++) sigma_tr[i] = sigma[i];

       /* return -> new stresses, dlam, kappa_t */
       mat_pl_epc_radi_ms(kappa_t,kappa_c,beta,k_bar,kappa_tu,kappa_cu,kappa_e,gamma3,  
                          ftm,fcm,G,K,alpha,sigma,dlam,ft[0],ft[1],n,1);    

       /*--------- calculate some needed prevalues for updated state --------------------*/
       mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                         *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

       /*check for cap region*/
       if (I1 < I1_0) /* cap region, make return again with trial state */
       {
          *yip = 5;

          /* return -> new stresses, dlam[1], kappa_c */
          mat_pl_epc_radi_cap(&kappa_c_tr,kappa_cu,kappa_e,gamma2,gamma3,fcm,I1_tr,dev_tr,
                              G,K,alpha[1],sigma_tr,&dlam[1]);

          /*--- write back to stored values -------*/
          *kappa_c = kappa_c_tr;
          for (i=0; i<6; i++) sigma[i] = sigma_tr[i];

          /*--------- calculate some needed prevalues for updated state --------------------*/
          mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                            *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

          /* algorithmic elasto plastic tangent material matrix */
          mat_pl_epc_mapl_cap(*kappa_c,kappa_cu,kappa_e,gamma2,gamma3,fcm,alpha[1],
                              dlam[1],Ec,vc,R,L,sigma,d);
       }
       else
       {
          *yip = 3;                    /* both yield surfaces active after return    */
          if (dlam[0] == 0.) *yip = 4; /* only 2nd yield surface active after return */
          if (dlam[1] == 0.) *yip = 2; /* only 1st yield surface active after return */
          
          if (*yip == 2 )       /* standard DP-Mapl for 1st yield surface */
          {
             /* algorithmic elasto plastic tangent material matrix */
             mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                 ftm,fcm,alpha[0],dlam[0],Ec,vc,d,dev,n,*yip);
          }
          else if (*yip == 4 )  /* standard DP-Mapl for 2nd yield surface */
          {
             /* algorithmic elasto plastic tangent material matrix */
             mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                 ftm,fcm,alpha[1],dlam[1],Ec,vc,d,dev,n,*yip);
          }
          else if (*yip == 3)   /* Mapl for coupled yield surfaces */
          {
             mat_pl_epc_mapl_ms(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                                 ftm,fcm,alpha,dlam,Ec,vc,d,dev,n,*yip);
          }
       }
    }
    /*-----------------------------------------------------------------------| 
    |       compression-region                           case IV             |
    |-----------------------------------------------------------------------*/ 
    else if ( ft[0]<=EPS10 && ft[1]>EPS10)
    {
       /*assume stresses have to be projected onto 2nd yield surface -> standard DP*/
       *yip = 4;
       
       /*save internal parameters and stresses if projection has to be done several times*/
       kappa_c_tr = *kappa_c;
       I1_tr      = I1;
       dev_tr     = dev;
       for (i=0; i<6; i++) sigma_tr[i] = sigma[i];
       
       /* return -> new stresses, dlam, kappa_t */
       mat_pl_epc_radi_dp(kappa_t,kappa_c,beta,k_bar,kappa_tu,kappa_cu,kappa_e,gamma3,
                          ftm,fcm,G,K,alpha[1],sigma,dlam,ft[1],n,*yip);

       /*--------- calculate some needed prevalues for updated state --------------------*/
       mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                         *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

       /*check for cap region */
       if (I1 < I1_0) /* cap region, make return again with trial state */
       {
          *yip = 5;

          /* return -> new stresses, dlam[1], kappa_c */
          mat_pl_epc_radi_cap(&kappa_c_tr,kappa_cu,kappa_e,gamma2,gamma3,fcm,I1_tr,dev_tr,
                              G,K,alpha[1],sigma_tr,&dlam[1]);

          /*--- write back to stored values -------*/
          *kappa_c = kappa_c_tr;
          for (i=0; i<6; i++) sigma[i] = sigma_tr[i];

          /*--------- calculate some needed prevalues for updated state --------------------*/
          mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                            *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

          /* algorithmic elasto plastic tangent material matrix */
          mat_pl_epc_mapl_cap(*kappa_c,kappa_cu,kappa_e,gamma2,gamma3,fcm,alpha[1],
                              dlam[1],Ec,vc,R,L,sigma,d);
       }
       else          /* standard DP-Mapl for 2nd yield surface */
       {
          /* algorithmic elasto plastic tangent material matrix */
          mat_pl_epc_mapl_dp(*kappa_t,*kappa_c,beta,kappa_tu,kappa_cu,kappa_e,gamma3,
                              ftm,fcm,alpha[1],dlam[1],Ec,vc,d,dev,n,*yip);
       }
    }
    /*-----------------------------------------------------------------------|
    |       cap-region                                   case V              |
    |-----------------------------------------------------------------------*/
    else if ( ft[0]<=EPS10 && ft[1]<=EPS10 && ft[3]>EPS10 )   /*&& I1 < I1_0*/
    {
       *yip = 5;

       /* return -> new stresses, dlam[1], kappa_c */
       mat_pl_epc_radi_cap(kappa_c,kappa_cu,kappa_e,gamma2,gamma3,fcm,I1,dev,
                           G,K,alpha[1],sigma,&dlam[1]);

       /*--------- calculate some needed prevalues for updated state --------------------*/
       mat_pl_epc_preval(sigma,&dev,n,&I1,gamma2,gamma3,ftm,fcm,beta,alpha,k_bar,
                         *kappa_t,kappa_tu,*kappa_c,kappa_e,kappa_cu,&L,&R,&I1_0);

       /* algorithmic elasto plastic tangent material matrix */
       mat_pl_epc_mapl_cap(*kappa_c,kappa_cu,kappa_e,gamma2,gamma3,fcm,alpha[1],
                           dlam[1],Ec,vc,R,L,sigma,d);
    }
    
    for (i=0; i<6; i++) stress[i] = sigma[i];

  }
/*----------------------------------------------------------------------*/
end:
/*Store values into Working-Array --> outside of this routine*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_main */


/*!----------------------------------------------------------------------
\brief return algorithm for concrete material model in the region of a
'standard' drucker prager cone

<pre>                                                            sh 10/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for a standard 'Drucker Prager' yield criterion 
with the hardenig law of the concrete material model  
</pre>
\param  DOUBLE   *kappa_t  (i/o) uniaxial equivalent strain for tension region -> WA
\param  DOUBLE   *kappa_c  (i/o) uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    beta[2]   (i)  beta values
\param  DOUBLE    k_bar[3]  (i)  k_bar = beta * sigym (of trial state)
\param  DOUBLE    kappa_tu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    ftm       (i)  tensile strength of concrete              
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    G         (i)  shear modulus  
\param  DOUBLE    K         (i)  bulk modulus  
\param  DOUBLE    alpha     (i)  coefficient of friction          
\param  DOUBLE    sigma[6] (i/o) trial stresses to be projected  
\param  DOUBLE    dlam[2]   (o)  plastic multiplier for tension and compression region
\param  DOUBLE    ft_tr     (i)  yield criterion of trial state 
\param  DOUBLE    n[6]      (i)  gradient n of trial state 
\param  INT       yip       (i)  flag to account for active yield surface 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()  [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_radi_dp(DOUBLE *kappa_t,
                        DOUBLE *kappa_c,
                        DOUBLE  beta[2],
                        DOUBLE  k_bar[3],
                        DOUBLE  kappa_tu,
                        DOUBLE  kappa_cu,
                        DOUBLE  kappa_e,
                        DOUBLE  gamma3,
                        DOUBLE  ftm,
                        DOUBLE  fcm,
                        DOUBLE  G,
                        DOUBLE  K,
                        DOUBLE  alpha,
                        DOUBLE  sigma[6],
                        DOUBLE  dlam[2],
                        DOUBLE  ft_tr,
                        DOUBLE  n[6],
                        INT     yip)
{
/*----------------------------------------------------------------------*/
INT    i;
DOUBLE ro23,q13;
DOUBLE dlamda,epst, epstn;
DOUBLE f,dfdl;
DOUBLE hard[2],dhdl[2];         /*evolution of hardening & their derivatives with respect to dlamda*/
DOUBLE sig_m,s[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_radi_dp");
#endif

/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;

/*--------------------------------------------- initialize variables ---*/
    i = 0;
    dlam[0] = 0.;
    dlam[1] = 0.;
    dlamda  = 0.;
 
if      (yip == 2) epstn = *kappa_t;   /* 1st yield surface */
else if (yip == 4) epstn = *kappa_c;   /* 2nd yield surface */
else     dserror("wrong parameter yip in 'mat_pl_epc_radi_dp' "); 


/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
      epst = epstn + dlamda;  
/*---- get evolution of hardening and their derivatives ----------------*/
      if      (yip == 2)   /* 1st yield surface */
      {
         *kappa_t = epst; 
          dlam[0] = dlamda;
      }  
      else if (yip == 4)   /* 2nd yield surface */
      {
         *kappa_c = epst; 
          dlam[1] = dlamda; 
      } 
      
      mat_pl_epc_hard(beta,kappa_tu,kappa_e,kappa_cu,*kappa_c,*kappa_t,ftm,fcm,gamma3,k_bar,hard,dhdl);

/*------------------------------------------------- new yield stress ---*/
      /*--- Zug-/Zug-Druck-Bereiche ---*/
	if (yip == 2) 
      {
         /*------------------------------ apply drucker prager yield criteria -----*/
	   f = ft_tr - (2.*G + 9.*K * alpha*alpha) * dlamda - hard[0];
         /*-- derivative of the yield criteria with respect to plastic increment --*/
	   dfdl = -(2.*G + 9.*K * alpha*alpha) - dhdl[0];
	   if (dfdl >= 0.)   dserror("CHECK THE SOFTENING PARAMETER for tension in 'mat_pl_epc_radi_dp' ");        
	}
      /*--- Druck-Bereich ---*/
	else if (yip == 4) 
      {
         /*------------------------------ apply drucker prager yield criteria -----*/
	   f = ft_tr - (2.*G + 9.*K * alpha*alpha) * dlamda - hard[1];
         /*-- derivative of the yield criteria with respect to plastic increment --*/
	   dfdl = -(2.*G + 9.*K * alpha*alpha) - dhdl[1];
	   if (dfdl >= 0.)   dserror("CHECK THE SOFTENING PARAMETER for compression in 'mat_pl_epc_radi_dp' ");        
	}

/*------------------------------------------- new plastic increment  ---*/
	dlamda -= f / dfdl;
/*------------------------------------------------ check convergence ---*/
	if (fabs(f) > EPS8) 
      {
	   if (i > 30) dserror("CONVERGENCE CHECK FAILED in 'mat_pl_epc_radi_dp' ");        
	   goto L600;
	}

/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/
if (i > 10) printf("Iter in 'epc_radi_dp' : %d \n",i);
/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/

/*--- update of uniaxial plastic strains and plastic multiplier --------*/
if      (yip == 2)   /* 1st yield surface */
{
   *kappa_t = epst; 
    dlam[0] = dlamda;
}  
else if (yip == 4)   /* 2nd yield surface */
{
   *kappa_c = epst; 
    dlam[1] = dlamda; 
} 

/*-- update of stresses and backstress   -------------------------------*/

/*--- hydrostatic part of sigma ------------------*/
sig_m = q13 * (sigma[0] + sigma[1] + sigma[2]);
sig_m = sig_m - 3.* K * dlamda * alpha;

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;

for (i=0; i<6; i++) s[i] = s[i] - 2.* G * dlamda * n[i];

sigma[0] = sig_m + s[0];
sigma[1] = sig_m + s[1];
sigma[2] = sig_m + s[2];
sigma[3] =         s[3];
sigma[4] =         s[4];
sigma[5] =         s[5];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_radi_dp*/


/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 standard 'Drucker Prager' within the Elasto-Plastic Concrete material
 model                                 

<pre>                                                            sh 10/03
This routine calculates the elasto-plastic consistent material tangent
for a standard 'Drucker Prager' cone
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    kappa_t   (i)  uniaxial equivalent strain for tension region -> WA
\param  DOUBLE    kappa_c   (i)  uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    beta[2]   (i)  beta values
\param  DOUBLE    kappa_tu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    ftm       (i)  tensile strength of concrete              
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    alpha     (i)  coefficient of friction of active yield surface         
\param  DOUBLE    dlam      (i)  plastic multiplier of active yield surface
\param  DOUBLE    e         (i)  youngs modulus  
\param  DOUBLE    vnu       (i)  poisson's ration  
\param  DOUBLE  **d         (o)  material matrix to be calculated   
                                 i: elastic material matric (not needed)
                                 o: konsistent elasto-plastic material matrix          
\param  DOUBLE    norm      (i)  norm of deviatoric stresses |s| (projected) 
\param  DOUBLE    n[6]      (i)  gradient n (s/|s|)  (projected)
\param  INT       yip       (i)  flag to account for active yield surface 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*-----------------------------------------------------------------------*/
void mat_pl_epc_mapl_dp(DOUBLE   kappa_t,
                        DOUBLE   kappa_c,
                        DOUBLE   beta[2],
                        DOUBLE   kappa_tu,
                        DOUBLE   kappa_cu,
                        DOUBLE   kappa_e,
                        DOUBLE   gamma3,
                        DOUBLE   ftm,
                        DOUBLE   fcm,
                        DOUBLE   alpha,
                        DOUBLE   dlam,
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6],
                        INT      yip)
{
/*----------------------------------------------------------------------*/
INT    i,j, cc, irc;
DOUBLE ro23,q13,q23;
DOUBLE fac;
DOUBLE df2dss[36];
DOUBLE dfds[6];
DOUBLE dkdk;
DOUBLE hm[6][6],CC[36], CI_vec[36];
DOUBLE nenner, zaehler1[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_mapl_dp");
#endif
/*----------------------------------------------------------------------*/

/*initialize needed matrices*/
nenner = 0.0;
for (i=0; i<6; i++)  zaehler1[i]     = 0.0;
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;
q23  = 2./3.;
/*----------------------------------------------------------------------*/

/*- get INVERSE of linear elastic isotropic material tensor --- Cel_-1 -*/
mat_el_iso_inv(e, vnu, CI_vec); /* This routine returns a vector [36] --*/

/*------  df2dss = 1/|n| * (P_dev - n dyad n) --------------------------*/
df2dss[ 0] = ( q23 - n[0] * n[0]); 
df2dss[ 1] = (-q13 - n[0] * n[1]); 
df2dss[ 2] = (-q13 - n[0] * n[2]); 
df2dss[ 3] = (     - n[0] * n[3]); 
df2dss[ 4] = (     - n[0] * n[4]); 
df2dss[ 5] = (     - n[0] * n[5]); 

df2dss[ 6] = (-q13 - n[1] * n[0]); 
df2dss[ 7] = ( q23 - n[1] * n[1]); 
df2dss[ 8] = (-q13 - n[1] * n[2]); 
df2dss[ 9] = (     - n[1] * n[3]); 
df2dss[10] = (     - n[1] * n[4]); 
df2dss[11] = (     - n[1] * n[5]); 

df2dss[12] = (-q13 - n[2] * n[0]); 
df2dss[13] = (-q13 - n[2] * n[1]); 
df2dss[14] = ( q23 - n[2] * n[2]); 
df2dss[15] = (     - n[2] * n[3]); 
df2dss[16] = (     - n[2] * n[4]); 
df2dss[17] = (     - n[2] * n[5]); 

df2dss[18] = (     - n[3] * n[0]); 
df2dss[19] = (     - n[3] * n[1]); 
df2dss[20] = (     - n[3] * n[2]); 
df2dss[21] = ( 2.  - n[3] * n[3]); 
df2dss[22] = (     - n[3] * n[4]); 
df2dss[23] = (     - n[3] * n[5]); 

df2dss[24] = (     - n[4] * n[0]); 
df2dss[25] = (     - n[4] * n[1]); 
df2dss[26] = (     - n[4] * n[2]); 
df2dss[27] = (     - n[4] * n[3]); 
df2dss[28] = ( 2.  - n[4] * n[4]); 
df2dss[29] = (     - n[4] * n[5]); 

df2dss[30] = (     - n[5] * n[0]); 
df2dss[31] = (     - n[5] * n[1]); 
df2dss[32] = (     - n[5] * n[2]); 
df2dss[33] = (     - n[5] * n[3]); 
df2dss[34] = (     - n[5] * n[4]); 
df2dss[35] = ( 2.  - n[5] * n[5]); 

/*------  Cel_-1 + dlam*df2dss -------------------------------------------*/
fac = dlam / norm;
for (i=0; i<36; i++) CI_vec[i] += fac * df2dss[i];

/*------ Theta(CC/hm) = ( Cel_-1 + dlam*df2dss )-1 -----------------------*/
c1inv6(CI_vec,CC,&irc); /*fortran routine, uses a vector*/

/*----------- write vector back to matrix -------------------*/
cc=0;
for (i=0; i<6; i++) for (j=0; j<6; j++) hm[i][j] = CC[cc++];

/*- dkdk: derivative of hardening law with respect to uniaxial equivalent strain (kappa_t, kappa_c) -*/
if (yip == 2)      /* tension */
{
   dkdk = - 1./kappa_tu * ro23 * beta[0] * ftm * exp(- kappa_t / kappa_tu);
}
else if (yip == 4) /* compression */
{
   if (kappa_c < kappa_e) 
   {
      dkdk = ro23 * beta[1] * fcm * ((2. - 2.*gamma3)/kappa_e * (1. - kappa_c/kappa_e));
   }
   else if (kappa_e <= kappa_c && kappa_c < kappa_cu)
   {
      dkdk = ro23 * beta[1] * fcm * (-2.) * (kappa_c - kappa_e)/((kappa_cu - kappa_e)*(kappa_cu - kappa_e));
   }
   else 
   {
      dkdk = ro23 * beta[1] * 1e-4;
   }
}
else  dserror("wrong parameter yip (2 or 4) in 'mat_pl_epc_mapl_dp' ");  



/*------ dfds[6] = n[6] + alpha*1 -------------------------------------------*/
dfds[0] = n[0] + alpha;
dfds[1] = n[1] + alpha;
dfds[2] = n[2] + alpha;
dfds[3] = n[3] * 2.;
dfds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
dfds[5] = n[5] * 2.;

/*----- nenner = dfds:hm:dfds + dkdk --------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) nenner += dfds[j]*hm[j][i]*dfds[i];
nenner = nenner + dkdk;

/*------- zaehler1 = hm * dfds --------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) zaehler1[i] += hm[i][j]*dfds[j];

/*---- Cep = hm - [(hm:dfds)dyad(hm:dfds)]/(dfds:hm:dfds - dfdk) ---------*/
fac = 1./nenner;

for (i=0; i<6; i++) for (j=0; j<6; j++) 
{
      d[i][j] = hm[i][j] - fac * zaehler1[i] * zaehler1[j]; 
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_mapl_dp */


/*!----------------------------------------------------------------------
\brief return algorithm for 'Drucker Prager' material model for the
inverted cone

<pre>                                                            sh 09/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for the 'Drucker Prager' yield criterion with 
combined linear iso/kin. hardenig law  
</pre>
\param  DOUBLE   *kappa_t  (i/o) uniaxial equivalent strain for tension region -> WA
\param  DOUBLE   *kappa_c  (i/o) uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    beta[2]   (i)  beta values
\param  DOUBLE    k_bar[3]  (i)  k_bar = beta * sigym (of trial state)
\param  DOUBLE    kappa_tu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    ftm       (i)  tensile strength of concrete              
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    G         (i)  shear modulus  
\param  DOUBLE    K         (i)  bulk modulus 
\param  DOUBLE    alpha[3]  (i)  coefficient of friction of 1st,2nd yield surface and inverted cone         
\param  DOUBLE    sigma[6] (i/o) trial stresses to be projected  
\param  DOUBLE    dlam[2]   (o)  plastic multiplier (dlam1 + dlam2 ; of main yield surface + inverted cone)
\param  DOUBLE    ft1       (i)  yilcr at trial state of 1st yield surface
\param  DOUBLE    ft2       (i)  yilcr at trial state of inverted cone or 2nd yield surface
\param  DOUBLE    n[6]      (i)  gradient n of trial state 
\param  INT       flag      (i)  flag if apex (flag==0) or normal ms (flag==1) 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()  [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_radi_ms(DOUBLE   *kappa_t, 
                        DOUBLE   *kappa_c, 
                        DOUBLE    beta[2], 
                        DOUBLE    k_bar[3],
                        DOUBLE    kappa_tu,
                        DOUBLE    kappa_cu,
                        DOUBLE    kappa_e, 
                        DOUBLE    gamma3,  
                        DOUBLE    ftm,     
                        DOUBLE    fcm,     
                        DOUBLE    G,       
                        DOUBLE    K,       
                        DOUBLE    alpha[3],
                        DOUBLE    sigma[6],
                        DOUBLE    dlam[2], 
                        DOUBLE    ft1,     
                        DOUBLE    ft2,     
                        DOUBLE    n[6],    
                        INT       flag)    
{
/*----------------------------------------------------------------------*/
INT    i,c1,c2,z1,z2,kflag;
DOUBLE ro23,q13;
DOUBLE epst1,epst2;
DOUBLE epstn1, epstn2;
DOUBLE hard2,dhdl2;
DOUBLE hard[2],dhdl[2];
DOUBLE alpha2;
DOUBLE f1,f2;
DOUBLE F1,F2;
DOUBLE j11,j12,j21,j22,det;
DOUBLE sig_m,s[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_radi_ms");
#endif

/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;
c1   = 1;
c2   = 1;
z1   = 1;
z2   = 1;
/*--------------------------------------------- initialize variables ---*/
i = 0;
dlam[0] = 0.;
dlam[1] = 0.;
 
epstn1 = *kappa_t;
epstn2 = *kappa_c;

if (flag==0) z2 = 0; 

/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
      epst1 = epstn1 + z1 * dlam[0];
      epst2 = epstn2 + z2 * dlam[1];
      
/*---- get evolution of hardening and their derivatives ----------------*/
      mat_pl_epc_hard(beta,kappa_tu,kappa_e,kappa_cu,epst2,epst1,ftm,fcm,gamma3,k_bar,hard,dhdl);

/*------------------------------------------------- new yield stress ---*/
      /*--- apex-Bereiche ---*/
	if (flag==0) 
      {
         hard2 = -1./(alpha[0]*alpha[0]) * hard[0];
         dhdl2 = -1./(alpha[0]*alpha[0]) * dhdl[0];
         /*------------------------------ apply drucker prager yield criteria -----*/
	   f1 = ft1 - (2.*G + 9.*K * alpha[0]*alpha[0]) * dlam[0] 
                  - (2.*G + 9.*K * alpha[0]*alpha[2]) * dlam[1]
                  - hard[0];
	   f2 = ft2 - (2.*G + 9.*K * alpha[2]*alpha[2]) * dlam[1] 
                  - (2.*G + 9.*K * alpha[0]*alpha[2]) * dlam[0]
                  - hard2;
         /*------------------------------ Jacobi-Matrix --------------------------*/
         j11 = -(2.*G + 9.*K * alpha[0]*alpha[0] + dhdl[0]);
         j12 = -(2.*G - 3.*K);
         j21 = -(2.*G - 3.*K + dhdl2);
         j22 = -(2.*G + K/(alpha[0]*alpha[0]));
	}
      /*--- normal multisurface region ---*/
	else if (flag==1) 
      {
         /*------------------------------ apply drucker prager yield criteria -----*/
	   f1 = ft1 - (2.*G + 9.*K * alpha[0]*alpha[0]) * dlam[0] 
                  - (2.*G + 9.*K * alpha[0]*alpha[1]) * dlam[1]
                  - hard[0];
	   f2 = ft2 - (2.*G + 9.*K * alpha[1]*alpha[1]) * dlam[1] 
                  - (2.*G + 9.*K * alpha[0]*alpha[1]) * dlam[0]
                  - hard[1];
         /*------------------------------ Jacobi-Matrix --------------------------*/
         j11 = -(2.*G + 9.*K * alpha[0]*alpha[0] + dhdl[0]);
         j12 = -(2.*G + 9.*K * alpha[0]*alpha[1]);
         j21 = -(2.*G + 9.*K * alpha[0]*alpha[1]);
         j22 = -(2.*G + 9.*K * alpha[1]*alpha[1] + dhdl[1]);
	}

/*---- multisurface ---*/
L700:
      kflag = 0;
      
      F1 = c1*f1 + (1 - c1) * dlam[0];
      F2 = c2*f2 + (1 - c2) * dlam[1];
      
      j11 = c1 * j11 + (1 - c1);
      j12 = c1 * j12;
      j21 = c2 * j21;
      j22 = c2 * j22 + (1 - c2);
      /*------- calculate new increment of dlam ---------*/
      det = j11*j22 - j12*j21;

      dlam[0] = dlam[0] - 1./det * ( j22*F1 - j12*F2);
      dlam[1] = dlam[1] - 1./det * (-j21*F1 + j11*F2);

      /*---- check for inactive yield surfaces ----------*/
      if (dlam[0] < 0.0) 
      {
         c1      = 0;
         kflag   = 1;
      }
      
      if (flag==0)       /*apex -> inverted Kuhn-Tucker conditions*/
      {
         if (dlam[1] > 0.0) 
         {       
             c2      = 0;
             kflag   = 1;
         }
      }
      else if (flag==1)  /* standard multisurface */
      {
         if (dlam[1] < 0.0) 
         {
             c2      = 0;
             kflag   = 1;
         }
      }

      if (kflag == 1)
      {
        dlam[0] = 0.;
        dlam[1] = 0.;
        f1      = ft1;
        f2      = ft2; 
	  goto L700;
      }
      
/*------------------------------------------------ check convergence ---*/
	if (fabs(c1*f1 + c2*f2) > EPS8) 
      {
	   if (i > 30) dserror("CONVERGENCE CHECK FAILED in 'mat_pl_epc_radi_ms' "); 
	   goto L600;
	}

/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/
if (i > 10) printf("Iter in 'epc_radi_ms' : %d \n",i);
/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/

/*-------- update of uniaxial plastic strain  ----*/
*kappa_t = epstn1 + z1 * c1 * dlam[0];
*kappa_c = epstn2 + z2 * c2 * dlam[1];


/*-- update of stresses -------------------------------*/
if      (flag==0) alpha2 = alpha[2];
else if (flag==1) alpha2 = alpha[1];

/*--- hydrostatic part of sigma ------------------*/
sig_m = q13 * (sigma[0] + sigma[1] + sigma[2]);
sig_m = sig_m - 3.* K * c1 * dlam[0] * alpha[0]
              - 3.* K * c2 * dlam[1] * alpha2;

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;

for (i=0; i<6; i++) s[i] = s[i] - 2.* G * c1 * dlam[0] * n[i]
                                - 2.* G * c2 * dlam[1] * n[i];

sigma[0] = sig_m + s[0];
sigma[1] = sig_m + s[1];
sigma[2] = sig_m + s[2];
sigma[3] =         s[3];
sigma[4] =         s[4];
sigma[5] =         s[5];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_radi_ms */


/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 multisurface 'Drucker Prager' within the Elasto-Plastic Concrete material
 model (yip=2: apex, yip=3: multisurface)                                

<pre>                                                            sh 10/03
This routine calculates the elasto-plastic consistent material tangent
for multisurface 'Drucker Prager'
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    kappa_t   (i)  uniaxial equivalent strain for tension region -> WA
\param  DOUBLE    kappa_c   (i)  uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    beta[2]   (i)  beta values
\param  DOUBLE    kappa_tu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    ftm       (i)  tensile strength of concrete              
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    alpha[3]  (i)  coefficient of friction          
\param  DOUBLE    dlam[2]   (i)  plastic multiplier (dlam1 & dlam2)
\param  DOUBLE    e         (i)  youngs modulus  
\param  DOUBLE    vnu       (i)  poisson's ration  
\param  DOUBLE  **d         (o)  material matrix to be calculated   
                                 i: elastic material matric (not needed)
                                 o: konsistent elasto-plastic material matrix          
\param  DOUBLE    norm      (i)  norm of deviatoric stresses |s| (projected) 
\param  DOUBLE    n[6]      (i)  gradient n (s/|s|)  (projected)
\param  INT       yip       (i)  flag to account for active yield surface 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*-----------------------------------------------------------------------*/
void mat_pl_epc_mapl_ms(DOUBLE   kappa_t,
                        DOUBLE   kappa_c,
                        DOUBLE   beta[2],
                        DOUBLE   kappa_tu,
                        DOUBLE   kappa_cu,
                        DOUBLE   kappa_e,
                        DOUBLE   gamma3,
                        DOUBLE   ftm,
                        DOUBLE   fcm,
                        DOUBLE   alpha[3],
                        DOUBLE   dlam[2],
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6],
                        INT      yip)
{
/*----------------------------------------------------------------------*/
INT    i,j,k,cc,irc;
DOUBLE ro23,q13,q23;
DOUBLE dlam12;
DOUBLE alpha2;
DOUBLE det;
DOUBLE dhdl[2];
DOUBLE df1ds[6], df2ds[6];
DOUBLE E[2][2],E_help[2][2],U[6][2];
DOUBLE hm[6][6];
DOUBLE H1[6][2], H2[6][2];
DOUBLE HM2[6][6];
DOUBLE CC[36], CI_vec[36];
DOUBLE fac;
DOUBLE df2dss[36];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_mapl_ms");
#endif
/*----------------------------------------------------------------------*/
dlam12 = dlam[0] + dlam[1];
ro23 = sqrt(2./3.);
q13  = 1./3.;
q23  = 2./3.;
/*----------------------------------------------------------------------*/

/*initialize needed matrices*/
for (i=0; i<2; i++) for (j=0; j<2; j++) E_help[i][j] = 0.0;
for (i=0; i<6; i++) for (j=0; j<2; j++) H1[i][j]     = 0.0;
for (i=0; i<6; i++) for (j=0; j<2; j++) H2[i][j]     = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++) HM2[i][j]    = 0.0;
/*----------------------------------------------------------------------*/

/*- get INVERSE of linear elastic isotropic material tensor --- Cel_-1 -*/
mat_el_iso_inv(e, vnu, CI_vec); /* This routine returns a vector [36] --*/

/*------  df2dss = 1/|n| * (P_dev - n dyad n) --------------------------*/
df2dss[ 0] = ( q23 - n[0] * n[0]); 
df2dss[ 1] = (-q13 - n[0] * n[1]); 
df2dss[ 2] = (-q13 - n[0] * n[2]); 
df2dss[ 3] = (     - n[0] * n[3]); 
df2dss[ 4] = (     - n[0] * n[4]); 
df2dss[ 5] = (     - n[0] * n[5]); 

df2dss[ 6] = (-q13 - n[1] * n[0]); 
df2dss[ 7] = ( q23 - n[1] * n[1]); 
df2dss[ 8] = (-q13 - n[1] * n[2]); 
df2dss[ 9] = (     - n[1] * n[3]); 
df2dss[10] = (     - n[1] * n[4]); 
df2dss[11] = (     - n[1] * n[5]); 

df2dss[12] = (-q13 - n[2] * n[0]); 
df2dss[13] = (-q13 - n[2] * n[1]); 
df2dss[14] = ( q23 - n[2] * n[2]); 
df2dss[15] = (     - n[2] * n[3]); 
df2dss[16] = (     - n[2] * n[4]); 
df2dss[17] = (     - n[2] * n[5]); 

df2dss[18] = (     - n[3] * n[0]); 
df2dss[19] = (     - n[3] * n[1]); 
df2dss[20] = (     - n[3] * n[2]); 
df2dss[21] = ( 2.  - n[3] * n[3]); 
df2dss[22] = (     - n[3] * n[4]); 
df2dss[23] = (     - n[3] * n[5]); 

df2dss[24] = (     - n[4] * n[0]); 
df2dss[25] = (     - n[4] * n[1]); 
df2dss[26] = (     - n[4] * n[2]); 
df2dss[27] = (     - n[4] * n[3]); 
df2dss[28] = ( 2.  - n[4] * n[4]); 
df2dss[29] = (     - n[4] * n[5]); 

df2dss[30] = (     - n[5] * n[0]); 
df2dss[31] = (     - n[5] * n[1]); 
df2dss[32] = (     - n[5] * n[2]); 
df2dss[33] = (     - n[5] * n[3]); 
df2dss[34] = (     - n[5] * n[4]); 
df2dss[35] = ( 2.  - n[5] * n[5]); 

/*------  Cel_-1 + dlam*df2dss -------------------------------------------*/
fac = dlam12 / norm;
for (i=0; i<36; i++) CI_vec[i] += fac * df2dss[i];

/*------ Theta(CC/hm) = ( Cel_-1 + dlam*df2dss )-1 -----------------------*/
c1inv6(CI_vec,CC,&irc); /*fortran routine, uses a vector*/
/*----------- write vector back to matrix -------------------*/
cc=0;
for (i=0; i<6; i++) for (j=0; j<6; j++) hm[i][j] = CC[cc++];


/*- derivative of hardening law with respect to uniaxial equivalent strain (kappa_t, kappa_c) -*/

/* tension */
dhdl[0] = - 1./kappa_tu *  ro23 * beta[0] * ftm * exp(- kappa_t / kappa_tu);

/* compression */
if (kappa_c < kappa_e) 
{
   dhdl[1] = ro23 * beta[1] * fcm * ((2. - 2.*gamma3)/kappa_e * (1. - kappa_c/kappa_e));
}
else if (kappa_e <= kappa_c && kappa_c < kappa_cu)
{
   dhdl[1] = ro23 * beta[1] * fcm * (-2.) * (kappa_c - kappa_e)/((kappa_cu - kappa_e)*(kappa_cu - kappa_e));
}
else 
{
   dhdl[1] = ro23 * beta[1] * 1e-4;
}
/* apex */
dhdl[2] = - 1./(alpha[0]*alpha[0] * dhdl[0]);

/*----- Matrix of plastic Modules --> see Menrath (A.8)  ---------------------*/
if (yip == 2)      /*apex region*/
{
   alpha2  = alpha[2]; 
   E[0][0] = dhdl[0];
   E[0][1] = 0.0;
   E[1][0] = dhdl[2];
   E[1][1] = 0.0;
}
else if (yip == 3) /*multisurface region*/
{
   alpha2  = alpha[1]; 
   E[0][0] = dhdl[0];
   E[0][1] = 0.0;
   E[1][0] = 0.0;
   E[1][1] = dhdl[1];
}
else dserror("wrong parameter yip ( != 2 or 3) -> 'mat_pl_epc_mapl_ms' ");

/*------ df1ds[6] = n[6] + alpha1*1 -------------------------------------------*/
df1ds[0] = n[0] + alpha[0];
df1ds[1] = n[1] + alpha[0];
df1ds[2] = n[2] + alpha[0];
df1ds[3] = n[3] * 2.;
df1ds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
df1ds[5] = n[5] * 2.;

/*------ df2ds[6] = n[6] + alpha2*1 -------------------------------------------*/
df2ds[0] = n[0] + alpha2;
df2ds[1] = n[1] + alpha2;
df2ds[2] = n[2] + alpha2;
df2ds[3] = n[3] * 2.;
df2ds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
df2ds[5] = n[5] * 2.;

/*----- E_help = UT*hm*U --------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[0][0] += df1ds[j] * hm[j][i] * df1ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[0][1] += df1ds[j] * hm[j][i] * df2ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[1][0] += df2ds[j] * hm[j][i] * df1ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[1][1] += df2ds[j] * hm[j][i] * df2ds[i];
              
/*------ E_help = E + UT*hm*U ------------------------------*/
E_help[0][0] += E[0][0];
E_help[0][1] += E[0][1];
E_help[1][0] += E[1][0];
E_help[1][1] += E[1][1];

/*------------ E = ( E + UT*hm*U )-1 -----------------------*/
det = E_help[0][0]*E_help[1][1] - E_help[0][1]*E_help[1][0];
det = 1./det;

E[0][0] =  det * E_help[1][1];
E[0][1] = -det * E_help[0][1];
E[1][0] = -det * E_help[1][0];
E[1][1] =  det * E_help[0][0];

/*-----------------------------------------------------------*/
U[0][0] = df1ds[0];
U[1][0] = df1ds[1];
U[2][0] = df1ds[2];
U[3][0] = df1ds[3];
U[4][0] = df1ds[4];
U[5][0] = df1ds[5];

U[0][1] = df2ds[0];
U[1][1] = df2ds[1];
U[2][1] = df2ds[2];
U[3][1] = df2ds[3];
U[4][1] = df2ds[4];
U[5][1] = df2ds[5];

/*--------- H1[6][2] = hm[6][6]*U[6][2] ---------------------*/
for (j=0; j<6; j++) for (k=0; k<6; k++) for (i=0; i<2; i++) H1[j][i] += hm[j][k] * U[k][i]; 

/*--------- H2[6][2] = hm[6][6]*U[6][2]*E[2][2] = H1[6][2]*E[2][2] -----*/
for (j=0; j<6; j++) for (k=0; k<2; k++) for (i=0; i<2; i++) H2[j][i] += H1[j][k] * E[k][i]; 

/*--------- HM2[6][6] = hm[6][6]*U[6][2]*E[2][2]*UT[2][6]*hm[6][6] = H2[6][2] * H1T[2][6] ----*/
for (j=0; j<6; j++) for (k=0; k<2; k++) for (i=0; i<6; i++) HM2[j][i] += H2[j][k] * H1[i][k]; 

/*--- Cep = hm - HM2 ---------------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) d[i][j] = hm[i][j] - HM2[i][j]; 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_mapl_ms */


/*!----------------------------------------------------------------------
\brief return algorithm for concrete material model in the cap region 

<pre>                                                            sh 10/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for the spherical cap region 
</pre>
\param  DOUBLE   *kappa_c  (i/o) uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma2    (i)  fitting parameter (constant) 
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    I1_tr     (i)  1st invariant of trial stresses             
\param  DOUBLE    dev_tr    (i)  norm of deviatoric trial stresses             
\param  DOUBLE    G         (i)  shear modulus  
\param  DOUBLE    K         (i)  bulk modulus  
\param  DOUBLE    alpha2    (i)  coefficient of friction of second DP cone         
\param  DOUBLE    sigma[6] (i/o) trial stresses to be projected  
\param  DOUBLE   *dlam      (o)  plastic multiplier for tension and compression region

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()  [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_radi_cap(DOUBLE *kappa_c,
                         DOUBLE  kappa_cu,
                         DOUBLE  kappa_e,
                         DOUBLE  gamma2,
                         DOUBLE  gamma3,
                         DOUBLE  fcm,
                         DOUBLE  I1_tr,
                         DOUBLE  dev_tr,
                         DOUBLE  G,
                         DOUBLE  K,
                         DOUBLE  alpha2,
                         DOUBLE  sigma[6],
                         DOUBLE *dlam)
{
/*----------------------------------------------------------------------*/
INT    i;
DOUBLE q13;
DOUBLE epst, epstn, depst;
DOUBLE f1,f2,f;
DOUBLE sig_m,s[6];
DOUBLE sig_bar;
DOUBLE R,L;
DOUBLE k1;		/* L/R remains constant */
DOUBLE rdevc,rdevn,rcomc,rcomn;
DOUBLE a,b;
DOUBLE root;
DOUBLE dSdk,dRdk,dLdk;  /* derivatives of sig_bar, R and L with respect to depst */
DOUBLE dadl,dadk,dbdl,dbdk;
DOUBLE drootdl, drootdk;
DOUBLE sol11,sol12,sol21,sol22,soli11,soli12,soli21,soli22;
DOUBLE det;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_radi_cap");
#endif

/*----------------------------------------------------------------------*/
q13  = 1./3.;

/*--------------------------------------------- initialize variables ---*/
    i = 0;
   *dlam  = 0.;
    depst = 0.;
 
    epstn = *kappa_c;   /* 2nd yield surface */

/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
     ++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
     epst = epstn + depst;  
     
/*------------------------------- calculate new sig_bar, L, R ----------*/
    *kappa_c = epst;
    
    if (*kappa_c < kappa_e) 
    {
       sig_bar  = fcm * (gamma3 + 2.*(1.-gamma3)* *kappa_c/(kappa_e) - (1.-gamma3)*(*kappa_c/kappa_e)*(*kappa_c/kappa_e));
       dSdk     = fcm * ((2. - 2.*gamma3)/kappa_e * (1. - *kappa_c/kappa_e));   
    }
    else if (kappa_e <= *kappa_c && *kappa_c < kappa_cu)
    {
       sig_bar  = fcm * (1. - ((*kappa_c - kappa_e)/(kappa_cu - kappa_e))*((*kappa_c - kappa_e)/(kappa_cu - kappa_e)));
       dSdk     = fcm * (-2.) * (*kappa_c - kappa_e)/((kappa_cu - kappa_e)*(kappa_cu - kappa_e));
    }
    else  
    {
       sig_bar  = EPS6;
       dSdk     = EPS6;
    }
    
    if (sig_bar < 0.03*fcm)  
    {
       sig_bar  = 0.03*fcm;
       dSdk     = EPS6;
    }

    /*------------ L and R for spherical cap -------------------------------*/
    L  = -(sqrt(54.)*alpha2 + 2.)*gamma2*sig_bar;
    R  = sqrt(2./3. + 6.*alpha2*alpha2)*gamma2*sig_bar;
    k1 = L/R;
    
    dLdk = L/sig_bar * dSdk;
    dRdk = R/sig_bar * dSdk;
    
/*---- cap yield criteria ...... -----------------*/
    rdevc = dev_tr;
    rdevn = R + 2.*G*(*dlam);
    rcomc = I1_tr - L;
    rcomn = R + K*(*dlam);
    
    a     = rdevc/rdevn;
    b     = rcomc/rcomn;
    root  = sqrt(a*a + b*b/9.);
    
    f1  = a*a + b*b/9. - 1.;
    f2  = depst - (*dlam)* (k1 * b / (root*9.) + 1.);

/*---- derivative of both equations with respect to dlam and kappa ----------*/
    dadl  = - rdevc * 2. * G / (rdevn*rdevn);
    dadk  = - rdevc * dRdk  / (rdevn*rdevn);
    dbdl  = - rcomc * K  / (rcomn*rcomn);
    dbdk  = - (rcomn*dLdk + rcomc*dRdk) / (rcomn*rcomn);
    
    drootdl = 1./(2.*root) * (2.*a*dadl + 2./9.*b*dbdl);
    drootdk = 1./(2.*root) * (2.*a*dadk + 2./9.*b*dbdk);
            
    sol11 = 2.*a*dadl + 2./9.*b*dbdl;
    sol12 = 2.*a*dadk + 2./9.*b*dbdk;

    sol21 = -(b*k1/(9.*root) +1. ) - (*dlam)*( k1/9. * (root*dbdl - b*drootdl)/(root*root) );
    sol22 = 1. - (*dlam)*(k1/9.* (root*dbdk - b*drootdk)/(root*root) );  
     
    if (sol11 >= 0.)  dserror("CHECK THE SOFTENING PARAMETER");
    
/*-------------- solution equation system ------------------------------*/
     det = sol11 *sol22 - sol21*sol12;
     soli11 =   sol22/det;
     soli21 = - sol21/det;
     soli12 = - sol12/det;
     soli22 =   sol11/det;
     
    *dlam  = *dlam  - (soli11*f1 + soli12*f2);
     depst =  depst - (soli21*f1 + soli22*f2);
     
     f = sqrt( (f1/EPS6)*(f1/EPS6) + (f2/EPS6)*(f2/EPS6) );

/*------------------------------------------------ check convergence ---*/
	if (fabs(f) > EPS6) 
      {
	   if (i > 60) dserror("CONVERGENCE CHECK FAILED in 'mat_pl_epc_radi_cap' ");     
	   goto L600;
	}

/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/
if (i > 10) printf("Iter in 'epc_radi_cap' : %d \n",i);
/*TEST**TEST**TEST**TEST**TEST**TEST**TEST*/

/*-- update of uniaxial plastic strain --*/
   *kappa_c = epst; 


/*-- update of stresses ------------------------------------------*/

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;

for (i=0; i<6; i++) s[i] = s[i]*(R/rdevn);;

/*--- hydrostatic part of sigma ------------------*/
sig_m = q13 * (sigma[0] + sigma[1] + sigma[2]);
sig_m = sig_m - 3.* K * (*dlam) * b/(9.*root);

sigma[0] = sig_m + s[0];
sigma[1] = sig_m + s[1];
sigma[2] = sig_m + s[2];
sigma[3] =         s[3];
sigma[4] =         s[4];
sigma[5] =         s[5];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_radi_cap*/


/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 standard 'Drucker Prager' within the Elasto-Plastic Concrete material
 model                                 

<pre>                                                            sh 10/03
This routine calculates the elasto-plastic consistent material tangent
for the spherical cap
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    kappa_c   (i)  uniaxial equivalent strain for compression region -> WA
\param  DOUBLE    kappa_cu  (i)  constant calculated in prevalues
\param  DOUBLE    kappa_e   (i)  constant calculated in prevalues
\param  DOUBLE    gamma2    (i)  fitting parameter (constant) 
\param  DOUBLE    gamma3    (i)  fitting parameter (constant) 
\param  DOUBLE    fcm       (i)  compressive strength of concrete              
\param  DOUBLE    alpha2    (i)  coefficient of friction of 2nd cone         
\param  DOUBLE    dlam      (i)  plastic multiplier of active yield surface
\param  DOUBLE    e         (i)  youngs modulus  
\param  DOUBLE    vnu       (i)  poisson's ration  
\param  DOUBLE    R         (i)  radius of cap  
\param  DOUBLE    L         (i)  center of cap  
\param  DOUBLE    sigma[6]  (i)  updated stresses  
\param  DOUBLE  **d         (o)  material matrix to be calculated   
                                 i: elastic material matric (not needed)
                                 o: konsistent elasto-plastic material matrix          

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*-----------------------------------------------------------------------*/
void mat_pl_epc_mapl_cap(DOUBLE   kappa_c,
                         DOUBLE   kappa_cu,
                         DOUBLE   kappa_e,
                         DOUBLE   gamma2,
                         DOUBLE   gamma3,
                         DOUBLE   fcm,
                         DOUBLE   alpha2,
                         DOUBLE   dlam,
                         DOUBLE   e,
                         DOUBLE   vnu,
                         DOUBLE   R,
                         DOUBLE   L,
                         DOUBLE   sigma[6],
                         DOUBLE **d)
{
/*----------------------------------------------------------------------*/
INT    i,j, cc, irc;
DOUBLE ro23,q19,q13,q23,fac;
DOUBLE I1,s[6];       /* 1st invariant of updated stresses and deviatoric stresses*/
DOUBLE k1;

DOUBLE z11,z12,z21,z22;
DOUBLE zi11,zi12,zi21,zi22;
DOUBLE det;

DOUBLE Hk,dHdk;     /* H(k) and dH(k)/dk */
DOUBLE dSdk,dS2dkk; /* dsig_bar/dk and d^2sig_bar/dk^2 */
DOUBLE IL_R;        /* (I1-L)/R */

DOUBLE dfds[6];
DOUBLE df2dss[36];
DOUBLE dfdk;
DOUBLE df2dsk[6];
DOUBLE df2dkk;

DOUBLE dkapp;

DOUBLE N1[6], N2[6];

DOUBLE hm[6][6],CC[36], CI_vec[36];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_mapl_cap");
#endif
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;
q23  = 2./3.;
q19  = 1./9.;

k1   = L/R;
 
/*----------------------------------------------------------------------*/
z11 = 0.;
z12 = 0.;
z21 = 0.;
z22 = 0.;

for (i=0; i<6; i++) N1[i] = 0.;
for (i=0; i<6; i++) N2[i] = 0.;

/*-------- I1 = trace(sigma) -------------------------------------------*/
I1 = sigma[0] + sigma[1] + sigma[2];

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;


/*- get INVERSE of linear elastic isotropic material tensor --- Cel_-1 -*/
mat_el_iso_inv(e, vnu, CI_vec); /* This routine returns a vector [36] --*/

/*----------------------------------------------------------------------*/

/*------ calculate useful prevalues and 1st order gradients ------------*/
IL_R  = (I1-L)/R;
dkapp = (q19*IL_R*k1 + 1.) * dlam;

if (kappa_c < kappa_e) 
{
   dSdk     = fcm * ((2. - 2.*gamma3)/kappa_e * (1. - kappa_c/kappa_e)); 
   dS2dkk   = -2.*fcm*(1.-gamma3)/(kappa_e*kappa_e);  
}
else if (kappa_e <= kappa_c && kappa_c < kappa_cu)
{
   dSdk     = fcm * (-2.) * (kappa_c - kappa_e)/((kappa_cu - kappa_e)*(kappa_cu - kappa_e));
   dS2dkk   = -2.*fcm/((kappa_cu*kappa_e)*(kappa_cu*kappa_e));
}
else  
{
   dSdk     = EPS6;
   dS2dkk   = 0.0;
}

Hk   = sqrt(q23 + 6.*alpha2*alpha2)*gamma2*dSdk;
dHdk = sqrt(q23 + 6.*alpha2*alpha2)*gamma2*dS2dkk;

dfds[0] = s[0]/R + q19*IL_R;
dfds[1] = s[1]/R + q19*IL_R;
dfds[2] = s[2]/R + q19*IL_R;
dfds[3] = s[3]/R;
dfds[4] = s[4]/R;
dfds[5] = s[5]/R;

dfdk = -(q19*IL_R*k1 + 1)*Hk;

/*------- calculate 2nd order gradients --------------------------------*/

/*------  df2dsk = - 1/(9*R) * (11 - (I1-L)/R * dfds) * k1 * Hk ----------*/
fac = q19/R;

df2dsk[0] = - ( fac * (1 - IL_R * dfds[0]) * k1 * Hk);
df2dsk[1] = - ( fac * (1 - IL_R * dfds[1]) * k1 * Hk);
df2dsk[2] = - ( fac * (1 - IL_R * dfds[2]) * k1 * Hk);
df2dsk[3] = - ( fac * (    IL_R * dfds[3]) * k1 * Hk);
df2dsk[4] = - ( fac * (    IL_R * dfds[4]) * k1 * Hk);
df2dsk[5] = - ( fac * (    IL_R * dfds[5]) * k1 * Hk);

/*--  df2dkk =  (k1^2 * Hk^2)/(9*R) * (1 - 1/9*((I1-L)/R)^2) + dfdk/Hk * dHdk ---*/
df2dkk = (k1*k1*Hk*Hk)/(9.*R) * ( 1.- q19*IL_R*IL_R) + dfdk/Hk * dHdk;

/*------  df2dss = 1/R * (P_dev - dfds dyad dfds + 1/9* 11 dyad 11) ----*/
fac = 1./R;

df2dss[ 0] = fac * ( q23 - dfds[0] * dfds[0] + q19); 
df2dss[ 1] = fac * (-q13 - dfds[0] * dfds[1] + q19); 
df2dss[ 2] = fac * (-q13 - dfds[0] * dfds[2] + q19); 
df2dss[ 3] = fac * (     - dfds[0] * dfds[3]); 
df2dss[ 4] = fac * (     - dfds[0] * dfds[4]); 
df2dss[ 5] = fac * (     - dfds[0] * dfds[5]); 

df2dss[ 6] = fac * (-q13 - dfds[1] * dfds[0] + q19); 
df2dss[ 7] = fac * ( q23 - dfds[1] * dfds[1] + q19); 
df2dss[ 8] = fac * (-q13 - dfds[1] * dfds[2] + q19); 
df2dss[ 9] = fac * (     - dfds[1] * dfds[3]); 
df2dss[10] = fac * (     - dfds[1] * dfds[4]); 
df2dss[11] = fac * (     - dfds[1] * dfds[5]); 

df2dss[12] = fac * (-q13 - dfds[2] * dfds[0] + q19); 
df2dss[13] = fac * (-q13 - dfds[2] * dfds[1] + q19); 
df2dss[14] = fac * ( q23 - dfds[2] * dfds[2] + q19); 
df2dss[15] = fac * (     - dfds[2] * dfds[3]); 
df2dss[16] = fac * (     - dfds[2] * dfds[4]); 
df2dss[17] = fac * (     - dfds[2] * dfds[5]); 

df2dss[18] = fac * (     - dfds[3] * dfds[0]); 
df2dss[19] = fac * (     - dfds[3] * dfds[1]); 
df2dss[20] = fac * (     - dfds[3] * dfds[2]); 
df2dss[21] = fac * ( 2.  - dfds[3] * dfds[3]); 
df2dss[22] = fac * (     - dfds[3] * dfds[4]); 
df2dss[23] = fac * (     - dfds[3] * dfds[5]); 

df2dss[24] = fac * (     - dfds[4] * dfds[0]); 
df2dss[25] = fac * (     - dfds[4] * dfds[1]); 
df2dss[26] = fac * (     - dfds[4] * dfds[2]); 
df2dss[27] = fac * (     - dfds[4] * dfds[3]); 
df2dss[28] = fac * ( 2.  - dfds[4] * dfds[4]); 
df2dss[29] = fac * (     - dfds[4] * dfds[5]); 

df2dss[30] = fac * (     - dfds[5] * dfds[0]); 
df2dss[31] = fac * (     - dfds[5] * dfds[1]); 
df2dss[32] = fac * (     - dfds[5] * dfds[2]); 
df2dss[33] = fac * (     - dfds[5] * dfds[3]); 
df2dss[34] = fac * (     - dfds[5] * dfds[4]); 
df2dss[35] = fac * ( 2.  - dfds[5] * dfds[5]); 

/*------  Cel_-1 + dlam*df2dss -------------------------------------------*/
for (i=0; i<36; i++) CI_vec[i] += dlam * df2dss[i];

/*------ Theta(CC/hm) = ( Cel_-1 + dlam*df2dss )-1 -----------------------*/
c1inv6(CI_vec,CC,&irc); /*fortran routine, uses a vector*/

/*----------- write vector back to matrix -------------------*/
cc=0;
for (i=0; i<6; i++) for (j=0; j<6; j++) hm[i][j] = CC[cc++];

/*------ dfds[6] -> write as vector ---------------------------------------*/
dfds[3] = dfds[3] * 2.;
dfds[4] = dfds[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
dfds[5] = dfds[5] * 2.;

/*------ df2dsk[6] -> write as vector ---------------------------------------*/
df2dsk[3] = df2dsk[3] * 2.;
df2dsk[4] = df2dsk[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
df2dsk[5] = df2dsk[5] * 2.;


if (fabs(dlam) < EPS10)  /*  in here are the limit values if dlam -> 0  */
{ 
   /*----- zi11 = (Hk - dHdk) / ((dfds:hm:dfds) * (Hk - dHdk) - (dfdk*dfdk)) ---*/
   for (i=0; i<6; i++) for (j=0; j<6; j++) zi11 += dfds[j]*hm[j][i]*dfds[i];
   zi11 = (Hk - dHdk) / (zi11 * (Hk - dHdk) - (dfdk * dfdk));

   /*----- zi12 = zi21 = zi22 = 0.0 --------------------------------------------*/
   zi12 = zi21 = zi22 = 0.0;
}
else   /* -- normal way of calculating the entries of the Z-Matrix -----*/
{
   /*----- z11 = dfds:hm:dfds  ------------------------------------------------*/
   for (i=0; i<6; i++) for (j=0; j<6; j++) z11 += dfds[j]*hm[j][i]*dfds[i];

   /*----- z12 = dfds:hm:df2dsk - 1/dlam*dfdk = z21 ----------------------------*/
   for (i=0; i<6; i++) for (j=0; j<6; j++) z12 += dfds[j]*hm[j][i]*df2dsk[i];
   z12 = z12 - 1./dlam * dfdk;

   z21 = z12;

   /*----- z22 = df2dsk:hm:df2dsk - 1/dlam*df2dkk - 1/(dlam^2)*(Hk - dHdk)  ---*/
   for (i=0; i<6; i++) for (j=0; j<6; j++) z22 += df2dsk[j]*hm[j][i]*df2dsk[i];
   z22 = z22 - 1./dlam * df2dkk - 1./(dlam*dlam) * (Hk + dkapp*dHdk);

   /*-------------- solution equation system ------------------------------*/
   det = z11 *z22 - z21*z12;
   zi11 =   z22/det;
   zi21 = - z21/det;
   zi12 = - z12/det;
   zi22 =   z11/det;
}

/*------- N1 = hm * dfds --------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) N1[i] += hm[i][j]*dfds[j];

/*------- N2 = hm * df2dsk --------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) N2[i] += hm[i][j]*df2dsk[j];

/*                 _2_  _2_                                            */
/*---- Cep = hm -  \   \     zi_ij * N_i dyad N_j              --------*/
/*                 /__ /__                                             */
/*                 i=1 i=1                                             */

for (i=0; i<6; i++) for (j=0; j<6; j++) 
{
      d[i][j] = hm[i][j] - (  zi11 * N1[i] * N1[j] 
                            + zi12 * N1[i] * N2[j]
                            + zi21 * N2[i] * N1[j] 
                            + zi22 * N2[i] * N2[j]); 
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_mapl_cap */


/*!----------------------------------------------------------------------
\brief checks if dia is in the allowable limits

<pre>                                                            sh 10/03

</pre>
\param  DOUBLE    Ec    (i)    youngs modulus of concrete
\param  DOUBLE    gt    (i)    tensile fracture energy of concrete
\param  DOUBLE    gc    (i)    compressive fracture energy of concrete
\param  DOUBLE    ftm   (i)    tensile strength of concrete
\param  DOUBLE    fcm   (i)    compressive strength of concrete
\param  DOUBLE   *dia  (i/o)   equivalent element length

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_dia(DOUBLE    Ec, 
                    DOUBLE    gt, 
                    DOUBLE    gc, 
                    DOUBLE    ftm,
                    DOUBLE    fcm,
                    DOUBLE   *dia)
{
/*----------------------------------------------------------------------*/
DOUBLE dia_t, dia_c, dia_h;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_dia");
#endif
/*----------------------------------------------------------------------*/
dia_t =         (Ec * gt)/(ftm*ftm);
dia_c = 3./4. * (Ec * gc)/(fcm*fcm);

if (dia_t < dia_c) dia_h = dia_t;
else               dia_h = dia_c;

if (dia_h < *dia)  *dia = dia_h;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_dia */



/*!----------------------------------------------------------------------
\brief calculates the evolution of the hardening and their derivatives
with respect to dlam

<pre>                                                            sh 10/03

</pre>
\param  DOUBLE    beta2]    (i)   beta values for 1st and 2nd yield surface
\param  DOUBLE    kappa_tu  (i)   
\param  DOUBLE    kappa_e   (i)   
\param  DOUBLE    kappa_cu  (i)   
\param  DOUBLE    kappa_c   (i)   actual uniaxial equivalent strain for compression region (kappa_c + dlam)
\param  DOUBLE    kappa_t   (i)   actual uniaxial equivalent strain for tension region (kappa_t + dlam)
\param  DOUBLE    ftm       (i)   tensile strength of concrete
\param  DOUBLE    fcm       (i)   compressive strength of concrete
\param  DOUBLE    gamma3    (i)   fitting parameter
\param  DOUBLE    k_bar[3]  (i)   uniaxial equivalent yield stress (trial)
\param  DOUBLE    hard[2]   (o)   incremental hardening for tension and compression
\param  DOUBLE    dhdl[2]   (o)   derivative of hard with respect to dlam

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_hard(DOUBLE    beta[2], 
                     DOUBLE    kappa_tu,
                     DOUBLE    kappa_e, 
                     DOUBLE    kappa_cu,
                     DOUBLE    kappa_c, 
                     DOUBLE    kappa_t,
                     DOUBLE    ftm,
                     DOUBLE    fcm,
                     DOUBLE    gamma3,
                     DOUBLE    k_bar[3],
                     DOUBLE    hard[2],
                     DOUBLE    dhdl[2])
{
/*----------------------------------------------------------------------*/
DOUBLE ro23;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_hard");
#endif
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);

/* Evolution for tension */
hard[0] = ro23 * beta[0] * ftm * exp(- kappa_t / kappa_tu);
dhdl[0] = - 1./kappa_tu * hard[0];                           /*derivative -> kappa_t + dlam */


/* Evolution for compression */
if (kappa_c < kappa_e) 
{
   hard[1] = ro23 * beta[1] * fcm * (gamma3 + 2.*(1.-gamma3)*kappa_c/(kappa_e) - (1.-gamma3)*(kappa_c/(kappa_e))*(kappa_c/(kappa_e)));
   dhdl[1] = ro23 * beta[1] * fcm * ((2. - 2.*gamma3)/kappa_e * (1. - kappa_c/kappa_e));
}
else if (kappa_e <= kappa_c && kappa_c < kappa_cu)
{
   hard[1] = ro23 * beta[1] * fcm * (1. - ((kappa_c - kappa_e)/(kappa_cu - kappa_e))*((kappa_c - kappa_e)/(kappa_cu - kappa_e)));
   dhdl[1] = ro23 * beta[1] * fcm * (-2.) * (kappa_c - kappa_e)/((kappa_cu - kappa_e)*(kappa_cu - kappa_e));
}
else 
{
   hard[1] = ro23 * beta[1] * 0.03*fcm;
   dhdl[1] = ro23 * beta[1] * 1e-4;
}

if (hard[1] < ro23*beta[1]*0.03*fcm)
{
   hard[1] = ro23 * beta[1] * 0.03*fcm;
   dhdl[1] = ro23 * beta[1] * 1e-4;
}

/*--- get difference of actual step (kappa_c+dlam) and last konverged step (kappa_c of trial state -> k_bar)*/
/* DELTA = (kappa_t + dlam) - kappa_t = actual - trial (converged from last step) */
hard[0] = hard[0] - ro23 * k_bar[0];    
hard[1] = hard[1] - ro23 * k_bar[1];  

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_hard */


/*!----------------------------------------------------------------------
\brief calculates some needed values, that change with update of stresses                            

<pre>                                                            sh 10/03
This routine calculates the 
- dev = norm of deviatoric stresses |s|
- the gradient n (dfds) (s/|s|)
- trace of stresses

with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    sig[6]   (i)   stresses  
\param  DOUBLE   *dev      (o)   norm of deviatoric stresses |s|
\param  DOUBLE    n[6]     (o)   gradient n of deviatoric stresses (s/|s|)
\param  DOUBLE   *I1       (o)   trace of stresses
\param  DOUBLE    gamma2   (i)   fitting parameter
\param  DOUBLE    gamma3   (i)   fitting parameter
\param  DOUBLE    ftm      (i)   tensile strength of concrete             
\param  DOUBLE    fcm      (i)   compressive strength of concrete             
\param  DOUBLE    beta[2]  (i)   beta values
\param  DOUBLE    alpha[3] (i)   alpha values
\param  DOUBLE    k_bar[3] (o)   k_bar = beta * sigym
\param  DOUBLE    kappa_t  (i)   internal variable for evolution in tension
\param  DOUBLE    kappa_tu (i)   value for evolution in tension
\param  DOUBLE    kappa_c  (i)   internal variable for evolution in compression
\param  DOUBLE    kappa_e  (i)   value for evolution in compression
\param  DOUBLE    kappa_cu (i)   value for evolution in compression
\param  DOUBLE   *L        (o)   Midpoint of spherical cap
\param  DOUBLE   *R        (o)   Radius of spherical cap
\param  DOUBLE   *I1_0     (o)   position on hyd. axis for changeover from 2nd DP-cone to CAP

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_epc_main()     [mat_pl_epc.c]

*----------------------------------------------------------------------*/
void mat_pl_epc_preval(DOUBLE    sig[6],
                       DOUBLE   *dev,
                       DOUBLE    n[6],
                       DOUBLE   *I1,
                       DOUBLE    gamma2,  
                       DOUBLE    gamma3, 
                       DOUBLE    ftm,    
                       DOUBLE    fcm,    
                       DOUBLE    beta[2],
                       DOUBLE    alpha[3],
                       DOUBLE    k_bar[3],
                       DOUBLE    kappa_t, 
                       DOUBLE    kappa_tu,
                       DOUBLE    kappa_c, 
                       DOUBLE    kappa_e, 
                       DOUBLE    kappa_cu,
                       DOUBLE   *L,       
                       DOUBLE   *R,       
                       DOUBLE   *I1_0)
{
/*----------------------------------------------------------------------*/
DOUBLE q13;
DOUBLE sigym[2];
DOUBLE xsi33;    /*helping variable for calculation of dev and n if |s| < EPS8*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_epc_preval");
#endif
/*----------------------------------------------------------------------*/
q13   = 1./3.;
xsi33 = sig[2];
/*-------- check for numerical problems  --> |s| = ZERO ----------------*/
if (FABS(sig[0]-sig[1]) < EPS8  &&    
    FABS(sig[0]-sig[2]) < EPS8  &&     /*sig[0]=sig[1]=sig[2]*/
    FABS(sig[1]-sig[2]) < EPS8)
{
   /*-- modify the stress vector slightly to avoid numerical problems --*/
   if(FABS(sig[2]) < EPS8)  xsi33 = EPS5;
   else                     xsi33 = 1.001*sig[2];
}    

/*-------- |s| ---------------------------------------------------------*/
*dev = sqrt(2./3.*(sig[0]*sig[0] + sig[1]*sig[1] + xsi33 *xsi33 -
                   sig[0]*sig[1] - sig[0]*xsi33  - sig[1]*xsi33 )
             + 2.*(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]));

/*-------- components of gradient n ------------------------------------*/
n[0]= q13 * (2. * sig[0] - sig[1] - xsi33 ) / *dev;
n[1]= q13 * (2. * sig[1] - sig[0] - xsi33 ) / *dev;
n[2]= q13 * (2. * xsi33  - sig[0] - sig[1]) / *dev;
n[3]= sig[3] / *dev;
n[4]= sig[4] / *dev;
n[5]= sig[5] / *dev;

/*-------- I1 = trace(sig) ---------------------------------------------*/
*I1 = sig[0] + sig[1] + sig[2];

/*-------- values for hardening/softening ------------------------------*/
/*tension*/
sigym[0] = ftm * exp(- kappa_t / kappa_tu);

/*compression*/
if (kappa_c < kappa_e) 
{
   sigym[1] = fcm * (gamma3 + 2.*(1.-gamma3)*kappa_c/(kappa_e) - (1.-gamma3)*(kappa_c/kappa_e)*(kappa_c/kappa_e));
}
else if (kappa_e <= kappa_c && kappa_c < kappa_cu)
{
   sigym[1] = fcm * (1. - ((kappa_c - kappa_e)/(kappa_cu - kappa_e))*((kappa_c - kappa_e)/(kappa_cu - kappa_e)));
}
else sigym[1] = 0.03*fcm;

if (sigym[1] < 0.03*fcm) sigym[1] = 0.03*fcm;

/*------------ k_bar values --------------------------------------------*/
k_bar[0] = beta[0] * sigym[0];
k_bar[1] = beta[1] * sigym[1];
k_bar[2] = -1./(3.*alpha[0]*alpha[0]) * k_bar[0];

/*------------ L and R for spherical cap -------------------------------*/
*L = -(sqrt(54.)*alpha[1] + 2.)*gamma2*sigym[1];
*R = sqrt(2./3. + 6.*alpha[1]*alpha[1])*gamma2*sigym[1];

/*----- calculate I1_0 for updated case --------------------------------*/
*I1_0 = -2.*gamma2*sigym[1];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_epc_preval */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/

