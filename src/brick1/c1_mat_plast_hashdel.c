/*!----------------------------------------------------------------------
\file
\brief contains the routine ' ' which calclate displacement
       derivatives for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief routines for hashin delamination 

<pre>                                                              al 06/02
This routine ... for an 3D-hex-element.

</pre>
\param **sigy  DOUBLE  (i)   input value

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1rad6(double *sig,
            double *sigy,
            double   a[6][6],     
            double   c_orth[6][6],     
            double   sig13,
            double   sig23,
            double   sig33,
            double   sinf33,
            double   gc,
            double   gamma,
            double   *kappa,
            double   *dlam, /* increment of plastic multiplier         */
            double   *zstrich,
            double   *mu,
            double   c1hgt,     /* element thickness ?!         */
            double   c1layhgt,  /* element layer thickness ?!   */
            double *hard,    /* hardening modulus                       */
            double *hard2);   /* hardening modulus                       */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void c1elpl6(
             double c_orth[6][6],
             double a[6][6],     
             double *sig,
             double *hard,  
             double dlam,    /* increment of plastic multiplier         */
             double sigy,
             double sinf33, 
             int   *imod,
             int  ivisco);    
/*----------------------------------------------------------------------*
 |    topic : hashin delamination crit                    al    9/01    |
 |                                                                      |
 |                           f = sqrt(s*p*s) - sigy                     |
 *----------------------------------------------------------------------*/

double c1yield_hashin_delam(
               double   *sig ,
               double   a[6][6],
               double    sigy ,
               double   *dsas 
             );  
/*----------------------------------------------------------------------*
 |    topic : form viscoplastic regularization            al    9/01    |
 *----------------------------------------------------------------------*/
/*    sig      stress to be regularized
      kappa    equiv strain to be regularized
      c        orthotropic materialtensor to be regularized */
void c1visco6(
               double   delta_t,
               double   eta_start,
               double   *eta_n ,
               double   *sig ,
               double   *sigl ,
               double   *kappa ,
               double   kappal ,
               double   c_orth[6][6],
               double   c_elas[6][6],
               double   *delsig 
             ) ; 
/*----------------------------------------------------------------------*
 |    topic : calculate the orthotropic material tensor c al    9/01    |
 *----------------------------------------------------------------------*/
void c1corth(
                   double   sinf33,
                   double   sigy,
                   int      ivisco,
                   double   c_orth[6][6],
                   double   emod1 ,
                   double   emod2 ,
                   double   emod3 ,
                   double   gmod12,
                   double   gmod13,
                   double   gmod23,
                   double   xnue12,
                   double   xnue13,
                   double   xnue23
                  ) ; 
/*----------------------------------------------------------------------*
 |    topic : calculate some needed values                al    9/01    |
 |                                                                      |
 |    hashin delamination                                               |
 *----------------------------------------------------------------------*/
void c1prevalue_hd(
                   double  *sigy_b,
                   double   sig13,
                   double   sig23,
                   double   sig33,
                   double   sinf33,
                   double   gc,
                   double   gamma,
                   double   kappa,
                   double   dlam,
                   double   a[6][6],
                   double  *zstrich_b,
                   double  *mu_b,
                   double   c1hgt,    /* element thickness ?!         */
                   double   c1layhgt  /* element layer thickness ?!   */
                  );  
/*----------------------------------------------------------------------*
 | - delamination hashin plasticity based - 3D                   al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_hashdel
                        (
                        MATERIAL  *mat, /* material data                */
                        ELEMENT   *ele, /* actual element               */
                        double **bop,   /* derivative operator          */
                        int       ip,   /* integration point            */
                        double *stress, /*ele stress (-resultant) vector*/      
                        double **c,     /* material matrix              */
                        double  *disd,  /* displacement derivatives     */
                        double g[6][6], /* transformation matrix        */
                        double gi[6][6],/* inverse of g                 */
                        int istore,     /* controls storing of stresses */
                        int newval)     /* controls eval. of stresses   */
{
/**/      
/*----------------------------------------------------------------------*
 |                                                                      |
 |    6-------18-------2           6----------------2                   |
 |    |\               |\          |\               |\                  |
 |    | \              | \         | \        S     | \                 |
 |    |  13            |  9        |  \       |     |  \                |
 |    |   \            |   \       |   \      |     |   \               |
 |   14    \           10   \      |    \  \  |     10   \              |
 |    |     5-------17-------1     |     5----------------1             |
 |    |     |          |     |     |     |   \|     |     |             |
 |    |     |          |     |     | T---|----o--------   |             |
 |    |     |          |     |     |     |    |\    |     |             |
 |    7-----|-19-------3     |     7-----|----|-\---3     |             |
 |     \    12          \    8      \    |    |  \   \    |             |
 |      \   |            \   |       \   |    |   R   \   |             |
 |       15 |             11 |        \  |    |        \  |             |
 |        \ |              \ |         \ |              \ |             |
 |         \|               \|          \|               \|             |
 |          4-------16-------0           4----------------0             |
 |                                                                      |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int i,j,k;
int yip;
int isoft;
int iupd;
double e1, e2, e3, a1, b1, c1, sum, epstn, ft;
double sig[6];
double eps[6];
double strain[6];
double delsig[6];
double deleps[6];
double tau[6];
double tol = 1.0E-10;
double dlam;
/* new: */
int    imod;
int    ivisco;
double sinf33;
double sigl[6];
double epsl[6];
double kappa, kappal;
double gamma;
double delta_t;
double eta_start;
double eta_n;
double sigy;
double mu;
double sig13, sig23, sig33, gc, zstrich;
double a[6][6];
double c_orth[6][6];
double c_help[6][6];
double c_elas[6][6];
double c_plas[6][6];
double eps_help[6];
double eps_sort[6];
double sig_orth[6];
double eps_orth[6];
double str[6];
double dsas;


double emod1 ,emod2 ,emod3 ;
double gmod12,gmod13,gmod23;
double xnue12,xnue13,xnue23;
/*--------------------------------------------------------*/
double c1hgt;    /* element thickness ?!         */
/*--------------------------------------------------------*/
double c1layhgt; /* element layer thickness ?!   */
/*--------------------------------------------------------*/
int is=1;         /* ????????????????????????? */
int ic=5;         /* ????????????????????????? */


double yld, sm, sx, sy, sz, sxy, sxz, syz, sig2, hard, hard2;
double expo  = 0.;
double alpha = 0.;
#ifdef DEBUG 
dstrc_enter("c1_mat_plast_hashdel");
#endif
/*------------------------------------------------------- initialize ---*/          
  iupd=0;
  dlam = 0.;
  sinf33 = mat->m.pl_hash->sinf33;
/*------------------------------------- get actual strains -> strain ---*/
 /*???????????????????????*/
  c1_eps (disd,eps,1);
/*--- orthotropic material tensor ---*/    
  c1_mat_elorth(mat->m.el_orth->emod1,
                mat->m.el_orth->emod2,
                mat->m.el_orth->emod3,
                mat->m.el_orth->xnue23,
                mat->m.el_orth->xnue13,
                mat->m.el_orth->xnue12,
                mat->m.el_orth->gmod12,
                mat->m.el_orth->gmod23,
                mat->m.el_orth->gmod13,
                c);
  /*c1gld (c,g);/* transform local to global material matrix */
/*--- stresses ---*/    
  c1mefm(strain, c, stress);
/*---------- get current state and variables of integration point ip ---*/         
  ivisco = mat->m.pl_hash->ivisco;

  for (i=0; i<6; i++)
  {
    sigl[i] = ele->e.c1->elewa[0].ipwa[ip].sig[i];
    epsl[i] = ele->e.c1->elewa[0].ipwa[ip].eps[i];
  }
  kappa  = ele->e.c1->elewa[0].ipwa[ip].kappa;
  imod   = ele->e.c1->elewa[0].ipwa[ip].imod;
  yip    = ele->e.c1->elewa[0].ipwa[ip].yip;
  kappal = kappa;
/*----------------------------------------------------------------------*/
/*---------------------------------- calculate some needed prevalues ---*/
  sig13 = mat->m.pl_hash->s13   ;
  sig23 = mat->m.pl_hash->s23   ;
  sig33 = mat->m.pl_hash->s33   ;
  gc    = mat->m.pl_hash->gc    ;
  gamma = mat->m.pl_hash->gamma ;
  emod1 = mat->m.pl_hash->emod1  ;
  emod2 = mat->m.pl_hash->emod2  ;
  emod3 = mat->m.pl_hash->emod3  ;
  xnue23= mat->m.pl_hash->xnue23 ;
  xnue13= mat->m.pl_hash->xnue13 ;
  xnue12= mat->m.pl_hash->xnue12 ;
  gmod12= mat->m.pl_hash->gmod12 ;
  gmod23= mat->m.pl_hash->gmod23 ;
  gmod13= mat->m.pl_hash->gmod13 ;
                          
  c1hgt    = mat->m.pl_hash->c1hgt   ;
  c1layhgt = mat->m.pl_hash->c1layhgt;
  
  
  c1prevalue_hd (
                 &sigy,
                 sig13,
                 sig23,
                 sig33,
                 sinf33,
                 gc,
                 gamma,
                 kappa,
                 dlam,
                 a,
                 &zstrich,
                 &mu,
                 c1hgt,           /* hoehe des elements !?! */
                 c1layhgt         /* hoehe des elements !?! */
                 );

/*-------------------------- determine orthotropic material tensor c ---*/

  c1corth(
           sinf33,
           sigy,
           ivisco,
           c_orth,
           emod1 ,
           emod2 ,
           emod3 ,
           gmod12,
           gmod13,
           gmod23,
           xnue12,
           xnue13,
           xnue23
                  );  
  for (i=0; i<6; i++) for (j=0; j<6; j++) c_help[i][j] = c_orth[i][j]; 
  for (i=0; i<6; i++) for (j=0; j<6; j++) c_elas[i][j] = c_orth[i][j]; 


/*-------------------------------- transform of epsilon to kartesian ---*/

  /*for (i=0; i<6; i++) eps_help[i] = eps[i]; 

  /*s9trans(eps_help,matrix,'e','cuor',gkonr,layer,nummat,
     &             numlay,lay) 
/*----------------------------------------------- sorting of epsilon ---*/
/*
                      ***_SORT=(***11,***22,***33,***12,***23,***13)
                           ***=(***11,***12,***22,***13,***23,***33) */

  for (i=0; i<6; i++) eps_help[i] = eps[i]; 
  /*
      eps_orth[0] = eps_help[0];
      eps_orth[1] = eps_help[2];
      eps_orth[2] = eps_help[5];
      eps_orth[3] = eps_help[1];
      eps_orth[4] = eps_help[4];
      eps_orth[5] = eps_help[3];
  */
  for (i=0; i<6; i++) eps_orth[i] = eps[i]; 

/*------------------------ transform local to global material matrix ---*/
/*  c1gld (d,g); */
/*----------------------------------------------------------------------*/
/*  if(newval==1)
  {
    for (i=0; i<6; i++)  stress[i] = sig[i];
    goto end;
  }
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++)
  {
    stress[i] = 0.0;
    tau[i]    = 0.0;
  }
/*-----------------------------------------------------------------------|
|     yip > 0  stresses are available from last update                   |
|         = 1  e l a s t i c                                             |
|         = 2  p l a s t i c                                             |
|     update flag must set to store change of parameter yip              |
|     no changes have been made on stress state                          |
|-----------------------------------------------------------------------*/
  if(yip>0)
  {
    if(yip==1)
    {
      yip=-yip;
    }
    else
    {
/*================================== praedictor step ===== yip = 2 =====*/
      if (ic>0) 
      {
        dlam=0.;
/*------------------------------------ elastoplast. materialtangente ---*/

        hard = -zstrich;

        c1elpl6(
                c_orth,
                a,     
                sig,
                &hard,  
                dlam,    /* increment of plastic multiplier         */
                sigy,
                sinf33, 
                &imod,
               ivisco) ;   
/*------------------------------------ viscoplastisic regularization ---*/


        if (ivisco>0 && imod!=1) 
        {
            delta_t   = mat->m.pl_hash->deltat;
            eta_start = mat->m.pl_hash->eta_i ; 
            for (i=0; i<6; i++) delsig[i] = 0.;

            c1visco6( delta_t,
                      eta_start,
                      &eta_n ,
                      sig ,
                      sigl ,
                      &kappa ,
                       kappal ,
                      c_orth,
                      c_elas,
                      delsig 
                    ) ; 
        }
/*c!--------------------------------transformation von c_orth nach krumml.
/*------------------------ transform local to global material matrix ---*/
  for (i=0; i<6; i++) for (j=0; j<6; j++) c_help[i][j] = c_orth[i][j];
  /*c1gld (c_orth,g);!!!!!!!!!!!!!!!!!!!1111*/ 
  
/*
               call mxcr8  (c_orth,6,6,c_help)
               call s9trans(vector,c_orth,'c','orcu',gkonr,layer,nummat,
     &                      numlay,lay) 
               call mxcr8  (c_orth,6,6,c)
*/        
      }
      yip=-yip;
    }
    for (i=0; i<6; i++) sig[i] = sigl[i];
    if(is>0) for (i=0; i<6; i++) sig_orth[i] = sigl[i];

/*------------------------------------ sorting of sigma according to ---*/
/*               shell9:        str=(str11,str12,str22,str13,str23,str33)
                 brick1: sig_orth=(sig11,sig22,sig33,sig12,sig23,sig13) */

    /* str(1) = sig_orth(1)
     str(2) = sig_orth(4)
     str(3) = sig_orth(2)
     str(4) = sig_orth(6)
     str(5) = sig_orth(5)
     str(6) = sig_orth(3)*/
    for (i=0; i<6; i++) stress[i] = sig_orth[i];

    c1trss2global (stress,g);

/* completed softening on curvilinear components!!! */
    if (imod==1) 
    {
       /*call mxir8 (str,6,1,1.0d-6)
       call mxir8 (c,6,6,1.0d-12)*/
       for (i=0; i<6; i++) stress[i] = 1.0E-6; 
       for (i=0; i<6; i++) for (j=0; j<6; j++)  c[i][j] = 1.0E-12;
       for (i=0; i<6; i++) c[i][i] = 1.0E-6; 
    }
    iupd=1;
    goto end;
  }
/*-----------------------------------------------------------------------|
|   1. calculate incremental strains     deleps                          |
|   2. calculate stress increment assuming elastic behaviour             |
|   3. calculate total stress                                            |
|   4. check stress deviator against current yield surface               |
|-----------------------------------------------------------------------*/

/*---------------------------------------------------------- STEP 1. ---*/
  for (i=0; i<6; i++) deleps[i] = eps_orth[i] - epsl[i];

/*---------------------------------------------------------- STEP 2. ---*/
  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += c_orth[i][j]*deleps[j];

/*---------------------------------------------------------- STEP 3. ---*/
  for (i=0; i<6; i++) sig[i] = sigl[i] + delsig[i];

/*---------------------------------------------------------- STEP 4. ---*/
/*------------------------------------------------ verify plasticity ---*/
  ft=c1yield_hashin_delam(sig,a,sigy,&dsas);
/* IF COMPLETE SOFTENING IS ALREADY REACHED, CONSEQUENTLY NO NEED TO
   CALCULATE SIG AND C VIA THE RADIAL RETURN PROCEDURE
   IMOD IS SET FOR THE FIRST TIME IN S9ELPL6 IF SIGY HAS REACHED SINF33*/
  if (imod==1) 
  {
    yip=2;

    for (i=0; i<6; i++) sig[i] = 1.0E-6; 
    for (i=0; i<6; i++) stress[i] = 1.0E-6; 
    for (i=0; i<6; i++) for (j=0; j<6; j++)  c[i][j] = 1.0E-12;
    for (i=0; i<6; i++) c[i][i] = 1.0E-6; 

    for (i=0; i<6; i++) epsl[i]     = eps_orth[i]; 
    for (i=0; i<6; i++) eps_sort[i] = eps_orth[i]; 
    goto end;
  }

/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<=tol) 
  {
    yip = 1;
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    yip = 2;
    
/*-- plastisches multiplikatorinkrement ---*/
    c1rad6(sig,
           &sigy,
           a,     
           c_orth,     
           sig13,
           sig23,
           sig33,
           sinf33,
           gc,
           gamma,
           &kappa,
           &dlam, /* increment of plastic multiplier         */
           &zstrich,
           &mu,
           c1hgt,     /* element thickness ?!         */
           c1layhgt,  /* element layer thickness ?!   */
           &hard,    /* hardening modulus                       */
           &hard2);   /* hardening modulus                       */
    
    if (ic>0) 
    {
      c1elpl6(
                c_orth,
                a,     
                sig,
                &hard2,  
                dlam,    /* increment of plastic multiplier         */
                sigy,
                sinf33, 
                &imod,
                ivisco) ;   
    
      for (i=0; i<6; i++) for (j=0; j<6; j++) c_plas[i][j] = c_orth[i][j];
      for (i=0; i<6; i++) for (j=0; j<6; j++) c[i][j] = c_orth[i][j];
      /*c1gld (c,g);/* transform local to global material matrix */
    } 
    
    
    
    /*transform stresses to local coordinate system */
    /* c1trss2local (tau,gi);
    /* c1matp1 (ym, fhard, uniax, pv, sig2,stress,epstn,dlam,d); */
    /*c1gld (c,g);/* transform local to global material matrix */
    /*transform stresses to global coordinate system */
    /*c1trss2global (tau,g);*/
    /**/
  }
/*----------------------------------------------------------------------*/
  /*transform stresses to global coordinate system */
  /*c1trss2global (tau,g);*/

  for (i=0; i<6; i++) sig_orth[i] = sig[i];
  for (i=0; i<6; i++)   stress[i] = sig[i];

/*-----------------------------------------------completed softening ---*/
/* completed softening on curvilinear components!!! */
  if (imod==1)
  {
    for (i=0; i<6; i++) stress[i] = 1.0E-6;

    for (i=0; i<6; i++)
    {
      for (j=0; j<6; j++)
      {
        c_orth[i][j] = 1.0E-12;
        if (i==j) c_orth[i][j] = 1.0E-6;
      }
    }
  }

  for (i=0; i<6; i++) epsl[i] = eps_orth[i];
  for (i=0; i<6; i++)  eps[i] = eps_orth[i];









/*----------------------------------------------------------------------*/
end:
/*------------------------------- update wa and visco regularization ---*/
  if(istore==1 || iupd==1)
  {
    for (i=0; i<6; i++)
    {
      ele->e.c1->elewa[0].ipwa[ip].sig[i] = sig[i];
      ele->e.c1->elewa[0].ipwa[ip].eps[i] = epsl[i];
    }
    ele->e.c1->elewa[0].ipwa[ip].kappa  = kappa;
    ele->e.c1->elewa[0].ipwa[ip].imod   = imod ;
    ele->e.c1->elewa[0].ipwa[ip].yip    = yip  ;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_mat_plast_hashdel */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1rad6(double *sig,
            double *sigy,
            double   a[6][6],     
            double   c_orth[6][6],     
            double   sig13,
            double   sig23,
            double   sig33,
            double   sinf33,
            double   gc,
            double   gamma,
            double   *kappa,
            double   *dlam, /* increment of plastic multiplier         */
            double   *zstrich,
            double   *mu,
            double   c1hgt,     /* element thickness ?!         */
            double   c1layhgt,  /* element layer thickness ?!   */
            double *hard,    /* hardening modulus                       */
            double *hard2)    /* hardening modulus                       */
{
/*----------------------------------------------------------------------*/
int i, j, cc, irc;
int ione = 1;
int isix = 6;
double dzero = 0.;
int iter;
double sig_old, fac, fact, f, dsas, npn, dfdl;
double signew[6], as[6], n[6], pn[6];
double  ma66h[36];
double fact_a[36];
double     ci[36], pci[36], p[36], pi[36];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1rad6");
#endif
/*-------------------------------------------- initialized variables ---*/
  *dlam = 0.;
/*----------------------------------------- check compression stress ---*/
  sig_old = 0.;
  if (sig[2]<0.) 
  {
     sig_old = sig[2];
     sig[2]  = 0.;
  }
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
  iter=0;
/* for non-linear hardening rule */
L500:
  ++iter;
/*--------------------------------------update of internal variables ---*/
  c1prevalue_hd (
                 sigy,
                 sig13,
                 sig23,
                 sig33,
                 sinf33,
                 gc,
                 gamma,
                 (*kappa),
                 (*dlam),
                 a,
                 zstrich,
                 mu,
                 c1hgt,           /* hoehe des elements !?! */
                 c1layhgt         /* hoehe des elements !?! */
                 );

/*----------------------------------------------- update of stresses ---*/
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = c_orth[i][j];  
  c1inv6(ma66h,ci,&irc);
  fact=(*dlam)/(*sigy);
  
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) fact_a[cc++] =  fact * a[i][j];  
  
  for (i=0; i<36; i++)      p[i] = ci[i] + fact_a[i];
  c1inv6(p,pi,&irc);
  
  c1ab (pi,ci,pci,&isix,&isix,&isix,&dzero);
  c1ab (pci,sig,signew,&isix,&isix,&ione,&dzero);

/*---------------------------------------------------------------- f ---*/
  f=c1yield_hashin_delam(signew,a,(*sigy),&dsas);
/*--------------------------------------------------- calculate dfdl ---*/
  fac = dsas/(*sigy)*(1.-(*dlam)*(-(*zstrich))/(*sigy));
  /* a is an array !!! */
  /* c1ab(a,signew,as,&isix,&isix,&ione,&dzero); */
  for (i=0; i<6; i++) n[i] =  1./dsas * as[i];

  c1ab (pi,n,pn,&isix,&isix,&ione,&dzero);
  mxmatb (n,pn,&npn,&isix,&ione,&ione,&dzero);
/* */     
  dfdl = -(fac*npn-(*zstrich));

/*--------------------------------------------------------- new dlam ---*/
/* new plastic increment */
/* dlam|(k+1) = dlam|k - phi/(dphi/dlambda)|k */

  (*dlam) = (*dlam) - f/dfdl;

/* check convergence */
  if (fabs(f) > 1e-8) {
     if (iter > 200) {
              dserror("local iteration exceeds limit");
          }
          goto L500; }
/*--------------------------------------update of internal variables ---*/
  c1prevalue_hd (
                 sigy,
                 sig13,
                 sig23,
                 sig33,
                 sinf33,
                 gc,
                 gamma,
                 (*kappa),
                 (*dlam),
                 a,
                 zstrich,
                 mu,
                 c1hgt,           /* hoehe des elements !?! */
                 c1layhgt         /* hoehe des elements !?! */
                 );
/*----------------------------------------------------------------------*/

  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = c_orth[i][j];  
  c1inv6(ma66h,ci,&irc);
  fact=(*dlam)/(*sigy);
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) fact_a[cc++] =  fact * a[i][j];  
  for (i=0; i<36; i++)      p[i] = ci[i] + fact_a[i];
  c1inv6(p,pi,&irc);
  
  c1ab (pi,ci,pci,&isix,&isix,&isix,&dzero);
  c1ab (pci,sig,signew,&isix,&isix,&ione,&dzero);

  for (i=0; i<6; i++) sig[i] = signew[i];  
  if (sig_old<0.) sig[2] = sig_old;
/*---------------------------------------------------------------- f ---*/
  f=c1yield_hashin_delam(sig,a,(*sigy),&dsas);
/*-------------------------------------------------- update of kappa ---*/
  (*kappa) = (*kappa) + (*dlam);
/*------------------------------------ calculate hardening parameter ---*/
  (*hard)  = (-(*zstrich))/(1.-(*dlam)*(-(*zstrich))/(*sigy));
  (*hard2) = -(*zstrich);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1rad6 */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void c1elpl6(
             double c[6][6],
             double a[6][6],     
             double *sig,
             double *hard,  
             double dlam,    /* increment of plastic multiplier         */
             double sigy,
             double sinf33, 
             int   *imod,
             int  ivisco)    
{
/*----------------------------------------------------------------------*/
int i, j, k, l, cc, irc;
int isix=6, ione=1;
double dzero = 0.;
double fact, sas, dsas;
double dum;
double sig_old, ntpn;
double tol = 1.0E-6;
double pn[6], pn1[6];
double pia[6][6];
double gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};

double    fact_a[36], ma66h[36], pci[36], ci[36], p[36], pi[36];
double    as[6], n[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1elpl6");
#endif
/*----------------------------------------- compression stress sig33 ---*/                                                                     
  sig_old = 0.;
  
  if (sig[2]<0.) 
  {
     sig_old = sig[2];
     sig[2] = 0. ;
  }
/*------ completed softening? set remainig stiffness values and imod ---*/

  if (abs(sigy-sinf33)<tol && (ivisco==0 || ivisco==2)) 
  {
    for (i=0; i<6; i++)
    {
      for (j=0; j<6; j++)
      {
       c[i][j] = 1.0E-12;
       if (i==j) c[i][j] = 1.0E-6;
      }
      sig [i] = 1.0E-12;
    }    
    (*hard)    = 1.0E-12;
    (*imod) = 1;
    goto end;
  }
/*-------------------------------------------------calculate matrix p---*/
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = c[i][j];  
  c1inv6(ma66h,ci,&irc);

  fact=dlam/sigy;
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) fact_a[cc++] =  fact * a[i][j];  
  for (i=0; i<36; i++)      p[i] = ci[i] + fact_a[i];
  c1inv6(p,pi,&irc);
/*
  c1ab (pi,ci,pci,&isix,&isix,&isix,0);
  c1ab (pci,sig,signew,&isix,&isix,&ione,&dzero);
/*-------------------------------------------------calculate vector n---*/
  /* a is an array !!! */
  /*c1ab (a,sig,as,&isix,&isix,&ione,&dzero);*/
  sas = 0.;
  for (i=0; i<6; i++) sas += sig[i] * as[i];
  dsas = sqrt(sas);
  for (i=0; i<6; i++) n[i] = as[i] * 1./dsas;
/*------------------------------------------------------calculate p*n---*/
  c1ab (pi,n,pn,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) pn1[i] = pn[i];
  mxmatb (n,pn,&ntpn,&isix,&ione,&ione,&dzero);
/*------------------------------------------------------ calculate c ---*/
  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) pia[i][j] = pi[cc++];
  
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
        c[i][j] = pia[i][j] - pn[i] * pn1[j] / (ntpn + (*hard));
    }
  }

  if (sig_old<0.) sig[2] = sig_old;
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif  
return; 
} /* end of c1elpl6 */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
/*----------------------------------------------------------------------*
 |    topic : calculate some needed values                al    9/01    |
 |                                                                      |
 |    hashin delamination                                               |
 *----------------------------------------------------------------------*/
void c1prevalue_hd(
                   double  *sigy_b,
                   double   sig13,
                   double   sig23,
                   double   sig33,
                   double   sinf33,
                   double   gc,
                   double   gamma,
                   double   kappa,
                   double   dlam,
                   double   a[6][6],
                   double  *zstrich_b,
                   double  *mu_b,
                   double   c1hgt,    /* element thickness ?!         */
                   double   c1layhgt  /* element layer thickness ?!   */
                  )  
{
/*----------------------------------------------------------------------*/
int i, j;
double mu        = 0.;
double z_strich  = 0.;
double alpha     = 0.;
double sigy      = 0.;
double tp        = 0.;
double wert      = 0.;
double dum       = 0.;
double tol       = -60.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1prevalue_hd");
#endif
/*-------------------------------- berechnung der vergleichsspannung ---*/                  
  mu = sig33*c1layhgt*c1hgt/(100.*2.*gc);      
/*     mu = sig33/(2.*gc)  */    
  tp = 1./mu-1./(mu*gamma)-sinf33/(mu*sig33);
  
  alpha = kappa + dlam;
  if (alpha>=0 && alpha<=tp) 
  {
    sigy = sig33*(1.-mu*alpha);
    z_strich = sig33*mu;
  } else if (alpha>tp) {
    wert = gamma-1.- log(gamma)-alpha*gamma*mu
                   + log(sig33/(sig33-sinf33))-gamma*sinf33/sig33;
    if (wert<tol) 
    {
      sigy     = sinf33;
      z_strich = 1.0E-12;
    } else {
      sigy = (sig33-sinf33)*exp(wert)+sinf33;
      z_strich = (sig33-sinf33)*gamma*mu*exp(wert);
    }
  }
/*----------------------------------- aufstellen der kopplungsmatrix ---*/                  
  for (i=0; i<6; i++) for (j=0; j<6; j++) a[i][j] = 0.; 

  a[2][2] = 1./(sig33*sig33)*(sig33*sig33);
  a[4][4] = 1./(sig23*sig23)*(sig33*sig33);
  a[5][5] = 1./(sig13*sig13)*(sig33*sig33);
/*----------------------------------------------------------------------*/
  *sigy_b     = sigy;
  *mu_b       = mu;
  *zstrich_b = z_strich;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1prevalue_hd */
/*----------------------------------------------------------------------*
 |    topic : calculate the orthotropic material tensor c al    9/01    |
 *----------------------------------------------------------------------*/
void c1corth(
                   double   sinf33,
                   double   sigy,
                   int      ivisco,
                   double   c_orth[6][6],
                   double   emod1 ,
                   double   emod2 ,
                   double   emod3 ,
                   double   gmod12,
                   double   gmod13,
                   double   gmod23,
                   double   xnue12,
                   double   xnue13,
                   double   xnue23
                  )  
{
/*----------------------------------------------------------------------*/
int i, j;
double tol=1.0E-6;
double xnue31        = 0.;
double xnue32        = 0.;
double xnue21        = 0.;
double emod          = 0.;
double delta         = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1corth");
#endif
/*-----------------------------------------------completed softening ---*/
  if (fabs(sigy-sinf33)<tol && (ivisco==0||ivisco==2))
  {
    for (i=0; i<6; i++)
    {
      for (j=0; j<6; j++)
      {
        c_orth[i][j] = 1.0E-12;
        if (i==j) c_orth[i][j] = 1.0E-6;
      }
    }
    goto end;
  }
/*---------------------------------- 3d-materialtensor (orthotropic) ---*/
  for (i=0; i<6; i++) for (j=0; j<6; j++) c_orth[i][j] = 0.; 

  xnue31=xnue13*emod3/emod1;
  xnue32=xnue23*emod3/emod2;
  xnue21=xnue12*emod2/emod1;
  emod=emod1*emod2*emod3;
  delta=1.-xnue13*xnue31-xnue23*xnue32-xnue12*xnue21;
  delta=(delta-2.*xnue31*xnue23*xnue12)/emod;
  c_orth[0][0]=(1.-xnue23*xnue32)/(emod2*emod3*delta);
  c_orth[1][1]=(1.-xnue13*xnue31)/(emod1*emod3*delta);
  c_orth[2][2]=(1.-xnue12*xnue21)/(emod1*emod2*delta);
  c_orth[1][0]=(xnue12+xnue13*xnue32)/(emod1*emod3*delta);
  c_orth[2][0]=(xnue13+xnue12*xnue23)/(emod1*emod2*delta);
  c_orth[2][1]=(xnue23+xnue21*xnue13)/(emod1*emod2*delta);
  c_orth[3][3]=gmod12;
  c_orth[5][5]=gmod13;
  c_orth[4][4]=gmod23;
  c_orth[0][1]=c_orth[1][0];
  c_orth[0][2]=c_orth[2][0];
  c_orth[1][2]=c_orth[2][1];
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1corth */
/*----------------------------------------------------------------------*
 |    topic : form viscoplastic regularization            al    9/01    |
 *----------------------------------------------------------------------*/
/*    sig      stress to be regularized
      kappa    equiv strain to be regularized
      c        orthotropic materialtensor to be regularized */
void c1visco6(
               double   delta_t,
               double   eta_start,
               double   *eta_n ,
               double   *sig ,
               double   *sigl ,
               double   *kappa ,
               double   kappal ,
               double   c_orth[6][6],
               double   c_elas[6][6],
               double   *delsig 
             )  
{
/*----------------------------------------------------------------------*/
int i, j;
double tol=1.0E-6;
double beta           = 0.;
double delta          = 0.;
double delta_e        = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1visco6");
#endif
/*-------------------------------------------- viscoplastic stresses ---*/
  if ((*eta_n)==0.) (*eta_n) = eta_start;
  delta = delta_t / (*eta_n);
  if (delta>70.) 
  {
    delta_e = 70.;
  } else {
    delta_e = delta;
  }
  beta  = exp(-delta_e);
  for (i=0; i<6; i++) 
  {
    sig[i] = beta*sigl[i] + (1.-beta)*(sig[i]+1./delta*delsig[i]);
  }
/*----------------------------------------------- viscoplastic kappa ---*/
  (*kappa) = beta * kappal + (1.-beta) * (*kappa);
/*---------------------------------------------- viscoplastic c_orth ---*/
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      c_orth[i][j] = (1.-beta)*(1./delta*c_elas[i][j]+c_orth[i][j]);
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1visco6 */
/*----------------------------------------------------------------------*
 |    topic : hashin delamination crit                    al    9/01    |
 |                                                                      |
 |                           f = sqrt(s*p*s) - sigy                     |
 *----------------------------------------------------------------------*/

double c1yield_hashin_delam(
               double   *sig ,
               double   a[6][6],
               double    sigy ,
               double   *dsas 
             )  
{
/*----------------------------------------------------------------------*/
int i, j;
double tol=1.0E-6;
double c1yield_hashin_delam_value  = 0.;
double sas                         = 0.;
double as[6];
double sig_old                     = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1yield_hashin_delam");
#endif
/*-------------------------------------------- viscoplastic stresses ---*/
  sig_old = 0.;
  if (sig[2]<0.) 
  {
     sig_old = sig[2];
     sig[2]  = 0.;
  }      
  /*call s9ab (a,sig,as,6,6,1,0)*/
  for (i=0; i<6; i++) as[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) as[i] += a[i][j]*sig[j];
  /*call s9skpr (sas,sig,as,6)*/
  sas=0.;                                                               
  for (i=0; i<6; i++) sas += sig[i]*as[i];
  
  c1yield_hashin_delam_value = sqrt(sas) - sigy;
  (*dsas) = sqrt(sas);

  if (sig_old<0.) sig[2] = sig_old;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return c1yield_hashin_delam_value;
} /* end of c1yield_hashin_delam */
#endif
/*! @} (documentation module close)*/
  
        
        
        
        
        
        
        
        
