/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_mat_plast_epc' which calculates the 
       constitutive matrix - forces - elastoplastic concrete - 2D 
       (planes stress, plane strain)

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | const. matrix - forces - elastoplastic concrete - 2D     al  9/01    |
 | plane stress, plane strain                                           |
 *----------------------------------------------------------------------*/
void w1_mat_plast_epc(  double dens      ,
                        double emod      ,
                        double vnu       ,
                        double alfat     ,
                        double xsi       ,
                        double sigyc     ,
                        double ftm       ,
                        double fcm       ,
                        double gt        ,
                        double gc        ,
                        double gamma1    ,
                        double gamma2    ,
                        int    nstiff    ,
                        int    maxreb    ,
                        int    *rebar     ,
                        double *reb_area  ,
                        double *reb_ang   ,
                        double *reb_so    ,
                        double *reb_ds    ,
                        double *reb_rgamma,
                        double *reb_dens  ,
                        double *reb_alfat ,
                        double *reb_emod  ,
                        double *reb_rebnue,
                        double *reb_sigy  ,
                        double *reb_hard  , 
                        ELEMENT   *ele,
                        WALL_TYPE wtype,
                        double **bop,
                        double  *gop,
                        double  *alph,
                        double **xjm,
                        int ip,
                        double *stressc,       
                        double **d,
                        int istore,/* controls storing of new stresses to wa */
                        int newval)/* controls evaluation of new stresses    */
{
/*----------------------------------------------------------------------*/
int i,j,k;
int iupd, yipold, yip, yip2;
double q23, betah; 
double e1, e2, e3, a1, b1, c1, sum, epstn, acrs;
double alpha[4];
double hards[4];
double alphac[2];
double hardsc[2];
double sigym[4];
double gmod, com, dfac, cappaet, cappaut, cappae, cappauc;
double sig3, fbd, hydv; 
double epsi[4];
double dfac1;
int    iflag;
double tau3[4];
double devsig[4];
double sm[4];
double dev;
double hyd;
double dfaci;
double epst2;
double devsigt[4],smt[4],dev3,hyd3,facmin,dlam12;
double angle, thick;

double disd[5];
double sig[4];
double ft[4];
double sigc[4];
double sigi[4];
double sigy[3];
double stress[4];
double di[4];
double dnc[2][4];
double gradc[2][4];
double dn[3][4];
double dcom[3][4];
double tauc[4];
double gradi[4];
double grad[3][4];
double eps[4];
double epst[2];
double strain[4];
double delsig[4];
double deleps[4];
double tau[4];
double qn[4];
double tol = 1.0E-10;
double tol2= 1.0E-5 ;
double dlam[2];
double dlamc;
double rad;
double hard = 0.;
double phi = 0.;
double dia = 0.;
double pr  = 0.;
int    isoft = 0;
WALL_TYPE local_wtype;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_epc");
#endif
/*----------------------------------------------------------------------*/
  rad   = atan(1.)/45.;
/*  phi   = phi * rad;          
/*------------ original global elastic matrix for current point -> D ---*/
  /* look at horst's work (pp.30): in case of plane stress ->
     switch to plane strain and condense stress and material 
     tensor at the end */
  local_wtype=plane_strain;
  
  w1_mat_linel(emod, vnu, local_wtype, d);
  
  if(wtype==plane_stress) pr = 0.;
  else                    pr = vnu;
/*--------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,gop,alph,wtype,disd) ;                  
/*------------------------------------- get actual strains -> strain ---*/
  w1_eps (disd,local_wtype,strain);
/*----------------------------- get old values -> sig, eps,epstn,yip ---*/
  for (i=0; i<4; i++)
  {
    sig[i]   = ele->e.w1->elewa[0].ipwa[ip].sig[ i];
    eps[i]   = ele->e.w1->elewa[0].ipwa[ip].eps[ i];
    sigc[i]  = ele->e.w1->elewa[0].ipwa[ip].sigc[i];
    gradi[i] = ele->e.w1->elewa[0].ipwa[ip].grad[i];
  }
  epst[0]  = ele->e.w1->elewa[0].ipwa[ip].dlam[0];
  epst[1]  = ele->e.w1->elewa[0].ipwa[ip].dlam[1];
  yip      = ele->e.w1->elewa[0].ipwa[ip].yip;
/*----------------------------------------------------------------------*/
  if(newval==1)
  {
    for (i=0; i<4; i++)  stressc[i] = sig[i];
    goto end;
  }
/*--------------------------------------------------- initialization ---*/
  q23    = 2./3.; 
  betah  = 1.;
  sum    = 0.;
  iupd   = 0;
  yipold = yip; 

  for (i=0; i<4; i++)
  {
    stress[ i] = 0.;
    stressc[i] = 0.;
    tau[ i]    = 0.;
    tauc[i]    = 0.;
    sigi[i]    = 0.;
    qn[i]      = 0.;
    di[i]      = d[i][3];
  }
  for (i=0; i<3; i++)
  {
    for (j=0; j<4; j++)
    {
      dn[i][j]      = 0.;
      dcom[i][j]    = 0.;
    }
  }
  for (i=0; i<3; i++)
  {
    ft[i]         = 0.;
    sigy[i]       = 0.;
  }
  dlam[0]  = 0.;   
  dlam[1]  = 0.;   
/*-------------------------------------------- average crack spacing ---*/
      if(maxreb>0)
      {
        w1acrs (ele->e.w1->thick,maxreb,sig,&angle,
                reb_area, reb_ang, reb_so, reb_ds, reb_rgamma,
                &thick, ele->e.w1->elewa[0].dia,xjm,&acrs);
      }
      else
      {
        acrs = ele->e.w1->elewa[0].dia;
      }
/*---------------------------------------------- material parameters ---*/
  w1cpreva (epst,
            sigym,alpha,hards,&emod,&gmod,&vnu,&com,
            &fcm,&gc,&ftm,&gt,&gamma1,&gamma2,&dfac,
            &ele->e.w1->elewa[0].dia,&acrs,
            &cappaet,&cappaut,&cappae,&cappauc,&sig3,&fbd,
            d,sig,eps);
  hydv = - 2. * gamma2 * sigym[2];
/*----------------- store material parameters for tension stiffening ---*/
      if (maxreb>0 && nstiff==1) 
      {
         ele->e.w1->elewa[0].ipwa[ip].stcappae  = cappaet;
         ele->e.w1->elewa[0].ipwa[ip].stcappaut = cappaut;
         ele->e.w1->elewa[0].ipwa[ip].stalpha   = angle;
         ele->e.w1->elewa[0].ipwa[ip].stthick   = thick;
         ele->e.w1->elewa[0].ipwa[ip].stfbd     = fbd;  
      }
/*-----------------------------------------------------------------------|
|     yip > 0  stresses are available from last update                   |
|         = 1  e l a s t i c                                             |
|         = 2  p l a s t i c                                             |
|     update flag must set to store change of parameter yip              |
|     no changes have been made on stress state                          |
|-----------------------------------------------------------------------*/
  if (yip>0) 
  {
    for (i=0; i<4; i++)
    {
      stress[ i] = sig[i];
      stressc[i] = sig[i];
      tau[ i]    = sig[i];
    }
    if (yip==1) 
    {/*yip==1*/
      yip=-yip ; 
      if(wtype==plane_stress) 
      { 
        for (i=0; i<4; i++)
        {
          ele->e.w1->elewa[0].ipwa[ip].sigi[i] = sig[i];
          ele->e.w1->elewa[0].ipwa[ip].epsi[i] = eps[i];
          ele->e.w1->elewa[0].ipwa[ip].di[  i] = d[i][3];
        }
        for (i=0; i<4; i++)
          for (j=0; j<4; j++) d[i][j] *= dfac;
        
        w1concep (d);

        for (i=0; i<4; i++)
        {
          strain[i] = eps[i];
        }
      }
/*----------------------------------------------------------------------*/
    }/*yip==1*/
    else
    {/*yip!=1*/
/*-------------------------------------------------- praedictor step ---*/  
      if (yip == 2 || yip == 3) 
      {
        w1pres(tau, devsig, sm, &dev, &hyd);

        w1preds(tau,&alpha[0],devsig,&dev,dn[0],dcom[0],grad[0]);

        for (i=0; i<4; i++)
        {
          dn[0][i]   = gradi[i]-dcom[0][i];
          grad[0][i] = gradi[i]; 
        }
	w1mapl2 (tau,d,&dlam[0],wtype,&alpha[0],&emod,&gmod,
                 &com,&betah,&hards[0],dn[0],grad[0],&dev);
      }
      else if (yip == 4) 
      {
        w1pres(tau, devsig, sm, &dev, &hyd);

        w1preds(tau,&alpha[2],devsig,&dev,dn[2],dcom[2],grad[2]);

        for (i=0; i<4; i++)
        {
          dn[2][i]   = gradi[i]-dcom[2][i];
          grad[2][i] = gradi[i]; 
        }
	w1mapl2 (tau,d,&dlam[1],wtype,&alpha[2],&emod,&gmod,
                 &com,&betah,&hards[2],dn[2],grad[2],&dev);
      }
      else if (yip == 5) 
      {
        w1pres(tau, devsig, sm, &dev, &hyd);
        w1cpreds (sigym,alpha,devsig,hyd,grad[2]);
        
        for (i=0; i<4; i++)
        {
          dn[2][i]   = gradi[i]-dcom[2][i];
          grad[2][i] = gradi[i]; 
        }
          w1mplcap(tau,d,&dlam[1],wtype,alpha,&emod,&vnu,hards,sigym,
                   grad[2],&cappae,&cappauc,&epst[1],&sig3);
      }
/*------------ condensed constitutive tensor in case of plane stress ---*/
      if(wtype==plane_stress) 
      {   
        for (i=0; i<4; i++)
        {
          ele->e.w1->elewa[0].ipwa[ip].sigi[i] = sig[i];
          ele->e.w1->elewa[0].ipwa[ip].epsi[i] = eps[i];
          ele->e.w1->elewa[0].ipwa[ip].di[  i] = d[i][3];
        }
        w1concep (d);
        for (i=0; i<4; i++) strain[i] = eps[i];
      }
/*----------------------------------------------------------------------*/
      yip  = -yip ;
/*----------------------------------------------------------------------*/
    }/*yip!=1*/
/*----------------------------------------------------------------------*/
      iupd = 1;
      goto end;
  }
/*-----------------------------------------------------------------------|
|  (1) get the values from the last iteration step and calculate the     |
|      total strain e33 ( = strain[3] ) in case of plane stress          |
|  (2) calculate incremental strains "deleps"                            |
|  (3) check unloading behaviour for damage                              |
|  (4) calculate stress increment assuming orthotropic elastic behaviour |
|  (5) calculate total stress                                            |
|  (6) check stress deviator against current yield surface               |
|-----------------------------------------------------------------------*/
  if(wtype==plane_stress)                                      /*  (1)  */
  {  
    for (i=0; i<4; i++)
    {
      sigi[i]   = ele->e.w1->elewa[0].ipwa[ip].sigi[i];
      epsi[i]   = ele->e.w1->elewa[0].ipwa[ip].epsi[i];
      di[  i]   = ele->e.w1->elewa[0].ipwa[ip].di[  i];
    }
    if (fabs(di[3]) - 0.0001 < 0.)
    {
      w1iwadi (emod, vnu, di);
    }
    w1de33 (sigi,epsi,di,strain);
  } 
  
  for (i=0; i<4; i++) deleps[i] = strain[i] - eps[i];          /*  (2)  */
  
  for (i=0; i<4; i++) delsig[i] = 0.0;                         /*  (3)  */
  for (i=0; i<4; i++) for (j=0; j<4; j++) delsig[i] += d[i][j]*deleps[j];
  
  dfac1 = 0.;
  for (i=0; i<4; i++) dfac1 += gradi[i]*delsig[i];

  iflag = 0;                                                   /*  (4)  */
  if(dfac1<0. || yip==1 || yip==-1) 
  {
    for (i=0; i<4; i++)
    {
      for (j=0; j<4; j++)
      {
         d[i][j] *= dfac;
      }
    }
    for (i=0; i<4; i++) delsig[i] = 0.0;
    for (i=0; i<4; i++) for (j=0; j<4; j++) delsig[i] += d[i][j]*deleps[j];
    iflag = 1;
  }
  for (i=0; i<4; i++)                                          /*  (5)  */
  {
    tau[ i] = sig[i] + delsig[i];
    tau3[i] = tau[i];
    tauc[i] = tau[i];
  }
                                               /* yield conditions (6)  */
  w1pres(tau,devsig,sm,&dev,&hyd)          ;
  w1yicsr(dev, hyd, sigym[0], alpha[0], &ft[0]);
  w1yicsr(dev, hyd, sigym[1], alpha[1], &ft[1]);
  w1yicsr(dev, hyd, sigym[2], alpha[2], &ft[2]);
  w1yiccap(dev,hyd, sigym   , alpha[3], &ft[3]);

/*------------- state of stress within yield surface - E L A S T I C ---*/

  if ( (ft[0]<=tol)&&(ft[2]<=tol)&&((ft[3]<=tol)||(hyd>=hydv)) ) 
  {   
    yip=1;  
    dlam[0] = 0.;
    dlam[1] = 0.;
        	
    if(wtype==plane_stress)
    { 
      w1consig (d,tau,tauc); 
          for (i=0; i<4; i++) di[i] = d[i][3];
          w1concep (d);
          for (i=0; i<4; i++)
          {
            ele->e.w1->elewa[0].ipwa[ip].sigi[i] = tau[i];
            ele->e.w1->elewa[0].ipwa[ip].epsi[i] = strain[i];
            ele->e.w1->elewa[0].ipwa[ip].di[  i] = di[i];
          }
    }

    for (i=0; i<4; i++)                                         
    {
      stress[ i] = tau[ i];
      stressc[i] = tauc[i];
    }
  }
/*------------ state of stress outside yield surface - P L A S T I C ----|
|                                                                        |
|       transform stresses to local coordinate system                    |
|       evaluate new stresses with stress projection                     |
|       evaluate new material matrix if requested                        |
|-----------------------------------------------------------------------*/
  else 
  {
    if(iflag==1) 
    {
      dfaci = 1. / dfac;
      dfac  = 1.;
      for (i=0; i<4; i++) for (j=0; j<4; j++) d[i][j] *= dfaci;

      for (i=0; i<4; i++) delsig[i] = 0.0;
      for (i=0; i<4; i++) for (j=0;j<4;j++) delsig[i]+=d[i][j]*deleps[j];
      
      for (i=0; i<4; i++)                   
      {
        tau[ i] = sig[i] + delsig[i];
        tau3[i] = tau[i];
        tauc[i] = tau[i];
      }
      w1pres(tau,devsig,sm,&dev,&hyd);
      w1yicsr(dev, hyd, sigym[0], alpha[0], &ft[0]);
      w1yicsr(dev, hyd, sigym[1], alpha[1], &ft[1]);
      w1yicsr(dev, hyd, sigym[2], alpha[2], &ft[2]);
      w1yiccap(dev,hyd, sigym   , alpha[3], &ft[3]);
     }
    
      w1preds(tau,&alpha[0],devsig,&dev,dn[0],dcom[0],grad[0]);
      w1preds(tau,&alpha[1],devsig,&dev,dn[1],dcom[1],grad[1]);
      w1preds(tau,&alpha[2],devsig,&dev,dn[2],dcom[2],grad[2]);

/*-----------------------------------------------------------------------| 
|       tension- / tension-compression-region                            |
|-----------------------------------------------------------------------*/ 
 
    if ((ft[0]>tol && ft[2]<=tol) || ft[1]<0.) 
    {
      yip=2;
      if (ft[1]>=0.) 
      {

        w1cradi (tau,&epst[0],&dlam[0],wtype,yip,
                     &alpha[0],&ft[0],&emod,&gmod,&com,&sigym[0],&hards[0],
                     &sigy[0],dn[0],dcom[0],grad[0],devsig,sm,
                     &fcm,&gc,&ftm,&gt,&gamma1,&gamma2,
                     &ele->e.w1->elewa[0].dia,&acrs); 

	w1mapl2 (tauc,d,&dlam[0],wtype,&alpha[0],&emod,&gmod,
                 &com,&betah,&hards[0],dn[0],grad[0],&dev);
      } 
      else
      {
        /* drucker-pragers apex (inverted cone) */
        w1cradms (tau ,epst,dlam,wtype,yip,
                      alpha,ft,emod,gmod,com,sigym,hards,
                      dn,dcom,devsig,sm,
                      fcm,gc,ftm,gt,gamma1,gamma2,
                      ele->e.w1->elewa[0].dia,acrs); 
        if (dlam[0]>=0.) 
        {
	  w1mapl2 (tauc,d,&dlam[0],wtype,&alpha[0],&emod,&gmod,
                   &com,&betah,&hards[0],dn[0],grad[0],&dev);
        }
        else
        {
	  dlam12 = dlam[0] + dlam[1];
          w1maplg (tauc,d,dlam12,wtype,yip,emod,gmod,com,hards,dnc,grad);
        } 
      } 
    }
/*-----------------------------------------------------------------------| 
|       multisurface-region (compression/tension)                        |
|-----------------------------------------------------------------------*/ 
    else if (ft[2]>tol && ft[0]>tol) 
    {
      yip2=3;


      w1cradms (tau3,epst,dlam,wtype,yip2,
                alpha,ft,emod,gmod,com,sigym,hards,
                dn,dcom,devsig,sm,
                fcm,gc,ftm,gt,gamma1,gamma2,
                ele->e.w1->elewa[0].dia,acrs); 
      w1conver(dlam,alpha,hards,dn,grad,&dlamc,alphac,hardsc,dnc,gradc);
      
      w1pres(tau3,devsigt,smt,&dev3,&hyd3);
      hydv = - 2. * gamma2 * sigym[2];



      if (hyd3>=hydv || (yip!=5&&yip!=-5)) 
      {
        yip=3;
        if (dlam[0] == 0.) yip=4;
        if (dlam[1] == 0.) yip=2;

	if (yip==2) w1mapl2 (tauc,d,&dlam[0],wtype,&alpha[0],&emod,&gmod,
                             &com,&betah,&hards[0],dn[2],grad[0],&dev);
                   
        else if (yip==3) w1maplg (tauc,d,dlamc,wtype,yip,
                         emod,gmod,com,hardsc,dnc,gradc);
        
        else if (yip==4) w1mapl2 (tauc,d,&dlam[1],wtype,&alpha[2],&emod,&gmod,
                                  &com,&betah,&hards[2],dn[2],grad[2],&dev);
        
        for (i=0; i<4; i++) tau[i] = tau3[i];   
      }
      else
      {  
            yip=5;
            epst[1] = epst[1]-dlam[1];
      
      w1radcap (tau ,&epst[1],&dlam[1],wtype,
                    alpha,&emod,&gmod,&com,sigym,hards,
                    grad[2],devsig,sm,&dev,&hyd,&hyd3,
                    &fcm,&gc,&ftm,&gt,&gamma2,&ele->e.w1->elewa[0].dia);

          w1mplcap(tau,d,&dlam[1],wtype,alpha,&emod,&vnu,hards,sigym,
                   grad[2],&cappae,&cappauc,&epst[1],&sig3);
      }
    }
/*-----------------------------------------------------------------------| 
|       compression-region                                               |
|-----------------------------------------------------------------------*/ 
    else if (ft[2]>tol && ft[0]<=tol) 
    {
      /* standard drucker-prager yield function */
      yip2=4;
      epst2 = epst[1];
      w1cradi (tau3,&epst2,&dlam[1],wtype,yip2,
                    &alpha[2],&ft[2],&emod,&gmod,&com,&sigym[2],&hards[2],
                    &sigy[2],dn[2],dcom[2],grad[2],devsig,sm,
                    &fcm,&gc,&ftm,&gt,&gamma1,&gamma2,
                    &ele->e.w1->elewa[0].dia,&acrs); 
      w1pres(tau3,devsigt,smt,&dev3,&hyd3);
      hydv = - 2. * gamma2 * sigy[2];

      /* check cap-region */
      if (hyd3>=hydv || (yip!=5&&yip!=-5)) 
      {
        yip=4;
        w1mapl2 (tauc,d,&dlam[1],wtype,&alpha[2],&emod,&gmod,
                 &com,&betah,&hards[2],dn[2],grad[2],&dev);
        epst[1] = epst2;
        for (i = 0; i < 4; i++) tau[i]  = tau3[i]; 
      }
      else 
      { 
        yip=5;
          w1radcap (tau ,&epst[1],&dlam[1],wtype,
                    alpha,&emod,&gmod,&com,sigym,hards,
                    grad[2],devsig,sm,&dev,&hyd,&hyd3,
                    &fcm,&gc,&ftm,&gt,&gamma2,&ele->e.w1->elewa[0].dia);
          w1mplcap(tau,d,&dlam[1],wtype,alpha,&emod,&vnu,hards,sigym,
                   grad[2],&cappae,&cappauc,&epst[1],&sig3);
      } 
    }
/*-----------------------------------------------------------------------| 
|       cap-region                                                       |
|-----------------------------------------------------------------------*/ 
    else if (ft[3]>tol && ft[0]<=tol && ft[2]<=tol) 
    {
      hydv = - 2. * gamma2 * sigym[2];
        if (hyd>=hydv) 
        {
          yip=1;
          dlam[0] = 0.;
          dlam[1] = 0.;
          for (i=0;i<4;i++) tau[i] = tauc[i];	
        }
        else
        {
          yip=5;
          w1radcap (tau3 ,&epst[1],&dlam[1],wtype,
                    alpha,&emod,&gmod,&com,sigym,hards,
                    grad[2],devsig,sm,&dev,&hyd,&hyd3,
                    &fcm,&gc,&ftm,&gt,&gamma2,&ele->e.w1->elewa[0].dia);
          w1mplcap(tau3,d,&dlam[1],wtype,alpha,&emod,&vnu,hards,sigym,
                   grad[2],&cappae,&cappauc,&epst[1],&sig3);
          for (i=0;i<4;i++) tau[i] = tau3[i];	
        }
    } 
/*----------------------------------------------------------------------*/
        if (yip==2)
        {
          for (i=0;i<4;i++) gradi[i] = grad[0][i];	
        }
        else if (yip==3)
        {
          if (dlam[0]>=dlam[1]) for (i=0;i<4;i++) gradi[i] = grad[0][i];	
          else                  for (i=0;i<4;i++) gradi[i] = grad[2][i];	
        }
        else if (yip>=4)
        {
            for (i=0;i<4;i++) gradi[i] = grad[2][i];	
        }
/*----------------------------------------------------------------------*/
    if(wtype==plane_stress)
    { 
      w1consig (d,tau,tauc); 
          for (i=0; i<4; i++) di[i] = d[i][3];
          w1concep (d);
    }
    else
    {
          for (i=0; i<4; i++) tauc[i] = tau[i];
    }

    for (i=0; i<4; i++)
    {
      stress[ i] = tau[i];
      stressc[i] = tauc[i];
    }
/*---------------------------------- initialization at singularities ---*/
    facmin = tol2 * emod;
    for (i=0; i<4; i++)
    {
      for (j=0; j<4; j++)
      {
        if(fabs(d[i][j])>facmin) facmin = fabs(d[i][j]);
      }
    }
    
    if(facmin < tol2*emod)
    {
      facmin = tol2 * emod;
      for (i=0; i<4; i++) for (j=0; j<4; j++) d[i][j] *= facmin; 
    }

    if(wtype==plane_stress)
    { 
      for (i=0; i<4; i++)
      {
        ele->e.w1->elewa[0].ipwa[ip].sigi[i] = stress[i];
        ele->e.w1->elewa[0].ipwa[ip].epsi[i] = strain[i];
        ele->e.w1->elewa[0].ipwa[ip].di[  i] = di[i];
      }
    }
 }/* outside yield surface */
/*----------------------------------------------------------------------*/
end:
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
  if(istore==1||iupd==1)
  {
    for (i=0; i<4; i++)
    {
      ele->e.w1->elewa[0].ipwa[ip].sig[ i] = stress[i]  ;
      ele->e.w1->elewa[0].ipwa[ip].eps[ i] = strain[i]  ;
      ele->e.w1->elewa[0].ipwa[ip].sigc[i] = stressc[i] ;
      ele->e.w1->elewa[0].ipwa[ip].grad[i] = gradi[i] ;
    }
    ele->e.w1->elewa[0].ipwa[ip].dlam[0] = epst[0] ;
    ele->e.w1->elewa[0].ipwa[ip].dlam[1] = epst[1] ;
    ele->e.w1->elewa[0].ipwa[ip].yip     = yip     ; 
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_epc */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
