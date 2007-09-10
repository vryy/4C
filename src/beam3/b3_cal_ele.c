/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_cal_ele'

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0711 - 685-6574
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief performs all calculation for beam elements

<pre>                                                              fh 09/02
This routine performs all calculation for beam elements (stiffness matrix,
element load vector, internal forces)

</pre>
\param *ele            ELEMENT    (i/o) actual element
\param *data           B3_DATA     (o)  data set for gauss points
\param *mat            MATERIAL    (o)  actual material
\param **estif_global  ARRAY       (o)  global element stiffness matrix
\param *force          DOUBLE      (o)  global force vector
\param *action         CALC_ACTION (i)  action to do
\param init            INT         (i)  initialization (1) or calculation (2,3)


\warning There is nothing special in this routine
\return void
\sa calling:   b3_cal_sec() , b3_cal_lst() , b3_con_dof() , b3_cal_trn() ,
               b3_trans_stf() , b3_funct_deriv() , b3_jaco() , b3_boplin() ,
	       b3_boplin3D() , b3_call_mat() , b3_keku() , b3_load() ,
	       b3_loadlin() , b3_cal_stress() , b3_cal_stresslin(), b3_edisp(),
	       b3_exforce() , b3_cal_eps(), b3_con_warp()
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_cal_ele(ELEMENT     *ele,
                B3_DATA     *data,
                MATERIAL    *mat,
                ARRAY       *estif_global,
                DOUBLE      *force,
                CALC_ACTION *action,
		INT          init)
{
INT                 i,j;              /* some loopers */
INT                 ip=0;             /* code for actual integration point */
INT                 nir;              /* num GP in r direction */
INT                 lr, ls, lt;       /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 ndof = 6;         /* number of nodal dofs */
INT                 nist=3;           /* number of divisions in rs-direction (Simpson-Integration) */
INT                 numeledof;        /* number of element dofs */
INT                 calcstep = 0;     /* number of calculation step */
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 0;/* controls evaluation of new stresses    */
const INT           numgpc = 3;/* number of gausss point coordinates */
INT                 numeps = 6;/* number of nodal strains */
INT                 ike       ;/* controls calculation of linear elastic element stiffness matrix (12x12) */
const INT           max = (MAXDOFPERNODE+2)*MAXNOD_BEAM3;

DOUBLE              a,b;              /* a = element height, b = element width */
DOUBLE              e1,e2,e3;         /* GP-coords */
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE              fac;              /* integration factor of 1 point r,s,t   */
DOUBLE              l2;	              /* dx/dr for integration */
DOUBLE              det;              /* determinant of jacobian matrix */
DOUBLE              detrs;            /* da = dydz = detrs * drds        */
DOUBLE              pv;               /* possions ratio */
DOUBLE              gs;               /* inverse of shear correction factor */

static ARRAY    intforce_a;   /* local internal force vector */
static DOUBLE  *intforce;
static ARRAY    eload_a;  /* static element load vector */
static DOUBLE  *eload;
static ARRAY    KLA_a;    /* Kl*a for calculation of internal forces */
static DOUBLE **KLA;
static ARRAY    K_a;      /* local / global element stiffness matrix */
static DOUBLE **K;
static ARRAY    A_a;      /* transformation matrix loc->glob*/
static DOUBLE **A;
static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    HC_a;     /* hinge code vector */
static INT     *HC;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    ijm_a;    /* inverse of jacobian matrix */
static DOUBLE **ijm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    edisp_a;  /* element displacements */
static DOUBLE  *edisp;
static ARRAY    dummy_a;  /* dummy vector */
static DOUBLE  *dummy;
static ARRAY    eps_a;    /* strain vector at integration point */
static DOUBLE  *eps;
static DOUBLE **estif;    /* element stiffness matrix ke */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("b3_cal_ele");
#endif
/*------------------------------------------------- some working arrays */
istore = 0;

if (init==1)
{
funct     = amdef("funct"   ,&funct_a,MAXNOD_BEAM3,1              ,"DV");
HC        = amdef("HC    "  ,&HC_a  ,max+1,1                      ,"IV");
deriv     = amdef("deriv"   ,&deriv_a,2,MAXNOD_BEAM3              ,"DA");
D         = amdef("D"       ,&D_a   ,numeps,numeps                ,"DA");
K         = amdef("K"       ,&K_a   ,max,max                      ,"DA");
A         = amdef("A"       ,&A_a   ,max,max                      ,"DA");
xjm       = amdef("xjm"     ,&xjm_a ,numgpc,numgpc                ,"DA");
ijm       = amdef("ijm"     ,&ijm_a ,numgpc,numgpc                ,"DA");
bop       = amdef("bop"     ,&bop_a ,numeps,((numeps+2)*MAXNOD_BEAM3) ,"DA");
KLA       = amdef("KLA"     ,&KLA_a ,max,max                      ,"DA");
intforce  = amdef("intforce",&intforce_a,max,1                    ,"DV");
eload     = amdef("eload"   ,&eload_a,max,1                       ,"DV");
edisp     = amdef("edisp"   ,&edisp_a,max,1		          ,"DV");
dummy     = amdef("dummy"   ,&dummy_a,max,1			  ,"DV");
eps       = amdef("eps"     ,&eps_a,max,1		          ,"DV");
goto end;
}
/*------------------------------------------- integration parameters ---*/
b3_intg(ele,data,1);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
amzero(&K_a);
amzero(&A_a);
amzero(&HC_a);
amzero(&KLA_a);
amzero(&intforce_a);
amzero(&eload_a);
amzero(&edisp_a);
amzero(&dummy_a);
amzero(&eps_a);

estif   = estif_global->a.da;
a       = ele->e.b3->height;
b       = ele->e.b3->width;
ike     = ele->e.b3->ike;
nir     = ele->e.b3->nGP[0];
iel     = ele->numnp;
numeledof = ndof * iel;
detrs   = 0.25*a*b;
gs      = ele->e.b3->gs;

switch (mat->mattyp)
{
   case m_stvenant: pv = mat->m.stvenant->possionratio;
   break;
   case m_pl_mises: pv = mat->m.pl_mises->possionratio;
   break;
   default:
     dserror(" unknown type of material law");
   break;
}
/*----------------------------------------------------------------------*
|  calcstep = 1:	initialization of beam3 element                 |
|  calcstep = 2: 	calculation of linear elastic stiffness matrix  |
|  calcstep = 3:        calculation of global element load vector       |
|  calcstep = 4:        calculation of local internal force vector      |
|-----------------------------------------------------------------------*/
switch (*action)
{
/*-------- calculate structural linear element stiffness matrix --------*/
case calc_struct_linstiff:
  calcstep=2;
/*----------------------------------------------------------------------*
|       ike = 1:	Euler-Bernoulli-Beam element (Direct Stiffness) |
|       ike = 2: 	Timoshenko-Beam element (Finite Element method) |
|       ike = 3:        spatial Beam element according to Bathe         |
|-----------------------------------------------------------------------*/
  if (ike==1)
  {
     b3_cal_sec(ele);
     b3_cal_lst(ele,mat,K,HC);
     if (HC[0] != 0) b3_con_dof(K,HC,numeledof);
     b3_cal_trn(ele,A);
     b3_trans_stf(K,A,estif,numeledof,calcstep);
  }
  else if (ike==2)
  {
     b3_cal_sec(ele);
     b3_cal_lst(ele,mat,K,HC);
     b3_cal_trn(ele,A);
     amzero(&K_a);
     l2=ele->e.b3->length/2.;
     for (lr=0; lr<nir; lr++)
     {
        /*=============================== gaussian point and weight at it =*/
        e1   = data->xgrr[lr];
        facr = data->wgtr[lr]*l2;
        /*------------------------- shape functions and their derivatives */
        b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
        /*--------------------------------------- calculate operator B ---*/
        amzero(&bop_a);
	b3_boplin(bop,deriv,funct,iel,l2);
        /*------------------------------------------ call material law ---*/
        amzero(&D_a);
	b3_call_mat(ele,mat,eps,bop,D,NULL,lr,0,0);
        /*----------------------------------------------------------------*/
        b3_keku(K,bop,D,facr,numeledof,numeps);
     }
     if (HC[0] != 0) b3_con_dof(K,HC,numeledof);
     b3_trans_stf(K,A,estif,numeledof,calcstep);
  }
  else if (ike==3)
  {
     b3_cal_lst(ele,mat,K,HC);
     b3_cal_sec(ele);
     ele->e.b3->alpha=0.;
     b3_cal_trn(ele,A);
     amzero(&K_a);
     for (lr=0; lr<nir; lr++)
     {
        /*=============================== gaussian point and weight at it =*/
	e1   = data->xgrr[lr];
	facr = data->wgtr[lr];
	/*------------------------- shape functions and their derivatives */
	b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
	for (ls=0; ls<nist; ls++)
	{
	   /*=========================== lobatto point s and weight at it =*/
	   e2   = data->xlst[ls];
	   facs = data->wlst[ls];
	   for (lt=0; lt<nist; lt++)
	   {
	      /*======================== lobatto point t and weight at it =*/
	      e3   = data->xlst[lt];
	      fact = data->wlst[lt];
	      /*-----------Jacobian matrix at point r,s,t -----------------*/
              b3_jaco(funct,deriv,xjm,ijm,A,&det,e2,e3,ele,iel);
              /*-----------calculate operator B ---------------------------*/
              amzero(&bop_a);
	      b3_boplin3D(bop,deriv,funct,xjm,ijm,A,e2,e3,b,a,iel,init);
              /*---------- call material law ------------------------------*/
              amzero(&D_a);
	      b3_call_mat(ele,mat,eps,bop,D,NULL,lr,0,0);
	      /*------factor for integration wxr*wxs*wxt*det---------------*/
	      fac=facr*facs*fact*det;
	      numeps=3;
	      b3_keku(K,bop,D,fac,numeledof+2*iel,numeps);
	   }
	}
     }

     /*   local element stiffness matrix corresponds to global dofs,      */
     /*   so in case of hinges transformation to local dofs is needed.    */
     math_matmattrndense(KLA,K,A,numeledof,numeledof,numeledof,0,1.);
     math_matmatdense(K,A,KLA,numeledof,numeledof,numeledof,0,1.);
     b3_con_dof(K,HC,numeledof);
     /*----Reomve shear correction factor from torsional terms of estif------*/
     for (i=3; i<numeledof; i=i+6)
     {
        for (j=3; j<numeledof; j=j+6) K[i][j]=K[i][j]*gs;
     }
     for (i=numeledof; i<numeledof+2*iel; i++)
     {
        for (j=0; j<numeledof+2*iel; j++)
	{
	   K[i][j]=K[i][j]*gs;
	   if (i!=j) K[j][i]=K[i][j];
	}
     }
     /*----statical condensation of warping dofs-----------------------------*/
     b3_con_warp(K,iel,numeledof+2*iel);
     math_mattrnmatdense(KLA,A,K,numeledof,numeledof,numeledof,0,1.);
     math_matmatdense(estif,KLA,A,numeledof,numeledof,numeledof,0,1.);

  }
break;

/*-------- calculate element load vector ----------------------------------*/
case calc_struct_eleload:
   calcstep=3;
   if (ike==1)
   {
      b3_cal_lst(ele,mat,K,HC);
      b3_load(ele,mat,K,eload,HC,calcstep);
      for (i=0; i<12; i++) force[i]=eload[i];
   }
   else if (ike==2)
   {
      l2=ele->e.b3->length/2.;
      b3_cal_lst(ele,mat,K,HC);
      b3_cal_trn(ele,A);
      for (lr=0; lr<nir; lr++)
      /*=============================== gaussian point and weight at it ===*/
      {
          e1   = data->xgrr[lr];
          facr = data->wgtr[lr];
	  b3_funct_deriv(funct,deriv,e1,ele->distyp,0);
          b3_loadlin(ele,mat,funct,e1,facr,intforce,calcstep);
      }
      if (HC[0] != 0) b3_con_loadlin(intforce,HC,iel);
      math_mattrnvecdense(eload,A,intforce,numeledof,numeledof,0,1.);
      for (i=0; i<numeledof; i++)
      {
	 force[i] = eload[i];
      }
   }
   else
   {
      b3_cal_sec(ele); /* only for nonlinear computation */
      l2=ele->e.b3->length/2.;
      b3_cal_lst(ele,mat,K,HC);
      b3_cal_trn(ele,A);
      for (lr=0; lr<nir; lr++)
      /*=============================== gaussian point and weight at it ===*/
      {
          e1   = data->xgrr[lr];
          facr = data->wgtr[lr];
	  b3_funct_deriv(funct,deriv,e1,ele->distyp,0);
          b3_loadlin(ele,mat,funct,e1,facr,intforce,calcstep);
      }
      if (HC[0] != 0) b3_con_loadlin(intforce,HC,iel);
      math_mattrnvecdense(eload,A,intforce,numeledof,numeledof,0,1.);
      for (i=0; i<numeledof; i++)
      {
	 force[i] = eload[i];
      }
   }
break;

/*-------- calculate element internal force vector ------------------------*/
case calc_struct_stress:
   calcstep=4;
   am4zero(&(ele->e.b3->force_GP));
   if (ike==1)
   {
      b3_cal_lst(ele,mat,K,HC);
      b3_load(ele,mat,K,intforce,HC,calcstep);
      b3_cal_trn(ele,A);
      b3_trans_stf(K,A,KLA,numeledof,calcstep);
      b3_cal_force(ele,KLA,intforce,calcstep);
   }
   else if (ike==2)
   {
      l2=ele->e.b3->length/2.;
      b3_cal_lst(ele,mat,K,HC);
      b3_cal_trn(ele,A);
      b3_edisp(ele,edisp);
      /*--------------calculate local displacements--------------------------*/
      math_matvecdense(dummy,A,edisp,numeledof,numeledof,0,1.);
      /*---calculate internal forces at Gauss points ----------*/
      for (lr=0; lr<nir; lr++)
      {
         e1   = data->xgrr[lr];
	 b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
	 amzero(&bop_a);
	 b3_boplin(bop,deriv,funct,iel,l2);
	 amzero(&eps_a);
	 b3_cal_eps(eps,pv,dummy,bop,ndof,numeledof);
	 amzero(&D_a);
	 b3_call_mat(ele,mat,eps,bop,D,intforce,lr,0,1);
	 b3_cal_forcelin(ele,intforce,lr,0);
      }
      /* extrapolate internal forces from integration points to element nodes */
      b3_exforce(ele,data);
   }
   else
   {
      b3_cal_lst(ele,mat,K,HC);
      b3_cal_trn(ele,A);
      b3_edisp(ele,edisp);
      /*---calculate internal forces at Gauss points ----------*/
      for (lr=0; lr<nir; lr++)
      {
	 e1   = data->xgrr[lr];
	 b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
	 for (ls=0; ls<nist; ls++)
	 {
         /*=========================== lobatto point s and weight at it ===*/
	   e2   = data->xlst[ls];
	   facs = data->wlst[ls];
	   for (lt=0; lt<nist; lt++)
	   {
	      /*====================== lobatto point t and weight at it ===*/
	      e3   = data->xlst[lt];
	      fact = data->wlst[lt];
	      /*-----------Jacobian matrix at point r,s,t -----------------*/
              b3_jaco(funct,deriv,xjm,ijm,A,&det,e2,e3,ele,iel);
	      /*-----------calculate operator B ---------------------------*/
	      amzero(&bop_a);
	      b3_boplin3D(bop,deriv,funct,xjm,ijm,A,e2,e3,b,a,iel,init);
	      /*------factor for integration da=dr*ds----------------------*/
	      fac=detrs*facs*fact;
	      amzero(&eps_a);
	      numeps=3;
	      b3_cal_eps(eps,pv,edisp,bop,numeps,numeledof);
	      amzero(&D_a);
	      b3_call_mat(ele,mat,eps,bop,D,intforce,ip,0,1);
	      b3_cal_forcelin3D(ele,intforce,fac,lr,e2,e3,0);
	      ip+=1;
	   }
	 }
      }
      /* extrapolate internal forces from integration points to element nodes */
      b3_exforce(ele,data);
   }
break;

/*-------- calculate structural nonlinear element stiffness matrix -----*/
case calc_struct_nlnstiff:
   b3_cal_lst(ele,mat,K,HC);
   b3_cal_sec(ele);
   ele->e.b3->alpha=0.;
   b3_cal_trn(ele,A);
   b3_edisp(ele,edisp);
   amzero(&K_a);
   for (lr=0; lr<nir; lr++)
   {
      /*=============================== gaussian point and weight at it ===*/
      e1   = data->xgrr[lr];
      facr = data->wgtr[lr];
      /*------------------------- shape functions and their derivatives */
      b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
      for (ls=0; ls<nist; ls++)
      {
  	 /*=========================== lobatto point s and weight at it ===*/
  	 e2   = data->xlst[ls];
  	 facs = data->wlst[ls];
  	 for (lt=0; lt<nist; lt++)
  	 {
  	    /*======================== lobatto point t and weight at it ===*/
  	    e3	= data->xlst[lt];
  	    fact = data->wlst[lt];
  	    /*-----------Jacobian matrix at point r,s,t -------------------*/
  	    b3_jaco(funct,deriv,xjm,ijm,A,&det,e2,e3,ele,iel);
  	    /*-----------calculate operator B -----------------------------*/
  	    amzero(&bop_a);
  	    b3_boplin3D(bop,deriv,funct,xjm,ijm,A,e2,e3,b,a,iel,init);
  	    /*---------- call material law --------------------------------*/
  	    amzero(&D_a);
	    numeps=3;
	    b3_cal_eps(eps,pv,edisp,bop,numeps,numeledof);
	    b3_call_mat(ele,mat,eps,bop,D,intforce,ip,0,0);
  	    ip+=1;
	    fac=facr*facs*fact*det;
	    /*------calculate Fictitious nodal forces----------------------*/
	    math_mattrnvecdense(force,bop,intforce,numeledof,3,1,fac);
	    /*------factor for integration wxr*wxs*wxt*det-----------------*/
  	    fac=facr*facs*fact*det;
  	    b3_keku(K,bop,D,fac,numeledof+2*iel,numeps);
  	 }
      }
   }
   /*----statical condensation of warping dofs-----------------------------*/
   b3_con_warp(K,iel,numeledof+2*iel);
   if (HC[0] != 0)
   {
      /*   local element stiffness matrix corresponds to global dofs,      */
      /*   so in case of hinges transformation to local dofs is needed.    */
      math_matmatdense(KLA,A,K,numeledof,numeledof,numeledof,0,1.);
      b3_con_dof(KLA,HC,numeledof);
      math_mattrnmatdense(estif,A,KLA,numeledof,numeledof,numeledof,0,1.);
   }
   else
   {
     for (i=0; i<numeledof; i++)
     {
       for (j=0; j<numeledof; j++)
       {
         estif[i][j] = K[i][j];
       }
     }
   }
break;

/*-------- calculate structural update of actual step ---------------------*/
case calc_struct_update_istep:
   b3_cal_lst(ele,mat,K,HC);
   b3_cal_sec(ele);
   ele->e.b3->alpha=0.;
   b3_cal_trn(ele,A);
   b3_edisp(ele,edisp);
   for (lr=0; lr<nir; lr++)
   {
      /*=============================== gaussian point and weight at it ===*/
      e1   = data->xgrr[lr];
      facr = data->wgtr[lr];
      /*------------------------- shape functions and their derivatives */
      b3_funct_deriv(funct,deriv,e1,ele->distyp,1);
      for (ls=0; ls<nist; ls++)
      {
  	 /*=========================== lobatto point s and weight at it ===*/
  	 e2   = data->xlst[ls];
  	 facs = data->wlst[ls];
  	 for (lt=0; lt<nist; lt++)
  	 {
  	    /*======================== lobatto point t and weight at it ===*/
  	    e3	= data->xlst[lt];
  	    fact = data->wlst[lt];
  	    /*-----------Jacobian matrix at point r,s,t -------------------*/
  	    b3_jaco(funct,deriv,xjm,ijm,A,&det,e2,e3,ele,iel);
  	    /*-----------calculate operator B -----------------------------*/
  	    amzero(&bop_a);
  	    b3_boplin3D(bop,deriv,funct,xjm,ijm,A,e2,e3,b,a,iel,init);
  	    /*---------- call material law --------------------------------*/
  	    numeps=3;
	    amzero(&D_a);
            b3_cal_eps(eps,pv,edisp,bop,numeps,numeledof);
	    b3_call_mat(ele,mat,eps,bop,D,intforce,ip,1,0);
	    ip+=1;
  	 }
      }
   }
break;

default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_cal_ele */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
#endif
