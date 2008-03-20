/*!----------------------------------------------------------------------
\file wall1_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_element.H"
#include "Epetra_SerialDenseSolver.h"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra


/*----------------------------------------------------------------------*
 |                                                        mgit 03/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Wall1::ActionType act = Wall1::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Wall1::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Wall1::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Wall1::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Wall1::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Wall1::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Wall1::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Wall1::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Wall1::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Wall1::calc_struct_update_istep;
  else if (action=="calc_struct_update_genalpha_imrlike")  act = Wall1::calc_struct_update_genalpha_imrlike;
  else dserror("Unknown type of action for Wall1");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);
  switch(act)
  {
    case Wall1::calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat);
    }
    break;
    case Wall1::calc_struct_nlnstiffmass:
    case Wall1::calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat);
    }
    break;
    case calc_struct_update_istep:
    {
      // do something with internal EAS, etc parameters
      if (iseas_)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;
        blas.COPY((*alphao).M()*(*alphao).N(), (*alpha).A(), (*alphao).A());  // alphao := alpha
      }
    }
    break;
    case calc_struct_update_genalpha_imrlike:
    {
      // do something with internal EAS, etc parameters
      // this depends on the applied solution technique (static, generalised-alpha, 
      // or other time integrators)
      if (iseas_)
      {
        double alphaf = params.get<double>("alpha f", 0.0);  // generalised-alpha TIS parameter alpha_f
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1-alphaf}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;
        blas.SCAL((*alphao).M()*(*alphao).N(), -alphaf/(1.0-alphaf), (*alphao).A());  // alphao *= -alphaf/(1.0-alphaf)
        blas.AXPY((*alphao).M()*(*alphao).N(), 1.0/(1.0-alphaf), (*alpha).A(), (*alphao).A());  // alphao += 1.0/(1.0-alphaf) * alpha
        blas.COPY((*alpha).M()*(*alpha).N(), (*alphao).A(), (*alpha).A());  // alpha := alphao
      }
    }
    break;
    case calc_struct_stress:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      const DRT::UTILS::IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule_);
      Epetra_SerialDenseMatrix stress(intpoints.nquad,NUMSTR_W1);
      Epetra_SerialDenseMatrix strain(intpoints.nquad,NUMSTR_W1);
      bool cauchy = params.get<bool>("cauchy", false);
      w1_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,cauchy);
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;
    default:
      dserror("Unknown type of action for Wall1 %d", act);
  }
  return 0;

}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mgit 05/07|
 *----------------------------------------------------------------------*/

int DRT::ELEMENTS::Wall1::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0]; 
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

 // general arrays
  int       ngauss  = 0;
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  double det;

  // no. of nodes on this surface
  const int iel = NumNode();

  const DiscretizationType distype = this->Shape();

  const int numdf = 2;

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule_);

  //  vector<double>* thick = data_.Get<vector<double> >("thick");
  //  if (!thick) dserror("Cannot find vector of nodal thickness");

  Epetra_SerialDenseVector      funct(iel);
  Epetra_SerialDenseMatrix deriv(2,iel);

  double xrefe[2][MAXNOD_WALL1];
  double xcure[2][MAXNOD_WALL1];


  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {

    xrefe[0][k] = Nodes()[k]->X()[0];
    xrefe[1][k] = Nodes()[k]->X()[1];

    xcure[0][k] = xrefe[0][k] + mydisp[k*numdf+0];
    xcure[1][k] = xrefe[1][k] + mydisp[k*numdf+1];
  }


  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val");

  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */


	const double e1 = intpoints.qxg[ip][0];
	const double e2 = intpoints.qxg[ip][1];
	const double wgt = intpoints.qwgt[ip];

    /*-------------------- shape functions at gp e1,e2 on mid surface */
    //w1_shapefunctions(funct,deriv,e1,e2,iel,1);

    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    /*--------------------------------------- compute jacobian Matrix */
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);
  
    /*------------------------------------ integration factor  -------*/
    double fac=0;
    fac = wgt * det;

    // load vector ar
    double ar[2];
    // loop the dofs of a node
    // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]
    for (int i=0; i<2; ++i)
    {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    }

    // add load components
    for (int node=0; node<NumNode(); ++node)
      for (int dof=0; dof<2; ++dof)
         elevec1[node*2+dof] += funct[node] *ar[dof];

      ngauss++;
  } // for (int ip=0; ip<totngp; ++ip)


return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_nlnstiffmass(vector<int>&               lm,
                                           vector<double>&           disp,
                                           vector<double>&           residual,
                                           Epetra_SerialDenseMatrix* stiffmatrix,
                                           Epetra_SerialDenseMatrix* massmatrix,
                                           Epetra_SerialDenseVector* force,
                                           Epetra_SerialDenseMatrix* elestress,
                                           Epetra_SerialDenseMatrix* elestrain,
                                           struct _MATERIAL*         material,
                                           const bool                cauchy)
{
  const int numnode = NumNode();
  const int numdf   = 2;
  int       ngauss  = 0;
  const int nd      = numnode*numdf;


   // general arrays
  Epetra_SerialDenseVector      funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix boplin;
  boplin.Shape(4,2*numnode);
  Epetra_SerialDenseVector F;
  F.Size(4);
  Epetra_SerialDenseVector strain;
  strain.Size(4);
  double det;
  double xrefe[2][MAXNOD_WALL1];
  double xcure[2][MAXNOD_WALL1];
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);

   // for eas, in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas
  Epetra_SerialDenseMatrix F_enh;   // EAS matrix F_enh
  F_enh.Shape(4,1);
  Epetra_SerialDenseMatrix F_tot;   // EAS vector F_tot 
  F_tot.Shape(4,3);
  Epetra_SerialDenseMatrix p_stress;// first piola-kirchhoff stress vector 
  p_stress.Shape(4,1);
  Epetra_SerialDenseMatrix xjm0;    // Jacobian Matrix (origin)
  xjm0.Shape(2,2);
  Epetra_SerialDenseVector F0;      // Deformation Gradient (origin)
  F0.Size(4);
  Epetra_SerialDenseMatrix boplin0; // B operator (origin)
  boplin0.Shape(4,2*numnode);   
  Epetra_SerialDenseMatrix W0;      // W operator (origin)
  W0.Shape(4,2*numnode);   
  Epetra_SerialDenseMatrix G;       // G operator
  G.Shape(4,neas);
  Epetra_SerialDenseMatrix Z;        // Z operator
  Z.Shape(2*numnode,neas);
  Epetra_SerialDenseMatrix FCF;     // FCF^T  
  FCF.Shape(4,4);
  Epetra_SerialDenseMatrix Kda;     // EAS matrix Kda
  Kda.Shape(2*numnode,neas);
  Epetra_SerialDenseMatrix Kaa;     // EAS matrix Kaa
  Kaa.Shape(neas,neas);
  Epetra_SerialDenseVector feas;    // EAS portion of internal forces
  feas.Size(neas); 
  double detJ0;                     // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history
  
  // ------------------------------------ check calculation of mass matrix
  int imass=0;
  double density=0.0;
  if (massmatrix)
  {
    imass=1;
    /*------------------------------------------------ switch material type */
    switch (material->mattyp)
    {
    case m_stvenant:/*-------------------------------------- linear elastic */
      density = material->m.stvenant->density;
      break;
    case m_neohooke:/*------------------------------ kompressible neo-hooke */
      density = material->m.neohooke->density;
      break;
    case m_stvenpor:/*------------------------ porous linear elastic ---*/
      density = material->m.stvenpor->density;
      break;
    case m_pl_mises:/*--------------------------- von mises material law ---*/
      dserror("Ilegal typ of material for this element");
      break;
    case m_pl_mises_3D:/*-------------Stefan's von mises 3D material law ---*/
      dserror("Ilegal typ of material for this element");
      break;
    case m_pl_dp:/*------------------------- drucker prager material law ---*/
      dserror("Ilegal typ of material for this element");
      break;
    default:
      dserror("Ilegal typ of material for this element");
      break;
    }
  }

  /*------- get integraton data ---------------------------------------- */

  const int iel = numnode;

  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule_);

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {

    xrefe[0][k] = Nodes()[k]->X()[0];
    xrefe[1][k] = Nodes()[k]->X()[1];

    xcure[0][k] = xrefe[0][k] + disp[k*numdf+0];
    xcure[1][k] = xrefe[1][k] + disp[k*numdf+1];

  }

  if (iseas_ == true) 
  { 
	  
    /*
  	** EAS Update of alphas:
  	** the current alphas are (re-)evaluated out of
  	** Kaa and Kda of previous step to avoid additional element call.
  	** This corresponds to the (innermost) element update loop
  	** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
  	*/
  	alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get alpha of previous iteration

    // get stored EAS history
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");
    
    // we need the (residual) displacement at the previous step
    Epetra_SerialDenseVector res_d(2*numnode);
    for (int i = 0; i < (2*numnode); ++i) {
      res_d(i) = residual[i];
    }
    
    // add Kda . res_d to feas
    (*oldfeas).Multiply('T','N',1.0,(*oldKda),res_d,1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);
    /* end of EAS Update ******************/
    
      /* evaluation of EAS variables (which are constant for the following):
      ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
      ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
      ** -> T0^{-T}
      */
      w1_eassetup(boplin0,F0,xjm0,detJ0,xrefe,xcure,distype);
    }
//    else {cout << "Warning: Wall_Element without EAS" << endl;}


  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */

	
	const double e1 = intpoints.qxg[ip][0];
	const double e2 = intpoints.qxg[ip][1];
	const double wgt = intpoints.qwgt[ip];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    /*--------------------------------------- compute jacobian Matrix */
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);

    /*------------------------------------ integration factor  -------*/
    double fac = wgt * det * thickness_;
    

    /*------------------------------compute mass matrix if imass-----*/
    if (imass)
    {
     double facm = fac * density;
     for (int a=0; a<iel; a++)
     {
      for (int b=0; b<iel; b++)
      {
       (*massmatrix)(2*a,2*b)     += facm * funct(a) * funct(b); /* a,b even */
       (*massmatrix)(2*a+1,2*b+1) += facm * funct(a) * funct(b); /* a,b odd  */
      }
     }
    }
    /*----------------------------------- calculate operator Blin  ---*/
    w1_boplin(boplin,deriv,xjm,det,iel);
    //cout.precision(16);
    /*------------ calculate defgrad F^u, Green-Lagrange-strain E^u --*/
    w1_defgrad(F,strain,xrefe,xcure,boplin,iel);

    /*-calculate defgrad F in matrix notation and Blin in current conf.*/
    w1_boplin_cure(b_cure,boplin,F,numeps,nd);

    // EAS technology: "enhance the deformation gradient"  ---- --- EAS  
    if (iseas_ == true)
    {
      if (elestress != NULL)
    	  dserror("stress output for eas not yet impelemented"); 
      /*-----calculate the enhanced deformation gradient and--------------------
      -----alsoe the operators G, W0 and Z------------------------------------*/
    
      w1_call_defgrad_enh(F_enh,xjm0,xjm,detJ0,det,F0,*alpha,e1,e2,G,W0,boplin0,Z);
      
      /*-----total deformation gradient, Green-Lagrange-strain E^F -----------*/
      w1_call_defgrad_tot(F_enh,F_tot,F,strain);
      /* call material law----------------------------------------------------*/ 
      w1_call_matgeononl(strain,stress,C,numeps,material);

       // return gp strains (only in case of stress/strain output)
       if (elestress != NULL)
       {
         for (int i = 0; i < NUMSTR_W1; ++i)
           (*elestrain)(ip,i) = strain(i);
       }
       
       // return gp stresses (only in case of stress/strain output)
       if (elestress != NULL)
       {
         if (!cauchy)
         {
           (*elestress)(ip,0) = stress(0,0);
           (*elestress)(ip,1) = stress(1,1);
           (*elestress)(ip,2) = stress(0,2);
         }
         else
         {
           double detf = F[0]*F[1]-F[2]*F[3];
           Epetra_SerialDenseMatrix defgrad(2,2);
           Epetra_SerialDenseMatrix pkstress(2,2);
           defgrad(0,0) = F[0];
           defgrad(0,1) = F[2];
           defgrad(1,0) = F[3];
           defgrad(1,1) = F[1];
           pkstress(0,0)= stress(0,0);
           pkstress(0,1)= stress(0,2);
           pkstress(1,0)= stress(0,2);
           pkstress(1,1)= stress(1,1);

           Epetra_SerialDenseMatrix temp(2,2);
           Epetra_SerialDenseMatrix cauchy(2,2);
           temp.Multiply('N','T',1.0,pkstress,defgrad,0.0);
           cauchy.Multiply('N','N',1/detf,defgrad,temp,0.0);
           (*elestress)(ip,0) = cauchy(0,0);
           (*elestress)(ip,1) = cauchy(1,1);
           (*elestress)(ip,2) = cauchy(0,1);
         }  
       }
      /*-----first piola-kirchhoff stress vector------------------------------*/
      w1_stress_eas(stress,F_tot,p_stress);
      
      /*-----stiffness matrix kdd---------------------------------------------*/  
      w1_kdd(boplin,W0,F_tot,C,stress,FCF,*stiffmatrix,fac);
      /*-----matrix kda-------------------------------------------------------*/  
      w1_kda(FCF,W0,boplin,stress,G,Z,Kda,p_stress,fac);	
      /*-----matrix kaa-------------------------------------------------------*/  
      w1_kaa(FCF,stress,G,Kaa,fac);	
      /*-----nodal forces ----------------------------------------------------*/  
      if (force) w1_fint_eas(W0,boplin,G,p_stress,*force,feas,fac);

   }  
   else
   {	   
     w1_call_matgeononl(strain,stress,C,numeps,material);
     
     // return gp strains (only in case of stress/strain output)
     if (elestress != NULL)
     {
       for (int i = 0; i < NUMSTR_W1; ++i)
         (*elestrain)(ip,i) = strain(i);
     }
     
     // return gp stresses (only in case of stress/strain output)
     if (elestress != NULL)
     {
       if (!cauchy)
       {
         (*elestress)(ip,0) = stress(0,0);
         (*elestress)(ip,1) = stress(1,1);
         (*elestress)(ip,2) = stress(0,2);
       }
       else
       {
         double detf = F[0]*F[1]-F[2]*F[3];
         Epetra_SerialDenseMatrix defgrad(2,2);
         Epetra_SerialDenseMatrix pkstress(2,2);
         defgrad(0,0) = F[0];
         defgrad(0,1) = F[2];
         defgrad(1,0) = F[3];
         defgrad(1,1) = F[1];
         pkstress(0,0)= stress(0,0);
         pkstress(0,1)= stress(0,2);
         pkstress(1,0)= stress(0,2);
         pkstress(1,1)= stress(1,1);

         Epetra_SerialDenseMatrix temp(2,2);
         Epetra_SerialDenseMatrix cauchy(2,2);
         temp.Multiply('N','T',1.0,pkstress,defgrad,0.0);
         cauchy.Multiply('N','N',1/detf,defgrad,temp,0.0);
         (*elestress)(ip,0) = cauchy(0,0);
         (*elestress)(ip,1) = cauchy(1,1);
         (*elestress)(ip,2) = cauchy(0,1);
       }  
     }
         
     /*---------------------- geometric part of stiffness matrix kg ---*/
     if (stiffmatrix) w1_kg(*stiffmatrix,boplin,stress,fac,nd,numeps);
     
     /*------------------ elastic+displacement stiffness matrix keu ---*/
     if (stiffmatrix) w1_keu(*stiffmatrix,b_cure,C,fac,nd,numeps);
 
     /*--------------- nodal forces fi from integration of stresses ---*/
     if (force) w1_fint(stress,b_cure,*force,fac,nd);
   }

     ngauss++;
   
  } // for (int ip=0; ip<totngp; ++ip)

  
   // EAS technology: ------------------------------------------------------ EAS
   // subtract EAS matrices from disp-based Kdd to "soften" element
    
    if (iseas_ == true) 
    {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();
      
      
      Epetra_SerialDenseMatrix KdaKaa(2*NumNode(),neas); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('N', 'N', 1.0, Kda, Kaa, 1.0);

       
      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
      (*stiffmatrix).Multiply('N', 'T', -1.0, KdaKaa, Kda, 1.0);

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, feas, 1.0);

      // store current EAS data in history
      for (int i=0; i<neas; ++i) 
        for (int j=0; j<neas; ++j)
          (*oldKaainv)(i,j) = Kaa(i,j);
      
      for (int i=0; i<(2*NumNode()); ++i) 
        for (int j=0; j<neas; ++j)
        {
          (*oldKda)(i,j) = Kda(i,j);
          (*oldfeas)(j,0) = feas(j);
        }
   
    }
   
    // -------------------------------------------------------------------- EAS
  
  return;
}

/*----------------------------------------------------------------------*
 |  jacobian matrix (private)                                  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_jacobianmatrix(double xrefe[2][MAXNOD_WALL1],
                          const Epetra_SerialDenseMatrix& deriv,
                          Epetra_SerialDenseMatrix& xjm,
			              double* det,
                          const int iel)
{

   memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double));

   for (int k=0; k<iel; k++)
   {
        xjm(0,0) += deriv(0,k) * xrefe[0][k];
        xjm(0,1) += deriv(0,k) * xrefe[1][k];
        xjm(1,0) += deriv(1,k) * xrefe[0][k];
        xjm(1,1) += deriv(1,k) * xrefe[1][k];
   }

/*------------------------------------------ determinant of jacobian ---*/
     *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

      if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
/*----------------------------------------------------------------------*/

   return;
} // DRT::ELEMENTS::Wall1::w1_jacobianmatrix

/*----------------------------------------------------------------------*
 |  Matrix boplin in reference configuration (private)         mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_boplin(Epetra_SerialDenseMatrix& boplin,
                           Epetra_SerialDenseMatrix& deriv,
                           Epetra_SerialDenseMatrix& xjm,
                           double& det,
                           const int iel)
{

  double dum;
  double xji[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  dum = 1.0/det;
  xji[0][0] = xjm(1,1)* dum;
  xji[0][1] =-xjm(0,1)* dum;
  xji[1][0] =-xjm(1,0)* dum;
  xji[1][1] = xjm(0,0)* dum;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,x    0   |
       |   0    Nk,y |
       | Nk,y    0   |
       |  0     Nk,x |
  */
  for (int inode=0; inode<iel; inode++)
  {
    int dnode = inode*2;

    boplin(0,dnode+0) = deriv(0,inode)*xji[0][0] + deriv(1,inode)*xji[0][1];
    boplin(1,dnode+1) = deriv(0,inode)*xji[1][0] + deriv(1,inode)*xji[1][1];
    boplin(2,dnode+0) = boplin(1,dnode+1);
    boplin(3,dnode+1) = boplin(0,dnode+0);
  } /* end of loop over nodes */
  return;
}

/* DRT::ELEMENTS::Wall1::w1_boplin */

/*----------------------------------------------------------------------*
 | Deformation gradient F and Green-Langrange strain (private)  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_defgrad(Epetra_SerialDenseVector& F,
                           Epetra_SerialDenseVector& strain,
                           const double xrefe[][MAXNOD_WALL1],
                           const double xcure[][MAXNOD_WALL1],
                           Epetra_SerialDenseMatrix& boplin,
                           const int iel)
{
  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */

  memset(F.A(),0,F.N()*F.M()*sizeof(double));

  F[0] = 1;
  F[1] = 1;
  for (int inode=0; inode<iel; inode++)
  {
     F[0] += boplin(0,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);  // F_11
     F[1] += boplin(1,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);  // F_22
     F[2] += boplin(2,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);  // F_12
     F[3] += boplin(3,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);  // F_21
  } /* end of loop over nodes */

  /*-----------------------calculate Green-Lagrange strain E -------------*/
  strain[0] = 0.5 * (F[0] * F[0] + F[3] * F[3] - 1.0);  // E_11
  strain[1] = 0.5 * (F[2] * F[2] + F[1] * F[1] - 1.0);  // E_22
  strain[2] = 0.5 * (F[0] * F[2] + F[3] * F[1]);        // E_12
  strain[3] = strain[2];                                // E_21

  return;
}

/* DRT::ELEMENTS::Wall1::w1_defgrad */


/*----------------------------------------------------------------------*
 | Deformation gradient F in matrix notation and B in
 reference configuration (private)                             mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_boplin_cure(Epetra_SerialDenseMatrix& b_cure,
                                          Epetra_SerialDenseMatrix& boplin,
                                          Epetra_SerialDenseVector& F,
                                          int numeps,
                                          int nd)
{


     Epetra_SerialDenseMatrix Fmatrix;
     Fmatrix.Shape(4,4);


  /*---------------------------write Vector F as a matrix Fmatrix*/

     Fmatrix(0,0) = F[0];
     Fmatrix(0,2) = 0.5 * F[2];
     Fmatrix(0,3) = 0.5 * F[2];
     Fmatrix(1,1) = F[1];
     Fmatrix(1,2) = 0.5 * F[3];
     Fmatrix(1,3) = 0.5 * F[3];
     Fmatrix(2,1) = F[2];
     Fmatrix(2,2) = 0.5 * F[0];
     Fmatrix(2,3) = 0.5 * F[0];
     Fmatrix(3,0) = F[3];
     Fmatrix(3,2) = 0.5 * F[1];
     Fmatrix(3,3) = 0.5 * F[1];

    /*-------------------------------------------------int_b_cure operator*/
      memset(b_cure.A(),0,b_cure.N()*b_cure.M()*sizeof(double));
      for(int i=0; i<numeps; i++)
        for(int j=0; j<nd; j++)
          for(int k=0; k<numeps; k++)
            b_cure(i,j) += Fmatrix(k,i)*boplin(k,j);
    /*----------------------------------------------------------------*/

  return;
}

/* DRT::ELEMENTS::Wall1::w1_boplin_cure */

/*----------------------------------------------------------------------*
 | Constitutive matrix C and stresses (private)                mgit 05/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_call_matgeononl(Epetra_SerialDenseVector& strain,
                                              Epetra_SerialDenseMatrix& stress,
                                              Epetra_SerialDenseMatrix& C,
                                              const int numeps,
                                              struct _MATERIAL* material)

{
  /*--------------------------- call material law -> get tangent modulus--*/
    switch(material->mattyp)
    {
    case m_stvenant:/*--------------------------------- linear elastic ---*/
    {
      double ym = material->m.stvenant->youngs;
      double pv = material->m.stvenant->possionratio;


  /*-------------- some comments, so that even fluid people are able to
     understand this quickly :-)
     the "strain" vector looks like:

         | EPS_xx |
         | EPS_yy |
         | EPS_xy |
         | EPS_yx |

  */
  /*---------------------------------material-tangente-- plane stress ---*/
    switch(wtype_)
    {
    case plane_stress:
      {
      double e1=ym/(1. - pv*pv);
      double e2=pv*e1;
      double e3=e1*(1. - pv)/2.;

      C(0,0)=e1;
      C(0,1)=e2;
      C(0,2)=0.;
      C(0,3)=0.;

      C(1,0)=e2;
      C(1,1)=e1;
      C(1,2)=0.;
      C(1,3)=0.;

      C(2,0)=0.;
      C(2,1)=0.;
      C(2,2)=e3;
      C(2,3)=e3;

      C(3,0)=0.;
      C(3,1)=0.;
      C(3,2)=e3;
      C(3,3)=e3;
      }
      break;
    case plane_strain:
     /*----------- material-tangente - plane strain, rotational symmetry ---*/
      {
      double c1=ym/(1.0+pv);
      double b1=c1*pv/(1.0-2.0*pv);
      double a1=b1+c1;

      C(0,0)=a1;
      C(0,1)=b1;
      C(0,2)=0.;
      C(0,3)=0.;

      C(1,0)=b1;
      C(1,1)=a1;
      C(1,2)=0.;
      C(1,3)=0.;

      C(2,0)=0.;
      C(2,1)=0.;
      C(2,2)=c1/2.;
      C(2,3)=c1/2.;

      C(3,0)=0.;
      C(3,1)=0.;
      C(3,2)=c1/2;
      C(3,3)=c1/2;
      }
      break;
      default:
	dserror("nonsense");
  }
  /*-------------------------- evaluate 2.PK-stresses -------------------*/
  /*------------------ Summenschleife -> += (2.PK stored as vecor) ------*/

  Epetra_SerialDenseVector svector;
  svector.Size(3);

  for (int k=0; k<3; k++)
  {
    for (int i=0; i<numeps; i++)
    {
      svector(k) += C(k,i) * strain(i);
    }
  }
  /*------------------ 2.PK stored as matrix -----------------------------*/
  stress(0,0)=svector(0);
  stress(0,2)=svector(2);
  stress(1,1)=svector(1);
  stress(1,3)=svector(2);
  stress(2,0)=svector(2);
  stress(2,2)=svector(1);
  stress(3,1)=svector(2);
  stress(3,3)=svector(0);
 
  }

    break;
    default:
    dserror(" unknown type of material law");
    }

  return;
}

/* DRT::ELEMENTS::Wall1::w1_call_matgeononl */


/*----------------------------------------------------------------------*
| geometric stiffness part (total lagrange)                   mgit 05/07|
*----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_kg(Epetra_SerialDenseMatrix& estif,
                                 Epetra_SerialDenseMatrix& boplin,
                                 Epetra_SerialDenseMatrix& stress,
                                 double fac,
                                 int nd,
                                 int numeps)
{
  /*---------------------------------------------- perform B^T * SIGMA * B*/
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
      for(int r=0; r<numeps; r++)
         for(int m=0; m<numeps; m++)
            estif(i,j) += boplin(r,i)*stress(r,m)*boplin(m,j)*fac;

  return;
   
}

/* DRT::ELEMENTS::Wall1::w1_kg */

/*----------------------------------------------------------------------*
| elastic and initial displacement stiffness (total lagrange)  mgit 05/07               
*----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_keu(Epetra_SerialDenseMatrix& estif,
                                  Epetra_SerialDenseMatrix& b_cure,
                                  Epetra_SerialDenseMatrix& C,
                                  double fac,
                                  int nd,
                                  int numeps)
{
	
  /*------------- perform B_cure^T * D * B_cure, whereas B_cure = F^T * B */
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
        for(int k=0; k<numeps; k++)
           for(int m=0; m<numeps; m++)
             estif(i,j) +=  b_cure(k,i)*C(k,m)*b_cure(m,j)*fac;

  return;
}

/* DRT::ELEMENTS::Wall1::w1_keu */

/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) mgit 05/07  |
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_fint(Epetra_SerialDenseMatrix& stress,
                                   Epetra_SerialDenseMatrix& b_cure,
                                   Epetra_SerialDenseVector& intforce,
                                   double fac,
                                   int nd)

{
  Epetra_SerialDenseVector st;
  st.Size(4);

  st[0] = fac * stress(0,0);
  st[1] = fac * stress(1,1);
  st[2] = fac * stress(0,2);
  st[3] = fac * stress(0,2);

  for(int i=0; i<nd; i++)
    for(int j=0; j<4; j++)
      intforce[i] += b_cure(j,i)*st[j];

  return;
}

/* DRT::ELEMENTS::Wall1::w1_fint */

/*----------------------------------------------------------------------*
 |  setup of constant EAS data (private)                       mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_eassetup(
		  Epetra_SerialDenseMatrix& boplin0,
		  Epetra_SerialDenseVector& F0,       // deformation gradient at origin
		  Epetra_SerialDenseMatrix& xjm0,     // jacobian matrix at origin
		  double& detJ0,                      // det of Jacobian at origin
          double xrefe[2][MAXNOD_WALL1],      // material element coords
          double xcure[2][MAXNOD_WALL1],      // current element coords
          const DRT::Element::DiscretizationType& distype)

{
  // derivatives at origin
   Epetra_SerialDenseMatrix deriv0;
   deriv0.Shape(2,NumNode());

   DRT::UTILS::shape_function_2D_deriv1(deriv0,0.0,0.0,distype);

  // compute jacobian matrix at origin
  memset(xjm0.A(),0,xjm0.N()*xjm0.M()*sizeof(double));
  for (int k=0; k<NumNode(); k++)
  {
	  xjm0(0,0) += deriv0(0,k) * xrefe[0][k];  // X,r += (X,r)^k
	  xjm0(0,1) += deriv0(0,k) * xrefe[1][k];  // Y,r += (Y,r)^k
	  xjm0(1,0) += deriv0(1,k) * xrefe[0][k];  // X,s += (X,s)^k
	  xjm0(1,1) += deriv0(1,k) * xrefe[1][k];  // Y,s += (Y,s)^k
  }

  /*------------------------------------------ determinant of jacobian ---*/
  detJ0 = xjm0[0][0] * xjm0[1][1] - xjm0[1][0] * xjm0[0][1];

  if (detJ0<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");


  //compute boplin at origin (boplin0)

  int inode0;
  int dnode0;
  double xji0[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  xji0[0][0] = xjm0(1,1) / detJ0;
  xji0[0][1] =-xjm0(0,1) / detJ0;
  xji0[1][0] =-xjm0(1,0) / detJ0;
  xji0[1][1] = xjm0(0,0) / detJ0;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,X    0   |
       |   0    Nk,Y |
       | Nk,Y    0   |
       |  0     Nk,X |
  */
  for (inode0=0; inode0<NumNode(); inode0++)
  {
    dnode0 = inode0*2;

    boplin0(0,dnode0+0) = deriv0(0,inode0)*xji0[0][0] + deriv0(1,inode0)*xji0[0][1];
    boplin0(1,dnode0+1) = deriv0(0,inode0)*xji0[1][0] + deriv0(1,inode0)*xji0[1][1];
    boplin0(2,dnode0+0) = boplin0(1,dnode0+1);
    boplin0(3,dnode0+1) = boplin0(0,dnode0+0);
  }


  //compute deformation gradient at origin (F0)

  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */

  memset(F0.A(),0,F0.N()*F0.M()*sizeof(double));

  F0[0]=1;
  F0[1]=1;
  for (int inode=0; inode<NumNode(); inode++)
  {
     F0[0] += boplin0(0,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);
     F0[1] += boplin0(1,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);
     F0[2] += boplin0(2,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);
     F0[3] += boplin0(3,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);
  } /* end of loop over nodes */


  return;
} // end of w1_eassetup


/*----------------------------------------------------------------------*
 |get the enhanced deformation gradient and
 also the operators G,W0 and Z (private) mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_defgrad_enh(
		  Epetra_SerialDenseMatrix& F_enh,
		  Epetra_SerialDenseMatrix xjm0,
		  Epetra_SerialDenseMatrix xjm,
		  double detJ0,
		  double det,
		  Epetra_SerialDenseVector F0,
		  Epetra_SerialDenseMatrix alpha,
		  double e1,
		  double e2,
		  Epetra_SerialDenseMatrix& G,
		  Epetra_SerialDenseMatrix& W0,
		  Epetra_SerialDenseMatrix boplin0,
		  Epetra_SerialDenseMatrix& Z)     
{
   
  // EAS 
	
  /* eas is the EAS interpolation of 4 modes, based on
  **
  **     M1 = r 0    M2 = 0 s    M3 = 0 0    M4 = 0 0
  **          0 0         0 0         r 0         0 s
  **
  **
  **     M = M1*alpha1 + M2*alpha2 + M3*alpha3 + M4*alpha4
  **
  */

  Epetra_SerialDenseMatrix M;
  M.Shape(2,2);

  Epetra_SerialDenseMatrix M_temp;
  M_temp.Shape(2,2);

  // fill up 4 EAS matrices at each gp

  M(0,0) = e1*alpha(0,0);
  M(0,1) = e2*alpha(1,0);
  M(1,0) = e1*alpha(2,0);
  M(1,1) = e2*alpha(3,0);

  // inverse of jacobian matrix at element origin
  Epetra_SerialDenseMatrix xjm_inv0;
  xjm_inv0.Shape(2,2);

  xjm_inv0(0,0) = xjm0(1,1) / det;
  xjm_inv0(0,1) =-xjm0(0,1) / det;
  xjm_inv0(1,0) =-xjm0(1,0) / det;
  xjm_inv0(1,1) = xjm0(0,0) / det;

  Epetra_SerialDenseMatrix A;       // A operator  
  A.Shape(2,2);
   
  M_temp.Multiply('N','T',1.0,M,xjm_inv0,0.0);
  A.Multiply('T','N',detJ0/det,xjm0,M_temp,0.0);
  
  // enhanced deformation gradient at origin (four rows, one column)
  
  F_enh(0,0) = A(0,0)*F0(0) + A(1,0)*F0(2);
  F_enh(1,0) = A(1,1)*F0(1) + A(0,1)*F0(3);
  F_enh(2,0) = A(0,1)*F0(0) + A(1,1)*F0(2);
  F_enh(3,0) = A(1,0)*F0(1) + A(0,0)*F0(3);
  
  // get operator W0 
  
  // write matrix A in a different way (matrix 4x4)
  
  Epetra_SerialDenseMatrix A_big;
  A_big.Shape(4,4);  // entries are zero

  A_big(0,0) = A(0,0);
  A_big(0,2) = A(1,0);
  A_big(1,1) = A(1,1);
  A_big(1,3) = A(0,1);
  A_big(2,0) = A(0,1);
  A_big(2,2) = A(1,1);
  A_big(3,1) = A(1,0);
  A_big(3,3) = A(0,0);
  
  // multiplication A_big x boplin0
  
  W0.Multiply('N','N',1.0,A_big,boplin0,0.0);

  // calculation operators G and Z, therfore matrices A are needed 
  // without alphas
  
  // vector M_ges, includes the matrices M1 to M4
  vector <Epetra_SerialDenseMatrix> M_ges(neas);
  
  // vector A_ges, includes the matrices A1 to A4
    vector <Epetra_SerialDenseMatrix> A_ges(4);  
    
   for (int ieas=0; ieas<neas; ieas++)   
   {
     (M_ges[ieas]).Shape(2,2);
     
     memset((M_ges[ieas]).A(),0,(M_ges[ieas]).N()*(M_ges[ieas]).M()*sizeof(double));
   
     (A_ges[ieas]).Shape(2,2);
   }
  
  // fill M-Matrixes, not including eas-parameters alpha
  
  (M_ges[0])(0,0) = e1;  //(M_ges[0])(0,0)=e1;
  (M_ges[1])(0,1) = e2;  //(M_ges[1])(0,1)=e2;
  (M_ges[2])(1,0) = e1;  //(M_ges[2])(1,0)=e1;  
  (M_ges[3])(1,1) = e2;  //(M_ges[3])(1,1)=e2;
  
  // declaration of matrix (4x4) without eas-parameter alpha
  
  Epetra_SerialDenseMatrix Awa_big;
  Awa_big.Shape(4,4);
  
  // declaration of matrix WO without eas-parameter alpha
  
  Epetra_SerialDenseMatrix W0wa;
  W0wa.Shape(4,2*NumNode());
 

  for (int i=0; i<neas; i++)  // loop over eas-parameter
  {  
    M_temp.Multiply('N','T',1.0,M_ges[i],xjm_inv0,0.0);
    A_ges[i].Multiply('T','N',detJ0/det,xjm0,M_temp,0.0);
    
    // fill G-operator
  
    G(0,i) = A_ges[i](0,0)*F0(0) + A_ges[i](1,0)*F0(2);
    G(1,i) = A_ges[i](1,1)*F0(1) + A_ges[i](0,1)*F0(3);
    G(2,i) = A_ges[i](0,1)*F0(0) + A_ges[i](1,1)*F0(2);
    G(3,i) = A_ges[i](1,0)*F0(1) + A_ges[i](0,0)*F0(3);

    
    memset(Awa_big.A(),0,Awa_big.N()*Awa_big.M()*sizeof(double));
        
    Awa_big(0,0) = A_ges[i](0,0);
    Awa_big(0,2) = A_ges[i](1,0);
    Awa_big(1,1) = A_ges[i](1,1);
    Awa_big(1,3) = A_ges[i](0,1);
    Awa_big(2,0) = A_ges[i](0,1);
    Awa_big(2,2) = A_ges[i](1,1);
    Awa_big(3,1) = A_ges[i](1,0);
    Awa_big(3,3) = A_ges[i](0,0);
      
    // calculate operator W0wa without eas-parameters alpha
    
    W0wa.Multiply('N','N',1.0,Awa_big,boplin0,0.0);

    // fill Z-operator
    
    for (int indof; indof<NumNode(); indof=indof+2)
    {
      Z(indof,i) = W0wa(0,indof);
      Z(indof+1,i) = W0wa(1,indof+1);
    }
    
     
  }  // end loop over eas parameter

  return;
}  // end of w1_call_fenh

/*----------------------------------------------------------------------*
 |total deformation gradient and green lagrange strain (private)mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_defgrad_tot(
		  Epetra_SerialDenseMatrix F_enh,
		  Epetra_SerialDenseMatrix& F_tot,
		  Epetra_SerialDenseVector F,
		  Epetra_SerialDenseVector& strain)
{

  // total deformation gradient in matrix notation
  F_tot(0,0) = F(0) + F_enh(0,0);
  F_tot(0,2) = F(2) + F_enh(2,0);
  F_tot(1,1) = F(1) + F_enh(1,0);
  F_tot(1,2) = F(3) + F_enh(3,0);
  F_tot(2,1) = F(2) + F_enh(2,0);
  F_tot(2,2) = F(0) + F_enh(0,0);
  F_tot(3,0) = F(3) + F_enh(3,0);
  F_tot(3,2) = F(1) + F_enh(1,0);
  

  /*-----------------------calculate Green-Lagrange strain ------------*/
  strain[0] = 0.5 * (F_tot(0,0) * F_tot(0,0) + F_tot(3,0) * F_tot(3,0) - 1.0);
  strain[1] = 0.5 * (F_tot(0,2) * F_tot(0,2) + F_tot(3,2) * F_tot(3,2) - 1.0);
  strain[2] = 0.5 * (F_tot(0,0) * F_tot(0,2) + F_tot(3,0) * F_tot(3,2));
  strain[3] = strain[2];

  return;
} // end of w1_call_defgrad_tot

/*-----------------------------------------------------------------------------*
 |first piola-kirchhoff stress vector                       (private)mgit 02/08|
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_stress_eas(
		  Epetra_SerialDenseMatrix stress,
		  Epetra_SerialDenseMatrix F_tot,
		  Epetra_SerialDenseMatrix& p_stress)     
{

  /*-------------reduce stress matrix-----------------------------------------*/

  Epetra_SerialDenseMatrix stress_red;  // 2. piola-krichhoff (vector (3-dim))
  stress_red.Shape(3,1);
	
  stress_red(0,0) = stress(0,0);     // S_11
  stress_red(1,0) = stress(1,1);     // S_22
  stress_red(2,0) = stress(0,2);     // S_12 (=S_21)

  /*-first piola-kirchhoff stress vector--------------------------------------*/
  p_stress.Multiply('N','N',1.0,F_tot,stress_red,0.0);
  
  return;
} // end of w1_p_stress


/*-----------------------------------------------------------------------------*
| calculate stiffness matrix kdd                                     mgit 03/08|               
*-----------------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_kdd(Epetra_SerialDenseMatrix boplin,
                                  Epetra_SerialDenseMatrix W0,
                                  Epetra_SerialDenseMatrix F_tot,
                                  Epetra_SerialDenseMatrix C,
                                  Epetra_SerialDenseMatrix stress,
                                  Epetra_SerialDenseMatrix& FCF,
                                  Epetra_SerialDenseMatrix& estif,
                                  double fac)
{
  
  // contitutive matrix (3x3)
  Epetra_SerialDenseMatrix C_red;
  C_red.Shape(3,3);
  
  C_red(0,0) = C(0,0); 
  C_red(0,1) = C(0,1);
  C_red(0,2) = C(0,2);
  C_red(1,0) = C(1,0);
  C_red(1,1) = C(1,1);
  C_red(1,2) = C(1,2);
  C_red(2,0) = C(2,0);
  C_red(2,1) = C(2,1);
  C_red(2,2) = C(2,2);
  
  Epetra_SerialDenseMatrix BplusW;
  BplusW.Shape(4,2*NumNode());
    
  Epetra_SerialDenseMatrix Temp;
  Temp.Shape(4,3);
  
  Epetra_SerialDenseMatrix Temp1;
  Temp1.Shape(4,2*NumNode());
    
  Epetra_SerialDenseMatrix Temp2;
  Temp2.Shape(2*NumNode(),2*NumNode());
  
  Epetra_SerialDenseMatrix Temp3;
  Temp3.Shape(4,2*NumNode());
  
  Epetra_SerialDenseMatrix Temp4;
  Temp4.Shape(2*NumNode(),2*NumNode());
    
  // BplusW = B+W0  
  for (int i=0; i<4; i++)
    for (int j=0; j<2*NumNode(); j++)
      BplusW(i,j) = boplin(i,j) + W0(i,j);
  
  // Temp (4x3) = F*C
  Temp.Multiply('N','N',1.0,F_tot,C_red,0.0);
  
  // FCF^T (4x4) = Temp*F^T
  FCF.Multiply('N','T',1.0,Temp,F_tot,0.0);
  
  // Temp1 (4x8) = FCF^T * (B+W0)
  Temp1.Multiply('N','N',1.0,FCF,BplusW,0.0);

  // Temp2 (8x8) = fac * (B+W0)^T * FCF^T * (B+W0)
  Temp2.Multiply('T','N',fac,BplusW,Temp1,0.0);

  // Temp3 (4x8) = S*(B+W0)
  Temp3.Multiply('N','N',1.0,stress,BplusW,0.0);

  // Temp4 = fac * (B+W0)^T * S * (B+W0)
  Temp4.Multiply('T','N',fac,BplusW,Temp3,0.0);

  // Kdd = (B+W0)^T*FCF^T*(B+W0) + (B+W0)^T*S*(B+W0)   
  for (int i=0; i<2*NumNode(); i++)
    for (int j=0; j<2*NumNode(); j++)
      estif(i,j) += Temp2(i,j) + Temp4(i,j);

  return;
}

/* DRT::ELEMENTS::Wall1::w1_kdd */

/*-----------------------------------------------------------------------------*
| calculate matrix kda                                               mgit 03/08|               
*-----------------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_kda(Epetra_SerialDenseMatrix FCF,
		                          Epetra_SerialDenseMatrix W0,
		                          Epetra_SerialDenseMatrix boplin,   
		                          Epetra_SerialDenseMatrix stress,
                                  Epetra_SerialDenseMatrix G,
                                  Epetra_SerialDenseMatrix Z,
                                  Epetra_SerialDenseMatrix& Kda,
                                  Epetra_SerialDenseMatrix p_stress,
                                  double fac)
{
  
  Epetra_SerialDenseMatrix BplusW;
  BplusW.Shape(4,2*NumNode());
  
  // BplusW = B+W0  
  for(int i=0; i<4; i++)
    for(int j=0; j<2*NumNode(); j++)
      BplusW (i,j)=boplin(i,j)+W0(i,j);
  
  Epetra_SerialDenseMatrix Temp1;
  Temp1.Shape(4,neas);
    
  Epetra_SerialDenseMatrix Temp2;
  Temp2.Shape(2*NumNode(),neas);
  
  Epetra_SerialDenseMatrix Temp3;
  Temp3.Shape(4,neas);
  
  Epetra_SerialDenseMatrix Temp4;
  Temp4.Shape(2*NumNode(),neas);
  
  Epetra_SerialDenseMatrix Temp5;
  Temp5.Shape(2*NumNode(),neas);
    
  // Temp1 = FCF^T*G
  Temp1.Multiply('N','N',1.0,FCF,G,0.0);
  
  // Temp2 = fac * (B+W0)^T*FCF^T*G
  Temp2.Multiply('T','N',fac,BplusW,Temp1,0.0);

  // Temp3 = S*(G)
  Temp3.Multiply('N','N',1.0,stress,G,0.0);

  // Temp4 = fac * (B+W0)^T*S*G
  Temp4.Multiply('T','N',fac,BplusW,Temp3,0.0);

  // Temp5 = fac * P*Z
  for (int i=0; i<NumNode(); i++)
  {
    for (int ieas=0; ieas<neas; ieas++)
    {  
	Temp5(i*2,ieas)= (p_stress(0,0)*Z(i*2,ieas)+p_stress(2,0)*Z(i*2+1,ieas))*fac;
	Temp5(i*2+1,ieas)= (p_stress(3,0)*Z(i*2,ieas)+p_stress(1,0)*Z(i*2+1,ieas))*fac;
    }
  }	

  // Kda (8x4) = (B+W0)^T*FCF^T*G) + (B+W0)^T*S*G + PZ   
  for (int i=0; i<2*NumNode(); i++)
    for (int j=0; j<4; j++)
      Kda(i,j) += Temp2(i,j) + Temp4(i,j) + Temp5(i,j);

  return;
}

/* DRT::ELEMENTS::Wall1::w1_kda */


/*-----------------------------------------------------------------------------*
| calculate matrix kaa                                               mgit 03/08|               
*-----------------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_kaa(Epetra_SerialDenseMatrix FCF,
		                          Epetra_SerialDenseMatrix stress,
                                  Epetra_SerialDenseMatrix G,
                                  Epetra_SerialDenseMatrix& Kaa,
                                  double fac)
{
  Epetra_SerialDenseMatrix Temp1;
  Temp1.Shape(4,neas);
    
  Epetra_SerialDenseMatrix Temp2;
  Temp2.Shape(neas,neas);
  
  Epetra_SerialDenseMatrix Temp3;
  Temp3.Shape(4,neas);
  
  Epetra_SerialDenseMatrix Temp4;
  Temp4.Shape(neas,neas);
  
  // Temp1 = FCF*G
  Temp1.Multiply('N','N',1.0,FCF,G,0.0);
  
  // Temp2 = fac * G^T*FCF^T*G
  Temp2.Multiply('T','N',fac,G,Temp1,0.0);

  // Temp3 = S*G
  Temp3.Multiply('N','N',1.0,stress,G,0.0);

  // Temp4 = fac * G^T*S*G
  Temp4.Multiply('T','N',fac,G,Temp3,0.0);

  // Kaa = G^T*FCF^T*G + G^T*S*G    
  for (int i=0; i<neas; i++)
    for (int j=0; j<neas; j++)
      Kaa(i,j) += Temp2(i,j) + Temp4(i,j);
     
  return;
}

/* DRT::ELEMENTS::Wall1::w1_kaa */


/*-----------------------------------------------------------------------------*
| calculate internal forces fint(displacements u) and feas           mgit 03/08|               
*-----------------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_fint_eas(Epetra_SerialDenseMatrix W0,
		                          Epetra_SerialDenseMatrix boplin,
                                  Epetra_SerialDenseMatrix G,
                                  Epetra_SerialDenseMatrix p_stress,
                                  Epetra_SerialDenseVector& intforce,
                                  Epetra_SerialDenseVector& feas,
                                  double fac)

{
  Epetra_SerialDenseMatrix BplusW;
  BplusW.Shape(4,2*NumNode());
		  
  // BplusW = B+W0  
  for(int i=0; i<4; i++)
    for(int j=0; j<2*NumNode(); j++)
	  BplusW(i,j) = boplin(i,j) + W0(i,j);
	
	
  Epetra_SerialDenseMatrix Temp1;
  Temp1.Shape(2*NumNode(),1);

  Epetra_SerialDenseMatrix Temp2;
  Temp2.Shape(neas,1);

   // Temp1 (8x1) = (BL+W0)^T*p_stress
  Temp1.Multiply('T','N',1.0,BplusW,p_stress,0.0);
  
  for (int i=0; i<2*NumNode(); i++)
     intforce(i) += fac*Temp1(i,0);
  
  // Temp2 = G^T*p_stress
  Temp2.Multiply('T','N',1.0,G,p_stress,0.0);

  for (int i=0; i<neas; i++)
     feas(i) += fac*Temp2(i,0);

  return;
}

/* DRT::ELEMENTS::Wall1::w1_fint_eas */






#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
