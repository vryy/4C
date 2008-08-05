/*======================================================================*/
/*!
\file wall1_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*/
// macros
#ifdef D_WALL1
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

/*----------------------------------------------------------------------*/
// headers
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/stvenantkirchhoff.H"

/*----------------------------------------------------------------------*/
// namespaces
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
int DRT::ELEMENTS::Wall1::Evaluate(ParameterList&            params,
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
  else if (action=="calc_struct_nlnstifflmass") act = Wall1::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_nlnstiff_gemm") act = Wall1::calc_struct_nlnstiff_gemm;
  else if (action=="calc_struct_stress")        act = Wall1::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Wall1::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Wall1::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Wall1::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Wall1::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = Wall1::calc_struct_reset_istep;
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
    case Wall1::calc_struct_nlnstifflmass:
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
      if (act==calc_struct_nlnstifflmass) w1_lumpmass(&elemat2);
    }
    break;
    // NULL-pointer for mass matrix in case of calculating only stiff matrix 
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
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;
    case Wall1::calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix (initialised to zero)
      // This matrix is not utterly useless. It is used to apply EAS-stuff in a linearised manner
      // onto the internal force vector.
      Epetra_SerialDenseMatrix myemat(lm.size(),lm.size());
      w1_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;
    case Wall1::calc_struct_nlnstiff_gemm:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> dispo = discretization.GetState("old displacement");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (dispo==Teuchos::null or disp==Teuchos::null or res==Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydispo(lm.size());
      DRT::UTILS::ExtractMyValues(*dispo,mydispo,lm);
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      FintStiffMassGEMM(params,lm,mydispo,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;
    case calc_struct_update_istep:
    {
      // do something with internal EAS, etc parameters
      if (iseas_)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;  // BLAS front-end dummy
        blas.COPY((*alphao).M()*(*alphao).N(), (*alpha).A(), (*alphao).A());  // alphao := alpha
      }
    }
    break;
    case calc_struct_update_imrlike:
    {
      // do something with internal EAS, etc parameters
      // this depends on the applied solution technique (static, generalised-alpha,
      // or other time integrators)
      if (iseas_)
      {
        double alphaf = params.get<double>("alpha f", 0.0);  // generalised-alpha TIS parameter alpha_f
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1-alphaf}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;  // BLAS front-end dummy
        blas.SCAL((*alphao).M()*(*alphao).N(), -alphaf/(1.0-alphaf), (*alphao).A());  // alphao *= -alphaf/(1.0-alphaf)
        blas.AXPY((*alphao).M()*(*alphao).N(), 1.0/(1.0-alphaf), (*alpha).A(), (*alphao).A());  // alphao += 1.0/(1.0-alphaf) * alpha
        blas.COPY((*alpha).M()*(*alpha).N(), (*alphao).A(), (*alpha).A());  // alpha := alphao
      }
    }
    break;
    case calc_struct_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (iseas_)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;
        blas.COPY((*alphao).M()*(*alphao).N(), (*alphao).A(), (*alpha).A());  // alpha := alphao
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
      Epetra_SerialDenseMatrix stress(intpoints.nquad,Wall1::numstr_);
      Epetra_SerialDenseMatrix strain(intpoints.nquad,Wall1::numstr_);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      if (iostrain!="euler_almansi") w1_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,cauchy);
      else dserror("requested strain output option not yet available for wall1");
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;
    case Wall1::calc_struct_eleload:
    {
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
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
  double curvefac = 1.0;  // default time curve factor
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // general arrays
  int ngauss = 0;  // total number of Gauss points
  Epetra_SerialDenseMatrix xjm(2,2);  // iso-parametric Jacobian
  double det;  // determinant of iso-parametric Jacobian

  // no. of nodes on this surface
  const int iel = NumNode();

  // quad, tri, etc
  const DiscretizationType distype = Shape();

  // number of DOFs at each element node
  const int numdf = 2;

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule_);

  //  vector<double>* thick = data_.Get<vector<double> >("thick");
  //  if (!thick) dserror("Cannot find vector of nodal thickness");

  // shape functions
  Epetra_SerialDenseVector funct(iel);
  // natural derivatives of shape funcions
  Epetra_SerialDenseMatrix deriv(2,iel);

  // reference co-ordinates of element nodes
  Epetra_SerialDenseMatrix xrefe(2,iel);
  // current co-ordinates of element nodes
  Epetra_SerialDenseMatrix xcure(2,iel);


  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];

    xcure(0,k) = xrefe(0,k) + mydisp[k*numdf+0];
    xcure(1,k) = xrefe(1,k) + mydisp[k*numdf+1];
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
      ar[i] = fac * (*onoff)[i] * (*val)[i] * curvefac;
    }

    // add load components
    for (int node=0; node<NumNode(); ++node)
      for (int dof=0; dof<2; ++dof)
         elevec1[node*2+dof] += funct[node] * ar[dof];

    ngauss++;
  } // for (int ip=0; ip<totngp; ++ip)

  // finished
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_nlnstiffmass(const vector<int>&        lm,
                                           const vector<double>&     disp,
                                           const vector<double>&     residual,
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
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);

  // for EAS, in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas
  Epetra_SerialDenseMatrix* F_enh;  // EAS matrix F_enh
  Epetra_SerialDenseMatrix* F_tot;  // EAS vector F_tot
  Epetra_SerialDenseMatrix* p_stress;  // first piola-kirchhoff stress vector
  Epetra_SerialDenseMatrix* xjm0;  // Jacobian Matrix (origin)
  Epetra_SerialDenseVector* F0;  // Deformation Gradient (origin)
  Epetra_SerialDenseMatrix* boplin0; // B operator (origin)
  Epetra_SerialDenseMatrix* W0;  // W operator (origin)
  Epetra_SerialDenseMatrix* G;  // G operator
  Epetra_SerialDenseMatrix* Z;  // Z operator
  Epetra_SerialDenseMatrix* FCF;  // FCF^T
  Epetra_SerialDenseMatrix* Kda;  // EAS matrix Kda
  Epetra_SerialDenseMatrix* Kaa;  // EAS matrix Kaa
  Epetra_SerialDenseVector* feas; // EAS portion of internal forces
  double detJ0;  // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history

/*
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
  G.Shape(4,Wall1::neas_);
  Epetra_SerialDenseMatrix Z;        // Z operator
  Z.Shape(2*numnode,Wall1::neas_);
  Epetra_SerialDenseMatrix FCF;     // FCF^T
  FCF.Shape(4,4);
  Epetra_SerialDenseMatrix Kda;     // EAS matrix Kda
  Kda.Shape(2*numnode,Wall1::neas_);
  Epetra_SerialDenseMatrix Kaa;     // EAS matrix Kaa
  Kaa.Shape(Wall1::neas_,Wall1::neas_);
  Epetra_SerialDenseVector feas;    // EAS portion of internal forces
  feas.Size(Wall1::neas_);
  double detJ0;                     // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history
*/

  // ------------------------------------ check calculation of mass matrix
  double density = 0.0;
  if (massmatrix) density = Density(material);

  /*------- get integraton data ---------------------------------------- */
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule_);

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
    xcure(0,k) = xrefe(0,k) + disp[k*numdf+0];
    xcure(1,k) = xrefe(1,k) + disp[k*numdf+1];
  }

  if (iseas_ == true)
  {
    // allocate EAS quantities
    F_enh = new Epetra_SerialDenseMatrix(4,1);
    F_tot = new Epetra_SerialDenseMatrix(4,3);
    p_stress = new Epetra_SerialDenseMatrix(4,1);
    xjm0 = new Epetra_SerialDenseMatrix(2,2);
    F0 = new Epetra_SerialDenseVector(4);
    boplin0 = new Epetra_SerialDenseMatrix(4,2*numnode);
    W0 = new Epetra_SerialDenseMatrix(4,2*numnode);
    G = new Epetra_SerialDenseMatrix(4,Wall1::neas_);
    Z = new Epetra_SerialDenseMatrix(2*numnode,Wall1::neas_);
    FCF = new Epetra_SerialDenseMatrix(4,4);
    Kda = new Epetra_SerialDenseMatrix(2*numnode,Wall1::neas_);
    Kaa = new Epetra_SerialDenseMatrix(Wall1::neas_,Wall1::neas_);
    feas = new Epetra_SerialDenseVector(Wall1::neas_);

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
    w1_eassetup(*boplin0,*F0,*xjm0,detJ0,xrefe,xcure,distype);
  }

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
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,numnode);

    /*------------------------------------ integration factor  -------*/
    double fac = wgt * det * thickness_;

    /*------------------------------compute mass matrix if imass-----*/
    if (massmatrix)
    {
      double facm = fac * density;
      for (int a=0; a<numnode; a++)
      {
        for (int b=0; b<numnode; b++)
        {
          (*massmatrix)(2*a,2*b)     += facm * funct(a) * funct(b); /* a,b even */
          (*massmatrix)(2*a+1,2*b+1) += facm * funct(a) * funct(b); /* a,b odd  */
        }
      }
    }

    /*----------------------------------- calculate operator Blin  ---*/
    w1_boplin(boplin,deriv,xjm,det,numnode);
    //cout.precision(16);
    /*------------ calculate defgrad F^u, Green-Lagrange-strain E^u --*/
    w1_defgrad(F,strain,xrefe,xcure,boplin,numnode);

    /*-calculate defgrad F in matrix notation and Blin in current conf.*/
    w1_boplin_cure(b_cure,boplin,F,numeps,nd);

    // EAS technology: "enhance the deformation gradient"  ---- --- EAS
    if (iseas_ == true)
    {
      /*-----calculate the enhanced deformation gradient and--------------------
      -----alsoe the operators G, W0 and Z------------------------------------*/

      w1_call_defgrad_enh(*F_enh,*xjm0,xjm,detJ0,det,*F0,*alpha,e1,e2,*G,*W0,*boplin0,*Z);

      /*-----total deformation gradient, Green-Lagrange-strain E^F -----------*/
      w1_call_defgrad_tot(*F_enh,*F_tot,F,strain);
      /* call material law----------------------------------------------------*/
      w1_call_matgeononl(strain,stress,C,numeps,material);

       // return gp strains (only in case of stress/strain output)
       if (elestrain != NULL)
       {
         for (int i = 0; i < Wall1::numstr_; ++i)
           (*elestrain)(ip,i) = strain(i);
       }

       // return gp stresses (only in case of stress/strain output)
       if (elestress != NULL)
       {
         if (cauchy)
         {
           StressCauchy(ip, (*F_tot)(0,0), (*F_tot)(1,1), (*F_tot)(0,2), (*F_tot)(1,2), stress, elestress);
         }
         else
         {
           (*elestress)(ip,0) = stress(0,0);
           (*elestress)(ip,1) = stress(1,1);
           (*elestress)(ip,2) = stress(0,2);
         }
       }

      /*-----first piola-kirchhoff stress vector------------------------------*/
      w1_stress_eas(stress,(*F_tot),(*p_stress));

      /*-----stiffness matrix kdd---------------------------------------------*/
      if (stiffmatrix) w1_kdd(boplin,*W0,*F_tot,C,stress,*FCF,*stiffmatrix,fac);
      /*-----matrix kda-------------------------------------------------------*/
      w1_kda(*FCF,*W0,boplin,stress,*G,*Z,*Kda,*p_stress,fac);
      /*-----matrix kaa-------------------------------------------------------*/
      w1_kaa(*FCF,stress,(*G),(*Kaa),fac);
      /*-----nodal forces ----------------------------------------------------*/
      if (force) w1_fint_eas(*W0,boplin,*G,*p_stress,*force,*feas,fac);

   }
   else
   {
     w1_call_matgeononl(strain,stress,C,numeps,material);

     // return gp strains (only in case of stress/strain output)
     if (elestress != NULL)
     {
       for (int i = 0; i < Wall1::numstr_; ++i)
         (*elestrain)(ip,i) = strain(i);
     }

     // return gp stresses (only in case of stress/strain output)
     if (elestress != NULL)
     {
       if (cauchy)
       {
         StressCauchy(ip, F[0], F[1], F[2], F[3], stress, elestress);
       }
       else
       {
         (*elestress)(ip,0) = stress(0,0);
         (*elestress)(ip,1) = stress(1,1);
         (*elestress)(ip,2) = stress(0,2);
       }
     }

     /*---------------------- geometric part of stiffness matrix kg ---*/
     if (stiffmatrix) w1_kg(*stiffmatrix,boplin,stress,fac,nd,numeps);

     /*------------------ elastic+displacement stiffness matrix keu ---*/
     if (stiffmatrix) w1_keu(*stiffmatrix,b_cure,C,fac,nd,numeps);

     /*--------------- nodal forces fi from integration of stresses ---*/
     if (force) w1_fint(stress,b_cure,*force,fac,nd);
   }

  } // for (int ip=0; ip<totngp; ++ip)


  // EAS technology: ------------------------------------------------------ EAS
  // subtract EAS matrices from disp-based Kdd to "soften" element

  if (force != NULL && stiffmatrix != NULL) 
  {
    if (iseas_ == true)
    {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(*Kaa);
      solve_for_inverseKaa.Invert();


      Epetra_SerialDenseMatrix KdaKaa(2*NumNode(),Wall1::neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('N', 'N', 1.0, (*Kda), (*Kaa), 1.0);


      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
      if (stiffmatrix) (*stiffmatrix).Multiply('N', 'T', -1.0, KdaKaa, *Kda, 1.0);

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, (*feas), 1.0);

      // store current EAS data in history
      for (int i=0; i<Wall1::neas_; ++i)
        for (int j=0; j<Wall1::neas_; ++j)
          (*oldKaainv)(i,j) = (*Kaa)(i,j);

      for (int i=0; i<(2*NumNode()); ++i)
        for (int j=0; j<Wall1::neas_; ++j)
        {
          (*oldKda)(i,j) = (*Kda)(i,j);
          (*oldfeas)(j,0) = (*feas)(j);
        }

    }
  }
  // -------------------------------------------------------------------- EAS

  // clean EAS data
  if (iseas_)
  {
    delete F_enh;
    delete F_tot;
    delete p_stress;
    delete xjm0;
    delete F0;
    delete boplin0;
    delete W0;
    delete G;
    delete Z;
    delete FCF;
    delete Kda;
    delete Kaa;
    delete feas;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  jacobian matrix (private)                                  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_jacobianmatrix(
  const Epetra_SerialDenseMatrix& xrefe,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix& xjm,
  double* det,
  const int iel
)
{

   memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double));

   for (int k=0; k<iel; k++)
   {
        xjm(0,0) += deriv(0,k) * xrefe(0,k);
        xjm(0,1) += deriv(0,k) * xrefe(1,k);
        xjm(1,0) += deriv(1,k) * xrefe(0,k);
        xjm(1,1) += deriv(1,k) * xrefe(1,k);
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
                           const Epetra_SerialDenseMatrix& xrefe,
                           const Epetra_SerialDenseMatrix& xcure,
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
     F[0] += boplin(0,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_11
     F[1] += boplin(1,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_22
     F[2] += boplin(2,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_12
     F[3] += boplin(3,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_21
  } /* end of loop over nodes */

  /*-----------------------calculate Green-Lagrange strain E -------------*/
  strain[0] = 0.5 * (F[0] * F[0] + F[3] * F[3] - 1.0);  // E_11
  strain[1] = 0.5 * (F[2] * F[2] + F[1] * F[1] - 1.0);  // E_22
  strain[2] = 0.5 * (F[0] * F[2] + F[3] * F[1]);        // E_12
  strain[3] = strain[2];                                // E_21

  /*-----------------------linear engineering strain eps -----------------*/
  /* (choose 2PK stresses for stress output, when using linear strains!)  */
  //strain[0] = 0.5 * (F[0] + F[0]) - 1.0;
  //strain[1] = 0.5 * (F[1] + F[1]) - 1.0;
  //strain[2] = 0.5 * (F[2] + F[3]);
  //strain[3] = strain[2];

  return;
}

/* DRT::ELEMENTS::Wall1::w1_defgrad */


/*----------------------------------------------------------------------*
 | Deformation gradient F in matrix notation and B in
 reference configuration (private)                             mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_boplin_cure(Epetra_SerialDenseMatrix& b_cure,
                                          const Epetra_SerialDenseMatrix& boplin,
                                          const Epetra_SerialDenseVector& F,
                                          const int numeps,
                                          const int nd)
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


//{
//  RefCountPtr<MAT::Material> mat = Material();
//  Epetra_SerialDenseMatrix cmat;
//  
//  switch(material->mattyp)
//  {
//    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
//    {
//      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
//
//      stvk->Evaluate(glstrain,cmat,stress);
//
//      *density = stvk->Density();
//
//      break;
//    }
//    default:
//      dserror("Illegal type %d of material for wall1 element ", mat->MaterialType());
//      break;
//  }
//
//  /*--------------------------------------------------------------------*/
//  return;
//}  // of w1_mat_sel

/*----------------------------------------------------------------------*
| geometric stiffness part (total lagrange)                   mgit 05/07|
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_kg(Epetra_SerialDenseMatrix& estif,
                                 const Epetra_SerialDenseMatrix& boplin,
                                 const Epetra_SerialDenseMatrix& stress,
                                 const double fac,
                                 const int nd,
                                 const int numeps)
{
  /*---------------------------------------------- perform B^T * SIGMA * B*/
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
      for(int r=0; r<numeps; r++)
         for(int m=0; m<numeps; m++)
            estif(i,j) += boplin(r,i)*stress(r,m)*boplin(m,j)*fac;

  return;

}  // DRT::ELEMENTS::Wall1::w1_kg

/*----------------------------------------------------------------------*
| elastic and initial displacement stiffness (total lagrange)  mgit 05/07
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_keu(Epetra_SerialDenseMatrix& estif,
                                  const Epetra_SerialDenseMatrix& b_cure,
                                  const Epetra_SerialDenseMatrix& C,
                                  const double fac,
                                  const int nd,
                                  const int numeps)
{

  /*------------- perform B_cure^T * D * B_cure, whereas B_cure = F^T * B */
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
        for(int k=0; k<numeps; k++)
           for(int m=0; m<numeps; m++)
             estif(i,j) +=  b_cure(k,i)*C(k,m)*b_cure(m,j)*fac;

  return;
}  // DRT::ELEMENTS::Wall1::w1_keu


/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) mgit 05/07  |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_fint(const Epetra_SerialDenseMatrix& stress,
                                   const Epetra_SerialDenseMatrix& b_cure,
                                   Epetra_SerialDenseVector& intforce,
                                   const double fac,
                                   const int nd)

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
}  // DRT::ELEMENTS::Wall1::w1_fint


/*-----------------------------------------------------------------------------*
| lump mass matrix                                                  bborn 07/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;  
      for (int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}  // w1_lumpmass


/*-----------------------------------------------------------------------------*
| deliver density                                                   bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1::Density(
  const struct _MATERIAL* material
)
{
  // switch material type
  switch (material->mattyp)
  {
  case m_stvenant :  // linear elastic
    return material->m.stvenant->density;
    break;
  case m_neohooke : // kompressible neo-Hooke
    return material->m.neohooke->density;
    break;
  case m_stvenpor :  //porous linear elastic
    return material->m.stvenpor->density;
    break;
  case m_pl_mises: // von Mises material law
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  case m_pl_mises_3D: // Stefan's von mises 3D material law (certainly not Stefan Lenz's law)
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  case m_pl_dp :  // Drucker-Prager material law
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  default:
    dserror("Illegal typ of material for this element");
    return 0;
    break;
  }
}  // Density

/*-----------------------------------------------------------------------------*
| deliver Cauchy stress                                             bborn 08/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::StressCauchy(
  const int ip,
  const double& F11,
  const double& F22,
  const double& F12,
  const double& F21,
  const Epetra_SerialDenseMatrix& stress,
  Epetra_SerialDenseMatrix* elestress
)
{
  // Question: Is this true for plane stress and/or plane strain mode?

  double detf = F11*F22 - F12*F21;
  // Def.grad. tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix defgrad(2,2);
  defgrad(0,0) = F11;
  defgrad(0,1) = F12;
  defgrad(1,0) = F21;
  defgrad(1,1) = F22;
  // PK2 stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix pk2stress(2,2);
  pk2stress(0,0) = stress(0,0);
  pk2stress(0,1) = stress(0,2);
  pk2stress(1,0) = stress(0,2);
  pk2stress(1,1) = stress(1,1);

  // PK1 stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix pk1stress(2,2);
  pk1stress.Multiply('N','T',1.0/detf,pk2stress,defgrad,0.0);

  // Cauchy stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix cauchystress(2,2);
  cauchystress.Multiply('N','N',1.0,defgrad,pk1stress,0.0);

  // copy results to array for output
  (*elestress)(ip,0) = cauchystress(0,0);
  (*elestress)(ip,1) = cauchystress(1,1);
  (*elestress)(ip,2) = cauchystress(0,1);
}  // StressCauchy


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
