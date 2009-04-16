/*!----------------------------------------------------------------------
\file so_nurbs27_evaluate.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nurbs27.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::Evaluate(
  ParameterList&            params        ,
  DRT::Discretization&      discretization,
  vector<int>&              lm            ,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<81,81> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<81,81> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<81, 1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<81, 1> elevec2(elevec2_epetra.A(),true);

  // start with "none"
  DRT::ELEMENTS::NURBS::So_nurbs27::ActionType act = So_nurbs27::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff"     ) act = So_nurbs27::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff"     ) act = So_nurbs27::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_nurbs27::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass" ) act = So_nurbs27::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass" ) act = So_nurbs27::calc_struct_nlnstiffmass;
  else dserror("Unknown type of action for So_nurbs27");
  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      sonurbs27_nlnstiffmass(lm            ,
			     discretization,
			     mydisp        ,
			     myres         ,
			     &elemat1      ,
			     NULL          ,
			     &elevec1      ,
			     params        );
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<81,81>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      sonurbs27_nlnstiffmass(lm            ,
			     discretization,
			     mydisp        ,
			     myres         ,
			     matptr        ,
			     NULL          ,
			     &elevec1      ,
			     params        );
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<81,81> myemat(true);
      sonurbs27_nlnstiffmass(lm            ,
			     discretization,
			     mydisp        ,
			     myres         ,
			     &myemat       ,
			     NULL          ,
			     &elevec1      ,
			     params        );
    }
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      
      sonurbs27_nlnstiffmass(lm            ,
			     discretization,
			     mydisp        ,
			     myres         ,
			     &elemat1      ,
			     &elemat2      ,
			     &elevec1      ,
			     params        );
    }
    break;

    default:
      dserror("Unknown type of action for So_nurbs27");
  }
  return 0;
} // DRT::ELEMENTS::So_nurbs27::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::EvaluateNeumann(
  ParameterList&            params        ,
  DRT::Discretization&      discretization,
  DRT::Condition&           condition     ,
  vector<int>&              lm            ,
  Epetra_SerialDenseVector& elevec1       )
{
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if(nurbsdis==NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele=false;
  (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());
  
  // there is nothing to be done for zero sized elements in knotspan
  if(zero_ele)
  {
    return(0);
  }

  LINALG::Matrix<27,1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode=0; inode<27; inode++)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
    
    weights(inode) = cp->W();
  }

  /*------------------------------------------------------------------*/
  /*                   update element geometry                        */
  /*------------------------------------------------------------------*/

  // material coord. of element
  LINALG::Matrix<27,3> xrefe;  
  for (int i=0; i<27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  /*------------------------------------------------------------------*/
  /*                 TIME CURVE BUSINESS                              */
  /*------------------------------------------------------------------*/
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

  /*------------------------------------------------------------------*/
  /*               for nurbs_27 with 27 GAUSS POINTS                  */
  /*                  o CONST SHAPE FUNCTIONS                         */
  /*                  o CONST DERIVATIVES                             */
  /*                  o CONST GAUSSWEIGHTS                            */
  /*------------------------------------------------------------------*/
  const static vector<LINALG::Matrix<27, 1> > shapefcts 
    = sonurbs27_shapefcts(myknots,weights);
  const static vector<LINALG::Matrix< 3,27> > derivs    
    = sonurbs27_derivs   (myknots,weights);
  const static vector<double>                 gpweights = sonurbs27_gpweights();

  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp=27;

  for (int gp=0; gp<numgp; ++gp) {

    // compute the Jacobian matrix
    LINALG::Matrix<3,3> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();

    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // integration factor
    double fac = gpweights[gp] * curvefac * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<3; dim++) 
    {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<27; ++nodid) 
      {
	elevec1[nodid*3+dim] += shapefcts[gp](nodid) * dim_fac;
      }
    }

  }/* end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_nurbs27::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::InitJacobianMapping(DRT::Discretization& dis)
{

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

  if(nurbsdis==NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele=false;
  (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());
  
  // there is nothing to be done for zero sized elements in knotspan
  if(zero_ele)
  {
    return;
  }

  LINALG::Matrix<27,1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode=0; inode<27; inode++)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
    
    weights(inode) = cp->W();
  }

  const static vector<LINALG::Matrix<3,27> > derivs 
    = sonurbs27_derivs(myknots,weights);
  LINALG::Matrix<27,3>                       xrefe;
  for (int i=0; i<27; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }

  const int numgp=27;

  invJ_.resize(numgp);
  detJ_.resize(numgp);
  for (int gp=0; gp<numgp; ++gp)
  {
    invJ_[gp].Multiply(derivs[gp],xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0) 
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0) 
      dserror("NEGATIVE JACOBIAN DETERMINANT");

  }
  return;
} // DRT::ELEMENTS::So_nurbs27::InitJacobianMapping()

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_nlnstiffmass(
  vector<int>&           lm            , // location matrix
  DRT::Discretization&   discretization, // discretisation to extract knot vector
  vector<double>&        disp          , // current displacements
  vector<double>&        residual      , // current residual displ
  LINALG::Matrix<81,81>* stiffmatrix   , // element stiffness matrix
  LINALG::Matrix<81,81>* massmatrix    , // element mass matrix
  LINALG::Matrix<81, 1>* force         , // element internal force vector
  ParameterList&         params        ) // strain output option
{

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if(nurbsdis==NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele=false;
  (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());
  
  // there is nothing to be done for zero sized elements in knotspan
  if(zero_ele)
  {
    return;
  }

  LINALG::Matrix<27,1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode=0; inode<27; inode++)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
    
    weights(inode) = cp->W();
  }

  /*------------------------------------------------------------------*/
  /*               for nurbs_27 with 27 GAUSS POINTS                  */
  /*                  o CONST SHAPE FUNCTIONS                         */
  /*                  o CONST DERIVATIVES                             */
  /*                  o CONST GAUSSWEIGHTS                            */
  /*------------------------------------------------------------------*/
  const static vector<LINALG::Matrix<27,1> > shapefcts
    = sonurbs27_shapefcts(myknots,weights);
  const static vector<LINALG::Matrix<3,27> > derivs    
    = sonurbs27_derivs   (myknots,weights);
  const static vector<double>                gpweights = sonurbs27_gpweights();

  // update element geometry
  LINALG::Matrix<27,3> xrefe;  // material coord. of element
  LINALG::Matrix<27,3> xcurr;  // current  coord. of element
  for (int i=0; i<27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*3  ];
    xcurr(i,1) = xrefe(i,1) + disp[i*3+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*3+2];

  }

  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp=27;

  LINALG::Matrix<3,27> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<3,3> defgrd(false);
  for (int gp=0; gp<numgp; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]);
    double detJ = detJ_[gp];

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<3,3> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(6);
    LINALG::Matrix<6,1> glstrain(glstrain_epetra.A(),true);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<6,81> bop;
    for (int i=0; i<27; ++i)
    {
      bop(0,3*i  ) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,3*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,3*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,3*i  ) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,3*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,3*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,3*i  ) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,3*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,3*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,3*i  ) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,3*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,3*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,3*i  ) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,3*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,3*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,3*i  ) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,3*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,3*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law 
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<6,6> cmat  (true);
    LINALG::Matrix<6,1> stress(true);
    sonurbs27_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law 

    double detJ_w = detJ*gpweights[gp];
    if (force != NULL && stiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,81> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ_w); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<27; ++inod) 
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<27; ++jnod) 
	{
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<3; ++idim)
	  {
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
	  }

          (*stiffmatrix)(3*inod  ,3*jnod  ) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness
    }

    if (massmatrix != NULL) // evaluate mass matrix
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<27; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<27; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(3*inod  ,3*jnod  ) += massfactor;
          (*massmatrix)(3*inod+1,3*jnod+1) += massfactor;
          (*massmatrix)(3*inod+2,3*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix 

  }/* end of Loop over GP */

  return;
} // DRT::ELEMENTS::So_nurbs7::sonurbs27_nlnstiffmass

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fcts at all 27 Gauss Points                     |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<27,1> > DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_shapefcts(
  const std::vector<Epetra_SerialDenseVector> & myknots,
  const LINALG::Matrix<27,1>                  & weights
  )
{
  const int numgp=27;

  vector<LINALG::Matrix<27,1> > shapefcts(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up nodal f at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp) 
  {
    LINALG::Matrix<3,1> gp;
    gp(0)=intpoints.qxg[igp][0];
    gp(1)=intpoints.qxg[igp][1];
    gp(2)=intpoints.qxg[igp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct
      (shapefcts[igp]        ,
       gp                    ,
       myknots               ,
       weights               ,
       DRT::Element::nurbs27);
  }
  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fct derivs at all 27 Gauss Points              |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<3,27> > DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_derivs(
  const std::vector<Epetra_SerialDenseVector> & myknots,
  const LINALG::Matrix<27,1>                  & weights
)
{
  const int numgp=27;

  vector<LINALG::Matrix<3,27> > derivs(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp) 
  {
    LINALG::Matrix<3,1> gp;
    gp(0)=intpoints.qxg[igp][0];
    gp(1)=intpoints.qxg[igp][1];
    gp(2)=intpoints.qxg[igp][2];

    LINALG::Matrix<27,1> dummyfct;

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
      (dummyfct              ,
       derivs[igp]           ,
       gp                    ,
       myknots               ,
       weights               ,
       DRT::Element::nurbs27);
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Weights at all 27 Gauss Points                     |         
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_gpweights()
{
  const int numgp=27;

  vector<double> gpweights(numgp);
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int i = 0; i < numgp; ++i) 
  {
    gpweights[i] = intpoints.qwgt[i];
  }
  return gpweights;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::Sonurbs27Register::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_nurbs27) continue;
    DRT::ELEMENTS::NURBS::So_nurbs27* actele = dynamic_cast<DRT::ELEMENTS::NURBS::So_nurbs27*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_nurbs27* failed");
    actele->InitJacobianMapping(dis);
  }
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

