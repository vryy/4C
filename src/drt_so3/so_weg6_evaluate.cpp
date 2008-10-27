/*!----------------------------------------------------------------------
\file so_weg6_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_weg6.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//#include "Epetra_SerialDenseSolver.h"

// inverse design object
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#include "inversedesign.H"
#endif

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_weg6::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6> elemat1(elemat1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6> elemat2(elemat2_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,          1> elevec1(elevec1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,          1> elevec2(elevec2_epetra.A(),true);
  // start with "none"
  DRT::ELEMENTS::So_weg6::ActionType act = So_weg6::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_weg6::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_weg6::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_weg6::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_weg6::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_weg6::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_weg6::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_weg6::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_weg6::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_weg6::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = So_weg6::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = So_weg6::calc_struct_reset_istep;
  else if (action=="postprocess_stress")        act = So_weg6::postprocess_stress;
#ifdef PRESTRESS
  else if (action=="calc_struct_prestress_update") act = So_weg6::prestress_update;
#endif
#ifdef INVERSEDESIGNCREATE
  else if (action=="calc_struct_inversedesign_update")            act = So_weg6::inversedesign_update;
#endif
  else dserror("Unknown type of action for So_weg6");

  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      sow6_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
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
#ifndef INVERSEDESIGNCREATE
      sow6_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
#else
      invdesign_->sow6_nlnstiffmass(this,lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
#endif
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
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6> myemat(true); // set to zero
      sow6_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
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
#ifndef INVERSEDESIGNCREATE
      sow6_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
#else
      invdesign_->sow6_nlnstiffmass(this,lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
#endif
    }
    break;

    // evaluate stresses and strains at gauss points
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
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6> stress;
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6> strain;
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      bool ea = (iostrain == "euler_almansi");
#ifndef INVERSEDESIGNCREATE
      sow6_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,ea);
#else
      invdesign_->sow6_nlnstiffmass(this,lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,ea);
#endif
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {

      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMSTR_WEG6> nodalstresses;
        soweg6_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_WEG6);

        for (int i=0;i<NUMNOD_WEG6;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz")
      {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_WEG6; ++i) {
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_WEG6; ++j) {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_WEG6;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMSTR_WEG6> nodalstresses;
        soweg6_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_WEG6);

        for (int i=0;i<NUMNOD_WEG6;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_WEG6; ++i) {
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_WEG6; ++j) {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_WEG6;
          }
        }
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_update_imrlike:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_reset_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

#ifdef PRESTRESS
    // in case of prestressing, make a snapshot of the current green-Lagrange strains and add them to
    // the previously stored GL strains in an incremental manner
    case prestress_update:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // build def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_WEG6,9);
      DefGradient(mydisp,gpdefgrd,*prestress_);

      // update deformation gradient and put back to storage
      LINALG::FixedSizeSerialDenseMatrix<3,3> deltaF;
      LINALG::FixedSizeSerialDenseMatrix<3,3> Fhist;
      LINALG::FixedSizeSerialDenseMatrix<3,3> Fnew;
      for (int gp=0; gp<NUMGPT_WEG6; ++gp)
      {
        prestress_->StoragetoMatrix(gp,deltaF,gpdefgrd);
        prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());
        Fnew.Multiply(deltaF,Fhist);
        prestress_->MatrixtoStorage(gp,Fnew,prestress_->FHistory());
      }

      // push-forward invJ for every gaussian point
      UpdateJacobianMapping(mydisp,*prestress_);
    }
    break;
#endif

#ifdef INVERSEDESIGNCREATE
    case inversedesign_update:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      invdesign_->sow6_StoreMaterialConfiguration(this,mydisp);
      invdesign_->IsInit() = true; // this is to make the restart work
    }
    break;
#endif

    default:
      dserror("Unknown type of action for Solid3");
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_weg6::EvaluateNeumann(ParameterList&           params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  dserror("Body force of wedge6 not implemented");
  return 0;
}

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::InitJacobianMapping()
{
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMGPT_WEG6>* shapefct;
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6*NUMDIM_WEG6,NUMNOD_WEG6>* deriv;
/* pointer to (static) weight factors at each gp */
  LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,1>* weights;
  sow6_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate

  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xrefe;
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_WEG6);
  detJ_.resize(NUMGPT_WEG6);
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMGPT_WEG6> deriv_gp;
  for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  {
    // get submatrix of deriv at actual gp
    for (int m=0; m<NUMDIM_WEG6; ++m)
      for (int n=0; n<NUMGPT_WEG6; ++n)
        deriv_gp(m,n)=(*deriv)(NUMDIM_WEG6*gp+m,n);

    invJ_[gp].Multiply(deriv_gp,xrefe);
    detJ_[gp] = invJ_[gp].Invert();

#ifdef PRESTRESS
    if (!(prestress_->IsInit()))
      prestress_->MatrixtoStorage(gp,invJ_[gp],prestress_->JHistory());
#endif
#ifdef INVERSEDESIGNUSE
    if (!(invdesign_->IsInit()))
    {
      invdesign_->MatrixtoStorage(gp,invJ_[gp],invdesign_->JHistory());
      invdesign_->DetJHistory()[gp] = detJ_[gp];
    }
#endif
  }
#ifdef PRESTRESS
  prestress_->IsInit() = true;
#endif
#ifdef INVERSEDESIGNUSE
  invdesign_->IsInit() = true;
#endif
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* stiffmatrix,    // element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* massmatrix,     // element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,1>* force,          // element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestress,      // element stresses
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestrain,      // strains at GP
      ParameterList&            params,         // algorithmic parameters e.g. time
      const bool                cauchy,         // stress output option
      const bool                euler_almansi)  // strain output option
{

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> > shapefcts = sow6_shapefcts();
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> > derivs = sow6_derivs();
  const static vector<double> gpweights = sow6_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xrefe;  // material coord. of element
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xcurr;  // current  coord. of element
#if defined(PRESTRESS) || defined(POSTSTRESS)
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xdisp;
#endif
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_WEG6+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_WEG6+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_WEG6+2];

#if defined(PRESTRESS) || defined(POSTSTRESS)
    xdisp(i,0) = disp[i*NODDOF_WEG6+0];
    xdisp(i,1) = disp[i*NODDOF_WEG6+1];
    xdisp(i,2) = disp[i*NODDOF_WEG6+2];
#endif
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> N_XYZ;
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]);
    double detJ = detJ_[gp];

    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> defgrd(false);
#if defined(PRESTRESS) || defined(POSTSTRESS)
    {
      // get Jacobian mapping wrt to the stored configuration
      LINALG::FixedSizeSerialDenseMatrix<3,3> invJdef;
      prestress_->StoragetoMatrix(gp,invJdef,prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> N_xyz;
      N_xyz.Multiply(invJdef,derivs[gp]);

      // build multiplicative incremental defgrd
      defgrd.MultiplyTT(xdisp,N_xyz);
      defgrd(0,0) += 1.0;
      defgrd(1,1) += 1.0;
      defgrd(2,2) += 1.0;

      // get stored old incremental F
      LINALG::FixedSizeSerialDenseMatrix<3,3> Fhist;
      prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());

      // build total defgrd = delta F * F_old
      LINALG::FixedSizeSerialDenseMatrix<3,3> Fnew;
      Fnew.Multiply(defgrd,Fhist);
      defgrd = Fnew;
    }
#else
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);
#endif

#ifdef INVERSEDESIGNUSE
    {
      // make the multiplicative update so that defgrd refers to
      // the reference configuration that resulted from the inverse
      // design analysis
      LINALG::FixedSizeSerialDenseMatrix<3,3> Fhist;
      invdesign_->StoragetoMatrix(gp,Fhist,invdesign_->FHistory());
      LINALG::FixedSizeSerialDenseMatrix<3,3> tmp3x3;
      tmp3x3.Multiply(defgrd,Fhist);
      defgrd = tmp3x3;  // defgrd is still a view to defgrd_epetra

      // make detJ and invJ refer to the ref. configuration that resulted from
      // the inverse design analysis
      detJ = invdesign_->DetJHistory()[gp];
      invdesign_->StoragetoMatrix(gp,tmp3x3,invdesign_->JHistory());
      N_XYZ.Multiply(tmp3x3,derivs[gp]);
    }
#endif

    // Right Cauchy-Green tensor = F^T * F
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::FixedSizeSerialDenseMatrix<6,1> glstrain(false);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    // return gp strains (only in case of stress/strain output)
    if (elestrain != NULL)
    {
      if (!euler_almansi)
      {
        for (int i = 0; i < 3; ++i)
          (*elestrain)(gp,i) = glstrain(i);
        for (int i = 3; i < 6; ++i)
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> gl;
        gl(0,0) = glstrain(0);
        gl(0,1) = 0.5*glstrain(3);
        gl(0,2) = 0.5*glstrain(5);
        gl(1,0) = gl(0,1);
        gl(1,1) = glstrain(1);
        gl(1,2) = 0.5*glstrain(4);
        gl(2,0) = gl(0,2);
        gl(2,1) = gl(1,2);
        gl(2,2) = glstrain(2);

        // inverse of deformation gradient
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> invdefgrd; // make a copy here otherwise defgrd is destroyed!
        invdefgrd.Invert(defgrd);

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> temp;
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> euler_almansi;
        temp.Multiply(gl,invdefgrd);
        euler_almansi.MultiplyTN(invdefgrd,temp);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }

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
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,NUMDOF_WEG6> bop;
    for (int i=0; i<NUMNOD_WEG6; ++i)
    {
      bop(0,NODDOF_WEG6*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_WEG6*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_WEG6*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_WEG6*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_WEG6*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_WEG6*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_WEG6*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_WEG6*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_WEG6*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_WEG6*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_WEG6*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_WEG6*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_WEG6*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_WEG6*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_WEG6*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_WEG6*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_WEG6*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_WEG6*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,NUMSTR_WEG6> cmat(true);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,1> stress(true);
    double density = 0.0;
    sow6_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL)
    {
      if (!cauchy)
        for (int i = 0; i < NUMSTR_WEG6; ++i)
          (*elestress)(gp,i) = stress(i);
      else
      {                               // return Cauchy stresses
        const double detF = defgrd.Determinant();

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> pkstress;
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> temp;
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> cauchystress;
        temp.Multiply(1.0/detF,defgrd,pkstress,0.);
        cauchystress.MultiplyNT(temp,defgrd);

        (*elestress)(gp,0) = cauchystress(0,0);
        (*elestress)(gp,1) = cauchystress(1,1);
        (*elestress)(gp,2) = cauchystress(2,2);
        (*elestress)(gp,3) = cauchystress(0,1);
        (*elestress)(gp,4) = cauchystress(1,2);
        (*elestress)(gp,5) = cauchystress(0,2);
      }
    }

    const double detJ_w = detJ * gpweights[gp];
    if (force != NULL && stiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w,bop,stress,1.0);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,NUMDOF_WEG6> cb;
      cb.Multiply(cmat,bop);          // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,1> sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ_w);     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(NUMDIM_WEG6);     // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_WEG6; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0,inod) + sfac(3) * N_XYZ(1,inod) + sfac(5) * N_XYZ(2,inod);
        SmB_L[1] = sfac(3) * N_XYZ(0,inod) + sfac(1) * N_XYZ(1,inod) + sfac(4) * N_XYZ(2,inod);
        SmB_L[2] = sfac(5) * N_XYZ(0,inod) + sfac(4) * N_XYZ(1,inod) + sfac(2) * N_XYZ(2,inod);
        for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod)
        {
          double bopstrbop = 0.0;            // intermediate value
          for (int idim=0; idim<NUMDIM_WEG6; ++idim)
            bopstrbop += N_XYZ(idim,jnod) * SmB_L[idim];
          (*stiffmatrix)(NUMDIM_WEG6*inod+0,NUMDIM_WEG6*jnod+0) += bopstrbop;
          (*stiffmatrix)(NUMDIM_WEG6*inod+1,NUMDIM_WEG6*jnod+1) += bopstrbop;
          (*stiffmatrix)(NUMDIM_WEG6*inod+2,NUMDIM_WEG6*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness ******************************
    }

    if (massmatrix != NULL)
    { // evaluate mass matrix +++++++++++++++++++++++++
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_WEG6; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod)
        {
          massfactor = ifactor * shapefcts[gp](jnod);     // intermediate factor
          (*massmatrix)(NUMDIM_WEG6*inod+0,NUMDIM_WEG6*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+1,NUMDIM_WEG6*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+2,NUMDIM_WEG6*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Wedge6 Shape fcts at all 6 Gauss Points           maf 09/08|
 *----------------------------------------------------------------------*/
const vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> > DRT::ELEMENTS::So_weg6::sow6_shapefcts()
{
  vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> > shapefcts(NUMGPT_WEG6);
  // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
  // fill up nodal f at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_wedge_6point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp) {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    DRT::UTILS::shape_function_3D(shapefcts[igp], r, s, t, wedge6);
  }
  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Wedge6 Shape fct-derivs at all 6 Gauss Points     maf 09/08|
 *----------------------------------------------------------------------*/
const vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> > DRT::ELEMENTS::So_weg6::sow6_derivs()
{
  vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> > derivs(NUMGPT_WEG6);
  // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_wedge_6point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp) {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    DRT::UTILS::shape_function_3D_deriv1(derivs[igp], r, s, t, wedge6);
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Wedge6 Weights at all 6 Gauss Points              maf 09/08|
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_weg6::sow6_weights()
{
  vector<double> weights(NUMGPT_WEG6);
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_wedge_6point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int i = 0; i < NUMGPT_WEG6; ++i) {
    weights[i] = intpoints.qwgt[i];
  }
  return weights;
}


/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex8                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_shapederiv(
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMGPT_WEG6>** shapefct,  // pointer to pointer of shapefct
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMNOD_WEG6>** deriv,     // pointer to pointer of derivs
  LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,1>** weights)   // pointer to pointer of weights
{
  // static matrix objects, kept in memory
  static LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMGPT_WEG6>  f;  // shape functions
  static LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMNOD_WEG6> df;  // derivatives
  static LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,1> weightfactors;   // weights for each gp
  static bool fdf_eval;                      // flag for re-evaluate everything


  if (fdf_eval==true) { // if true f,df already evaluated
    *shapefct = &f; // return adress of static object to target of pointer
    *deriv = &df; // return adress of static object to target of pointer
    *weights = &weightfactors; // return adress of static object to target of pointer
    return;
  }
  else {
    // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
    // fill up nodal f at each gp
    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_wedge_6point;
    const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule_);
    for (int igp = 0; igp < intpoints.nquad; ++igp) {
      const double r = intpoints.qxg[igp][0];
      const double s = intpoints.qxg[igp][1];
      const double t = intpoints.qxg[igp][2];

      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> funct;
      LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6, NUMNOD_WEG6> deriv;
      DRT::UTILS::shape_function_3D(funct, r, s, t, wedge6);
      DRT::UTILS::shape_function_3D_deriv1(deriv, r, s, t, wedge6);
      for (int inode = 0; inode < NUMNOD_WEG6; ++inode) {
        f(inode, igp) = funct(inode);
        df(igp*NUMDIM_WEG6+0, inode) = deriv(0, inode);
        df(igp*NUMDIM_WEG6+1, inode) = deriv(1, inode);
        df(igp*NUMDIM_WEG6+2, inode) = deriv(2, inode);
        weightfactors(igp) = intpoints.qwgt[igp];
      }
    }
    // return adresses of just evaluated matrices
    *shapefct = &f; // return adress of static object to target of pointer
    *deriv = &df; // return adress of static object to target of pointer
    *weights = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = true; // now all arrays are filled statically
  }
  return;
}  // of sow6_shapederiv

/*----------------------------------------------------------------------*
 |  lump mass matrix                                         bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_lumpmass(LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 04/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Sow6Register::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_weg6) continue;
    DRT::ELEMENTS::So_weg6* actele = dynamic_cast<DRT::ELEMENTS::So_weg6*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_weg6* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

#if defined(PRESTRESS) || defined(POSTSTRESS)
/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::DefGradient(const vector<double>& disp,
                                         Epetra_SerialDenseMatrix& gpdefgrd,
                                         DRT::ELEMENTS::PreStress& prestress)
{
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> > shapefcts = sow6_shapefcts();
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> > derivs = sow6_derivs();

  // update element geometry
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xdisp;  // current  coord. of element
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_WEG6+0];
    xdisp(i,1) = disp[i*NODDOF_WEG6+1];
    xdisp(i,2) = disp[i*NODDOF_WEG6+2];
  }

  for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  {
#if 0
    // get submatrix of deriv at actual gp
    LINALG::SerialDenseMatrix deriv_gp(NUMDIM_WEG6,NUMGPT_WEG6);
    for (int m=0; m<NUMDIM_WEG6; ++m)
      for (int n=0; n<NUMGPT_WEG6; ++n)
        deriv_gp(m,n)=(*deriv)(NUMDIM_WEG6*gp+m,n);
#endif

    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::FixedSizeSerialDenseMatrix<3,3> invJdef;
    prestress.StoragetoMatrix(gp,invJdef,prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> N_xyz;
    N_xyz.Multiply(invJdef,derivs[gp]);

    // build defgrd (independent of xrefe!)
    LINALG::FixedSizeSerialDenseMatrix<3,3> defgrd;
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;

    prestress.MatrixtoStorage(gp,defgrd,gpdefgrd);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::UpdateJacobianMapping(
                                            const vector<double>& disp,
                                            DRT::ELEMENTS::PreStress& prestress)
{
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,1> > shapefcts = sow6_shapefcts();
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> > derivs = sow6_derivs();

  // get incremental disp
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6> xdisp;
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_WEG6+0];
    xdisp(i,1) = disp[i*NODDOF_WEG6+1];
    xdisp(i,2) = disp[i*NODDOF_WEG6+2];
  }

  LINALG::FixedSizeSerialDenseMatrix<3,3> invJhist;
  LINALG::FixedSizeSerialDenseMatrix<3,3> invJ;
  LINALG::FixedSizeSerialDenseMatrix<3,3> defgrd;
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> N_xyz;
  LINALG::FixedSizeSerialDenseMatrix<3,3> invJnew;
  for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp,invJhist,prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist,derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd,invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp,invJnew,prestress.JHistory());
  } // for (int gp=0; gp<NUMGPT_WEG6; ++gp)

  return;
}
#endif // #if defined(PRESTRESS) || defined(POSTSTRESS)

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WEG6
