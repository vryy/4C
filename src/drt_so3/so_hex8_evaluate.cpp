/*!----------------------------------------------------------------------
\file so_hex8_evaluate.cpp
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

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_weg6.H"
#include "so_hex8.H"
#include "so_tet10.H"
#include "so_tet4.H"
#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/visconeohooke.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_hex8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_hex8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_hex8::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_hex8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                         act = So_hex8::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")                    act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")                  act = So_hex8::calc_struct_update_imrlike;
  else if (action=="eas_init_multi")                              act = So_hex8::eas_init_multi;
  else if (action=="eas_set_multi")                               act = So_hex8::eas_set_multi;
  else if (action=="calc_homog_stressdens")                       act = So_hex8::calc_homog_stressdens;
  else if (action=="postprocess_stress")                          act = So_hex8::postprocess_stress;
  else if (action=="multi_readrestart")                           act = So_hex8::multi_readrestart;
#ifdef PRESTRESS
  else if (action=="calc_struct_prestress_update_green_lagrange") act = So_hex8::update_gl;
#endif
  else dserror("Unknown type of action for So_hex8");

//  const double time = params.get("total time",-1.0);
//  const double dt = params.get("delta time",-1.0);

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
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
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
      Epetra_SerialDenseMatrix* matptr = NULL;
      if (elemat1.N()) matptr = &elemat1;
      soh8_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params);
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
      Epetra_SerialDenseMatrix myemat(lm.size(),lm.size());
      soh8_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
      if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
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
      Epetra_SerialDenseMatrix stress(NUMGPT_SOH8,NUMSTR_SOH8);
      Epetra_SerialDenseMatrix strain(NUMGPT_SOH8,NUMSTR_SOH8);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      bool ea = (iostrain == "euler_almansi");
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,ea);
      AddtoPack(*stressdata, stress);
#if defined(PRESTRESS) || defined(POSTSTRESS)
      {
        RCP<Epetra_SerialDenseMatrix>& gl = PreStrains();
        if (gl==null)
          dserror("Cannot output prestrains");
        if (gl->M() != strain.M() || gl->N() != strain.N())
          dserror("Mismatch in dimension");
        // the element outputs 0.5* strains[3-5], but we have the computational quantity here
        Epetra_SerialDenseMatrix tmp(*gl);
        for (int i=0; i<NUMGPT_SOH8; ++i)
          for (int j=3; j<6; ++j)
            tmp(i,j) *= 0.5;
        strain += tmp;
      }
#endif
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];

      if (stresstype=="ndxyz") 
      {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOH8,NUMSTR_SOH8);
        soh8_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOH8);

        for (int i=0;i<NUMNOD_SOH8;++i)
        {
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOH8;++i)
        {
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOH8;++i)
        {
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
        const Epetra_BlockMap& elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) 
        {
          for (int i = 0; i < NUMSTR_SOH8; ++i) 
          {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j) 
            {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOH8 * (*gpstress)(j,i);
            }
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") 
      {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOH8,NUMSTR_SOH8);
        soh8_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOH8);

        for (int i=0;i<NUMNOD_SOH8;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOH8;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOH8;++i){
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
          for (int i = 0; i < NUMSTR_SOH8; ++i) 
          {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j) 
            {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOH8 * (*gpstress)(j,i);
            }
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
      // do something with internal EAS, etc parameters
      if (eastype_ != soh8_easnone) 
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;
        blas.COPY((*alphao).M()*(*alphao).N(), (*alpha).A(), (*alphao).A());  // alphao := alpha
      }
      // Update of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_struct_update_imrlike: 
    {
      // do something with internal EAS, etc parameters
      // this depends on the applied solution technique (static, generalised-alpha,
      // or other time integrators)
      if (eastype_ != soh8_easnone) 
      {
        double alphaf = params.get<double>("alpha f", 0.0);  // generalised-alpha TIS parameter alpha_f
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1-alphaf}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS::Epetra_BLAS blas;
        blas.SCAL((*alphao).M()*(*alphao).N(), -alphaf/(1.0-alphaf), (*alphao).A());  // alphao *= -alphaf/(1.0-alphaf)
        blas.AXPY((*alphao).M()*(*alphao).N(), 1.0/(1.0-alphaf), (*alpha).A(), (*alphao).A());  // alphao += 1.0/(1.0-alphaf) * alpha
        blas.COPY((*alpha).M()*(*alpha).N(), (*alphao).A(), (*alpha).A());  // alpha := alphao
      }
      // Update of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_homog_stressdens: 
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8_homog(params, mydisp, myres);
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // allocation and initializiation of these data arrays can only be
    // done in the elements that know the number of EAS parameters
    case eas_init_multi: 
    {
      if (eastype_ != soh8_easnone) 
      {
        soh8_eas_init_multi(params);
      }
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // before any microscale simulation, EAS internal data has to be
    // set accordingly
    case eas_set_multi: 
    {
      if (eastype_ != soh8_easnone) 
      {
        soh8_set_eas_multi(params);
      }
    }
    break;

    // read restart of microscale
    case multi_readrestart: 
    {
      RefCountPtr<MAT::Material> mat = Material();

      if (mat->MaterialType()==m_struct_multiscale)
        soh8_read_restart_multi(params);
    }
    break;

#ifdef PRESTRESS
    // in case of prestressing, make a snapshot of the current green-Lagrange strains and add them to
    // the previously stored GL strains in an incremental manner
    case update_gl:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Epetra_SerialDenseMatrix strain(NUMGPT_SOH8,NUMSTR_SOH8);
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,NULL,&strain,params,false,false);
      // the element outputs 0.5* strains[3-5], but we want the computational quantity here
      for (int i=0; i<NUMGPT_SOH8; ++i)
        for (int j=3; j<6; ++j) strain(i,j) *= 2.0;
      RCP<Epetra_SerialDenseMatrix>& gl = PreStrains();
      if (gl==null) dserror("Prestress array not initialized");
      if (gl->M() != strain.M() || gl->N() != strain.N())
        dserror("Prestress arrauy not initialized");
      (*gl) += strain;
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
int DRT::ELEMENTS::So_hex8::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
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
  // **

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }

  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // compute the Jacobian matrix
    Epetra_SerialDenseMatrix jac(NUMDIM_SOH8,NUMDIM_SOH8);
    jac.Multiply('N','N',1.0,int_hex8.deriv_gp[gp],xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = int_hex8.weights(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
        elevec1[nodid*NUMDIM_SOH8+dim] += int_hex8.shapefct_gp[gp](nodid) * (*onoff)[dim] *\
        								  (*val)[dim] * fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex8::EvaluateNeumann

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::InitJacobianMapping()
{
  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply('N','N',1.0,int_hex8.deriv_gp[gp],xrefe,0.0);
    detJ_[gp] = LINALG::NonsymInverse3x3(invJ_[gp]);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass( 
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      ParameterList&            params,         // algorithmic parameters e.g. time
      const bool                cauchy,         // stress output option
      const bool                euler_almansi)  // strain output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
/* ============================================================================*/

  // update element geometry
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha = NULL;  // EAS alphas
  vector<Epetra_SerialDenseMatrix>* M_GP = NULL;   // EAS matrix M at all GPs
  LINALG::SerialDenseMatrix M;      // EAS matrix M at current GP
  Epetra_SerialDenseVector feas;    // EAS portion of internal forces
  Epetra_SerialDenseMatrix Kaa;     // EAS matrix Kaa
  Epetra_SerialDenseMatrix Kda;     // EAS matrix Kda
  double detJ0;                     // detJ(origin)
  Epetra_SerialDenseMatrix T0invT;  // trafo matrix
  Epetra_SerialDenseMatrix* oldfeas = NULL;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv = NULL; // EAS history
  Epetra_SerialDenseMatrix* oldKda = NULL;    // EAS history
  if (eastype_ != soh8_easnone) {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get alpha of previous iteration
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
    LINALG::SerialDenseVector res_d(NUMDOF_SOH8);
    for (int i = 0; i < NUMDOF_SOH8; ++i) {
      res_d(i) = residual[i];
    }
    // add Kda . res_d to feas
    (*oldfeas).Multiply('N','N',1.0,(*oldKda),res_d,1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);
    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.Size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.Shape(neas_,neas_);

    // EAS matrix K_{d alpha}
    Kda.Shape(neas_,NUMDOF_SOH8);

    // transformation matrix T0, maps M-matrix evaluated at origin
    // between local element coords and global coords
    // here we already get the inverse transposed T0
    T0invT.Shape(NUMSTR_SOH8,NUMSTR_SOH8);

    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    soh8_eassetup(&M_GP,detJ0,T0invT,xrefe);
  } // -------------------------------------------------------------------- EAS


  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::SerialDenseMatrix N_XYZ(NUMDIM_SOH8,NUMNOD_SOH8);
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply('N','N',1.0,invJ_[gp],int_hex8.deriv_gp[gp],0.0);
    const double detJ = detJ_[gp];

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    LINALG::SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
    defgrd.Multiply('T','T',1.0,xcurr,N_XYZ,0.0);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8,NUMDIM_SOH8);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_SOH8);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      M.LightShape(NUMSTR_SOH8,neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      //Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      M.Multiply('N','N',detJ0/detJ,T0invT,M_GP->at(gp),0.0);
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      glstrain.Multiply('N','N',1.0,M,(*alpha),1.0);
    } // ------------------------------------------------------------------ EAS

    // return gp strains (only in case of stress/strain output)
    if (elestrain != NULL){
      if (!euler_almansi) {
        for (int i = 0; i < 3; ++i) {
          (*elestrain)(gp,i) = glstrain(i);
        }
        for (int i = 3; i < 6; ++i) {
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
        }
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::SerialDenseMatrix gl(NUMDIM_SOH8,NUMDIM_SOH8);
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
        Epetra_SerialDenseMatrix invdefgrd(defgrd); // make a copy here otherwise defgrd is destroyed!
        LINALG::NonsymInverse3x3(invdefgrd);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix euler_almansi(NUMDIM_SOH8,NUMDIM_SOH8);
        temp.Multiply('N','N',1.0,gl,invdefgrd,0.);
        euler_almansi.Multiply('T','N',1.0,invdefgrd,temp,0.);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }

#if defined(PRESTRESS) || defined(POSTSTRESS)
    {
      // note: must be AFTER strains are output above!
      RCP<Epetra_SerialDenseMatrix>& gl = PreStrains();
      if (gl==null) dserror("Prestress array not initialized");
      if (gl->M() != NUMGPT_SOH8 || gl->N() != NUMSTR_SOH8)
        dserror("Prestress array not initialized");
      for (int i=0; i<6; ++i)
        glstrain(i) += (*gl)(gp,i);
    }
#endif

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
    LINALG::SerialDenseMatrix bop(NUMSTR_SOH8,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) {
      bop(0,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    double density;
    soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL){                // return 2nd Piola-Kirchhoff stresses
      if (!cauchy) {
        for (int i = 0; i < NUMSTR_SOH8; ++i) {
          (*elestress)(gp,i) = stress(i);
        }
      }
      else {                               // return Cauchy stresses
        double detF = defgrd(0,0)*defgrd(1,1)*defgrd(2,2) +
                      defgrd(0,1)*defgrd(1,2)*defgrd(2,0) +
                      defgrd(0,2)*defgrd(1,0)*defgrd(2,1) -
                      defgrd(0,2)*defgrd(1,1)*defgrd(2,0) -
                      defgrd(0,0)*defgrd(1,2)*defgrd(2,1) -
                      defgrd(0,1)*defgrd(1,0)*defgrd(2,2);

        LINALG::SerialDenseMatrix pkstress(NUMDIM_SOH8,NUMDIM_SOH8);
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix cauchystress(NUMDIM_SOH8,NUMDIM_SOH8);
        temp.Multiply('N','N',1.0/detF,defgrd,pkstress,0.);
        cauchystress.Multiply('N','T',1.0,temp,defgrd,0.);

        (*elestress)(gp,0) = cauchystress(0,0);
        (*elestress)(gp,1) = cauchystress(1,1);
        (*elestress)(gp,2) = cauchystress(2,2);
        (*elestress)(gp,3) = cauchystress(0,1);
        (*elestress)(gp,4) = cauchystress(1,2);
        (*elestress)(gp,5) = cauchystress(0,2);
      }
    }

    if (force != NULL && stiffmatrix != NULL) {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T', 'N', detJ * int_hex8.weights(gp), bop, stress,
          1.0);

      LINALG::SerialDenseMatrix cb(NUMSTR_SOH8, NUMDOF_SOH8);
      cb.Multiply('N', 'N', 1.0, cmat, bop, 0.0); // temporary C . Bl

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      (*stiffmatrix).Multiply('T', 'N', detJ * int_hex8.weights(gp), bop,
          cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ * int_hex8.weights(gp)); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(NUMDIM_SOH8); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(NUMDIM_SOH8*inod+0, NUMDIM_SOH8*jnod+0) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOH8*inod+1, NUMDIM_SOH8*jnod+1) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOH8*inod+2, NUMDIM_SOH8*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone) {
        double integrationfactor = detJ * int_hex8.weights(gp);
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        LINALG::SerialDenseMatrix cM(NUMSTR_SOH8, neas_); // temporary c . M
        cM.Multiply('N', 'N', 1.0, cmat, M, 0.0);
        Kaa.Multiply('T', 'N', integrationfactor, M, cM, 1.0);

        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        Kda.Multiply('T', 'N', integrationfactor, M, cb, 1.0);

        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        feas.Multiply('T', 'N', integrationfactor, M, stress, 1.0);
      } // ---------------------------------------------------------------- EAS
    }

    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate consistent mass matrix
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double massfactor = (int_hex8.shapefct_gp[gp])(inod) * density * (int_hex8.shapefct_gp[gp])(jnod)
                            * detJ * int_hex8.weights(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  if (force != NULL && stiffmatrix != NULL) {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ != soh8_easnone) {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();

      LINALG::SerialDenseMatrix KdaKaa(NUMDOF_SOH8, neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('T', 'N', 1.0, Kda, Kaa, 0.0);

      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      (*stiffmatrix).Multiply('N', 'N', -1.0, KdaKaa, Kda, 1.0);

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      (*force).Multiply('N', 'N', -1.0, KdaKaa, feas, 1.0);

      // store current EAS data in history
      for (int i=0; i<neas_; ++i) {
        for (int j=0; j<neas_; ++j) (*oldKaainv)(i, j) = Kaa(i,j);
        for (int j=0; j<NUMDOF_SOH8; ++j) (*oldKda)(i, j) = Kda(i,j);
        (*oldfeas)(i, 0) = feas(i);
      }
    } // -------------------------------------------------------------------- EAS
  }
  return;
} // DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_lumpmass(Epetra_SerialDenseMatrix* emass)
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
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex8                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_shapederiv(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseMatrix** deriv,     // pointer to pointer of derivs
      Epetra_SerialDenseVector** weights)   // pointer to pointer of weights
{
  // static matrix objects, kept in memory
  static Epetra_SerialDenseMatrix  f(NUMNOD_SOH8,NUMGPT_SOH8);  // shape functions
  static Epetra_SerialDenseMatrix df(NUMDOF_SOH8,NUMNOD_SOH8);  // derivatives
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOH8);   // weights for each gp
  static bool fdf_eval;                      // flag for re-evaluate everything

  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double gpw      = 1.0;              // weight at every gp for linear fct

  if (fdf_eval==true){             // if true f,df already evaluated
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv    = &df;            // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    return;
  } else {
    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
    const double w[NUMGPT_SOH8] = {   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw};

    // fill up nodal f at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
      f(0,i) = (1.0-r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      f(1,i) = (1.0+r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      f(2,i) = (1.0+r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      f(3,i) = (1.0-r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      f(4,i) = (1.0-r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      f(5,i) = (1.0+r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      f(6,i) = (1.0+r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      f(7,i) = (1.0-r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      weightfactors[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
    }

    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
        // df wrt to r "+0" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;

        // df wrt to s "+1" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;

        // df wrt to t "+2" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }

    // return adresses of just evaluated matrices
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv = &df;               // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = true;               // now all arrays are filled statically
  }
  return;
}  // of soh8_shapederiv


/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 04/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Soh8Register::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_hex8) continue;
    DRT::ELEMENTS::So_hex8* actele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}




#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
