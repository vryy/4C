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

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/visconeohooke.H"

// inverse design object
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#include "inversedesign.H"
#endif

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

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
  else if (action=="calc_struct_reset_istep")                     act = So_hex8::calc_struct_reset_istep;
  else if (action=="eas_init_multi")                              act = So_hex8::eas_init_multi;
  else if (action=="eas_set_multi")                               act = So_hex8::eas_set_multi;
  else if (action=="calc_homog_dens")                             act = So_hex8::calc_homog_dens;
  else if (action=="postprocess_stress")                          act = So_hex8::postprocess_stress;
  else if (action=="multi_readrestart")                           act = So_hex8::multi_readrestart;
#ifdef PRESTRESS
  else if (action=="calc_struct_prestress_update")                act = So_hex8::prestress_update;
#endif
#ifdef INVERSEDESIGNCREATE
  else if (action=="calc_struct_inversedesign_update")            act = So_hex8::inversedesign_update;
#endif
  else dserror("Unknown type of action for So_hex8");
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

#ifndef INVERSEDESIGNCREATE
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;
      soh8_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params);
#else
      Epetra_SerialDenseMatrix* matptr = NULL;
      if (elemat1_epetra.N()) matptr = &elemat1_epetra;
      invdesign_->soh8_nlnstiffmass(this,lm,mydisp,myres,matptr,NULL,&elevec1_epetra,NULL,NULL,params);
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
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8> myemat(true);
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
#ifndef INVERSEDESIGNCREATE
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
#else
      invdesign_->soh8_nlnstiffmass(this,lm,mydisp,myres,&elemat1_epetra,&elemat2_epetra,&elevec1_epetra,NULL,NULL,params);
#endif
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
      if (stressdata==null) dserror("Cannot get 'stress' data");
      if (straindata==null) dserror("Cannot get 'strain' data");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      bool ea = (iostrain == "euler_almansi");
#ifndef INVERSEDESIGNCREATE
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,ea);
#else
      LINALG::SerialDenseMatrix stress(NUMGPT_SOH8,NUMSTR_SOH8);
      LINALG::SerialDenseMatrix strain(NUMGPT_SOH8,NUMSTR_SOH8);
      invdesign_->soh8_nlnstiffmass(this,lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,ea);
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
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMSTR_SOH8> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOH8);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOH8;++i)
        {
          DRT::Node* node = nodes[i];
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
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH8;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMSTR_SOH8> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOH8);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOH8;++i){
          DRT::Node* node=nodes[i];
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
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH8;
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
        switch(eastype_) {
        case DRT::ELEMENTS::So_hex8::soh8_easfull : LINALG::DENSEFUNCTIONS::update<soh8_easfull, 1>(*alphao,*alpha); break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild : LINALG::DENSEFUNCTIONS::update<soh8_easmild, 1>(*alphao,*alpha); break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8: LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alphao,*alpha); break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
        }
      }
      // Update of history for visco material
      RCP<MAT::Material> mat = Material();
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
        switch(eastype_) {
        case DRT::ELEMENTS::So_hex8::soh8_easfull:
          LINALG::DENSEFUNCTIONS::update<soh8_easfull,1>(-alphaf/(1.0-alphaf),*alphao,1.0/(1.0-alphaf),*alpha);
          LINALG::DENSEFUNCTIONS::update<soh8_easfull,1>(*alpha,*alphao);
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild:
          LINALG::DENSEFUNCTIONS::update<soh8_easmild,1>(-alphaf/(1.0-alphaf),*alphao,1.0/(1.0-alphaf),*alpha);
          LINALG::DENSEFUNCTIONS::update<soh8_easmild,1>(*alpha,*alphao);
          break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
          LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(-alphaf/(1.0-alphaf),*alphao,1.0/(1.0-alphaf),*alpha);
          LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alpha,*alphao);
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
        }
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

    case calc_struct_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (eastype_ != soh8_easnone)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        switch(eastype_) {
        case DRT::ELEMENTS::So_hex8::soh8_easfull : LINALG::DENSEFUNCTIONS::update<soh8_easfull,1> (*alpha, *alphao); break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild : LINALG::DENSEFUNCTIONS::update<soh8_easmild,1> (*alpha, *alphao); break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8: LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alpha, *alphao); break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
        }
      }
      // Update of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Reset();
      }
    }
    break;

    case calc_homog_dens:
    {
      soh8_homog(params);
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
    case prestress_update:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // build def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOH8,9);
      DefGradient(mydisp,gpdefgrd,*prestress_);

      // update deformation gradient and put back to storage
      LINALG::SerialDenseMatrix deltaF(3,3);
      LINALG::SerialDenseMatrix Fhist(3,3);
      LINALG::SerialDenseMatrix Fnew(3,3);
      for (int gp=0; gp<NUMGPT_SOH8; ++gp)
      {
        prestress_->StoragetoMatrix(gp,deltaF,gpdefgrd);
        prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());
        Fnew.Multiply('N','N',1.0,deltaF,Fhist,0.0);
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
      invdesign_->soh8_StoreMaterialConfiguration(this,mydisp);
      invdesign_->IsInit() = true; // this is to make the restart work
    }
    break;
#endif


    default:
      dserror("Unknown type of action for So_hex8");
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
   LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // compute the Jacobian matrix
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> jac;
    LINALG::DENSEFUNCTIONS::multiply<NUMDIM_SOH8,NUMNOD_SOH8,NUMDIM_SOH8>(jac.A(),int_hex8.deriv_gp[gp].A(),xrefe.A());

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = int_hex8.weights(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
        elevec1[nodid*NUMDIM_SOH8+dim] += int_hex8.shapefct_gp[gp](nodid) * dim_fac;
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
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;
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
    //invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    LINALG::DENSEFUNCTIONS::multiply<NUMDIM_SOH8,NUMNOD_SOH8,NUMDIM_SOH8>(invJ_[gp].A(),int_hex8.deriv_gp[gp].A(),xrefe.A());
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
void DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix, // element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,1>* force,                 // element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,   // stresses at GP
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,   // strains at GP
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
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
#if defined(PRESTRESS) || defined(POSTSTRESS)
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOH8,NUMDIM_SOH8);
#endif
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];

#if defined(PRESTRESS) || defined(POSTSTRESS)
    xdisp(i,0) = disp[i*NODDOF_SOH8+0];
    xdisp(i,1) = disp[i*NODDOF_SOH8+1];
    xdisp(i,2) = disp[i*NODDOF_SOH8+2];
#endif
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
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> T0invT;  // trafo matrix
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
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    switch(eastype_)
    {
    case DRT::ELEMENTS::So_hex8::soh8_easfull:
      LINALG::DENSEFUNCTIONS::multiply<soh8_easfull,NUMDOF_SOH8,1>(1.0, *oldfeas, 1.0, *oldKda, res_d);
      LINALG::DENSEFUNCTIONS::multiply<soh8_easfull,soh8_easfull,1>(1.0,*alpha,-1.0,*oldKaainv,*oldfeas);
      break;
    case DRT::ELEMENTS::So_hex8::soh8_easmild:
      LINALG::DENSEFUNCTIONS::multiply<soh8_easmild,NUMDOF_SOH8,1>(1.0, *oldfeas, 1.0, *oldKda, res_d);
      LINALG::DENSEFUNCTIONS::multiply<soh8_easmild,soh8_easmild,1>(1.0,*alpha,-1.0,*oldKaainv,*oldfeas);
      break;
    case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
      LINALG::DENSEFUNCTIONS::multiply<soh8_eassosh8,NUMDOF_SOH8,1>(1.0, *oldfeas, 1.0, *oldKda, res_d);
      LINALG::DENSEFUNCTIONS::multiply<soh8_eassosh8,soh8_eassosh8,1>(1.0,*alpha,-1.0,*oldKaainv,*oldfeas);
      break;
    case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
    }
    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.Size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.Shape(neas_,neas_);

    // EAS matrix K_{d alpha}
    Kda.Shape(neas_,NUMDOF_SOH8);


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
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  Epetra_SerialDenseMatrix defgrd_epetra(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(defgrd_epetra.A(),true);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    LINALG::DENSEFUNCTIONS::multiply<NUMDIM_SOH8,NUMDIM_SOH8,NUMNOD_SOH8>
      (N_XYZ.A(), invJ_[gp].A(), int_hex8.deriv_gp[gp].A());
    const double detJ = detJ_[gp];

#if defined(PRESTRESS) || defined(POSTSTRESS)
    {
      // get Jacobian mapping wrt to the stored configuration
      LINALG::SerialDenseMatrix invJdef(3,3);
      prestress_->StoragetoMatrix(gp,invJdef,prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      LINALG::SerialDenseMatrix N_xyz(NUMDIM_SOH8,NUMNOD_SOH8);
      N_xyz.Multiply('N','N',1.0,invJdef,int_hex8.deriv_gp[gp],0.0);

      // build multiplicative incremental defgrd
      defgrd.Multiply('T','T',1.0,xdisp,N_xyz,0.0);
      defgrd(0,0) += 1.0;
      defgrd(1,1) += 1.0;
      defgrd(2,2) += 1.0;

      // get stored old incremental F
      LINALG::SerialDenseMatrix Fhist(3,3);
      prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());

      // build total defgrd = delta F * F_old
      LINALG::SerialDenseMatrix Fnew(3,3);
      Fnew.Multiply('N','N',1.0,defgrd,Fhist,0.0);
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
      LINALG::SerialDenseMatrix Fhist(3,3);
      invdesign_->StoragetoMatrix(gp,Fhist,invdesign_->FHistory());
      LINALG::SerialDenseMatrix tmp3x3(3,3);
      tmp3x3.Multiply('N','N',1.0,defgrd_epetra,Fhist,0.0);
      // this copies into the existing memory, so that the fixed size
      // view still works.
      defgrd_epetra = tmp3x3;

      // make detJ and invJ refer to the ref. configuration that resulted from
      // the inverse design analysis
      detJ = invdesign_->DetJHistory()[gp];
      invdesign_->StoragetoMatrix(gp,tmp3x3,invdesign_->JHistory());
      N_XYZ.Multiply('N','N',1.0,tmp3x3,int_hex8.deriv_gp[gp],0.0);
    }
#endif

    // Right Cauchy-Green tensor = F^T * F
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
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
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch(eastype_)
      {
      case DRT::ELEMENTS::So_hex8::soh8_easfull:
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_easfull>(M.A(), detJ0/detJ, T0invT.A(), (M_GP->at(gp)).A());
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,soh8_easfull,1>(1.0,glstrain.A(),1.0,M.A(),alpha->A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easmild:
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_easmild>(M.A(), detJ0/detJ, T0invT.A(), (M_GP->at(gp)).A());
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,soh8_easmild,1>(1.0,glstrain.A(),1.0,M.A(),alpha->A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_eassosh8>(M.A(), detJ0/detJ, T0invT.A(), (M_GP->at(gp)).A());
        LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,soh8_eassosh8,1>(1.0,glstrain.A(),1.0,M.A(),alpha->A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
      }
    } // ------------------------------------------------------------------ EAS

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
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> gl;
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
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
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
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
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
    Epetra_SerialDenseMatrix cmat_epetra(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress_epetra(NUMSTR_SOH8);
    double density;
    soh8_mat_sel(&stress_epetra,&cmat_epetra,&density,&glstrain_epetra,&defgrd_epetra,gp,params);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(cmat_epetra.A(),true);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> stress(stress_epetra.A(),true);

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL) // return 2nd Piola-Kirchhoff stresses
    {
      if (!cauchy)
      {
        for (int i = 0; i < NUMSTR_SOH8; ++i)
          (*elestress)(gp,i) = stress(i);
      }
      else // return Cauchy stresses
      {
        const double detF = defgrd.Determinant();

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> pkstress;
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);

        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
        LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchystress;
        temp.Multiply(1.0/detF,defgrd,pkstress,0.0);
        cauchystress.MultiplyNT(temp,defgrd);

        (*elestress)(gp,0) = cauchystress(0,0);
        (*elestress)(gp,1) = cauchystress(1,1);
        (*elestress)(gp,2) = cauchystress(2,2);
        (*elestress)(gp,3) = cauchystress(0,1);
        (*elestress)(gp,4) = cauchystress(1,2);
        (*elestress)(gp,5) = cauchystress(0,2);
      }
    }

    double detJ_w = detJ*int_hex8.weights(gp);
    if (force != NULL && stiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
      LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8, NUMDOF_SOH8> cb;
      cb.Multiply(cmat,bop);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ_w); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
      if (eastype_ != soh8_easnone)
      {
        double integrationfactor = detJ * int_hex8.weights(gp);
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        LINALG::SerialDenseMatrix cM(NUMSTR_SOH8, neas_); // temporary c . M
        switch(eastype_)
        {
        case DRT::ELEMENTS::So_hex8::soh8_easfull:
          LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_easfull>(cM.A(), cmat.A(), M.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easfull,NUMSTR_SOH8,soh8_easfull>(1.0, Kaa, integrationfactor, M, cM);
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easfull,NUMSTR_SOH8,NUMDOF_SOH8>(1.0, Kda.A(), integrationfactor, M.A(), cb.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easfull,NUMSTR_SOH8,1>(1.0, feas.A(), integrationfactor, M.A(), stress.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild:
          LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_easmild>(cM.A(), cmat.A(), M.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easmild,NUMSTR_SOH8,soh8_easmild>(1.0, Kaa, integrationfactor, M, cM);
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easmild,NUMSTR_SOH8,NUMDOF_SOH8>(1.0, Kda.A(), integrationfactor, M.A(), cb.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_easmild,NUMSTR_SOH8,1>(1.0, feas.A(), integrationfactor, M.A(), stress.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
          LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_eassosh8>(cM.A(), cmat.A(), M.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,soh8_eassosh8>(1.0, Kaa, integrationfactor, M, cM);
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,NUMDOF_SOH8>(1.0, Kda.A(), integrationfactor, M.A(), cb.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,1>(1.0, feas.A(), integrationfactor, M.A(), stress.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
        }
      } // ---------------------------------------------------------------- EAS
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        ifactor = (int_hex8.shapefct_gp[gp])(inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          massfactor = (int_hex8.shapefct_gp[gp])(jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  if (force != NULL && stiffmatrix != NULL)
  {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ != soh8_easnone)
    {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();

      LINALG::SerialDenseMatrix KdaKaa(NUMDOF_SOH8, neas_); // temporary Kda.Kaa^{-1}
      switch(eastype_)
      {
      case DRT::ELEMENTS::So_hex8::soh8_easfull:
        LINALG::DENSEFUNCTIONS::multiplyTN<NUMDOF_SOH8,soh8_easfull,soh8_easfull>(KdaKaa, Kda, Kaa);
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easfull,NUMDOF_SOH8>(1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easfull,1>(1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easmild:
        LINALG::DENSEFUNCTIONS::multiplyTN<NUMDOF_SOH8,soh8_easmild,soh8_easmild>(KdaKaa, Kda, Kaa);
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easmild,NUMDOF_SOH8>(1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easmild,1>(1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
        LINALG::DENSEFUNCTIONS::multiplyTN<NUMDOF_SOH8,soh8_eassosh8,soh8_eassosh8>(KdaKaa, Kda, Kaa);
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_eassosh8,NUMDOF_SOH8>(1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
        LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_eassosh8,1>(1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easnone: break;
      }

      // store current EAS data in history
      for (int i=0; i<neas_; ++i)
      {
        for (int j=0; j<neas_; ++j)
          (*oldKaainv)(i,j) = Kaa(i,j);
        for (int j=0; j<NUMDOF_SOH8; ++j)
          (*oldKda)(i,j) = Kda(i,j);
        (*oldfeas)(i,0) = feas(i);
      }
    } // -------------------------------------------------------------------- EAS
  }
  return;
} // DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass








/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_lumpmass(LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r=0; r<(*emass).M(); ++r)  // parse rows
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
      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMGPT_SOH8>** shapefct,   // pointer to pointer of shapefct
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMNOD_SOH8>** deriv,     // pointer to pointer of derivs
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,1>** weights)   // pointer to pointer of weights
{
  // static matrix objects, kept in memory
  static LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMGPT_SOH8>  f;  // shape functions
  static LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMNOD_SOH8> df;  // derivatives
  static LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,1> weightfactors;   // weights for each gp
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
      weightfactors(i) = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
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

#if defined(PRESTRESS) || defined(POSTSTRESS)
/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::DefGradient(const vector<double>& disp,
                                         Epetra_SerialDenseMatrix& gpdefgrd,
                                         DRT::ELEMENTS::PreStress& prestress)
{

  // update element geometry
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOH8+0];
    xdisp(i,1) = disp[i*NODDOF_SOH8+1];
    xdisp(i,2) = disp[i*NODDOF_SOH8+2];
  }

  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::SerialDenseMatrix invJdef(3,3);
    prestress.StoragetoMatrix(gp,invJdef,prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::SerialDenseMatrix N_xyz(NUMDIM_SOH8,NUMNOD_SOH8);
    N_xyz.Multiply('N','N',1.0,invJdef,int_hex8.deriv_gp[gp],0.0);

    // build defgrd (independent of xrefe!)
    LINALG::SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
    defgrd.Multiply('T','T',1.0,xdisp,N_xyz,0.0);
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
void DRT::ELEMENTS::So_hex8::UpdateJacobianMapping(
                                            const vector<double>& disp,
                                            DRT::ELEMENTS::PreStress& prestress)
{
  // get incremental disp
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOH8,NUMDIM_SOH8);
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOH8+0];
    xdisp(i,1) = disp[i*NODDOF_SOH8+1];
    xdisp(i,2) = disp[i*NODDOF_SOH8+2];
  }

  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  LINALG::SerialDenseMatrix invJhist(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::SerialDenseMatrix invJ(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::SerialDenseMatrix N_xyz(NUMDIM_SOH8,NUMNOD_SOH8);
  LINALG::SerialDenseMatrix invJnew(NUMDIM_SOH8,NUMDIM_SOH8);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp,invJhist,prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply('N','N',1.0,invJhist,int_hex8.deriv_gp[gp],0.0);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.Multiply('T','T',1.0,xdisp,N_xyz,0.0);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    LINALG::NonsymInverse3x3(defgrd);
    // push-forward of Jinv
    invJnew.Multiply('T','N',1.0,defgrd,invJhist,0.0);
    // store new reference configuration
    prestress.MatrixtoStorage(gp,invJnew,prestress.JHistory());
  } // for (int gp=0; gp<NUMGPT_SOH8; ++gp)

  return;
}
#endif // #if defined(PRESTRESS) || defined(POSTSTRESS)

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

