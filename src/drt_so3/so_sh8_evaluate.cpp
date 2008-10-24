/*!----------------------------------------------------------------------
\file so_sh8_evaluate.cpp
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
#include "so_sh8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_io/io_gmsh.H"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8::Evaluate(ParameterList&            params,
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
  else if (action=="calc_struct_linstiff")        act = So_hex8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")        act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")   act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")    act = So_hex8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")    act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")   act = So_hex8::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")          act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_eleload")         act = So_hex8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")         act = So_hex8::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")    act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")  act = So_hex8::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")     act = So_hex8::calc_struct_reset_istep;
  else if (action=="postprocess_stress")          act = So_hex8::postprocess_stress;
  else if (action=="eas_init_multi")                              act = So_hex8::eas_init_multi;
  else if (action=="eas_set_multi")                               act = So_hex8::eas_set_multi;
  else if (action=="calc_homog_dens")                             act = So_hex8::calc_homog_dens;
  else if (action=="multi_readrestart")                           act = So_hex8::multi_readrestart;
  else dserror("Unknown type of action for So_hex8");

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (Type() == DRT::Element::element_sosh8){
        sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
      } else if (Type() == DRT::Element::element_so_hex8){
        soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
      }
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (Type() == DRT::Element::element_sosh8){
        sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
      } else if (Type() == DRT::Element::element_so_hex8){
        soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params);
      }
    }
    break;

    // internal force vector only
    case calc_struct_internalforce: {
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
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (Type() == DRT::Element::element_sosh8) {
        sosh8_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params);
      } else if (Type() == DRT::Element::element_so_hex8) {
        soh8_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params);
      }
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent/lumped mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (Type() == DRT::Element::element_sosh8){
        sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
      } else if (Type() == DRT::Element::element_so_hex8){
        soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params);
      }
      // lump mass
      if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:{
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
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;

      if (iostrain == "euler_almansi") dserror ("euler_almansi strains not implemented due to missing defgrd");

      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (Type() == DRT::Element::element_sosh8){
        sosh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy);
      } else if (Type() == DRT::Element::element_so_hex8){
        soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,cauchy,false);
      }
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:{

      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz") {
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
      else if (stresstype=="cxyz") {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_SOH8; ++i)
          {
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH8;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
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
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH8;
          }
        }
      }
      else{
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

    case calc_struct_update_istep: {
      // do something with internal EAS, etc parameters
      if (eastype_ == soh8_eassosh8) {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        // alphao := alpha
        LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alphao,*alpha);
      }
      // Update of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_struct_update_imrlike: {
      // do something with internal EAS, etc parameters
      // this depends on the applied solution technique (static, generalised-alpha,
      // or other time integrators)
      if (eastype_ == soh8_eassosh8) {
        double alphaf = params.get<double>("alpha f", 0.0);  // generalised-alpha TIS parameter alpha_f
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1-alphaf}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        // alphao = (-alphaf/(1.0-alphaf))*alphao  + 1.0/(1.0-alphaf) * alpha
        LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(-alphaf/(1.0-alphaf),*alphao,1.0/(1.0-alphaf),*alpha);
        LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alpha,*alphao); // alpha := alphao
      }
      // Update of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_struct_reset_istep: {
      // do something with internal EAS, etc parameters
      if (eastype_ == soh8_eassosh8) {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        // alpha := alphao
        LINALG::DENSEFUNCTIONS::update<soh8_eassosh8,1>(*alpha, *alphao);
      }
      // Reset of history for visco material
      RefCountPtr<MAT::Material> mat = Material();
      if (mat->MaterialType() == m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Reset();
      }
      else if (mat->MaterialType() == m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
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

    default:
      dserror("Unknown type of action for So_sh8");
  }
  return 0;
}




/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::sosh8_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix, // element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,1>* force,                 // element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,   // stresses at GP
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,   // strains at GP
      ParameterList&            params,         // algorithmic parameters e.g. time
      const bool                cauchy)         // stress output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
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
  }

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha = NULL;         // EAS alphas
  vector<Epetra_SerialDenseMatrix>* M_GP = NULL;  // EAS matrix M at all GPs
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,soh8_eassosh8> M; // EAS matrix M at current GP, fixed for sosh8
  Epetra_SerialDenseVector feas;                  // EAS portion of internal forces
  Epetra_SerialDenseMatrix Kaa;                   // EAS matrix Kaa
  Epetra_SerialDenseMatrix Kda;                   // EAS matrix Kda
  double detJ0;                                   // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas = NULL;       // EAS history
  Epetra_SerialDenseMatrix* oldKaainv = NULL;     // EAS history
  Epetra_SerialDenseMatrix* oldKda = NULL;        // EAS history

  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> T0invT;  // trafo matrix

  if (eastype_ == soh8_eassosh8) {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get old alpha
    // evaluate current (updated) EAS alphas (from history variables)
    // get stored EAS history
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
    LINALG::DENSEFUNCTIONS::multiply<soh8_eassosh8, NUMDOF_SOH8,1>(1.0, *oldfeas, 1.0, *oldKda, res_d);
    // "new" alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    LINALG::DENSEFUNCTIONS::multiply<soh8_eassosh8,soh8_eassosh8,1>(1.0,*alpha,-1.0,*oldKaainv,*oldfeas);
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
  } else if (eastype_ == soh8_easnone){
  //cout << "Warning: Solid-Shell8 without EAS" << endl;
  } else dserror("Solid-Shell8 only with eas_sosh8");// ------------------- EAS

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  //LINALG::FixedSizeSerialDenseMatrix<num_ans*num_sp,NUMDOF_SOH8> B_ans_loc(true); //set to 0
  LINALG::FixedSizeSerialDenseMatrix<num_ans*num_sp,NUMDOF_SOH8> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> > jac_sps(num_sp);
  // CURRENT Jacobian evaluated at all ANS sampling points
  vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> > jac_cur_sps(num_sp);
  // pointer to derivs evaluated at all sampling points
  vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> >* deriv_sp = NULL;   //derivs eval. at all sampling points
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe,xcurr,&deriv_sp,jac_sps,jac_cur_sps,B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
  const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> jac_cur;
    jac_cur.Multiply(derivs[gp],xcurr);

    // set up B-Operator in local(parameter) element space including ANS
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMDOF_SOH8> bop_loc;
    for (int inode = 0; inode < NUMNOD_SOH8; ++inode) {
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
        // B_loc_rr = N_r.X_r
        bop_loc(0,inode*3+dim) = derivs[gp](0,inode) * jac_cur(0,dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1,inode*3+dim) = derivs[gp](1,inode) * jac_cur(1,dim);
        // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
        //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
        //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
        bop_loc(2,inode*3+dim) = 0.25*(1-r[gp])*(1-s[gp]) * B_ans_loc(0+4*num_ans,inode*3+dim)
                                +0.25*(1+r[gp])*(1-s[gp]) * B_ans_loc(0+5*num_ans,inode*3+dim)
                                +0.25*(1+r[gp])*(1+s[gp]) * B_ans_loc(0+6*num_ans,inode*3+dim)
                                +0.25*(1-r[gp])*(1+s[gp]) * B_ans_loc(0+7*num_ans,inode*3+dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3,inode*3+dim) = derivs[gp](0,inode) * jac_cur(1,dim)
                                +derivs[gp](1,inode) * jac_cur(0,dim);
        // B_loc_st = interpolation along r of ANS B_loc_st
        //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
        bop_loc(4,inode*3+dim) = 0.5*(1.0+r[gp]) * B_ans_loc(1+1*num_ans,inode*3+dim)
                                +0.5*(1.0-r[gp]) * B_ans_loc(1+3*num_ans,inode*3+dim);
        // B_loc_rt = interpolation along s of ANS B_loc_rt
        //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
        bop_loc(5,inode*3+dim) = 0.5*(1.0-s[gp]) * B_ans_loc(2+0*num_ans,inode*3+dim)
                                +0.5*(1.0+s[gp]) * B_ans_loc(2+2*num_ans,inode*3+dim);
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> TinvT;
    sosh8_evaluateT(jac,TinvT);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
    bop.Multiply(TinvT,bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    // but with modified ANS strains E33, E23 and E13
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> lstrain;
    // evaluate glstrains in local(parameter) coords
    // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
    lstrain(0)= 0.5 * (
       +(jac_cur(0,0)*jac_cur(0,0) + jac_cur(0,1)*jac_cur(0,1) + jac_cur(0,2)*jac_cur(0,2))
       -(jac(0,0)*jac(0,0)         + jac(0,1)*jac(0,1)         + jac(0,2)*jac(0,2)));
    // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
    lstrain(1)= 0.5 * (
       +(jac_cur(1,0)*jac_cur(1,0) + jac_cur(1,1)*jac_cur(1,1) + jac_cur(1,2)*jac_cur(1,2))
       -(jac(1,0)*jac(1,0)         + jac(1,1)*jac(1,1)         + jac(1,2)*jac(1,2)));
    // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
    lstrain(3)= (
       +(jac_cur(0,0)*jac_cur(1,0) + jac_cur(0,1)*jac_cur(1,1) + jac_cur(0,2)*jac_cur(1,2))
       -(jac(0,0)*jac(1,0)         + jac(0,1)*jac(1,1)         + jac(0,2)*jac(1,2)));

    // ANS modification of strains ************************************** ANS
    double dydt_A = 0.0; double dYdt_A = 0.0;
    double dxdt_B = 0.0; double dXdt_B = 0.0;
    double dydt_C = 0.0; double dYdt_C = 0.0;
    double dxdt_D = 0.0; double dXdt_D = 0.0;
    double dzdt_E = 0.0; double dZdt_E = 0.0;
    double dzdt_F = 0.0; double dZdt_F = 0.0;
    double dzdt_G = 0.0; double dZdt_G = 0.0;
    double dzdt_H = 0.0; double dZdt_H = 0.0;

    // vector product of rows of jacobians at corresponding sampling point    cout << jac_cur_sps;
    for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
      dydt_A += jac_cur_sps[0](0,dim) * jac_cur_sps[0](2,dim);
      dYdt_A += jac_sps[0](0,dim)     * jac_sps[0](2,dim);
      dxdt_B += jac_cur_sps[1](1,dim) * jac_cur_sps[1](2,dim);
      dXdt_B += jac_sps[1](1,dim)     * jac_sps[1](2,dim);
      dydt_C += jac_cur_sps[2](0,dim) * jac_cur_sps[2](2,dim);
      dYdt_C += jac_sps[2](0,dim)     * jac_sps[2](2,dim);
      dxdt_D += jac_cur_sps[3](1,dim) * jac_cur_sps[3](2,dim);
      dXdt_D += jac_sps[3](1,dim)     * jac_sps[3](2,dim);

      dzdt_E += jac_cur_sps[4](2,dim) * jac_cur_sps[4](2,dim);
      dZdt_E += jac_sps[4](2,dim)     * jac_sps[4](2,dim);
      dzdt_F += jac_cur_sps[5](2,dim) * jac_cur_sps[5](2,dim);
      dZdt_F += jac_sps[5](2,dim)     * jac_sps[5](2,dim);
      dzdt_G += jac_cur_sps[6](2,dim) * jac_cur_sps[6](2,dim);
      dZdt_G += jac_sps[6](2,dim)     * jac_sps[6](2,dim);
      dzdt_H += jac_cur_sps[7](2,dim) * jac_cur_sps[7](2,dim);
      dZdt_H += jac_sps[7](2,dim)     * jac_sps[7](2,dim);
    }
    // E33: remedy of curvature thickness locking
    // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
    lstrain(2) = 0.5 * (
       0.25*(1-r[gp])*(1-s[gp]) * (dzdt_E - dZdt_E)
      +0.25*(1+r[gp])*(1-s[gp]) * (dzdt_F - dZdt_F)
      +0.25*(1+r[gp])*(1+s[gp]) * (dzdt_G - dZdt_G)
      +0.25*(1-r[gp])*(1+s[gp]) * (dzdt_H - dZdt_H));
    // E23: remedy of transverse shear locking
    // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
    lstrain(4) = 0.5*(1+r[gp]) * (dxdt_B - dXdt_B) + 0.5*(1-r[gp]) * (dxdt_D - dXdt_D);
    // E13: remedy of transverse shear locking
    // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
    lstrain(5) = 0.5*(1-s[gp]) * (dydt_A - dYdt_A) + 0.5*(1+s[gp]) * (dydt_C - dYdt_C);
    // ANS modification of strains ************************************** ANS

    // transformation of local glstrains 'back' to global(material) space
    LINALG::SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
    glstrain.Multiply(TinvT,lstrain);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone) {
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,NUMSTR_SOH8,soh8_eassosh8>(M.A(),detJ0/detJ,T0invT.A(),M_GP->at(gp).A());
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      LINALG::DENSEFUNCTIONS::multiply<NUMSTR_SOH8,soh8_eassosh8,1>(1.0,glstrain.A(),1.0,M.A(),(*alpha).A());
    } // ------------------------------------------------------------------ EAS

    // return gp strains (only in case of stress/strain output)
    if (elestrain != NULL){
      for (int i = 0; i < 3; ++i) {
        (*elestrain)(gp,i) = glstrain(i);
      }
      for (int i = 3; i < 6; ++i) {
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat_epetra(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress_epetra(NUMSTR_SOH8);
    double density;
    // Caution!! the defgrd can not be modified with ANS to remedy locking
    // therefore it is empty and passed only for compatibility reasons
    Epetra_SerialDenseMatrix defgrd_epetra; // Caution!! empty!!
//#define disp_based_F
#ifdef disp_based_F
    defgrd_epetra.Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> invJ;
    invJ.Multiply(derivs[gp],xrefe);
    invJ.Invert();
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ,derivs[gp]);
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd;
    defgrd.MultiplyTT(xcurr,N_XYZ);
    for (int i = 0; i < NUMDIM_SOH8; ++i) {
      for (int j = 0; j < NUMDIM_SOH8; ++j) {
        defgrd_epetra(i,j) = defgrd(i,j);
      }
    }
#endif
    soh8_mat_sel(&stress_epetra,&cmat_epetra,&density,&glstrain_epetra,&defgrd_epetra,gp,params);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(cmat_epetra.A(),true);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> stress(stress_epetra.A(),true);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL){
      if (!cauchy) {                 // return 2nd Piola-Kirchhoff stresses
        for (int i = 0; i < NUMSTR_SOH8; ++i) {
          (*elestress)(gp,i) = stress(i);
        }
      }
      else {                         // return Cauchy stresses
        sosh8_Cauchy(elestress,gp,derivs[gp],xrefe,xcurr,glstrain,stress);
      }
    }

    double detJ_w = detJ*gpweights[gp];
    if (force != NULL && stiffmatrix != NULL) {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8, NUMDOF_SOH8> cb;
      cb.Multiply(cmat,bop); // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // intergrate `geometric' stiffness matrix and add to keu *****************
      // here also the ANS interpolation comes into play
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> G_ij;
          G_ij(0) = derivs[gp](0, inod) * derivs[gp](0, jnod); // rr-dir
          G_ij(1) = derivs[gp](1, inod) * derivs[gp](1, jnod); // ss-dir
          G_ij(3) = derivs[gp](0, inod) * derivs[gp](1, jnod)
                  + derivs[gp](1, inod) * derivs[gp](0, jnod); // rs-dir
          // ANS modification in tt-dir
          G_ij(2) = 0.25*(1-r[gp])*(1-s[gp]) * (*deriv_sp)[4](2,inod) * (*deriv_sp)[4](2,jnod)
                   +0.25*(1+r[gp])*(1-s[gp]) * (*deriv_sp)[5](2,inod) * (*deriv_sp)[5](2,jnod)
                   +0.25*(1+r[gp])*(1+s[gp]) * (*deriv_sp)[6](2,inod) * (*deriv_sp)[6](2,jnod)
                   +0.25*(1-r[gp])*(1+s[gp]) * (*deriv_sp)[7](2,inod) * (*deriv_sp)[7](2,jnod);
          // ANS modification in st-dir
          G_ij(4) = 0.5*((1+r[gp]) * ((*deriv_sp)[1](1,inod) * (*deriv_sp)[1](2,jnod)
                                     +(*deriv_sp)[1](2,inod) * (*deriv_sp)[1](1,jnod))
                        +(1-r[gp]) * ((*deriv_sp)[3](1,inod) * (*deriv_sp)[3](2,jnod)
                                     +(*deriv_sp)[3](2,inod) * (*deriv_sp)[3](1,jnod)));
          // ANS modification in rt-dir
          G_ij(5) = 0.5*((1-s[gp]) * ((*deriv_sp)[0](0,inod) * (*deriv_sp)[0](2,jnod)
                                     +(*deriv_sp)[0](2,inod) * (*deriv_sp)[0](0,jnod))
                        +(1+s[gp]) * ((*deriv_sp)[2](0,inod) * (*deriv_sp)[2](2,jnod)
                                     +(*deriv_sp)[2](2,inod) * (*deriv_sp)[2](0,jnod)));
          // transformation of local(parameter) space 'back' to global(material) space
          LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_SOH8*inod+0, NUMDIM_SOH8*jnod+0) += Gij;
          (*stiffmatrix)(NUMDIM_SOH8*inod+1, NUMDIM_SOH8*jnod+1) += Gij;
          (*stiffmatrix)(NUMDIM_SOH8*inod+2, NUMDIM_SOH8*jnod+2) += Gij;
        }
      } // end of intergrate `geometric' stiffness ******************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone) {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,soh8_eassosh8> cM; // temporary c . M
        cM.Multiply(cmat, M);
        LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,soh8_eassosh8>(1.0, Kaa.A(), detJ_w, M.A(), cM.A());
        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,NUMDOF_SOH8>(1.0, Kda.A(), detJ_w, M.A(), cb.A());
        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        LINALG::DENSEFUNCTIONS::multiplyTN<soh8_eassosh8,NUMSTR_SOH8,1>(1.0, feas.A(), detJ_w, M.A(), stress.A());
      } // ------------------------------------------------------------------ EAS
    }

    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
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
    if (eastype_ == soh8_eassosh8) {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();

      LINALG::SerialDenseMatrix KdaTKaa(NUMDOF_SOH8, soh8_eassosh8); // temporary Kda^T.Kaa^{-1}
      LINALG::DENSEFUNCTIONS::multiplyTN<NUMDOF_SOH8,soh8_eassosh8,soh8_eassosh8>(KdaTKaa, Kda, Kaa);
      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_eassosh8,NUMDOF_SOH8>(1.0, stiffmatrix->A(), -1.0, KdaTKaa.A(), Kda.A());
      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_eassosh8,1>(1.0, force->A(), -1.0, KdaTKaa.A(), feas.A());

      // store current EAS data in history
      for (int i=0; i<soh8_eassosh8; ++i) {
        for (int j=0; j<soh8_eassosh8; ++j) (*oldKaainv)(i,j) = Kaa(i,j);
        for (int j=0; j<NUMDOF_SOH8; ++j) (*oldKda)(i, j) = Kda(i,j);
        (*oldfeas)(i, 0) = feas(i);
      }
    } // -------------------------------------------------------------------- EAS
  }

  return;
} // DRT::ELEMENTS::So_sh8::sosh8_nlnstiffmass



/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::sosh8_anssetup(
          const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8>& xrefe, // material element coords
          const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8>& xcurr, // current element coords
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> >** deriv_sp,   // derivs eval. at all sampling points
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> >& jac_sps,     // jac at all sampling points
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> >& jac_cur_sps, // current jac at all sampling points
          LINALG::FixedSizeSerialDenseMatrix<num_ans*num_sp,NUMDOF_SOH8>& B_ans_loc) // modified B
{
  // static matrix object of derivs at sampling points, kept in memory
  static vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> > df_sp(num_sp);
  static bool dfsp_eval;                      // flag for re-evaluate everything

  if (dfsp_eval!=0){             // if true f,df already evaluated
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
  } else {
  /*====================================================================*/
  /* 8-node hexhedra Solid-Shell node topology
   * and location of sampling points A to H                             */
  /*--------------------------------------------------------------------*/
  /*                      t
   *                      |
   *             4========|================7
   *          // |        |              //||
   *        //   |        |            //  ||
   *      //     |        |   D      //    ||
   *     5=======E=================6       H
   *    ||       |        |        ||      ||
   *    ||   A   |        o--------||-- C -------s
   *    ||       |       /         ||      ||
   *    F        0----- B ---------G ------3
   *    ||     //     /            ||    //
   *    ||   //     /              ||  //
   *    || //     r                ||//
   *     1=========================2
   *
   */
  /*====================================================================*/
    // (r,s,t) gp-locations of sampling points A,B,C,D,E,F,G,H
    // numsp = 8 here set explicitly to allow direct initializing
    double r[8] = { 0.0, 1.0, 0.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double s[8] = {-1.0, 0.0, 1.0, 0.0,-1.0,-1.0, 1.0, 1.0};
    double t[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i=0; i<num_sp; ++i) {
        // df wrt to r "+0" for each node(0..7) at each sp [i]
        df_sp[i](0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        df_sp[i](0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        df_sp[i](0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        df_sp[i](0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        df_sp[i](0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        df_sp[i](0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        df_sp[i](0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        df_sp[i](0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;

        // df wrt to s "+1" for each node(0..7) at each sp [i]
        df_sp[i](1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        df_sp[i](1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        df_sp[i](1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        df_sp[i](1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        df_sp[i](1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        df_sp[i](1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        df_sp[i](1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        df_sp[i](1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;

        // df wrt to t "+2" for each node(0..7) at each sp [i]
        df_sp[i](2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        df_sp[i](2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        df_sp[i](2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        df_sp[i](2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        df_sp[i](2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        df_sp[i](2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        df_sp[i](2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        df_sp[i](2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
    dfsp_eval = 1;               // now all arrays are filled statically
  }

  for (int sp=0; sp<num_sp; ++sp){
    // compute Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp],xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp],xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> jac_cur;
  for (int sp = 0; sp < num_sp; ++sp) {
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    jac_cur.Multiply(df_sp[sp],xcurr);

    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_SOH8; ++inode) {
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp*num_ans+0,inode*3+dim) = df_sp[sp](2,inode)*jac_cur(2,dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp*num_ans+1,inode*3+dim) = df_sp[sp](1,inode)*jac_cur(2,dim)
                                            +df_sp[sp](2,inode)*jac_cur(1,dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp*num_ans+2,inode*3+dim) = df_sp[sp](0,inode)*jac_cur(2,dim)
                                            +df_sp[sp](2,inode)*jac_cur(0,dim);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  evaluate 'T'-transformation matrix )                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::sosh8_evaluateT(const LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8>& jac,
                                                  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8>& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0,0) = jac(0,0) * jac(0,0);
  TinvT(1,0) = jac(1,0) * jac(1,0);
  TinvT(2,0) = jac(2,0) * jac(2,0);
  TinvT(3,0) = 2 * jac(0,0) * jac(1,0);
  TinvT(4,0) = 2 * jac(1,0) * jac(2,0);
  TinvT(5,0) = 2 * jac(0,0) * jac(2,0);

  TinvT(0,1) = jac(0,1) * jac(0,1);
  TinvT(1,1) = jac(1,1) * jac(1,1);
  TinvT(2,1) = jac(2,1) * jac(2,1);
  TinvT(3,1) = 2 * jac(0,1) * jac(1,1);
  TinvT(4,1) = 2 * jac(1,1) * jac(2,1);
  TinvT(5,1) = 2 * jac(0,1) * jac(2,1);

  TinvT(0,2) = jac(0,2) * jac(0,2);
  TinvT(1,2) = jac(1,2) * jac(1,2);
  TinvT(2,2) = jac(2,2) * jac(2,2);
  TinvT(3,2) = 2 * jac(0,2) * jac(1,2);
  TinvT(4,2) = 2 * jac(1,2) * jac(2,2);
  TinvT(5,2) = 2 * jac(0,2) * jac(2,2);

  TinvT(0,3) = jac(0,0) * jac(0,1);
  TinvT(1,3) = jac(1,0) * jac(1,1);
  TinvT(2,3) = jac(2,0) * jac(2,1);
  TinvT(3,3) = jac(0,0) * jac(1,1) + jac(1,0) * jac(0,1);
  TinvT(4,3) = jac(1,0) * jac(2,1) + jac(2,0) * jac(1,1);
  TinvT(5,3) = jac(0,0) * jac(2,1) + jac(2,0) * jac(0,1);


  TinvT(0,4) = jac(0,1) * jac(0,2);
  TinvT(1,4) = jac(1,1) * jac(1,2);
  TinvT(2,4) = jac(2,1) * jac(2,2);
  TinvT(3,4) = jac(0,1) * jac(1,2) + jac(1,1) * jac(0,2);
  TinvT(4,4) = jac(1,1) * jac(2,2) + jac(2,1) * jac(1,2);
  TinvT(5,4) = jac(0,1) * jac(2,2) + jac(2,1) * jac(0,2);

  TinvT(0,5) = jac(0,0) * jac(0,2);
  TinvT(1,5) = jac(1,0) * jac(1,2);
  TinvT(2,5) = jac(2,0) * jac(2,2);
  TinvT(3,5) = jac(0,0) * jac(1,2) + jac(1,0) * jac(0,2);
  TinvT(4,5) = jac(1,0) * jac(2,2) + jac(2,0) * jac(1,2);
  TinvT(5,5) = jac(0,0) * jac(2,2) + jac(2,0) * jac(0,2);

  // now evaluate T^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<NUMSTR_SOH8,NUMSTR_SOH8,1> solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Tinv (Jacobian) failed");
  return;
}

/*----------------------------------------------------------------------*
 |  return Cauchy stress at gp                                 maf 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::sosh8_Cauchy(LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,
                                         const int gp,
                                         const LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8>& deriv,
                                         const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8>& xrefe,
                                         const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8>& xcurr,
                                         const LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1>& glstrain,
                                         const LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1>& stress)
{
  // with ANS you do NOT have the correct (locking-free) F, so we
  // compute it here JUST for mapping of correct (locking-free) stresses
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> invJ;
  invJ.Multiply(deriv,xrefe);
  invJ.Invert();
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  N_XYZ.Multiply(invJ,deriv);
  // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd;
  defgrd.MultiplyTT(xcurr,N_XYZ);

# if consistent_F
  //double disp1 = defgrd.NormOne();
  //double dispinf = defgrd.NormInf();

  /* to get the consistent (locking-free) F^mod, we need two spectral
   * compositions. First, find R (rotation tensor) from F=RU,
   * then from E^mod = 1/2((U^mod)^2 - 1) find U^mod,
   * and finally F^mod = RU^mod */

  // polar decomposition of displacement based F
  LINALG::SerialDenseMatrix u(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::SerialDenseMatrix s(NUMDIM_SOH8,NUMDIM_SOH8);
  LINALG::SerialDenseMatrix v(NUMDIM_SOH8,NUMDIM_SOH8);
  SVD(defgrd,u,s,v); // Singular Value Decomposition
  LINALG::SerialDenseMatrix rot(NUMDIM_SOH8,NUMDIM_SOH8);
  rot.Multiply('N','N',1.0,u,v,0.0);
  //temp.Multiply('N','N',1.0,v,s,0.0);
  //LINALG::SerialDenseMatrix stretch_disp(NUMDIM_SOH8,NUMDIM_SOH8);
  //stretch_disp.Multiply('N','T',1.0,temp,v,0.0);
  //defgrd.Multiply('N','N',1.0,rot,stretch_disp,0.0);
  //cout << defgrd;

  // get modified squared stretch (U^mod)^2 from glstrain
  LINALG::SerialDenseMatrix Usq_mod(NUMDIM_SOH8,NUMDIM_SOH8);
  for (int i = 0; i < NUMDIM_SOH8; ++i) Usq_mod(i,i) = 2.0 * glstrain(i) + 1.0;
  // off-diagonal terms are already twice in the Voigt-GLstrain-vector
  Usq_mod(0,1) =  glstrain(3);  Usq_mod(1,0) =  glstrain(3);
  Usq_mod(1,2) =  glstrain(4);  Usq_mod(2,1) =  glstrain(4);
  Usq_mod(0,2) =  glstrain(5);  Usq_mod(2,0) =  glstrain(5);
  // polar decomposition of (U^mod)^2
  SVD(Usq_mod,u,s,v); // Singular Value Decomposition
  LINALG::SerialDenseMatrix U_mod(NUMDIM_SOH8,NUMDIM_SOH8);
  for (int i = 0; i < NUMDIM_SOH8; ++i) s(i,i) = sqrt(s(i,i));
  LINALG::SerialDenseMatrix temp2(NUMDIM_SOH8,NUMDIM_SOH8);
  temp2.Multiply('N','N',1.0,u,s,0.0);
  U_mod.Multiply('N','N',1.0,temp2,v,0.0);

  // F^mod = RU^mod
  LINALG::SerialDenseMatrix defgrd_consistent(NUMDIM_SOH8,NUMDIM_SOH8);
  defgrd_consistent.Multiply('N','N',1.0,rot,U_mod,0.0);
  defgrd.SetView(defgrd_consistent.A());

  /*
  double mod1 = defgrd.NormOne();
  double modinf = defgrd.NormInf();
  if(((mod1-disp1)/mod1 > 0.03) || ((modinf-dispinf)/modinf > 0.03)){
    cout << "difference in F! mod1= " << mod1 << " disp1= " << disp1 << " modinf= " << modinf << " dispinf= " << dispinf << endl;
    cout << "Fmod" << endl << defgrd;
  }
  */
#endif

  double detF = defgrd.Determinant();

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

  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchystress;
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
  temp.Multiply(1.0/detF,defgrd,pkstress);
  cauchystress.MultiplyNT(temp,defgrd);

  (*elestress)(gp,0) = cauchystress(0,0);
  (*elestress)(gp,1) = cauchystress(1,1);
  (*elestress)(gp,2) = cauchystress(2,2);
  (*elestress)(gp,3) = cauchystress(0,1);
  (*elestress)(gp,4) = cauchystress(1,2);
  (*elestress)(gp,5) = cauchystress(0,2);
}




/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Sosh8Register::Initialize(DRT::Discretization& dis)
{
  //sosh8_gmshplotdis(dis);

  int num_morphed_so_hex8 = 0;

  // Loop through all elements
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::ELEMENTS::So_sh8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");

    if (!actele->nodes_rearranged_) {
      // check for automatic definition of thickness direction
      if (actele->thickdir_ == DRT::ELEMENTS::So_sh8::autoj) {
        actele->thickdir_ = actele->sosh8_findthickdir();
      }

      int new_nodeids[NUMNOD_SOH8];

      switch (actele->thickdir_) {
      case DRT::ELEMENTS::So_sh8::autor:
      case DRT::ELEMENTS::So_sh8::globx: {
        // resorting of nodes to arrive at local t-dir for global x-dir
        new_nodeids[0] = actele->NodeIds()[7];
        new_nodeids[1] = actele->NodeIds()[4];
        new_nodeids[2] = actele->NodeIds()[0];
        new_nodeids[3] = actele->NodeIds()[3];
        new_nodeids[4] = actele->NodeIds()[6];
        new_nodeids[5] = actele->NodeIds()[5];
        new_nodeids[6] = actele->NodeIds()[1];
        new_nodeids[7] = actele->NodeIds()[2];
//        actele->sosh8_gmshplotlabeledelement(actele->NodeIds());
//        actele->sosh8_gmshplotlabeledelement(new_nodeids);
        actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::autos:
      case DRT::ELEMENTS::So_sh8::globy: {
        // resorting of nodes to arrive at local t-dir for global y-dir
        new_nodeids[0] = actele->NodeIds()[4];
        new_nodeids[1] = actele->NodeIds()[5];
        new_nodeids[2] = actele->NodeIds()[1];
        new_nodeids[3] = actele->NodeIds()[0];
        new_nodeids[4] = actele->NodeIds()[7];
        new_nodeids[5] = actele->NodeIds()[6];
        new_nodeids[6] = actele->NodeIds()[2];
        new_nodeids[7] = actele->NodeIds()[3];
        actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::autot:
      case DRT::ELEMENTS::So_sh8::globz: {
        // no resorting necessary
        for (int node = 0; node < 8; ++node) {
          new_nodeids[node] = actele->NodeIds()[node];
        }
        actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::undefined: {
        // here comes plan B: morph So_sh8 to So_hex8
        actele->SetType(DRT::Element::element_so_hex8);
        actele->soh8_reiniteas(DRT::ELEMENTS::So_hex8::soh8_easmild);
        actele->InitJacobianMapping();
        num_morphed_so_hex8++;
        break;
      }
      case DRT::ELEMENTS::So_sh8::none: break;
      default:
        dserror("no thickness direction for So_sh8");
      }
      //actele->sosh8_gmshplotlabeledelement(actele->NodeIds());
    }
  }

  if (num_morphed_so_hex8>0){
    cout << endl << num_morphed_so_hex8
    << " Sosh8-Elements have no clear 'thin' direction and have morphed to So_hex8 with eas_mild" << endl;
  }


  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false,false,false);

  // **************** debug printout ot gmesh **********************************
  //sosh8_gmshplotdis(dis);

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
