/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_evaluate.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef D_SOLID3
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_sh8p8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io_gmsh.H"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/yeoh.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            bborn 03/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8::Evaluate(
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  std::vector<int>& lm,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_,NUMDOF_> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_,NUMDOF_> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
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
  else if (action=="eas_init_multi")              act = So_hex8::eas_init_multi;
  else if (action=="eas_set_multi")               act = So_hex8::eas_set_multi;
  else if (action=="calc_homog_dens")             act = So_hex8::calc_homog_dens;
  else if (action=="multi_readrestart")           act = So_hex8::multi_readrestart;
  else dserror("Unknown type of action for So_hex8");

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need zero current displacement and residual forces
      LINALG::Matrix<NUMDISP_,1> mydisp(true);
      LINALG::Matrix<NUMPRES_,1> mypres(true);
      LINALG::Matrix<NUMDISP_,NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_,NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_,NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_,1> force(true);
      LINALG::Matrix<NUMPRES_,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
      BuildElementMatrix(&elemat1,&stiffmatrix,&gradmatrix,NULL,&stabmatrix);
      BuildElementVector(&elevec1,&force,&incomp);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mystat,lm);
      LINALG::Matrix<NUMDISP_,1> mydisp;
      LINALG::Matrix<NUMPRES_,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_,NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_,NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_,NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_,1> force(true);
      LINALG::Matrix<NUMPRES_,1> incomp(true);
      double volume = 0.0;
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,&volume,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
//      cout << " disp=" << mydisp.Norm2() << " pres=" << mypres.Norm2();
//      cout << " stiff=" << stiffmatrix.Norm2() << " grad=" << gradmatrix.Norm2() << " stab=" << stabmatrix.Norm2();
//      cout << " force=" << force.Norm2() << " incomp=" << incomp.Norm2() << endl;
      BuildElementMatrix(&elemat1,&stiffmatrix,&gradmatrix,NULL,&stabmatrix);
      BuildElementVector(&elevec1,&force,&incomp);
      AssembleVolume(params,volume);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce: {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mystat,lm);
      LINALG::Matrix<NUMDISP_,1> mydisp;
      LINALG::Matrix<NUMPRES_,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_,NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_,NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_,NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_,1> force(true);
      LINALG::Matrix<NUMPRES_,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
      BuildElementMatrix(&elemat1,&stiffmatrix,&gradmatrix,NULL,&stabmatrix);
      BuildElementVector(&elevec1,&force,&incomp);
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
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mystat,lm);
      LINALG::Matrix<NUMDISP_,1> mydisp;
      LINALG::Matrix<NUMPRES_,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_,NUMDISP_> massmatrix(true);
      LINALG::Matrix<NUMDISP_,NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_,NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_,NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_,1> force(true);
      LINALG::Matrix<NUMPRES_,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     &massmatrix,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
      // lump mass
      if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&massmatrix);
      // assemble displacement pressure parts
      BuildElementMatrix(&elemat2,&massmatrix,NULL,NULL,NULL);
      BuildElementMatrix(&elemat1,&stiffmatrix,&gradmatrix,NULL,&stabmatrix);
      BuildElementVector(&elevec1,&force,&incomp);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:{
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<std::vector<char> > stressdata
        = params.get<Teuchos::RCP<std::vector<char> > >("stress", Teuchos::null);
      Teuchos::RCP<std::vector<char> > straindata 
        = params.get<Teuchos::RCP<std::vector<char> > >("strain", Teuchos::null);
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==Teuchos::null) dserror("Cannot get stress 'data'");
      if (straindata==Teuchos::null) dserror("Cannot get strain 'data'");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mystat,lm);
      LINALG::Matrix<NUMDISP_,1> mydisp;
      LINALG::Matrix<NUMPRES_,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMGPT_,NUMSTR_> stress;
      LINALG::Matrix<NUMGPT_,NUMSTR_> strain;
      INPAR::STR::StressType iostress 
        = params.get<INPAR::STR::StressType>("iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain 
        = params.get<INPAR::STR::StrainType>("iostrain", INPAR::STR::strain_none);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,NULL,NULL,NULL,NULL,NULL,
                     &stress,&strain,NULL,params,iostress,iostrain);
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:{

      const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap
        = params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_,NUMSTR_> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_,NUMSTR_> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        std::vector<int> numadjele(NUMNOD_);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_;++i)
        {
          DRT::Node* node = nodes[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_;++i)
        {
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_;++i)
        {
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz") {
        Teuchos::RCP<Epetra_MultiVector> elestress
          = params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
        if (elestress==Teuchos::null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_; ++i)
          {
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_,NUMSTR_> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        std::vector<int> numadjele(NUMNOD_);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_;++i){
          DRT::Node* node=nodes[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        Teuchos::RCP<Epetra_MultiVector> elestress
          = params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
        if (elestress==Teuchos::null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_; ++i)
          {
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_;
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

      // Update of history for visco material
      Teuchos::RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast<MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast<MAT::ViscoAnisotropic*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_struct_update_imrlike: {
      // do something with internal EAS, etc parameters
      // this depends on the applied solution technique (static, generalised-alpha,
      // or other time integrators)

      // Update of history for visco material
      Teuchos::RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast<MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast<MAT::ViscoAnisotropic*>(mat.get());
        visco->Update();
      }
    }
    break;

    case calc_struct_reset_istep: {
      // do something with internal EAS, etc parameters

      // Reset of history for visco material
      Teuchos::RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast<MAT::ViscoNeoHooke*>(mat.get());
        visco->Reset();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast<MAT::ViscoAnisotropic*>(mat.get());
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

      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
        soh8_read_restart_multi(params);
    }
    break;

    default:
      dserror("Unknown type of action for So_sh8p8");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ForceStiffMass(
  const std::vector<int>& lm,             // location matrix
  const LINALG::Matrix<NUMDISP_,1>& disp,           // current displacements
  const LINALG::Matrix<NUMPRES_,1>& pres,       // current pressures
  LINALG::Matrix<NUMDISP_,NUMDISP_>* massmatrix,  // element mass matrix
  LINALG::Matrix<NUMDISP_,NUMDISP_>* stiffmatrix, // element stiffness matrix
  LINALG::Matrix<NUMDISP_,NUMPRES_>* gradmatrix, // element gradient matrix
  LINALG::Matrix<NUMPRES_,NUMPRES_>* stabmatrix,  // element stabilisation matrix
  LINALG::Matrix<NUMDISP_,1>* force,                 // element internal force vector
  LINALG::Matrix<NUMPRES_,1>* incomp,   // incompressibility residual
  LINALG::Matrix<NUMGPT_,NUMSTR_>* elestress,   // stresses at GP
  LINALG::Matrix<NUMGPT_,NUMSTR_>* elestrain,   // strains at GP
  double* volume,  // element volume
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress, // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_,1> > shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_,NUMNOD_> > derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_,NUMDIM_> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_,NUMDIM_> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp(i*NODDISP_+0,0);
    xcurr(i,1) = xrefe(i,1) + disp(i*NODDISP_+1,0);
    xcurr(i,2) = xrefe(i,2) + disp(i*NODDISP_+2,0);
  }

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  //LINALG::Matrix<NUMANS_*NUMSP_,NUMDOF_> B_ans_loc(true); //set to 0
  LINALG::Matrix<NUMANS_*NUMSP_,NUMDISP_> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_,NUMDIM_> > jac_sps(NUMSP_);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_,NUMDIM_> > jac_cur_sps(NUMSP_);
  // pointer to derivs evaluated at all sampling points
  std::vector<LINALG::Matrix<NUMDIM_,NUMNOD_> >* deriv_sp = NULL;
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe,xcurr,&deriv_sp,jac_sps,jac_cur_sps,B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double r[NUMGPT_] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
  const double s[NUMGPT_] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};

  // ---------------------------------------------------------------------
  // first loop over Gauss point
  // stabilisation matrices
  LINALG::Matrix<1,1> estabd(true); // element volume
  LINALG::Matrix<1,NUMNOD_> estabe(true); // integral of pressure shape functions across element
  LINALG::Matrix<NUMNOD_,NUMNOD_> estabm(true);  // mass-like matrix
  std::vector<LINALG::Matrix<1,1> > estaba(NUMGPT_);  // shape functions for projected Q0/constant pressure
  for (int gp=0; gp<NUMGPT_; ++gp)
  {
    // (transposed) material-to-parametric Jacobian J = (X_{,xi})^T
    LINALG::Matrix<NUMDIM_,NUMDIM_> Jac;
    Jac.Multiply(derivs[gp],xrefe);
    const double detJ = Jac.Determinant();

    estaba[gp].PutScalar(1.0);

    const double wdetJ = detJ*gpweights[gp];
    estabd.MultiplyTN(wdetJ,estaba[gp],estaba[gp],1.0);
    estabe.MultiplyTT(wdetJ,estaba[gp],shapefcts[gp],1.0);
    estabm.MultiplyNT(wdetJ,shapefcts[gp],shapefcts[gp],1.0);
  }

  // stabilisation matrix
  if (stabmatrix != NULL)
  {
    // shear modulus
    const double shearmod = ShearMod();
    // (-Cem) = -1./shearmod*( Mem - Eem'*inv(Dem)*Eem );
    stabmatrix->Update(estabm);
    stabmatrix->MultiplyTN(-1.0/estabd(0,0),estabe,estabe,1.0);
    stabmatrix->Scale(-1.0/shearmod);
  }
 

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_; ++gp)
  {

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<NUMDIM_,NUMDIM_> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ = jac.Determinant();
    if (fabs(detJ) <= EPS10) dserror("JACOBIAN DETERMINANT CLOSE TO ZERO");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    LINALG::Matrix<NUMDIM_,NUMDIM_> jac_cur;
    jac_cur.Multiply(derivs[gp],xcurr);

    // set up B-Operator in local(parameter) element space including ANS
    LINALG::Matrix<NUMSTR_,NUMDISP_> bop_loc;
    for (int inode = 0; inode < NUMNOD_; ++inode) {
      for (int dim = 0; dim < NUMDIM_; ++dim) {
        // B_loc_rr = N_r.X_r
        bop_loc(0,inode*3+dim) = derivs[gp](0,inode) * jac_cur(0,dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1,inode*3+dim) = derivs[gp](1,inode) * jac_cur(1,dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3,inode*3+dim) = derivs[gp](0,inode) * jac_cur(1,dim)
                               + derivs[gp](1,inode) * jac_cur(0,dim);
        if (ans_ == ans_none) {
          // B_loc_tt = N_t.X_t
          bop_loc(2,inode*3+dim) = derivs[gp](2,inode) * jac_cur(2,dim);
          // B_loc_st = N_s.X_t + N_t.X_s
          bop_loc(4,inode*3+dim) = derivs[gp](1,inode) * jac_cur(2,dim)
                                 + derivs[gp](2,inode) * jac_cur(1,dim);
          // B_loc_rt = N_r.X_t + N_t.X_r
          bop_loc(5,inode*3+dim) = derivs[gp](0,inode) * jac_cur(2,dim)
                                 + derivs[gp](2,inode) * jac_cur(0,dim);
        }
        else {
          // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
          //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
          //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
          bop_loc(2,inode*3+dim) = 0.25*(1-r[gp])*(1-s[gp]) * B_ans_loc(0+4*NUMANS_,inode*3+dim)
                                 + 0.25*(1+r[gp])*(1-s[gp]) * B_ans_loc(0+5*NUMANS_,inode*3+dim)
                                 + 0.25*(1+r[gp])*(1+s[gp]) * B_ans_loc(0+6*NUMANS_,inode*3+dim)
                                 + 0.25*(1-r[gp])*(1+s[gp]) * B_ans_loc(0+7*NUMANS_,inode*3+dim);
          // B_loc_st = interpolation along r of ANS B_loc_st
          //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
          bop_loc(4,inode*3+dim) = 0.5*(1.0+r[gp]) * B_ans_loc(1+1*NUMANS_,inode*3+dim)
                                 + 0.5*(1.0-r[gp]) * B_ans_loc(1+3*NUMANS_,inode*3+dim);
          // B_loc_rt = interpolation along s of ANS B_loc_rt
          //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
          bop_loc(5,inode*3+dim) = 0.5*(1.0-s[gp]) * B_ans_loc(2+0*NUMANS_,inode*3+dim)
                                 + 0.5*(1.0+s[gp]) * B_ans_loc(2+2*NUMANS_,inode*3+dim);
        }
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    LINALG::Matrix<NUMSTR_,NUMSTR_> TinvT;
    sosh8_evaluateT(jac,TinvT);
    LINALG::Matrix<NUMSTR_,NUMDISP_> bop;
    bop.Multiply(TinvT,bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    // but with modified ANS strains E33, E23 and E13
    LINALG::Matrix<NUMSTR_,1> lstrain;
    // evaluate glstrains in local(parameter) coords
    // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
    lstrain(0) = 0.5 * (
       +(jac_cur(0,0)*jac_cur(0,0) + jac_cur(0,1)*jac_cur(0,1) + jac_cur(0,2)*jac_cur(0,2))
       -(jac(0,0)*jac(0,0)         + jac(0,1)*jac(0,1)         + jac(0,2)*jac(0,2)));
    // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
    lstrain(1) = 0.5 * (
       +(jac_cur(1,0)*jac_cur(1,0) + jac_cur(1,1)*jac_cur(1,1) + jac_cur(1,2)*jac_cur(1,2))
       -(jac(1,0)*jac(1,0)         + jac(1,1)*jac(1,1)         + jac(1,2)*jac(1,2)));
    // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
    lstrain(3) = (
       +(jac_cur(0,0)*jac_cur(1,0) + jac_cur(0,1)*jac_cur(1,1) + jac_cur(0,2)*jac_cur(1,2))
       -(jac(0,0)*jac(1,0)         + jac(0,1)*jac(1,1)         + jac(0,2)*jac(1,2)));
    // remaining natural strains
    if (ans_ == ans_none) {
      // Ett = 0.5 * (dz/dt * dz/dt^T - dZ/dt * dZ/dt^T)
      lstrain(2) = 0.5 * (
        +(jac_cur(2,0)*jac_cur(2,0) + jac_cur(2,1)*jac_cur(2,1) + jac_cur(2,2)*jac_cur(2,2))
        -(jac(2,0)*jac(2,0)         + jac(2,1)*jac(2,1)         + jac(2,2)*jac(2,2)));
      // Est = (dx/dt * dy/ds^T - dX/dt * dY/ds^T)
      lstrain(4) = (
        +(jac_cur(1,0)*jac_cur(2,0) + jac_cur(1,1)*jac_cur(2,1) + jac_cur(1,2)*jac_cur(2,2))
        -(jac(1,0)*jac(2,0)         + jac(1,1)*jac(2,1)         + jac(1,2)*jac(2,2)));
      // Etr = (dx/dr * dy/dt^T - dX/dr * dY/dt^T)
      lstrain(5) = (
        +(jac_cur(2,0)*jac_cur(0,0) + jac_cur(2,1)*jac_cur(0,1) + jac_cur(2,2)*jac_cur(0,2))
        -(jac(2,0)*jac(0,0)         + jac(2,1)*jac(0,1)         + jac(2,2)*jac(0,2)));
    }
    // ANS modification of strains ************************************** ANS
    else {
      double dydt_A = 0.0; double dYdt_A = 0.0;
      double dxdt_B = 0.0; double dXdt_B = 0.0;
      double dydt_C = 0.0; double dYdt_C = 0.0;
      double dxdt_D = 0.0; double dXdt_D = 0.0;
      double dzdt_E = 0.0; double dZdt_E = 0.0;
      double dzdt_F = 0.0; double dZdt_F = 0.0;
      double dzdt_G = 0.0; double dZdt_G = 0.0;
      double dzdt_H = 0.0; double dZdt_H = 0.0;

      // vector product of rows of jacobians at corresponding sampling point    cout << jac_cur_sps;
      for (int dim = 0; dim < NUMDIM_; ++dim) {
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
        + 0.25*(1+r[gp])*(1-s[gp]) * (dzdt_F - dZdt_F)
        + 0.25*(1+r[gp])*(1+s[gp]) * (dzdt_G - dZdt_G)
        + 0.25*(1-r[gp])*(1+s[gp]) * (dzdt_H - dZdt_H));
      // E23: remedy of transverse shear locking
      // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
      lstrain(4) = 0.5*(1+r[gp]) * (dxdt_B - dXdt_B) + 0.5*(1-r[gp]) * (dxdt_D - dXdt_D);
      // E13: remedy of transverse shear locking
      // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
      lstrain(5) = 0.5*(1-s[gp]) * (dydt_A - dYdt_A) + 0.5*(1+s[gp]) * (dydt_C - dYdt_C);
    }
    // ANS modification of strains ************************************** ANS

    // push local/natural/parametric glstrains forward to global/material space
    LINALG::Matrix<NUMSTR_,1> glstrain;
    glstrain.Multiply(TinvT,lstrain);

    // recover deformation gradient incoperating assumed natural GL strain
    double detdefgrad;  // determinant of assumed def.grad.
    static LINALG::Matrix<NUMDIM_,NUMDIM_> defgrad;  // assumed def.grad.
    static LINALG::Matrix<NUMDIM_,NUMDIM_> invdefgrad; // inverse of deformation gradient and its derivative
    static LINALG::Matrix<NUMDIM_,NUMDIM_> rgtstr;  // assumed material stretch
    static LINALG::Matrix<NUMDIM_,NUMDIM_> defgradD;  // pure disp-based def.grad.
    static LINALG::Matrix<NUMDIM_,NUMDIM_> rgtstrD;  // pure disp.-based material stretch
    static LINALG::Matrix<NUMDIM_,NUMDIM_> invrgtstrD;  // inverse of pure disp-based material stretch tensor U^{d;-1}
    AssDefGrad(detdefgrad,defgrad,invdefgrad,rgtstr,defgradD,rgtstrD,invrgtstrD,invJ_[gp],jac,jac_cur,glstrain);
//    cout << defgrad << invdefgrad << rgtstr << defgradD << rgtstrD << invrgtstrD << endl << endl;
    
    // assumend right material stretch 6-Voigt vector
//    LINALG::Matrix<NUMSTR_,1> rgtstrv;
//    Matrix2TensorToVector6Voigt(rgtstrv,rgtstr,voigt6_strain);

    // inverse of deformation gradient and its derivative 
//    LINALG::Matrix<NUMDIM_,NUMDIM_> invdefgrad(defgrad);
//    double detdefgrad = invdefgrad.Invert();

    // return gp strains if necessary
    if (iostrain != INPAR::STR::strain_none)
      Strain(elestrain,iostrain,gp,detdefgrad,defgrad,invdefgrad,glstrain);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_,NUMSTR_> cmat(true);
    LINALG::Matrix<NUMSTR_,1> stress(true);
    soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrad,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // pressure at Gauss point
    const double pressure = (shapefcts[gp]).Dot(pres);

    // return Gauss point stresses if necessary
    if (iostress != INPAR::STR::stress_none)
      Stress(elestress,iostress,gp,detdefgrad,defgrad,glstrain,stress,pressure);

    // effective shape function of scalar pressure field at current Gauss point
    static LINALG::Matrix<NUMPRES_,1> prshfct;
    if (stab_ == stab_nonaffine)
      prshfct.MultiplyTN(1.0/estabd(0,0),estabe,estaba[gp]);
    else if (stab_ == stab_affine)
      prshfct.Update(shapefcts[gp]);
 
    // integration factor
    const double detJ_w = detJ*gpweights[gp];

    // internal force
    if (force != NULL)
    {
      // integrate internal force vector
      // fint := fint 
      //      + (B^T . sigma) * detJ * w(gp) 
      //      + (-G) . ep   // will be done _after_ Gauss point loop
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }
    // incompressiblity equation
    if (incomp != NULL)
    {
      // pint := pint 
      //       - He . (Fdet - 1.0) * detJ * wp(gp)
      //       + (-Ce) . ep   // will be done _after_ Gauss point loop
      incomp->Update(-(detdefgrad-1.0)*detJ_w,prshfct,1.0);
    }
    // stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<NUMSTR_,NUMDISP_> cb;
      cb.Multiply(cmat,bop); // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu
      // here also the ANS interpolation comes into play
      Teuchos::RCP<LINALG::Matrix<NUMSTR_,NUMNOD_*NUMNOD_> > bopbydisp = Teuchos::null;
      if (lin_ > lin_sixth)
        bopbydisp = Teuchos::rcp(new LINALG::Matrix<NUMSTR_,NUMNOD_*NUMNOD_>());
      for (int inod=0; inod<NUMNOD_; ++inod)
      {
        for (int jnod=0; jnod<NUMNOD_; ++jnod)
        {
          LINALG::Matrix<NUMSTR_,1> G_ij;
          G_ij(0) = derivs[gp](0, inod) * derivs[gp](0, jnod); // rr-dir
          G_ij(1) = derivs[gp](1, inod) * derivs[gp](1, jnod); // ss-dir
          G_ij(3) = derivs[gp](0, inod) * derivs[gp](1, jnod)
                  + derivs[gp](1, inod) * derivs[gp](0, jnod); // rs-dir
          if (ans_ == ans_none) {
            G_ij(2) = derivs[gp](2, inod) * derivs[gp](2, jnod); // tt-dir
            G_ij(4) = derivs[gp](1, inod) * derivs[gp](2, jnod)
                    + derivs[gp](2, inod) * derivs[gp](1, jnod); // st-dir
            G_ij(5) = derivs[gp](2, inod) * derivs[gp](0, jnod)
                    + derivs[gp](0, inod) * derivs[gp](2, jnod); // tr-dir
          }
          else {
            // ANS modification in tt-dir
            G_ij(2) = 0.25*(1-r[gp])*(1-s[gp]) * (*deriv_sp)[4](2,inod) * (*deriv_sp)[4](2,jnod)
                    + 0.25*(1+r[gp])*(1-s[gp]) * (*deriv_sp)[5](2,inod) * (*deriv_sp)[5](2,jnod)
                    + 0.25*(1+r[gp])*(1+s[gp]) * (*deriv_sp)[6](2,inod) * (*deriv_sp)[6](2,jnod)
                    + 0.25*(1-r[gp])*(1+s[gp]) * (*deriv_sp)[7](2,inod) * (*deriv_sp)[7](2,jnod);
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
          }

          // transformation of local(parameter) space 'back' to global(material) space
          LINALG::Matrix<NUMSTR_,1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);

          // store B_{aBd,k}
          if (lin_ > lin_sixth)
            for (int istr=0; istr<NUMSTR_; ++istr)
              (*bopbydisp)(istr,NUMNOD_*inod+jnod) = G_ij_glob(istr);
            
          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_*inod+0, NUMDIM_*jnod+0) += Gij;
          (*stiffmatrix)(NUMDIM_*inod+1, NUMDIM_*jnod+1) += Gij;
          (*stiffmatrix)(NUMDIM_*inod+2, NUMDIM_*jnod+2) += Gij;
        }
      } // end of integrate `geometric' stiffness

      // add (incomplete) derivative of pressure-proportional force w.r.t. displacements
      // Kp = (-Gm*ep')_,d
      //    = (-dFv'*fvT*pN * Fdet*detJ*wp(i))_,d * ep'
      // Kp = Kp + dFv'*fvT*(pN*ep')*fvT'*dFv * Fdet*detJ*wp(i)  // due to Fdet_,d = fvT'*dFv * Fdet
      //         + (pN*ep')*dFv'*WmT*dFv * Fdet*detJ*wp(i)       // due to fvT_,d = WmT*dFv
      //         + ddFv * Fdet*detJ*wp(i)                        // due to ddFv = dFv'_,d*fvT = (fv'*dFv_,d)'
      // Ke = Keu + Kg - Kp;
      {
        // effective pressure at Gauss point
        const double effpressure = prshfct.Dot(pres);  // pN*ep'
        // Voigt 9-vector of transposed & inverted deformation gradient fvT := F^{-T}
        LINALG::Matrix<NUMDFGR_,1> tinvdefgrad;
        Matrix2TensorToVector9Voigt(tinvdefgrad,invdefgrad,true);

        // derivative of WmT := F^{-T}_{,F} in Voigt vector notation
        LINALG::Matrix<NUMDFGR_,NUMDFGR_> WmT;
        InvVector9VoigtDiffByItself(WmT,invdefgrad,true);
        // WmT := WmT + fvT*fvT'
        WmT.MultiplyNT(1.0,tinvdefgrad,tinvdefgrad,1.0);

        // Voigt vector indices
        const int* voigt6row = NULL;
        const int* voigt6col = NULL;
        Indices6VoigtTo2Tensor(voigt6row,voigt6col);
        const int* voigt9row = NULL;
        const int* voigt9col = NULL;
        Indices9VoigtTo2Tensor(voigt9row,voigt9col);
        const int* voigt3x3 = NULL;
        Indices2TensorTo9Voigt(voigt3x3);  // access is via (i,j) -> 3*i+j
        const int* voigt3x3sym = NULL;
        Indices2TensorTo6Voigt(voigt3x3sym);  // access is via (i,j) -> 3*i+j
        
        // material derivatives of shape functions
        LINALG::Matrix<NUMDIM_,NUMNOD_> derivsmat;
        derivsmat.MultiplyNN(invJ_[gp],derivs[gp]);

        // linear B-op
        // derivative of displ-based def.grad with respect to nodal displacements
        // F^d_{aC,d}
        int iboplin[NUMDFGR_*NUMNOD_];  // index entries which are non-equal zero
        LINALG::Matrix<NUMDFGR_,NUMDISP_> boplin(true);
        for (int ij=0; ij<NUMDFGR_; ++ij) {
          const int i = voigt9row[ij];
          const int j = voigt9col[ij];
          for (int k=0; k<NUMNOD_; ++k) {
            const int K = k*NODDISP_ + i;
            iboplin[NUMDFGR_*ij+k] = K;
            boplin(ij,K) = derivsmat(j,k);
          }
        }

        // derivative of def.grad. w.r.t. displacements Fv_{,d}
        LINALG::Matrix<NUMDFGR_,NUMDISP_> defgradbydisp;
        if (ans_ == ans_none) {
          defgradbydisp.Update(boplin);
        }
        else {
//          // inverse of pure disp-based material stretch tensor
//          // U^{d;-1}
//          LINALG::Matrix<NUMDIM_,NUMDIM_> invrgtstrD(rgtstrD);
//          const double detrgtstrD = invrgtstrD.Invert();
//          if (detrgtstrD < 0.0) dserror("Trouble in inverting right stretch tensor");
//          LINALG::Matrix<NUMSTR_,1> invrgtstrDv;
//          Matrix2TensorToVector6Voigt(invrgtstrDv,invrgtstrD,voigt6_strain);
          
          // derivative of pure-disp. inverse material stretch tensor with respect to nodal displacements
          // U^{d;-1}_{,d} = U^{d;-1}_{,U} . (C^d_{,U^d})^{-1} . C^d_{,d}
          // on exit of this block the following variables are going to hold ...
          static LINALG::Matrix<NUMSTR_,NUMDISP_> invrgtstrDbydisp; // ...U^{d;-1}_{,d}
          static LINALG::Matrix<NUMSTR_,NUMSTR_> invrgtstrDbyrgtstrD; // ...U^{d;-1}_{,U}
          static LINALG::Matrix<NUMSTR_,NUMSTR_> rcgDbyrgtstrD; // ...(C^d_{,U^d})^{-1}
          static LINALG::Matrix<NUMSTR_,NUMDISP_> rgtstrDbydisp; // ...U^d_{,d}
          {
            // U^{d;-1}_{,U}
            //LINALG::Matrix<NUMSTR_,NUMSTR_> invrgtstrDbyrgtstrD;
            InvVector6VoigtDiffByItself(invrgtstrDbyrgtstrD,invrgtstrD);

            // C^d_{,U^d} = (U^d . U^d)_{,U^d}
            //LINALG::Matrix<NUMSTR_,NUMSTR_> rcgDbyrgtstrD;
            SqVector6VoigtDiffByItself(rcgDbyrgtstrD,rgtstrD);

            // deformation gradient as Voigt matrix
            LINALG::Matrix<NUMSTR_,NUMDFGR_> defgradm;
            Matrix2TensorToMatrix6x9Voigt(defgradm,defgrad,true);

            // C^d_{,d} = 2 * Fm * Boplin, 6x24
            LINALG::Matrix<NUMSTR_,NUMDISP_> rcgDbydisp;
            rcgDbydisp.MultiplyNN(2.0,defgradm,boplin);
            
            // U^d_{,d} = (C^d_{,U^d})^{-1} . C^d_{,d}
            //LINALG::Matrix<NUMSTR_,NUMDISP_> rgtstrDbydisp;
            {
              LINALG::FixedSizeSerialDenseSolver<NUMSTR_,NUMSTR_,NUMDISP_> rcgDbyrgtstrDsolver;
              rcgDbyrgtstrDsolver.SetMatrix(rcgDbyrgtstrD);  // LHS
              rcgDbyrgtstrDsolver.SetVectors(rgtstrDbydisp,rcgDbydisp);  // SOL, RHS
              const int err = rcgDbyrgtstrDsolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
              if (lin_ >= lin_one) {
                const int err = rcgDbyrgtstrDsolver.Invert();
                if (err != 0) dserror("Failed to invert, error=%d", err);
              }
            }

            // U^{d;-1}_{,d} = U^{d;-1}_{,U} . U^d_{,d}
            invrgtstrDbydisp.MultiplyNN(invrgtstrDbyrgtstrD,rgtstrDbydisp);
          }

          // derivative of ass. mat. stretch tensor with respect to nodal displacements
          // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
          {
            // derivative of ass. right Cauchy-Green with respect to ass. material stretch tensor
            // C^{ass}_{,U^{ass}}
            LINALG::Matrix<NUMSTR_,NUMSTR_> rcgbyrgtstr;
            SqVector6VoigtDiffByItself(rcgbyrgtstr,rgtstr);

            // C^{ass}_{,d} = 2 * bop
            // C^{ass}_{AB,k} = 2 B_{ABk}

            // derivative of ass. mat. stretch tensor with respect to nodal displacements
            // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
            LINALG::Matrix<NUMSTR_,NUMDISP_> rgtstrbydisp;  // ... U^{ass}_{,d}
            LINALG::FixedSizeSerialDenseSolver<NUMSTR_,NUMSTR_,NUMDISP_> rcgbyrgtstrsolver;
            {
              rcgbyrgtstrsolver.SetMatrix(rcgbyrgtstr);  // LHS
              rcgbyrgtstrsolver.SetVectors(rgtstrbydisp,bop);  // SOL, RHS
              const int err = rcgbyrgtstrsolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
            }
            rgtstrbydisp.Scale(2.0);

            // derivative of def.grad. with respect to k nodal displacements d^k
            // F_{aB,k} = F^d_{aC,k} . U^{d;-1}_{CD} . U^{ass}_{DB}
            //          + F^d_{aC} . U^{d;-1}_{CD,k} . U^{ass}_{DB}
            //          + F^d_{aC} . U^{d;-1}_{CD} . U^{ass}_{DB,k}
            for (int aB=0; aB<NUMDFGR_; ++aB) {
              for (int k=0; k<NUMDISP_; ++k) {
                double defgradbydisp_aBk = 0.0;
                const int a = voigt9row[aB];
                const int B = voigt9col[aB];
                for (int C=0; C<NUMDIM_; ++C) {
                  for (int D=0; D<NUMDIM_; ++D) {
                    const int aC = voigt3x3[NUMDIM_*a+C];
                    const int CD = voigt3x3sym[NUMDIM_*C+D];
                    const int DB = voigt3x3sym[NUMDIM_*D+B];
                    const double CDfact = (C==D) ? 1.0 : 0.5;
                    const double DBfact = (D==B) ? 1.0 : 0.5;
                    defgradbydisp_aBk 
                      += boplin(aC,k) * invrgtstrD(C,D) * rgtstr(D,B)
                      + defgradD(a,C) * CDfact*invrgtstrDbydisp(CD,k) * rgtstr(D,B)
                      + defgradD(a,C) * invrgtstrD(C,D) * DBfact*rgtstrbydisp(DB,k);
                  }
                }
                defgradbydisp(aB,k) = defgradbydisp_aBk;
              }
            }


            // ext(p)ensive computation to achieve full tangent
            if (lin_ > lin_sixth) {
              // on #rcgbyrgtstr is stored the inverse of C^{ass}_{,U^{ass}}
              {
                const int err = rcgbyrgtstrsolver.Invert();
                if (err != 0) dserror("Failed to invert, error=%d", err);
              }

              // second derivative of assumed right Cauchy-Green tensor 
              // w.r.t. to right stretch tensor
              // C^{ass}_{,U^{ass} U^{ass}} = const
              int ircgbybyrgtstr[NUMSTR_*6];  // for sparse access
              LINALG::Matrix<NUMSTR_,6> rcgbybyrgtstr;
              SqVector6VoigtTwiceDiffByItself(ircgbybyrgtstr,rcgbybyrgtstr);

              // second derivative of disp-based right Cauchy-Green tensor 
              // w.r.t. to right stretch tensor
              // C^{d}_{,U^{d} U^{d}} = const
              // MARK: an extra variable is not needed as same as for assumed right CG tensor (above)

              // second derivative of disp-based inverse right stretch tensor
              // w.r.t. disp-based right stretch tensor
              // U^{d-1}_{,U^d U^d}
              LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_> invrgtstrDbybyrgtstrD;
              InvVector6VoigtTwiceDiffByItself(invrgtstrDbybyrgtstrD,invrgtstrD);

              // second derivative of assumed right stretch tensor w.r.t. displacements
              // U^{ass}_{DB,dk} = (C^{ass}_{,U^{ass}})_{DBEF}^{-1} 
              //                 . ( C^{ass}_{EF,dk} - C^{ass}_{EF,GHIJ}  U^{ass}_{IJ,d}  U^{ass}_{GH,k} )
              // and
              // second derivative of pure-disp right stretch tensor w.r.t. displacements
              // U^{d}_{DB,dk} = (C^{d}_{,U^{d}})_{DBEF}^{-1} 
              //               . ( C^{d}_{EF,dk} - C^{d}_{EF,GHIJ}  U^{d}_{IJ,d}  U^{d}_{GH,k} )
              LINALG::Matrix<NUMSTR_,NUMDISP_*NUMDISP_> rgtstrbybydisp;
              LINALG::Matrix<NUMSTR_,NUMDISP_*NUMDISP_> rgtstrDbybydisp;
              for (int DB=0; DB<NUMSTR_; ++DB) {
                for (int dk=0; dk<NUMDISP_*NUMDISP_; ++dk) {
                  const int d = dk / NUMDISP_;
                  const int k = dk % NUMDISP_;
                  if (k < d)  // symmetric in d and k : only upper 'triangle' is computed
                    continue;
                  const int kd = NUMDISP_*k + d;
                  int ndnk = -1;
                  if (d%NODDISP_ == k%NODDISP_) {
                    const int nd = d / NODDISP_;
                    const int nk = k / NODDISP_;
                    ndnk = nd*NUMNOD_ + nk;
                  }
                  double rgtstrbybydisp_DBdk = 0.0;
                  double rgtstrDbybydisp_DBdk = 0.0;
                  for (int EF=0; EF<NUMSTR_; ++EF) {
                    const int E = voigt9row[EF];
                    const int F = voigt9col[EF];
                    // C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,d}
                    // and
                    // C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,d}
                    double temp_EFdk = 0.0;
                    double tempD_EFdk = 0.0;
                    for (int GHIJ=0; GHIJ<6; ++GHIJ) {
                      if (ircgbybyrgtstr[NUMSTR_*EF+GHIJ] != -1) {
                        const int GH = ircgbybyrgtstr[NUMSTR_*EF+GHIJ] / NUMSTR_;
                        const int IJ = ircgbybyrgtstr[NUMSTR_*EF+GHIJ] % NUMSTR_;
                        // C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,d}
                        temp_EFdk += rcgbybyrgtstr(EF,GHIJ)*rgtstrbydisp(IJ,d)*rgtstrbydisp(GH,k);
                        // C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,d}
                        // C^{d}_{,U^d U^d} = C^{ass}_{,U^ass U^ass} = const
                        tempD_EFdk += rcgbybyrgtstr(EF,GHIJ)*rgtstrDbydisp(IJ,d)*rgtstrDbydisp(GH,k);
                      }
                    }
                    // U^{ass}_{DB,dk}
                    if (ndnk != -1)
                      rgtstrbybydisp_DBdk += rcgbyrgtstr(DB,EF) * 2.0*(*bopbydisp)(EF,ndnk);
                    rgtstrbybydisp_DBdk -= rcgbyrgtstr(DB,EF) * temp_EFdk;
                    // C^{d}_{EF,dk}
                    double rcgDbybydisp_EFdk = 0.0;
                    for (int m=0; m<NUMDIM_; ++m) {
                      const int mE = voigt3x3[NUMDIM_*m+E];
                      const int mF = voigt3x3[NUMDIM_*m+F];
                      if (E == F)  // make strain-like 6-Voigt vector
                        rcgDbybydisp_EFdk 
                          += 2.0*boplin(mE,d)*boplin(mF,k);
                      else  // thus setting  V_EF + V_FE if E!=F
                        rcgDbybydisp_EFdk 
                          += 2.0*boplin(mE,d)*boplin(mF,k)
                          +  2.0*boplin(mF,d)*boplin(mE,k);
                    }
                    // (C^{d}_{,U^{d}})_{DBEF}^{-1}
                    const double rcgDbyrgtstrD_DBEF = rcgDbyrgtstrD(DB,EF);
                    // U^{d}_{DB,dk}
                    rgtstrDbybydisp_DBdk += rcgDbyrgtstrD_DBEF * ( rcgDbybydisp_EFdk - tempD_EFdk);
                  }
                  rgtstrbybydisp(DB,dk) = rgtstrbybydisp_DBdk;
                  if (k != d) rgtstrbybydisp(DB,kd) = rgtstrbybydisp_DBdk;
                  rgtstrDbybydisp(DB,dk) = rgtstrDbybydisp_DBdk;
                  if (k != d) rgtstrDbybydisp(DB,kd) = rgtstrDbybydisp_DBdk;
                }
              }

              // second derivative of pure-disp inverse right stretch tensor w.r.t. displacements
              // U^{d-1}_{CD,dk} = U^{d-1}_{CD,EFGH} U^{d}_{GH,k} U^{d}_{EF,d}
              //                 + U^{d-1}_{CD,EF} U^{d}_{EF,dk}
              Teuchos::RCP<LINALG::Matrix<NUMSTR_,NUMDISP_*NUMDISP_> > invrgtstrDbybydisp = Teuchos::null; // ... U^{d-1}_{DB,dk}
              if (lin_ >= lin_one) {
                invrgtstrDbybydisp = Teuchos::rcp(new LINALG::Matrix<NUMSTR_,NUMDISP_*NUMDISP_>());

                // compute U^{d-1}_{DB,dk}
                for (int dk=0; dk<NUMDISP_*NUMDISP_; ++dk) {
                  const int d = dk / NUMDISP_;
                  const int k = dk % NUMDISP_;
                  if (k < d)  // symmetric in d and k : only upper triangle is computed
                    continue;
                  const int kd = NUMDISP_*k + d;
                  for (int CD=0; CD<NUMSTR_; ++CD) {
                    double invrgtstrDbybydisp_CDdk = 0.0;
                    for (int EF=0; EF<NUMSTR_; ++EF) {
                      const double rgtstrDbybydisp_EFdk = rgtstrDbybydisp(EF,dk);
                      invrgtstrDbybydisp_CDdk 
                        += invrgtstrDbyrgtstrD(CD,EF) * rgtstrDbybydisp_EFdk;
                      for (int GH=0; GH<NUMSTR_; ++GH) {
                        const int EFGH = NUMSTR_*EF + GH;
                        invrgtstrDbybydisp_CDdk 
                          += invrgtstrDbybyrgtstrD(CD,EFGH)  // col are strain-like 6-Voigt
                          * rgtstrDbydisp(GH,k)  // row are strain-like 6-Voigt too
                          * rgtstrDbydisp(EF,d);  // row are strain-like 6-Voigt too
                      }
                    }
                    (*invrgtstrDbybydisp)(CD,dk) = invrgtstrDbybydisp_CDdk;
                    if (k != d) (*invrgtstrDbybydisp)(CD,kd) = invrgtstrDbybydisp_CDdk;
                  }
                }
              } // if (lin_ >= lin_one) else
        

              // inverse assumed deformation gradient times boblin
              // F^{-T} . B_L = F^{-T} . F^d_{,d}
              LINALG::Matrix<NUMSTR_,NUMDISP_> invdefgradtimesboplin(true);  // sparse, 1/3 non-zeros
              for (int n=0; n<NUMNOD_; ++n) {
                const int d = iboplin[n];
                for (int BC=0; BC<NUMSTR_; ++BC) {
                  const int B = voigt6row[BC];
                  const int C = voigt6col[BC];
                  double invdefgradtimesboplin_BCd = 0.0;
                  for (int a=0; a<NUMDIM_; ++a) {
                    const int aC = voigt3x3[NUMDIM_*a+C];
                    const int aB = voigt3x3[NUMDIM_*a+B];
                    if (B == C)  // make strain-like 6-Voigt vector
                      invdefgradtimesboplin_BCd 
                        += invdefgrad(a,B)*boplin(aC,d);
                    else  // thus setting  V_BC + V_CB if C!=B
                      invdefgradtimesboplin_BCd 
                        += invdefgrad(a,B)*boplin(aC,d)
                        +  invdefgrad(a,C)*boplin(aB,d);
                  }
                  invdefgradtimesboplin(BC,d) = invdefgradtimesboplin_BCd;
                }
              }

              // I^{assd}_{BC} = F^{-T}_{Ba} . F^{d}_{aC}
              LINALG::Matrix<NUMDIM_,NUMDIM_> invdefgradtimesdefgradD;
              invdefgradtimesdefgradD.MultiplyTN(invdefgrad,defgradD);


              // contribute stuff containing second derivatives in displacements
              // F^{-T}_{aB} F_{aB,dk}
              // = F^{-1}_{Ba}  (F^d_{aC} U^{d;-1}_{CD} U^{ass}_{DB})_{,dk} 
              // = F^{-1}_{Ba}  F^d_{aC,dk} U^{d;-1}_{CD} U^{ass}_{DB}         |  = 0
              // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,dk} U^{ass}_{DB}         |  # 0, very pricy
              // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD} U^{ass}_{DB,dk}         |  # 0, pricy
              // + F^{-1}_{Ba}  F^d_{aC,d} U^{d;-1}_{CD,k} U^{ass}_{DB}        |  # 0, okay
              // + F^{-1}_{Ba}  F^d_{aC,d} U^{d;-1}_{CD} U^{ass}_{DB,k}        |  # 0, okay
              // + F^{-1}_{Ba}  F^d_{aC,k} U^{d;-1}_{CD,d} U^{ass}_{DB}        |  # 0, okay
              // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,d} U^{ass}_{DB,k}        |  # 0, okay
              // + F^{-1}_{Ba}  F^d_{aC,k} U^{d;-1}_{CD} U^{ass}_{DB,d}        |  # 0, okay
              // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,k} U^{ass}_{DB,d}        |  # 0, okay
              for (int d=0; d<NUMDISP_; ++d) {
                const int n = d / NODDISP_;
                for (int k=d; k<NUMDISP_; ++k) {  // symmetric matrix : only upper right triangle is computed 
                  const int m = k / NODDISP_;
                  double defgradbybydisp_dk = 0.0;
                  for (int B=0; B<NUMDIM_; ++B) {
                    for (int D=0; D<NUMDIM_; ++D) {
                      const int DB = voigt3x3sym[NUMDIM_*D+B];
                      const double DBfact = (D==B) ? 1.0 : 0.5;
                      double rgtstrbybydisp_DBdk = 0.0;
                      if (lin_ >= lin_half)
                        rgtstrbybydisp_DBdk = rgtstrbybydisp(DB,NUMDISP_*d+k);
                      for (int C=0; C<NUMDIM_; ++C) {
                        const int CD = voigt3x3sym[NUMDIM_*C+D];
                        const double CDfact = (C==D) ? 1.0 : 0.5;
                        const int BC = voigt3x3sym[NUMDIM_*B+C];
                        const double BCfact = (B==C) ? 1.0 : 0.5;
                        if ( (lin_ >= lin_third) and (d == iboplin[n]) )
                          defgradbybydisp_dk
                            += BCfact*invdefgradtimesboplin(BC,d) * CDfact*invrgtstrDbydisp(CD,k) * rgtstr(D,B)
                            + BCfact*invdefgradtimesboplin(BC,d) * invrgtstrD(C,D) * DBfact*rgtstrbydisp(DB,k);
                        if ( (lin_ >= lin_third) and (k == iboplin[m]) )
                          defgradbybydisp_dk
                            += BCfact*invdefgradtimesboplin(BC,k) * CDfact*invrgtstrDbydisp(CD,d) * rgtstr(D,B)
                            + BCfact*invdefgradtimesboplin(BC,k) * invrgtstrD(C,D) * DBfact*rgtstrbydisp(DB,d);
                        if ( (lin_ >= lin_third) )
                          defgradbybydisp_dk
                            += invdefgradtimesdefgradD(B,C) * CDfact*invrgtstrDbydisp(CD,d) * DBfact*rgtstrbydisp(DB,k)
                            + invdefgradtimesdefgradD(B,C) * CDfact*invrgtstrDbydisp(CD,k) * DBfact*rgtstrbydisp(DB,d);
                        if (lin_ >= lin_half)
                          defgradbybydisp_dk
                            += invdefgradtimesdefgradD(B,C) * invrgtstrD(C,D) * DBfact*rgtstrbybydisp_DBdk;
                        if (lin_ >= lin_one)
                        {
                          const double invrgtstrDbybydisp_CDdk = (*invrgtstrDbybydisp)(CD,NUMDISP_*d+k);
                          defgradbybydisp_dk
                            += invdefgradtimesdefgradD(B,C) * CDfact*invrgtstrDbybydisp_CDdk * rgtstr(D,B);
                        }
                      }
                    }
                  }
                  (*stiffmatrix)(d,k) -= defgradbybydisp_dk * effpressure*detdefgrad*detJ_w;
                  if (k != d) (*stiffmatrix)(k,d) -= defgradbybydisp_dk * effpressure*detdefgrad*detJ_w;
                }
              }

            } //  if (lin_ > lin_sixth)
          } // end block

        }  // if (ans_ == ans_none) else

        // contribute stuff containing first derivatives in displacements
        if (stab_ != stab_puredisp) {
          // AUX = (WmT + fvT*fvT') * dFv
          LINALG::Matrix<NUMDFGR_,NUMDISP_> aux;
          aux.MultiplyNN(WmT,defgradbydisp);
          // K -= dFv' * AUX * detJ * w(gp)
          stiffmatrix->MultiplyTN(-effpressure*detdefgrad*detJ_w,defgradbydisp,aux,1.0);
        }

        // derivative of incompressibility residual with respect to displacements
        // (-G) = -dFv'*fvT*pN * Fdet*detJ*wp(i);
        if (gradmatrix != NULL) {
          // AUX = fvT*pN
          LINALG::Matrix<NUMDISP_,1> aux;
          aux.MultiplyTN(defgradbydisp,tinvdefgrad);
          // contribute to G-op
          gradmatrix->MultiplyNT(-detdefgrad*detJ_w,aux,prshfct,1.0);
        }
      }
    }

    // mass matrix
    if (massmatrix != NULL) {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      for (int inod=0; inod<NUMNOD_; ++inod) {
        const double ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_; ++jnod) {
          const double massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_*inod+0,NUMDIM_*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_*inod+1,NUMDIM_*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_*inod+2,NUMDIM_*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix

    // store volume
    if (volume != NULL)
      *volume += detdefgrad*detJ_w;

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  // finish with internal force etc
  // internal force
  if (force != NULL)
  {
    // integrate internal force vector
    // fint := fint 
    //      + (B^T . sigma) * detJ * w(gp)   // already done
    //      + (-G) . ep
    if (stab_ != stab_puredisp)
      force->MultiplyNN(1.0,*gradmatrix,pres,1.0);
  }
  // incompressiblity equation
  if (incomp != NULL)
  {
    // pint := pint 
    //       - H . (Fdet - 1.0) * detJ * wp(gp)  // already done
    //       + (-Ce) . ep
    if (stab_ != stab_puredisp)
      incomp->MultiplyNN(1.0,*stabmatrix,pres,1.0);
    else
      incomp->Clear();
  }

  // fake pure-disp based approach (ANS might be active)
  if (stab_ == stab_puredisp)
  {
    if (gradmatrix != NULL)
      gradmatrix->Clear();
    if (stabmatrix != NULL)
    {
      stabmatrix->Clear();
      for (int i=0; i<NUMPRES_; ++i) (*stabmatrix)(i,i) = 1.0;
    }
  }

  // get away from here
  return;
} // DRT::ELEMENTS::So_sh8p8::ForceStiffMass


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Stress(
  LINALG::Matrix<NUMGPT_,NUMSTR_>* elestress,
  const INPAR::STR::StressType iostress,
  const int gp,
  const double& detdefgrd,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrd,
  const LINALG::Matrix<NUMSTR_,1>& glstrain,
  const LINALG::Matrix<NUMSTR_,1>& stress,
  const double& pressure
)
{
  switch (iostress)
  {
  case INPAR::STR::stress_2pk:  // 2nd Piola-Kirchhoff stress
    {
      if (elestress == NULL) dserror("stress data not available");
      // determine stress
      if (stab_ == stab_puredisp)
      {
        // store stress
        for (int i=0; i<NUMSTR_; ++i)
          (*elestress)(gp,i) = stress(i);
      }
      else
      {
        // inverted right Cauchy-Green strain tensor
        LINALG::Matrix<NUMDIM_,NUMDIM_> invcg;
        invcg.MultiplyTN(defgrd,defgrd);
        invcg.Invert();
        LINALG::Matrix<NUMSTR_,1> invcgv;
        Matrix2TensorToVector6Voigt(invcgv,invcg);
        // store stress
        for (int i=0; i<NUMSTR_; ++i)
          (*elestress)(gp,i) = stress(i) - pressure*detdefgrd*invcgv(i);
      }
    }
    break;
  case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      // pull back
      LINALG::Matrix<NUMSTR_,NUMSTR_> defgraddefgradT;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(defgraddefgradT,defgrd,
                                                    true,voigt6_stress,voigt6_stress);
      // (deviatoric) Cauchy stress vector
      LINALG::Matrix<NUMSTR_,1> cauchyv;
      cauchyv.MultiplyNN(1.0/detdefgrd,defgraddefgradT,stress);


      // determine stress
      if (stab_ == stab_puredisp)
      {
        // store stress
        for (int i=0; i<NUMSTR_; ++i)
          (*elestress)(gp,i) = cauchyv(i);
      }
      else
      {
        // above computed #cauchyv is deviatoric true stress
        // isochoric Cauchy stress vector
        LINALG::Matrix<NUMSTR_,1> isocauchyv(true);
        for (int i=0; i<NUMDIM_; ++i) isocauchyv = -pressure;
        // store stress
        for (int i=0; i<NUMSTR_; ++i)
          (*elestress)(gp,i) = cauchyv(i) + isocauchyv(i);
      }
    }
    break;
  case INPAR::STR::stress_none:
    break;
  default:
    dserror("requested stress option not available");
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Strain(
  LINALG::Matrix<NUMGPT_,NUMSTR_>* elestrain,  ///< store the strain herein
  const INPAR::STR::StrainType iostrain,
  const int gp,  ///< Gauss point index
  const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrd,  ///< (assumed) deformation gradient
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& invdefgrd,  ///< (assumed) inverted deformation gradient
  const LINALG::Matrix<NUMSTR_,1>& glstrain  ///< Green-Lagrange strain vector
  )
{
  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // store
      for (int i=0; i<NUMDIM_; ++i)
        (*elestrain)(gp,i) = glstrain(i);
      for (int i=NUMDIM_; i<NUMSTR_; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
    }
    break;
  case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // create push forward 6x6 matrix
      LINALG::Matrix<NUMSTR_,NUMSTR_> invdefgradTdefgrad;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(invdefgradTdefgrad,invdefgrd,
                                                    false,voigt6_strain,voigt6_strain);
      // push forward
      LINALG::Matrix<NUMSTR_,1> eastrain;
      eastrain.MultiplyNN(invdefgradTdefgrad,glstrain);
      // store
      for (int i=0; i<NUMDIM_; ++i)
        (*elestrain)(gp,i) = eastrain(i);
      for (int i=NUMDIM_; i<NUMSTR_; ++i)
        (*elestrain)(gp,i) = 0.5 * eastrain(i);
    }
    break;
  case INPAR::STR::strain_none:
    break;
  default:
    dserror("requested strain option not available");
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AssDefGrad(
  double& detdefgrad,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrad,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& invdefgrad,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& rgtstr,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& defgradD,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& rgtstrD,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& invrgtstrD,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& Jinv,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& Jac,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& jac,
  const LINALG::Matrix<NUMSTR_,1>& glstrain
  )
{
/*
  // inverse material Jacobian (X_{,xi})^T
  LINALG::Matrix<NUMDIM_,NUMDIM_> JinvX(Jac);
  double Jdet = JinvX.Invert();  // (X_{,xi})^{-T}
  if (Jdet < 0.0) dserror("Trouble during inversion of Jacobian");
*/

  // pure displacement-based deformation gradient
  // F = x_{,X} = x_{,xi} . xi_{,X} = x_{,xi} . (X_{,xi})^{-1} = jac^T . Jinv^T
  defgradD.MultiplyTT(jac,Jinv);
  
  // pure displacement-based right Cauchy-Green strain
  LINALG::Matrix<NUMDIM_,NUMDIM_> cgD;
  cgD.MultiplyTN(defgradD,defgradD);

  // rotation matrix in pure displacement based deformation gradient
  // and pure disp-based material stretch tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> rot(true);
  {
#if 0
    LINALG::Matrix<NUMDIM_,NUMDIM_> nd(true);
    LINALG::Matrix<NUMDIM_,NUMDIM_> lamd(true);
    LINALG::Matrix<NUMDIM_,NUMDIM_> NdT(true);
    LINALG::SVD(defgradD,nd,lamd,NdT);
    rot.MultiplyNN(nd,NdT);
    // pure disp-based material stretch tensor
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyTN(NdT,lamd);
    rgtstrD.MultiplyNN(aux,NdT);
#else
    // spectral decomposition of disp-based right Cauchy-Green tensor
    LINALG::Matrix<NUMDIM_,NUMDIM_> NdT;
    LINALG::Matrix<NUMDIM_,NUMDIM_> lamd;
#if 0
    LINALG::Matrix<NUMDIM_,NUMDIM_> Nd(true);
    LINALG::SVD(cgD,NdT,lamd,Nd);
    // spectral composition of disp-based right stretch tensor
    for (int i=0; i<NUMDIM_; ++i) lamd(i,i) = sqrt(lamd(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NdT,lamd);
    rgtstrD.MultiplyNN(aux,Nd);
#else
    LINALG::SYEV(cgD,lamd,NdT);
    // spectral composition of disp-based right stretch tensor
    for (int i=0; i<NUMDIM_; ++i) lamd(i,i) = sqrt(lamd(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NdT,lamd);
    rgtstrD.MultiplyNT(aux,NdT);
#endif
    // inverse disp-based right stretch tensor
    invrgtstrD.Update(rgtstrD);
    const double detrgtstrD = invrgtstrD.Invert();
    if (detrgtstrD < 0.0) dserror("Trouble during inversion of right stretch tensor");
    // rotation matrix
    rot.MultiplyNN(defgradD,invrgtstrD);
#endif
  }

  // assumed material stretch tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> invrgtstr;
  {
    LINALG::Matrix<NUMDIM_,NUMDIM_> cga;
    for (int i=0; i<NUMDIM_; ++i) cga(i,i) = 2.0*glstrain(i,0) + 1.0;
    // off-diagonal terms are already twice in the Voigt-GLstrain-vector
    cga(0,1) = cga(1,0) = glstrain(3);
    cga(1,2) = cga(2,1) = glstrain(4);
    cga(0,2) = cga(2,0) = glstrain(5);
    LINALG::Matrix<NUMDIM_,NUMDIM_> lama(true);
    LINALG::Matrix<NUMDIM_,NUMDIM_> NaT(true);
#if 0
    LINALG::Matrix<NUMDIM_,NUMDIM_> Na(true);
    LINALG::SVD(cga,NaT,lama,Na);
    for (int i=0; i<NUMDIM_; ++i) lama(i,i) = sqrt(lama(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NaT,lama);
    rgtstr.MultiplyNN(aux,Na);
#else
    LINALG::SYEV(cga,lama,NaT);
    for (int i=0; i<NUMDIM_; ++i) lama(i,i) = sqrt(lama(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NaT,lama);
    rgtstr.MultiplyNT(aux,NaT);
#endif
#if 0
    invrgtstr.Update(rgtstr);
    detdefgrad = invrgtstr.Invert();
#else
    invrgtstr.Clear();
    detdefgrad = 1.0;
    for (int al=0; al<NUMDIM_; ++al)
    {
      detdefgrad *= lama(al,al);
      for (int j=0; j<NUMDIM_; ++j)
        for (int i=0; i<NUMDIM_; ++i)
          invrgtstr(i,j) += NaT(i,al)*NaT(j,al)/lama(al,al);
    }
#endif
  }

  // assumed deformation gradient
  defgrad.MultiplyNN(rot,rgtstr);

  // inverse of assumed deformation gradient
  invdefgrad.MultiplyNT(invrgtstr,rot);

  // done
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::ELEMENTS::So_sh8p8::ShearMod() const
{
  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
  case INPAR::MAT::m_stvenant: /*-------- st.venant-kirchhoff-material */
  {
    MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
    return stvk->ShearMod();
    break;
  }
  case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
  {
    MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
    return neo->ShearMod();
    break;
  }
  case INPAR::MAT::m_aaaneohooke: /*-----------AAA NeoHookean Material */
  {
    MAT::AAAneohooke* aaaneo = static_cast<MAT::AAAneohooke*>(mat.get());
    return aaaneo->ShearMod();
    break;
  }
  case INPAR::MAT::m_viscoanisotropic: /*-----------AAA NeoHookean Material */
  {
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    return visco->ShearMod();
    break;
  }
  case INPAR::MAT::m_yeoh: /*-----------AAA NeoHookean Material */
  {
    MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(Material().get());
    return yeoh->ShearMod();
    break;
  }
  case INPAR::MAT::m_visconeohooke: /*-----------AAA NeoHookean Material */
  {
    MAT::ViscoNeoHooke* visconeo = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    return visconeo->ShearMod();
    break;
  }
  default:
    dserror("Cannot ask material for shear modulus");
    break;
  } // switch (mat->MaterialType())

  return 0;
}




/*======================================================================*/
/*======================================================================*/
/*======================================================================*/
/*======================================================================*/


/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Sosh8p8Register::Initialize(DRT::Discretization& dis)
{
  //sosh8_gmshplotdis(dis);

  int num_morphed_so_hex8_easnone = 0;

  // Loop through all elements
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8p8) continue;
    DRT::ELEMENTS::So_sh8p8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8p8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8p8* failed");

    if (!actele->nodes_rearranged_) {
      // check for automatic definition of thickness direction
      if (actele->thickdir_ == DRT::ELEMENTS::So_sh8p8::autoj) {
        actele->thickdir_ = actele->sosh8_findthickdir();
      }

      int new_nodeids[DRT::ELEMENTS::So_sh8p8::NUMNOD_];

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
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
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
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::autot:
      case DRT::ELEMENTS::So_sh8::globz: {
        // no resorting necessary
        for (int node = 0; node < 8; ++node) {
          new_nodeids[node] = actele->NodeIds()[node];
        }
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::undefined: {
        // here comes plan B: morph So_sh8p8 to So_hex8
        actele->SetANS(DRT::ELEMENTS::So_sh8p8::ans_none);
        num_morphed_so_hex8_easnone++;
        break;
      }
      case DRT::ELEMENTS::So_sh8::none: break;
      default:
        dserror("no thickness direction for So_sh8p8");
      }
      //actele->sosh8p8_gmshplotlabeledelement(actele->NodeIds());
    }
  }

  if (num_morphed_so_hex8_easnone>0){
    std::cout << endl << num_morphed_so_hex8_easnone
              << " Sosh8p8-Elements have no clear 'thin' direction and ANS is disabled!"
              << endl;
  }

  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false,false,false);

  // loop again to init Jacobian for Sosh8p8's
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8p8) continue;
    DRT::ELEMENTS::So_sh8p8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8p8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8p8* failed");
    actele->InitJacobianMapping();  // this sets #invJ_ in So_hex8
  }

  // **************** debug printout ot gmesh **********************************
  //sosh8_gmshplotdis(dis);

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
