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
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
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
  LINALG::Matrix<NUMDOF_SOSH8P8,NUMDOF_SOSH8P8> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOSH8P8,NUMDOF_SOSH8P8> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOSH8P8,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOSH8P8,1> elevec2(elevec2_epetra.A(),true);
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
      LINALG::Matrix<NUMDISP_SOSH8P8,1> mydisp(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,1> mypres(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8> gradmatrix(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8> stabmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,1> force(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
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
      LINALG::Matrix<NUMDISP_SOSH8P8,1> mydisp;
      LINALG::Matrix<NUMPRES_SOSH8P8,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8> gradmatrix(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8> stabmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,1> force(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
      BuildElementMatrix(&elemat1,&stiffmatrix,&gradmatrix,NULL,&stabmatrix);
      BuildElementVector(&elevec1,&force,&incomp);
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
      LINALG::Matrix<NUMDISP_SOSH8P8,1> mydisp;
      LINALG::Matrix<NUMPRES_SOSH8P8,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8> gradmatrix(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8> stabmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,1> force(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
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
      LINALG::Matrix<NUMDISP_SOSH8P8,1> mydisp;
      LINALG::Matrix<NUMPRES_SOSH8P8,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8> massmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8> gradmatrix(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8> stabmatrix(true);
      LINALG::Matrix<NUMDISP_SOSH8P8,1> force(true);
      LINALG::Matrix<NUMPRES_SOSH8P8,1> incomp(true);
      ForceStiffMass(lm,mydisp,mypres,
                     &massmatrix,&stiffmatrix,&gradmatrix,&stabmatrix,&force,&incomp,
                     NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
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
      LINALG::Matrix<NUMDISP_SOSH8P8,1> mydisp;
      LINALG::Matrix<NUMPRES_SOSH8P8,1> mypres;
      ExtractDispAndPres(mystat,mydisp,mypres);
      LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8> stress;
      LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8> strain;
      INPAR::STR::StressType iostress 
        = params.get<INPAR::STR::StressType>("iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain 
        = params.get<INPAR::STR::StrainType>("iostrain", INPAR::STR::strain_none);
      ForceStiffMass(lm,mydisp,mypres,
                     NULL,NULL,NULL,NULL,NULL,NULL,
                     &stress,&strain,params,iostress,iostrain);
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:{

      const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_SOSH8P8,NUMSTR_SOSH8P8> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        std::vector<int> numadjele(NUMNOD_SOSH8P8);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOSH8P8;++i)
        {
          DRT::Node* node = nodes[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOSH8P8;++i)
        {
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOSH8P8;++i)
        {
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz") {
        Teuchos::RCP<Epetra_MultiVector> elestress=params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
        if (elestress==Teuchos::null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_SOSH8P8; ++i)
          {
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOSH8P8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOSH8P8;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_SOSH8P8,NUMSTR_SOSH8P8> nodalstresses;
        soh8_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        std::vector<int> numadjele(NUMNOD_SOSH8P8);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOSH8P8;++i){
          DRT::Node* node=nodes[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOSH8P8;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOSH8P8;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        Teuchos::RCP<Epetra_MultiVector> elestress=params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
        if (elestress==Teuchos::null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_SOSH8P8; ++i)
          {
            double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOSH8P8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOSH8P8;
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
  const LINALG::Matrix<NUMDISP_SOSH8P8,1>& disp,           // current displacements
  const LINALG::Matrix<NUMPRES_SOSH8P8,1>& pres,       // current pressures
  LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* massmatrix,  // element mass matrix
  LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* stiffmatrix, // element stiffness matrix
  LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* gradmatrix, // element gradient matrix
  LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8>* stabmatrix,  // element stabilisation matrix
  LINALG::Matrix<NUMDISP_SOSH8P8,1>* force,                 // element internal force vector
  LINALG::Matrix<NUMPRES_SOSH8P8,1>* incomp,   // incompressibility residual
  LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestress,   // stresses at GP
  LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestrain,   // strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress, // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOSH8P8,1> > shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOSH8P8,NUMNOD_SOSH8P8> > derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOSH8P8,NUMDIM_SOSH8P8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOSH8P8,NUMDIM_SOSH8P8> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOSH8P8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp(i*NODDISP_SOSH8P8+0,0);
    xcurr(i,1) = xrefe(i,1) + disp(i*NODDISP_SOSH8P8+1,0);
    xcurr(i,2) = xrefe(i,2) + disp(i*NODDISP_SOSH8P8+2,0);
  }

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  //LINALG::Matrix<NUMANS_SOSH8P8*NUMSP_SOSH8P8,NUMDOF_SOSH8P8> B_ans_loc(true); //set to 0
  LINALG::Matrix<NUMANS_SOSH8P8*NUMSP_SOSH8P8,NUMDISP_SOSH8P8> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> > jac_sps(NUMSP_SOSH8P8);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> > jac_cur_sps(NUMSP_SOSH8P8);
  // pointer to derivs evaluated at all sampling points
  std::vector<LINALG::Matrix<NUMDIM_SOSH8P8,NUMNOD_SOSH8P8> >* deriv_sp = NULL;
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe,xcurr,&deriv_sp,jac_sps,jac_cur_sps,B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double r[NUMGPT_SOSH8P8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
  const double s[NUMGPT_SOSH8P8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};

  // ---------------------------------------------------------------------
  // first loop over Gauss point
  // stabilisation matrices
  LINALG::Matrix<1,1> estabd(true); // element volume
  LINALG::Matrix<1,NUMNOD_SOSH8P8> estabe(true); // integral of pressure shape functions across element
  LINALG::Matrix<NUMNOD_SOSH8P8,NUMNOD_SOSH8P8> estabm(true);  // mass-like matrix
  std::vector<LINALG::Matrix<1,1> > estaba(NUMGPT_SOSH8P8);  // shape functions for projected Q0/constant pressure
  for (int gp=0; gp<NUMGPT_SOSH8P8; ++gp)
  {
    // (transposed) material-to-parametric Jacobian J = (X_{,xi})^T
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> Jac;
    Jac.Multiply(derivs[gp],xrefe);
    double detJ = Jac.Determinant();

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
    // Cem = 1./shearmod*( Mem - Eem'*inv(Dem)*Eem );
    stabmatrix->Update(estabm);
    stabmatrix->MultiplyTN(-1.0/estabd(0,0),estabe,estabe,1.0);
    stabmatrix->Scale(-1.0/shearmod);
  }
 

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOSH8P8; ++gp)
  {

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> jac;
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
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> jac_cur;
    jac_cur.Multiply(derivs[gp],xcurr);

    // set up B-Operator in local(parameter) element space including ANS
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> bop_loc;
    for (int inode = 0; inode < NUMNOD_SOSH8P8; ++inode) {
      for (int dim = 0; dim < NUMDIM_SOSH8P8; ++dim) {
        // B_loc_rr = N_r.X_r
        bop_loc(0,inode*3+dim) = derivs[gp](0,inode) * jac_cur(0,dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1,inode*3+dim) = derivs[gp](1,inode) * jac_cur(1,dim);
        if (ans_ == ans_none)
          // B_loc_tt = N_t.X_t
          bop_loc(2,inode*3+dim) = derivs[gp](2,inode) * jac_cur(2,dim);
        else
          // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
          //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
          //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
          bop_loc(2,inode*3+dim) = 0.25*(1-r[gp])*(1-s[gp]) * B_ans_loc(0+4*NUMANS_SOSH8P8,inode*3+dim)
                                 + 0.25*(1+r[gp])*(1-s[gp]) * B_ans_loc(0+5*NUMANS_SOSH8P8,inode*3+dim)
                                 + 0.25*(1+r[gp])*(1+s[gp]) * B_ans_loc(0+6*NUMANS_SOSH8P8,inode*3+dim)
                                 + 0.25*(1-r[gp])*(1+s[gp]) * B_ans_loc(0+7*NUMANS_SOSH8P8,inode*3+dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3,inode*3+dim) = derivs[gp](0,inode) * jac_cur(1,dim)
                                +derivs[gp](1,inode) * jac_cur(0,dim);
        if (ans_ == ans_none)
          // B_loc_st = N_s.X_t + N_t.X_s
          bop_loc(4,inode*3+dim) = derivs[gp](1,inode) * jac_cur(2,dim)
                                 + derivs[gp](2,inode) * jac_cur(1,dim);
        else
          // B_loc_st = interpolation along r of ANS B_loc_st
          //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
          bop_loc(4,inode*3+dim) = 0.5*(1.0+r[gp]) * B_ans_loc(1+1*NUMANS_SOSH8P8,inode*3+dim)
                                 + 0.5*(1.0-r[gp]) * B_ans_loc(1+3*NUMANS_SOSH8P8,inode*3+dim);
        if (ans_ == ans_none)
          // B_loc_rt = N_r.X_t + N_t.X_r
          bop_loc(5,inode*3+dim) = derivs[gp](0,inode) * jac_cur(2,dim)
                                 + derivs[gp](2,inode) * jac_cur(0,dim);
        else
          // B_loc_rt = interpolation along s of ANS B_loc_rt
          //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
          bop_loc(5,inode*3+dim) = 0.5*(1.0-s[gp]) * B_ans_loc(2+0*NUMANS_SOSH8P8,inode*3+dim)
                                 + 0.5*(1.0+s[gp]) * B_ans_loc(2+2*NUMANS_SOSH8P8,inode*3+dim);
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> TinvT;
    sosh8_evaluateT(jac,TinvT);
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> bop;
    bop.Multiply(TinvT,bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    // but with modified ANS strains E33, E23 and E13
    LINALG::Matrix<NUMSTR_SOSH8P8,1> lstrain;
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
      for (int dim = 0; dim < NUMDIM_SOSH8P8; ++dim) {
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
    LINALG::Matrix<NUMSTR_SOSH8P8,1> glstrain;
    glstrain.Multiply(TinvT,lstrain);

    // recover deformation gradient incoperating assumed natural GL strain
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> defgrad;  // assumed def.grad.
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> rgtstr;  // assumed material stretch
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> defgradD;  // pure disp-based def.grad.
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> rgtstrD;  // pure disp.-based material stretch
    AssDefGrad(defgrad,rgtstr,defgradD,rgtstrD,invJ_[gp],jac,jac_cur,glstrain);

    // inverse of deformation gradient and its derivative 
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> invdefgrad(defgrad);
    double detdefgrad = invdefgrad.Invert();

    // return gp strains if necessary
    if (iostrain != INPAR::STR::strain_none)
      Strain(elestrain,iostrain,gp,detdefgrad,defgrad,invdefgrad,glstrain);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> cmat(true);
    LINALG::Matrix<NUMSTR_SOSH8P8,1> stress(true);
    soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrad,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // pressure at Gauss point
    double pressure = (shapefcts[gp]).Dot(pres);

    // return Gauss point stresses if necessary
    if (iostress != INPAR::STR::stress_none)
      Stress(elestress,iostress,gp,detdefgrad,defgrad,glstrain,stress,pressure);

    // effective shape function of scalar pressure field at current Gauss point
    LINALG::Matrix<NUMPRES_SOSH8P8,1> prshfct;
    if (stab_ == stab_nonaffine)
      prshfct.MultiplyTN(1.0/estabd(0,0),estabe,estaba[gp]);  //???????????
    else if (stab_ == stab_affine)
      prshfct.Update(shapefcts[gp]);
 
    // integration factor
    double detJ_w = detJ*gpweights[gp];

    // internal force
    if (force != NULL)
    {
      // integrate internal force vector
      // fint := fint 
      //      + (B^T . sigma) * detJ * w(gp) 
      //      + G . ep   // will be done _after_ Gauss point loop
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
      //force->MultiplyNN(1.0, *gradmatrix, pres, 1.0)
    }
    // incompressiblity equation
    if (incomp != NULL)
    {
      // pint := pint 
      //       - He . (Fdet - 1.0) * detJ * wp(gp)
      //       + Ce . ep   // will be done _after_ Gauss point loop
      incomp->Update(-(detdefgrad-1.0)*detJ_w,prshfct,1.0);
    }
    // stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<NUMSTR_SOSH8P8, NUMDISP_SOSH8P8> cb;
      cb.Multiply(cmat,bop); // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // intergrate `geometric' stiffness matrix and add to keu
      // here also the ANS interpolation comes into play
      for (int inod=0; inod<NUMNOD_SOSH8P8; ++inod)
      {
        for (int jnod=0; jnod<NUMNOD_SOSH8P8; ++jnod)
        {
          LINALG::Matrix<NUMSTR_SOSH8P8,1> G_ij;
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
          LINALG::Matrix<NUMSTR_SOSH8P8,1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);
            
          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_SOSH8P8*inod+0, NUMDIM_SOSH8P8*jnod+0) += Gij;
          (*stiffmatrix)(NUMDIM_SOSH8P8*inod+1, NUMDIM_SOSH8P8*jnod+1) += Gij;
          (*stiffmatrix)(NUMDIM_SOSH8P8*inod+2, NUMDIM_SOSH8P8*jnod+2) += Gij;
        }
      } // end of intergrate `geometric' stiffness

      // add (incomplete) derivative of pressure-proportional force w.r.t. displacements
      // Kp = Kp + dFv'*fvT*(pN*ep')*fvT'*dFv * Fdet*detJ*wp(i) ...
      //          + (pN*ep')*dFv'*WmT*dFv * Fdet*detJ*wp(i)
      // Ke = Keu + Kg - Kp;
      {
        // effective pressure at Gauss point
        const double effpressure = prshfct.Dot(pres);  // pN*ep'
        // Voigt 9-vector of transposed & inverted deformation gradient fvT := F^{-T}
        LINALG::Matrix<NUMDFGR_SOSH8P8,1> tinvdefgrad;
        Matrix2TensorToVector9Voigt(tinvdefgrad,invdefgrad,true);
        // derivative of WmT := F^{-T}_{,F} in Voigt vector notation
        LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8> WmT;
        InvVector9VoigtDiffByItself(WmT,invdefgrad,true);
        // WmT := WmT + fvT*fvT'
        WmT.MultiplyNT(1.0,tinvdefgrad,tinvdefgrad,1.0);
        // derivative of def.grad. w.r.t. displacements Fv_{,d}
        LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDISP_SOSH8P8> defgradbydisp;
        {
          // Voigt 9-vector indices
          const int* voigtrow9 = NULL;
          const int* voigtcol9 = NULL;
          Indices9VoigtTo2Tensor(voigtrow9,voigtcol9);
          const int* voigt3x3 = NULL;
          Indices2TensorTo9Voigt(voigt3x3);  // access is via (i,j) -> 3*i+j
          const int* voigt3x3sym = NULL;
          Indices2TensorTo6Voigt(voigt3x3sym);  // access is via (i,j) -> 3*i+j

          // material derivatives of shape functions
          LINALG::Matrix<NUMDIM_SOSH8P8,NUMNOD_SOSH8P8> derivsmat;
          derivsmat.MultiplyNN(invJ_[gp],derivsmat);

          // linear B-op
          // derivative of displ-based def.grad with respect to nodal displacements
          // F^d_{aC,d}
          LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDISP_SOSH8P8> boplin(true);
          for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
          {
            const int i = voigtrow9[I];
            const int j = voigtcol9[I];
            for (int k=0; k<NUMNOD_SOSH8P8; ++k)
            {
              const int K = j + k*NUMDIM_SOSH8P8;
              boplin(I,K) = derivsmat(j,i);
            }
          }

          // inverse of pure disp-based material stretch tensor
          // U^{d;-1}
          LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> invrgtstrD(rgtstrD);
          double detrgtstrD = invrgtstrD.Invert();
          if (detrgtstrD < 0.0) dserror("Trouble in inverting right stretch tensor");
          
          // derivative of pure-disp. inverse material stretch tensor with respect to nodal displacements
          // U^{d;-1}_{,d} = U^{d;-1}_{,U} . (C^d_{,U^d})^{-1} . C^d_{,d}
          LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> invrgtstrDbydisp;
          {
             // U^{d;-1}_{,U}
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> invrgtstrDbyrgtstrD;
            InvVector6VoigtDiffByItself(invrgtstrDbyrgtstrD,invrgtstrD);

            // C^d_{,U^d} = (U^d . U^d)_{,U^d}
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> rcgDbyrgtstrD;
            SqVector6VoigtDiffByItself(rcgDbyrgtstrD,rgtstrD);

            // deformation gradient as Voigt matrix
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMDFGR_SOSH8P8> defgradm;
            Matrix2TensorToMatrix6x9Voigt(defgradm,defgrad);

            // C^d_{,d} = 2 * Fm * Boplin, 6x24
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> rcgDbydisp;
            rcgDbydisp.MultiplyNN(2.0,defgradm,boplin);
            
            // AUX = (C^d_{,U^d})^{-1} . C^d_{,d}
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> aux;
            {
              LINALG::FixedSizeSerialDenseSolver<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> asolver;
              asolver.SetMatrix(rcgDbyrgtstrD);
              asolver.SetVectors(aux,rcgDbydisp);
              const int err = asolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
            }

            // U^{d;-1}_{,d} = U^{d;-1}_{,U} . AUX
            invrgtstrDbydisp.Multiply(invrgtstrDbyrgtstrD,aux);
          }

          // derivative of ass. mat. stretch tensor with respect to nodal displacements
          // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
          LINALG::Matrix<NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> rgtstrbydisp;
          {
            // derivative of ass. right Cauchy-Green with respect to ass. material stretch tensor
            // C^{ass}_{,U^{ass}}
            LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> rcgbyrgtstr;
            SqVector6VoigtDiffByItself(rcgbyrgtstr,rgtstr);

            // C^{ass}_{,d} = 2 * bop

            // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
            {
              LINALG::FixedSizeSerialDenseSolver<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8,NUMDISP_SOSH8P8> asolver;
              asolver.SetMatrix(rcgbyrgtstr);
              asolver.SetVectors(rgtstrbydisp,bop);
              const int err = asolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
            }
            rgtstrbydisp.Scale(2.0);
          }
          
          // derivative of def.grad. with respect to k nodal displacements d^k
          // F_{aB,k} = F^d_{aC,k} . U^{d;-1}_{CD} . U^{ass}_{DB}
          //          + F^d_{aC} . U^{d;-1}_{CD,k} . U^{ass}_{DB}
          //          + F^d_{aC} . U^{d;-1}_{CD} . U^{ass}_{DB,k}
          for (int ab=0; ab<NUMDFGR_SOSH8P8; ++ab)
          {
            for (int k=0; k<NUMDISP_SOSH8P8; ++k)
            {
              double defgradbydispabk = 0.0;
              const int a = voigtrow9[ab];
              const int b = voigtcol9[ab];
              for (int c=0; c<NUMDIM_SOSH8P8; ++c) {
                for (int d=0; d<NUMDIM_SOSH8P8; ++d) {
                  const int ac = voigt3x3[NUMDIM_SOSH8P8*a+c];
                  const int cd = voigt3x3sym[NUMDIM_SOSH8P8*c+d];
                  const int db = voigt3x3sym[NUMDIM_SOSH8P8*d+b];
                  const double cdfact = (c!=d) ? 0.5 : 1.0;
                  const double dbfact = (d!=b) ? 0.5 : 1.0;
                  defgradbydispabk += boplin(ac,k) * invrgtstrD(c,d) * rgtstr(d,b)
                    + defgradD(a,c) * cdfact*invrgtstrDbydisp(cd,k) * rgtstr(d,b)
                    + defgradD(a,c) * invrgtstrD(c,d) * dbfact*rgtstrbydisp(db,k);
                }
              }
              defgradbydisp(ab,k) = defgradbydispabk;
            }
          }
        }

        // finally contribute
        if (stab_ != stab_puredisp) {
          // AUX = (WmT + fvT*fvT') * dFv
          LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDISP_SOSH8P8> aux(true);
          aux.MultiplyNN(1.0,WmT,defgradbydisp,0.0);
          // K -= dFv' * AUX * detJ * w(gp)
          stiffmatrix->MultiplyTN(-effpressure*detdefgrad*detJ_w,defgradbydisp,aux,1.0);
        }

        // derivative of incompressibility residual with respect to displacements
        // G = dFv'*fvT*pN * Fdet*detJ*wp(i);
        if (gradmatrix != NULL) {
          // AUX = fvT*pN
          LINALG::Matrix<NUMDISP_SOSH8P8,1> aux;
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
      for (int inod=0; inod<NUMNOD_SOSH8P8; ++inod) {
        const double ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOSH8P8; ++jnod) {
          const double massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOSH8P8*inod+0,NUMDIM_SOSH8P8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOSH8P8*inod+1,NUMDIM_SOSH8P8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOSH8P8*inod+2,NUMDIM_SOSH8P8*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix
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
    //      + G . ep
    if (stab_ != stab_puredisp)
      force->MultiplyNN(1.0,*gradmatrix,pres,1.0);
  }
  // incompressiblity equation
  if (incomp != NULL)
  {
    // pint := pint 
    //       - H . (Fdet - 1.0) * detJ * wp(gp)  // already done
    //       + Ce . ep
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
      for (int i=0; i<NUMPRES_SOSH8P8; ++i) (*stabmatrix)(i,i) = 1.0;
    }
  }


  return;
} // DRT::ELEMENTS::So_sh8p8::ForceStiffMass


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Stress(
  LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestress,
  const INPAR::STR::StressType iostress,
  const int gp,
  const double& detdefgrd,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrd,
  const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain,
  const LINALG::Matrix<NUMSTR_SOSH8P8,1>& stress,
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
        for (int i=0; i<NUMSTR_SOSH8P8; ++i)
          (*elestress)(gp,i) = stress(i);
      }
      else
      {
        // inverted right Cauchy-Green strain tensor
        LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> invcg;
        invcg.MultiplyTN(defgrd,defgrd);
        invcg.Invert();
        LINALG::Matrix<NUMSTR_SOSH8P8,1> invcgv;
        Matrix2TensorToVector6Voigt(invcgv,invcg);
        // store stress
        for (int i=0; i<NUMSTR_SOSH8P8; ++i)
          (*elestress)(gp,i) = stress(i) - pressure*detdefgrd*invcgv(i);
      }
    }
    break;
  case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      // pull back
      LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> defgraddefgradT;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(defgraddefgradT,defgrd,
                                                    true,voigt6_stress,voigt6_stress);
      // (deviatoric) Cauchy stress vector
      LINALG::Matrix<NUMSTR_SOSH8P8,1> cauchyv;
      cauchyv.MultiplyNN(1.0/detdefgrd,defgraddefgradT,stress);


      // determine stress
      if (stab_ == stab_puredisp)
      {
        // store stress
        for (int i=0; i<NUMSTR_SOSH8P8; ++i)
          (*elestress)(gp,i) = cauchyv(i);
      }
      else
      {
        // above computed #cauchyv is deviatoric true stress
        // isochoric Cauchy stress vector
        LINALG::Matrix<NUMSTR_SOSH8P8,1> isocauchyv(true);
        for (int i=0; i<NUMDIM_SOSH8P8; ++i) isocauchyv = -pressure;
        // store stress
        for (int i=0; i<NUMSTR_SOSH8P8; ++i)
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
  LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestrain,  ///< store the strain herein
  const INPAR::STR::StrainType iostrain,
  const int gp,  ///< Gauss point index
  const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrd,  ///< (assumed) deformation gradient
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invdefgrd,  ///< (assumed) inverted deformation gradient
  const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain  ///< Green-Lagrange strain vector
  )
{
  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // store
      for (int i=0; i<NUMDIM_SOSH8P8; ++i)
        (*elestrain)(gp,i) = glstrain(i);
      for (int i=NUMDIM_SOSH8P8; i<NUMSTR_SOSH8P8; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
    }
    break;
  case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // create push forward 6x6 matrix
      LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8> invdefgradTdefgrad;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(invdefgradTdefgrad,invdefgrd,
                                                    false,voigt6_strain,voigt6_strain);
      // push forward
      LINALG::Matrix<NUMSTR_SOSH8P8,1> eastrain;
      eastrain.MultiplyNN(invdefgradTdefgrad,glstrain);
      // store
      for (int i=0; i<NUMDIM_SOSH8P8; ++i)
        (*elestrain)(gp,i) = eastrain(i);
      for (int i=NUMDIM_SOSH8P8; i<NUMSTR_SOSH8P8; ++i)
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
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrad,
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& rgtstr,
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgradD,
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& rgtstrD,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& Jinv,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& Jac,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& jac,
  const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain
  )
{
  // inverse material Jacobian (X_{,xi})^T
//  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> Jinv(Jac);
//  double Jdet = Jinv.Invert();  // (X_{,xi})^{-T}

  // pure displacement-based deformation gradient
  // F = x_{,X} = x_{,xi} . xi_{,X} = x_{,xi} . (X_{,xi})^{-1} = jac^T . Jinv^T
  defgradD.MultiplyTT(jac,Jinv);
//cout << defgradD << endl;
  
  // pure displacement-based right Cauchy-Green strain
//  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> cgD;
//  cgD.MultiplyTN(dgd,dgd);

  // rotation matrix in pure displacement based deformation gradient
  // and pure disp-based material stretch tensor
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> rot;
  {
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> nd(true);
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> lamd(true);
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> NdT(true);
    LINALG::SVD(defgradD,nd,lamd,NdT);
    rot.MultiplyNN(nd,NdT);
    // pure disp-based material stretch tensor
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> aux;
    aux.MultiplyNN(NdT,lamd);
    rgtstrD.MultiplyNT(aux,NdT);
  }

  // assumed material stretch tensor
  {
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> cga;
    for (int i=0; i<NUMDIM_SOSH8P8; ++i) cga(i,i) = 2.0*glstrain(i,0) + 1.0;
    // off-diagonal terms are already twice in the Voigt-GLstrain-vector
    cga(0,1) = glstrain(3);  cga(1,0) = glstrain(3);
    cga(1,2) = glstrain(4);  cga(2,1) = glstrain(4);
    cga(0,2) = glstrain(5);  cga(2,0) = glstrain(5);
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> NaT(true);
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> lama(true);
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> Na(true);
    LINALG::SVD(cga, NaT, lama, Na);
    for (int i=0; i<NUMDIM_SOSH8P8; ++i) lama(i,i) = sqrt(lama(i,i));
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> aux;
    aux.MultiplyNN(NaT,lama);
    rgtstr.MultiplyNN(aux,Na);
  }

  // assumed deformation gradient
  defgrad.MultiplyNN(rot,rgtstr);

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
  default:
    dserror("Cannot ask material for shear modulus");
    break;
  } // switch (mat->MaterialType())

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices6VoigtTo2Tensor(
  const int*& voigtrow6,
  const int*& voigtcol6,
  const bool transpose
  )
{
  const int VoigtRow6[NUMSTR_SOSH8P8] = {0,1,2, 0,1,2};
  const int VoigtCol6[NUMSTR_SOSH8P8] = {0,1,2, 1,2,0};

  if (transpose)
  {
    voigtrow6 = &(VoigtCol6[0]);
    voigtcol6 = &(VoigtRow6[0]);
  }
  else
  {
    voigtrow6 = &(VoigtRow6[0]);
    voigtcol6 = &(VoigtCol6[0]);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices9VoigtTo2Tensor(
  const int*& voigtrow9,
  const int*& voigtcol9,
  const bool transpose
  )
{
  // 9-Voigt C-index                      0 1 2  3 4 5  6 7 8
  const int VoigtRow9[NUMDFGR_SOSH8P8] = {0,1,2, 0,1,2, 0,2,1};
  const int VoigtCol9[NUMDFGR_SOSH8P8] = {0,1,2, 1,2,0, 2,1,0};

  if (transpose)
  {
    voigtrow9 = &(VoigtCol9[0]);
    voigtcol9 = &(VoigtRow9[0]);
  }
  else
  {
    voigtrow9 = &(VoigtRow9[0]);
    voigtcol9 = &(VoigtCol9[0]);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices2TensorTo9Voigt(
  const int*& voigt3x3
  )
{
  // tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
  // C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
  // Access : 3*i+j
  // 9-Voigt C-indices    0   3   6   8   1   4   5   7   2
  const int Voigt3x3[NUMDFGR_SOSH8P8] = {0,3,6, 8,1,4, 5,7,2};

  voigt3x3 = &(Voigt3x3[0]);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices2TensorTo6Voigt(
  const int*& voigt3x3
  )
{
  // tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
  // C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
  // Access : 3*i+j
  // 9-Voigt C-indices    0   3   5   3   1   4   5   4   2
  const int Voigt3x3[NUMDFGR_SOSH8P8] = {0,3,5, 3,1,4, 5,4,2};

  voigt3x3 = &(Voigt3x3[0]);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector9Voigt(
  LINALG::Matrix<NUMDFGR_SOSH8P8,1>& fvct,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat,
  const bool transpose
  )
{
  const int* voigtrow9 = NULL;
  const int* voigtcol9 = NULL;
  Indices9VoigtTo2Tensor(voigtrow9,voigtcol9,transpose);
    
  for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
    fvct(I,0) = fmat(voigtrow9[I],voigtcol9[I]);  // F_ij

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector6Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,1>& fvct,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat
  )
{
  const int* voigtrow6 = NULL;
  const int* voigtcol6 = NULL;
  Indices6VoigtTo2Tensor(voigtrow6,voigtcol6);
    
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
    if (I < NUMDIM_SOSH8P8)
      fvct(I) = fmat(voigtrow6[I],voigtcol6[I]);  // F_ij
    else
      fvct(I) = fmat(voigtrow6[I],voigtcol6[I])
        + fmat(voigtcol6[I],voigtrow6[I]);  // F_ij+F_ji

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector9VoigtDiffByItself(
  LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8>& invfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat,
  const bool transpose
  )
{
  const int* voigtrow9 = NULL;
  const int* voigtcol9 = NULL;
  Indices9VoigtTo2Tensor(voigtrow9,voigtcol9,transpose);

  for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
  {
    const int i = voigtrow9[I];
    const int k = voigtcol9[I];
    for (int J=0; J<NUMDFGR_SOSH8P8; ++J)
    {
      const int j = voigtrow9[J];
      const int l = voigtcol9[J];
      invfderf(I,J) = -invfmat(i,k)*invfmat(j,l);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& invfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat
  )
{
  const int voigtrow6[NUMSTR_SOSH8P8] = {0,1,2, 0,1,2};
  const int voigtcol6[NUMSTR_SOSH8P8] = {0,1,2, 1,2,0};

  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigtrow6[I];
    const int j = voigtcol6[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigtrow6[K];
      const int l = voigtcol6[K];
      invfderf(I,K) = -0.5*(invfmat(i,k)*invfmat(l,j) + invfmat(i,l)*invfmat(k,j));
      if (I >= NUMDIM_SOSH8P8)
        invfderf(I,K) += -0.5*(invfmat(j,k)*invfmat(l,i) + invfmat(j,l)*invfmat(k,i));
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& sqfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat
  )
{
  const int* voigtrow6 = NULL;
  const int* voigtcol6 = NULL;
  Indices6VoigtTo2Tensor(voigtrow6,voigtcol6);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> id(true);
  for (int i=0; i<NUMDIM_SOSH8P8; ++i) id(i,i) = 1.0;

  // (F.F)_{,F} with F^T=F
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigtrow6[I];
    const int j = voigtcol6[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigtrow6[K];
      const int l = voigtcol6[K];
      sqfderf(I,K) = id(i,k)*fmat(j,l) + id(j,l)*fmat(i,k);
      if (I >= NUMDIM_SOSH8P8)
        sqfderf(I,K) += id(j,k)*fmat(i,l) + id(i,l)*fmat(j,k);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToMatrix6x9Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMDFGR_SOSH8P8>& bm,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt
  )
{
  const int* voigtrow6 = NULL;
  const int* voigtcol6 = NULL;
  Indices6VoigtTo2Tensor(voigtrow6,voigtcol6);

  const int* voigtrow9 = NULL;
  const int* voigtcol9 = NULL;
  Indices9VoigtTo2Tensor(voigtrow9,voigtcol9);

  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigtrow6[I];
    const int j = voigtcol6[I];
    for (int K=0; K<NUMDFGR_SOSH8P8; ++K)
    {
      const int k = voigtrow9[K];
      const int l = voigtcol9[K];
      if (j == l)
        bm(I,K) = bt(k,i);
      else if ( (I >= NUMDIM_SOSH8P8) and (i == l) )
        bm(I,K) = bt(k,j);
      else
        bm(I,K) = 0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToLeftRightProductMatrix6x6Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& bm,  ///< (out) 6x6 Voigt matrix
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt,  ///< (in) 3x3 matrix of 2-tensor
  const bool transpose, ///< 3x3 input matrix is transposed
  const VoigtType outvoigt6,  ///< 6-Voigt vector layout on rows of 6x6 matrix
  const VoigtType invoigt6  ///< 6-Voigt vector layout on columns of 6x6 matrix
  )
{
  const int* voigtrow6 = NULL;
  const int* voigtcol6 = NULL;
  Indices6VoigtTo2Tensor(voigtrow6,voigtcol6);

  for (int ab=0; ab<NUMSTR_SOSH8P8; ++ab)
  {
    const int a = voigtrow6[ab];
    const int b = voigtcol6[ab];
    for (int AB=0; AB<NUMSTR_SOSH8P8; ++AB)
    {
      const int A = voigtrow6[AB];
      const int B = voigtcol6[AB];
      if (transpose)
      {
        bm(AB,ab) = bt(A,a)*bt(B,b);
        if (ab >= NUMSTR_SOSH8P8) bm(AB,ab) += bt(A,b)*bt(B,a);
      }
      else
      {
        bm(AB,ab) = bt(a,A)*bt(b,B);
        if (ab >= NUMSTR_SOSH8P8) bm(AB,ab) += bt(b,A)*bt(a,B);
      }
      if ( (invoigt6 == voigt6_stress) and (ab >= NUMSTR_SOSH8P8) )
        bm(AB,ab) *= 2.0;
      if ( (outvoigt6 == voigt6_stress) and (AB >= NUMSTR_SOSH8P8) )
        bm(AB,ab) *= 0.5;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ExtractDispAndPres(
  std::vector<double>& mystat,
  LINALG::Matrix<NUMDISP_SOSH8P8,1>& mydisp,
  LINALG::Matrix<NUMPRES_SOSH8P8,1>& mypres
  )
{
  for (int inod=0; inod<NUMNOD_SOSH8P8; ++inod)
  {
    for (int idis=0; idis<NODDISP_SOSH8P8; ++idis)
      mydisp(idis+(inod*NODDISP_SOSH8P8),0) = mystat[idis+(inod*NODDOF_SOSH8P8)];
    mypres(inod,0) = mystat[3+(inod*NODDOF_SOSH8P8)];
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementMatrix(
  LINALG::Matrix<NUMDOF_SOSH8P8,NUMDOF_SOSH8P8>* mat,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* matdd,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matdp,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matpd,
  const LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8>* matpp
)
{
  const int d2dp[NUMDISP_SOSH8P8] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_SOSH8P8] = {      3,      7,       11,         15,         19,         23,         27,        31};
  mat->Clear();
  for (int i=0; i<NUMDISP_SOSH8P8; ++i)
  {
    const int I = d2dp[i];
    if (matdd != NULL)
    {
      for (int j=0; j<NUMDISP_SOSH8P8; ++j)
      {
        const int J = d2dp[j];
        (*mat)(I,J) = (*matdd)(i,j);
      }
    }
    if (matdp != NULL)
    {
      for (int l=0; l<NUMPRES_SOSH8P8; ++l)
      {
        const int L = p2dp[l];
        (*mat)(I,L) = (*matdp)(i,l);
      }
    }
  }
  for (int k=0; k<NUMPRES_SOSH8P8; ++k)
  {
    const int K = p2dp[k];
    if ( (matpd != NULL) or (matdp != NULL) )
    {
      for (int j=0; j<NUMDISP_SOSH8P8; ++j)
      {
        const int J = d2dp[j];
        if (matpd != NULL)
          (*mat)(K,J) = (*matpd)(k,j);
        else
          (*mat)(K,J) = (*matdp)(j,k);
      }
    }
    if (matpp != NULL)
    {
      for (int l=0; l<NUMPRES_SOSH8P8; ++l)
      {
        const int L = p2dp[l];
        (*mat)(K,L) = (*matpp)(k,l);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementVector(
  LINALG::Matrix<NUMDOF_SOSH8P8,1>* vct,
  const LINALG::Matrix<NUMDISP_SOSH8P8,1>* vctd,
  const LINALG::Matrix<NUMPRES_SOSH8P8,1>* vctp
)
{
  const int d2dp[NUMDISP_SOSH8P8] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_SOSH8P8] = {      3,      7,       11,         15,         19,         23,         27,        31};
  vct->Clear();
  if (vctd != NULL)
  {
    for (int i=0; i<NUMDISP_SOSH8P8; ++i)
    {
      const int I = d2dp[i];
      (*vct)(I,0) = (*vctd)(i,0);
    }
  }
  if (vctp != NULL)
  {
    for (int k=0; k<NUMPRES_SOSH8P8; ++k)
    {
      const int K = p2dp[k];
      (*vct)(K,0) = (*vctp)(k,0);
    }
  }
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

  int num_morphed_so_hex8_easmild = 0;
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

      int new_nodeids[DRT::ELEMENTS::So_sh8p8::NUMNOD_SOSH8P8];

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
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_SOSH8P8, new_nodeids);
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
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_SOSH8P8, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::autot:
      case DRT::ELEMENTS::So_sh8::globz: {
        // no resorting necessary
        for (int node = 0; node < 8; ++node) {
          new_nodeids[node] = actele->NodeIds()[node];
        }
        actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_SOSH8P8, new_nodeids);
        actele->nodes_rearranged_ = true;
        break;
      }
      case DRT::ELEMENTS::So_sh8::undefined: {
        // here comes plan B: morph So_sh8p8 to So_hex8
        actele->SetType(DRT::Element::element_so_hex8);
        actele->soh8_reiniteas(DRT::ELEMENTS::So_hex8::soh8_easnone);
        actele->InitJacobianMapping();
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

  if (num_morphed_so_hex8_easmild>0){
    cout << endl << num_morphed_so_hex8_easmild
    << " Sosh8p8-Elements have no clear 'thin' direction and have morphed to So_hex8 with eas_mild" << endl;
  } 
  if (num_morphed_so_hex8_easnone>0){
    cout << endl << num_morphed_so_hex8_easnone
    << " Sosh8p8-Elements have no clear 'thin' direction and have morphed to So_hex8 with eas_none" << endl;
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
    actele->InitJacobianMapping();
  }

  // **************** debug printout ot gmesh **********************************
  //sosh8_gmshplotdis(dis);

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
