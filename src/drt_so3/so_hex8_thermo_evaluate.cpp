/*!----------------------------------------------------------------------
\file so_hex8_thermo_evaluate.cpp
\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/micromaterial.H"
#include <iterator>

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_globalproblem.H"

//#include "Sacado.hpp"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 02/10 |
 | originally by maf 04/07                                              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{
  // call the simple Evaluate(lm) if there is only the structure discretization
  if (la.Size()==1)
  {
    // "normal" Evaluate
    return Evaluate(
      params,
      discretization,
      la[0].lm_, // location vector is build by the first column of la
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }

  // TSI: volume coupling stuff

  // type of kinematic of the problem
  // solve a geometric linear system
  if (kintype_ == DRT::ELEMENTS::So_hex8::soh8_geolin)
  {
    return LinEvaluate(
      params,
      discretization,
      la,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }
  // else: geometrically non-linear with Total Lagrangean approach

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_nlnstiff")              act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")          act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")                act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_update_istep")          act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")           act = So_hex8::calc_struct_reset_istep;  // needed for TangDis predictor
  else if (action=="postprocess_stress")                act = So_hex8::postprocess_stress;
  else if (action=="calc_struct_stifftemp")             act = So_hex8::calc_struct_stifftemp;
  else if (action=="calc_poroelast_nlnstiff")   		    act = So_hex8::calc_poroelast_nlnstiff;
  else if (action=="calc_poroelast_structurecoupling")  act = So_hex8::calc_poroelast_structurecoupling;
  else if (action=="calc_poroelast_internalforce")   	  act = So_hex8::calc_poroelast_internalforce;
  else if (action=="calc_poroelast_nlnstiffmass")   	  act = So_hex8::calc_poroelast_nlnstiffmass;
  else dserror("Unknown type of action for So_hex8: %s",action.c_str());
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // nonlinear stiffness and internal force vector
  case calc_struct_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // elemat2, elevec2+3 is not used anyway

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState(0,"residual displacement");
    if (disp==null || res==null)
      dserror("Cannot get state vectors 'displacement' and/or residual");
    vector<double> mydisp((la[0].lm_).size());
    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);  // global, local, lm
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;
    // call the well-known soh8_nlnstiffmass for the normal structure solution
    soh8_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,NULL,params,
      INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
       = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");
      // call the temperature discretization in the location array
      std::vector<double> mytempnp((la[1].lm_).size());

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
      // number of nodes per element
      const int nen_ = 8;
      if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);

      // calculate the THERMOmechanical solution
      soh8_nlnstifftemp(la,mydisp,myres,mytempnp,&elemat1,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }

  }
  break;

  //==================================================================================
  // nonlinear stiffness, internal force vector, and consistent mass matrix
  case calc_struct_nlnstiffmass:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // massmatrix
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // eleve2+elevec3 is not used anyway

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0, "displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0, "residual displacement");
    if (disp==Teuchos::null || res==Teuchos::null)
      dserror("Cannot get state vectors 'displacement' and/or residual");

    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // lm now contains only u-dofs
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm); // lm now contains only u-dofs
    // call the well-known soh8_nlnstiffmass for the normal structure solution
    soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,NULL,
      params,INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

    // need current temperature state,
    // call the temperature discretization: thermo equates 2nd dofset
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // call the temperature discretization in the location array
      std::vector<double> mytempnp((la[1].lm_).size());
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
       = discretization.GetState(1,"temperature");
      if (tempnp==Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_ = 1;
      // number of nodes per element
      const int nen_ = 8;

      if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
      // build the current temperature vector
      LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[1]),true);  // view only!
      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);
      // calculate the THERMOmechanical solution
      soh8_nlnstifftemp(la,mydisp,myres,mytempnp,&elemat1,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

    }
  }
  break;

  //==================================================================================
  // allowing the predictor TangDis in .dat --> can be decisive in compressible case!
  case calc_struct_reset_istep:
  {}
  break;

  //==================================================================================
  // evaluate stresses and strains at gauss points
  case calc_struct_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==Owner())
    {
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;
      Teuchos::RCP<const Epetra_Vector> disp
        = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res
        = discretization.GetState(0,"residual displacement");
      Teuchos::RCP<vector<char> > stressdata
        = params.get<Teuchos::RCP<vector<char> > >("stress", Teuchos::null);
      Teuchos::RCP<vector<char> > straindata
        = params.get<Teuchos::RCP<vector<char> > >("strain", Teuchos::null);
      Teuchos::RCP<vector<char> > plstraindata
        = params.get<Teuchos::RCP<vector<char> > >("plstrain", Teuchos::null);
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
      if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
      if (plstraindata==Teuchos::null) dserror("Cannot get 'plastic strain' data");
      vector<double> mydisp((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> plstrain;
      INPAR::STR::StressType iostress
        = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
      INPAR::STR::StrainType ioplstrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "ioplstrain", INPAR::STR::strain_none);
      // call the well-known soh8_nlnstiffmass for the normal structure solution
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,&plstrain,
        params,iostress,iostrain,ioplstrain);

      // need current temperature state,
      // call the temperature discretization: thermo equates 2nd dofset
      // disassemble temperature
      if (discretization.HasState(1,"temperature"))
      {
        // call the temperature discretization in the location array
        std::vector<double> mytempnp((la[1].lm_).size());
        // check if you can get the temperature state
        Teuchos::RCP<const Epetra_Vector> tempnp
          = discretization.GetState(1,"temperature");
        if (tempnp==Teuchos::null)
          dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_ = 1;
        // number of nodes per element
        const int nen_ = 8;

        if (la[1].Size() != nen_*numdofpernode_)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
        // get the temperature dependent stress
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stresstemp;
        // calculate the THERMOmechanical solution: temperature stresses
        soh8_nlnstifftemp(la,mydisp,myres,mytempnp,NULL,NULL,&stresstemp,NULL,
          params,iostress,INPAR::STR::strain_none);

        // 08.04.11
        // calculate the THERMOmechanical term for fint: temperature stresses
        soh8_finttemp(la,mydisp,myres,mytempnp,NULL,&stresstemp,NULL,params,
          iostress,INPAR::STR::strain_none);

        // total stress
        // add stresstemp to the mechanical stress
        stress.Update(1.0,stresstemp,1.0);
      }

      {
        DRT::PackBuffer data;
        AddtoPack(data, stress);
        data.StartPacking();
        AddtoPack(data, stress);
        std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
      }

      {
        DRT::PackBuffer data;
        AddtoPack(data, strain);
        data.StartPacking();
        AddtoPack(data, strain);
        std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
      }

      {
        DRT::PackBuffer data;
        AddtoPack(data, plstrain);
        data.StartPacking();
        AddtoPack(data, plstrain);
        std::copy(data().begin(),data().end(),std::back_inserter(*plstraindata));
      }
    }
  }
  break;

  //==================================================================================
  case calc_struct_update_istep:
  {
    // Update of history for visco material if they exist
    RefCountPtr<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Update();
    }
  }
  break;

  //==================================================================================
  // postprocess stresses/strains at gauss points

  // note that in the following, quantities are always referred to as
  // "stresses" etc. although they might also apply to strains
  // (depending on what this routine is called for from the post filter)
  case postprocess_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==Owner())
    {
      const Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap
        = params.get<Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);
      Teuchos::RCP<Epetra_MultiVector> poststress
        = params.get<RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
      if (poststress==Teuchos::null)
        dserror("No element stress/strain vector available");

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        soh8_expol(gpstress, *poststress);

      }
      else if (stresstype=="cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_SOH8; ++i)
          {
            double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
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
  }
  break;

//==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case calc_poroelast_nlnstiff:
    {
      // stiffness
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
      // internal force vector
      LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
      // elemat2,elevec2+3 are not used anyway

      // need current displacement, velocities and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState(0,"residual displacement");
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0,"velocity");
      if (disp==null )
      dserror("Cannot get state vector 'displacement' ");
      if (vel==null )
      dserror("Cannot get state vector 'velocity' ");
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;

      vector<double> mydisp((lm).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // global, local, lm

      vector<double> myres((lm).size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      vector<double> myvel((lm).size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;
      // call the well-known soh8_nlnstiffmass for the normal structure solution
      // soh8_linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
      //        INPAR::STR::stress_none,INPAR::STR::strain_none);
      soh8_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures

      if (discretization.HasState(1,"fluidvel"))
      {
        //  dof per node
        const int numdofpernode = NumDofPerNode(1,*(Nodes()[0]));
        // number of nodes per element
        const int nen = 8;

        LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> myfluidvel(true);
        LINALG::Matrix<NUMNOD_SOH8,1> myepreaf(true);

        if (la[1].Size() != nen*numdofpernode)
        dserror("Location vector length for velocities does not match!");

        // check if you can get the velocity state
        Teuchos::RCP<const Epetra_Vector> velnp
        = discretization.GetState(1,"fluidvel");
        //if there are no velocities or pressures
        if (velnp==Teuchos::null)
        {
          dserror("Cannot get state vector 'fluidvel' ");
        }
        else
        {
          // extract local values of the global vectors
          std::vector<double> mymatrix(la[1].lm_.size());
          DRT::UTILS::ExtractMyValues(*velnp,mymatrix,la[1].lm_);

          for (int inode=0; inode<NUMNOD_SOH8; ++inode) // number of nodes

          {
            for(int idim=0; idim<NUMDIM_SOH8; ++idim) // number of dimensions

            {
              (myfluidvel)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
            } // end for(idim)

            (myepreaf)(inode,0) = mymatrix[NUMDIM_SOH8+(inode*numdofpernode)];
          }
        }

        soh8_nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,matptr,&elemat2,&elevec1,NULL,NULL,NULL,params,
            INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness, mass matrix and internal force vector for poroelasticity
    case calc_poroelast_nlnstiffmass:
    {
      // stiffness
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
      // mass
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
      // internal force vector
      LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
      // elemat2,elevec2+3 are not used anyway

      // need current displacement, velocities and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState(0,"residual displacement");
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0,"velocity");
      if (disp==null )
      dserror("Cannot get state vector 'displacement' ");
      if (vel==null )
      dserror("Cannot get state vector 'velocity' ");
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;

      vector<double> mydisp((lm).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // global, local, lm

      vector<double> myres((lm).size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      vector<double> myvel((lm).size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;
      // call the well-known soh8_nlnstiffmass for the normal structure solution
      // soh8_linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
      //        INPAR::STR::stress_none,INPAR::STR::strain_none);
      soh8_nlnstiffmass(lm,mydisp,myres,matptr,&elemat2,&elevec1,NULL,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

      //scale mass matrix
      const double initporosity = params.get<double>("initporosity");
      elemat2.Scale(1-initporosity);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures

      if (discretization.HasState(1,"fluidvel"))
      {
        //  dof per node
        const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
        // number of nodes per element
        const int nen_ = 8;

        LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> myfluidvel(true);
        LINALG::Matrix<NUMNOD_SOH8,1> myepreaf(true);

        if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for velocities does not match!");

        // check if you can get the velocity state
        Teuchos::RCP<const Epetra_Vector> velnp
        = discretization.GetState(1,"fluidvel");
        //if there are no velocities or pressures, set them to zero
        if (velnp==Teuchos::null)
        {
          dserror("Cannot get state vector 'fluidvel' ");
        }
        else
        {
          // extract local values of the global vectors
          std::vector<double> mymatrix(la[1].lm_.size());
          DRT::UTILS::ExtractMyValues(*velnp,mymatrix,la[1].lm_);

          for (int inode=0; inode<NUMNOD_SOH8; ++inode) // number of nodes

          {
            for(int idim=0; idim<NUMDIM_SOH8; ++idim) // number of dimensions

            {
              (myfluidvel)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
            } // end for(idim)

            (myepreaf)(inode,0) = mymatrix[NUMDIM_SOH8+(inode*numdofpernode_)];
          }
        }

        soh8_nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,matptr,NULL,&elevec1,NULL,NULL,NULL,params,
            INPAR::STR::stress_none,INPAR::STR::strain_none);
      }

    }
    break;

    //==================================================================================
    // coupling terms in force-vector and stiffness matrix for poroelasticity
    case calc_poroelast_structurecoupling:
    {
      // stiffness
      LINALG::Matrix<NUMDOF_SOH8,(NUMDIM_SOH8+1)*NUMNOD_SOH8> elemat1(elemat1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,(NUMDIM_SOH8+1)*NUMNOD_SOH8> elemat2(elemat2_epetra.A(),true);

      // internal force vector
      //LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
      //LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);

      // elemat2,elevec2+3 are not used anyway

      // need current displacement, velocities and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0,"velocity");
      if (disp==null )
      dserror("Cannot get state vector 'displacement' ");
      if (vel==null )
      dserror("Cannot get state vector 'velocity' ");
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;

      vector<double> mydisp((lm).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // global, local, lm

      vector<double> myvel((lm).size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      LINALG::Matrix<NUMDOF_SOH8,(NUMDIM_SOH8+1)*NUMNOD_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1,"fluidvel"))
      {
        //  dof per node
        const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
        // number of nodes per element
        const int nen_ = 8;

        LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> myvelnp(true);
        LINALG::Matrix<NUMNOD_SOH8,1> myepreaf(true);

        // check if you can get the velocity state
        Teuchos::RCP<const Epetra_Vector> velnp
        = discretization.GetState(1,"fluidvel");
        if (velnp==Teuchos::null)
        dserror("Cannot get state vector 'fluidvel'");

        if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for fluid velocities and pressures does not match!");

        // extract the current velocitites and pressures of the global vectors
        std::vector<double> mymatrix(la[1].lm_.size());
        DRT::UTILS::ExtractMyValues(*velnp,mymatrix,la[1].lm_);

        for (int inode=0; inode<NUMNOD_SOH8; ++inode) // number of nodes

        {
          for(int idim=0; idim<NUMDIM_SOH8; ++idim) // number of dimensions

          {
            (myvelnp)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
          } // end for(idim)

          (myepreaf)(inode,0) = mymatrix[NUMDIM_SOH8+(inode*numdofpernode_)];
        }

        soh8_coupling_poroelast(lm,mydisp,myvel,myvelnp,myepreaf,matptr,NULL,NULL,NULL,params);
      }

    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case calc_poroelast_internalforce:
    {
      // stiffness
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
      // internal force vector
      LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
      LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
      // elemat2,elevec2+3 are not used anyway

      // need current displacement, velocities and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState(0,"residual displacement");
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0,"velocity");
      if (disp==null )
      dserror("Cannot get state vector 'displacement' ");
      if (vel==null )
      dserror("Cannot get state vector 'velocity' ");
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;

      vector<double> mydisp((lm).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // global, local, lm

      vector<double> myres((lm).size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      vector<double> myvel((lm).size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      //  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      //  if (elemat1.IsInitialized()) matptr = &elemat1;
      // call the well-known soh8_nlnstiffmass for the normal structure solution
      // soh8_linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
      //        INPAR::STR::stress_none,INPAR::STR::strain_none);
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,&elevec1,NULL,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures

      if (discretization.HasState(1,"fluidvel"))
      {
        //  dof per node
        const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
        // number of nodes per element
        const int nen_ = 8;

        LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> myfluidvel(true);
        LINALG::Matrix<NUMNOD_SOH8,1> myepreaf(true);

        if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for velocities does not match!");

        // check if you can get the velocity state
        Teuchos::RCP<const Epetra_Vector> velnp
        = discretization.GetState(1,"fluidvel");
        //if there are no velocities or pressures, set them to zero
        if (velnp==Teuchos::null)
        {
          dserror("Cannot get state vector 'fluidvel' ");
        }
        else
        {
          // extract local values of the global vectors
          std::vector<double> mymatrix(la[1].lm_.size());
          DRT::UTILS::ExtractMyValues(*velnp,mymatrix,la[1].lm_);

          for (int inode=0; inode<NUMNOD_SOH8; ++inode) // number of nodes

          {
            for(int idim=0; idim<NUMDIM_SOH8; ++idim) // number of dimensions

            {
              (myfluidvel)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
            } // end for(idim)

            (myepreaf)(inode,0) = mymatrix[NUMDIM_SOH8+(inode*numdofpernode_)];
          }
        }

        soh8_nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,NULL,NULL,&elevec1,NULL,NULL,NULL,params,
            INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
    }
    break;
    //==================================================================================
    // linear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_stifftemp:
    {
      // mechanical-thermal system matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMNOD_SOH8> stiffmatrixcoupl(elemat1_epetra.A(),true);
      // elemat2,elevec1-3 are not used anyway

      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
      if (disp==null)
      dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp((la[0].lm_).size());
      // build the location vector only for the structure field
      DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);

      // calculate the mechanical-thermal system matrix
      soh8_stifftemp(la,mydisp,&stiffmatrixcoupl);
    }
    break;

    //==================================================================================
    default:
    dserror("Unknown type of action for So_hex8");
    } // action

    return 0;
  } // Evaluate


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 03/10 |
 | originally by maf 04/07  (private)                                   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstifftemp(
  DRT::Element::LocationArray& la,  // location array
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  vector<double>& temp, // current temperature
  LINALG::Matrix<NUMDOF_SOH8,1>* tempstiffmatrix, // coupling stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  // strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry (8x3)
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // X, material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // x, current  coord. of element
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

  // vector of the current element temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> etemp(true);
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    etemp(i,0) = temp[i+0];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // TODO 08.04.11 defgrd(false) or true??? thermo_lin: true, hex8_eval false!
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(false);

  // identical shapefunctions for the displacements and the temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> shapetemp;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ); //  (6.17)

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
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
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

    // copy structural shape functions needed for the thermo field
    shapetemp.Update(shapefcts[gp]);

    // product of shapefunctions and element temperatures for stresstemp
    LINALG::Matrix<1,1> Ntemp(true);
    Ntemp.MultiplyTN(shapetemp,etemp);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    // calculate the stress part dependent on the temperature in the material
    LINALG::Matrix<NUMSTR_SOH8,1> ctemp(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stresstemp(true);
    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) stresstemp = C . Delta T
    soh8_mat_temp(&stresstemp,&ctemp,&density,&Ntemp,&defgrd);

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH8; ++i)
        (*elestress)(gp,i) = stresstemp(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstresstemp;
      pkstresstemp(0,0) = stresstemp(0);
      pkstresstemp(0,1) = stresstemp(3);
      pkstresstemp(0,2) = stresstemp(5);
      pkstresstemp(1,0) = pkstresstemp(0,1);
      pkstresstemp(1,1) = stresstemp(1);
      pkstresstemp(1,2) = stresstemp(4);
      pkstresstemp(2,0) = pkstresstemp(0,2);
      pkstresstemp(2,1) = pkstresstemp(1,2);
      pkstresstemp(2,2) = stresstemp(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystresstemp;
      temp.Multiply(1.0/detF,defgrd,pkstresstemp,0.0);
      cauchystresstemp.MultiplyNT(temp,defgrd);

      (*elestress)(gp,0) = cauchystresstemp(0,0);
      (*elestress)(gp,1) = cauchystresstemp(1,1);
      (*elestress)(gp,2) = cauchystresstemp(2,2);
      (*elestress)(gp,3) = cauchystresstemp(0,1);
      (*elestress)(gp,4) = cauchystresstemp(1,2);
      (*elestress)(gp,5) = cauchystresstemp(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
    }

    double detJ_w = detJ*gpweights[gp];
    if (force != NULL && tempstiffmatrix != NULL)
    {
      // integrate internal force vector
      // f = f + (B^T . sigma_temp) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stresstemp, 1.0);

      //-----------------------------------------------------------------------
      // integrate `initial-displacement-temperature' stiffness matrix  02/10
      // update the stiffness part depending on the temperature
      // k_theta += 1.0 * k_theta + detJ * w(gp) * (B^T . C_Theta . N_Theta)
      //-----------------------------------------------------------------------
      LINALG::Matrix<NUMSTR_SOH8,1> cn; // 6x1 = 6x1n
      // extract the shapefunctions on the gp_i and multiply it with cn
      for (int stresscomp=0; stresscomp<NUMSTR_SOH8; ++stresscomp)
      {
        cn.Multiply(ctemp,Ntemp); // (6x1) = (6x1)(1x1)
      }
      tempstiffmatrix->MultiplyTN(detJ_w,bop,cn,1.0); // [(3nx6)(6x1)=(3nx1)] (24x6)(6x1)=(24x1)


      //-----------------------------------------------------------------------
      // integrate `geometric' stiffness matrix and add to keu            04/10
      // kgeo += (B_L^T . sigma_Theta . B_L^T) * detJ * w(gp)  with B_L = Ni,Xj
      //-----------------------------------------------------------------------
      // \f {\mathbf S}:\Delta\delta{\mathbf E} \f
      // \f = \delta {\mathbf d}^T \,{\mathbf k}_{geo}\, \Delta{\mathbf d} \f
      // auxiliary integrated temperature stress
      LINALG::Matrix<6,1> sfac(stresstemp);

      // detJ * w(gp) * [S11,S22,S33,S12=S21,S23=S32,S13=S31]
      sfac.Scale(detJ_w);

      // intermediate Sm.B_L
      vector<double> SmB_L(3);

      // like in so_hex8_evaluate no changes for tempstress
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)
      // with B_L = Ni,Xj: LINEAR B-operator see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
                     + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
                     + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
                     + sfac(2) * N_XYZ(2, inod);
        SmB_L[0] = sfac(3) * N_XYZ(1, inod) + sfac(0) * N_XYZ(0, inod)
                     + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
                     + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
                     + sfac(2) * N_XYZ(2, inod);

        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          // <3n,1>
          // 3n: displacement:3(3FHG/node)*8(8nodes/element)+3(dim),
          //  1: temperature is a scalar
          (*tempstiffmatrix)(3*inod+0,0) += bopstrbop;
          (*tempstiffmatrix)(3*inod+1,0) += bopstrbop;
          (*tempstiffmatrix)(3*inod+2,0) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::soh8_nlnstifftemp


/*----------------------------------------------------------------------*
 | material law with temperature part for So_hex8            dano 05/10 |
 | originally by gee in so_material.cpp 10/08                           |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_temp(
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stresstemp,
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* ctemp,
  double* density,
  LINALG::Matrix<1,1>* Ntemp,  // temperature of element
  LINALG::Matrix<3,3>* defgrd //,
//  const int gp,
//  Teuchos::ParameterList& params
  )
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stresstemp) dserror("No stress vector supplied");
  if (!ctemp) dserror("No material tangent matrix supplied");
  if (!Ntemp) dserror("No temperature supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    // st.venant-kirchhoff-material with temperature
    case INPAR::MAT::m_thermostvenant:
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast <MAT::ThermoStVenantKirchhoff*>(mat.get());
      thrstvk->Evaluate(*Ntemp,*ctemp,*stresstemp);
      *density = thrstvk->Density();
      return;
      break;
    }
    // small strain von Mises thermoelastoplastic material
    case INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      thrpllinelast->Evaluate(*Ntemp,*ctemp,*stresstemp);
      *density = thrpllinelast->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of temperature dependent material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh8_mat_temp



/*----------------------------------------------------------------------*
 | get the constant temperature fraction for stresstemp      dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Stempconst(
  LINALG::Matrix<6,1>* ctemp,
  LINALG::Matrix<6,1>* stempconst
  )
{
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    /*-------- thermo st.venant-kirchhoff-material */
    case INPAR::MAT::m_thermostvenant:
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast<MAT::ThermoStVenantKirchhoff*>(mat.get());
       return thrstvk->Stempconst(*ctemp,*stempconst);
       break;
    }
    // small strain von Mises thermoelastoplastic material
    case INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      return thrpllinelast->Stempconst(*ctemp,*stempconst);
      break;
    }
    default:
      dserror("Cannot ask material for the temperature rhs");
      break;
  } // switch (mat->MaterialType())

  return;

} // So_hex8::Stempconst


/*----------------------------------------------------------------------*
 | get the constant temperature fraction for stresstemp      dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Ctemp(LINALG::Matrix<6,1>* ctemp)
{
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    /*-------- thermo st.venant-kirchhoff-material */
    case INPAR::MAT::m_thermostvenant:
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast<MAT::ThermoStVenantKirchhoff*>(mat.get());
       return thrstvk->SetupCthermo(*ctemp);
       break;
    }
    // small strain von Mises thermoelastoplastic material
    case INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      return thrpllinelast->SetupCthermo(*ctemp);
      break;
    }
    default:
      dserror("Cannot ask material for the temperature rhs");
      break;

  } // switch (mat->MaterialType())

}


/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (private)
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstiff_poroelast(
    vector<int>& lm, // location matrix
    vector<double>& disp, // current displacements
    vector<double>& vel, // current velocities
    // vector<double>&           residual,       // current residual displ
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> & evelnp, LINALG::Matrix<
        NUMNOD_SOH8, 1> & epreaf,
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix, // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* reamatrix, // element reactive matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force, // element internal force vector
    LINALG::Matrix<NUMDOF_SOH8, 1>* forcerea, // element reactive force vector
    LINALG::Matrix<NUMGPT_SOH8, NUMSTR_SOH8>* elestress, // stresses at GP
    LINALG::Matrix<NUMGPT_SOH8, NUMSTR_SOH8>* elestrain, // strains at GP
    ParameterList& params, // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress, // stress output option
    const INPAR::STR::StrainType iostrain) // strain output option
{
  /* ============================================================================*
   ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
   ** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8, 1> > shapefcts =
      soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> > derivs =
      soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // get global id of the structure element
  int id = Id();
  //access fluid discretization
  RCP<DRT::Discretization> fluiddis = null;
  fluiddis = DRT::Problem::Instance()->Dis(genprob.numff, 0);
  //get corresponding fluid element (it has the same global ID as the structure element)
  DRT::Element* fluidele = fluiddis->gElement(id);
  if (fluidele == NULL)
    dserror("Fluid element %i not on local processor", id);

  //get fluid material
  const MAT::FluidPoro* fluidmat = static_cast<const MAT::FluidPoro*>((fluidele->Material()).get());
  if(fluidmat->MaterialType() != INPAR::MAT::m_fluidporo)
    dserror("invalid fluid material for poroelasticity");

  //get structure material
  MAT::StructPoro* structmat = static_cast<MAT::StructPoro*>((Material()).get());
  if(structmat->MaterialType() != INPAR::MAT::m_structporo)
    dserror("invalid structure material for poroelasticity");

  double reacoeff = fluidmat->ComputeReactionCoeff();
  const double bulkmodulus = structmat->Bulkmodulus();
  const double penalty = structmat->Penaltyparameter();
  const double initporosity = params.get<double>("initporosity");
  double dt = params.get<double>("delta time");

  // update element geometry
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xrefe; // material coord. of element
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xcurr; // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;

  // (r,s,t) gp-locations of fully integrated linear 8-node Hex
  const double gploc = 1.0/sqrt(3.0); // gp sampling point value for linear fct
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xgp; // coord. of guass point in reference coordinates
  xgp(0,0)=-gploc;
  xgp(1,0)=-gploc;
  xgp(2,0)=-gploc;

  xgp(0,1)= gploc;
  xgp(1,1)=-gploc;
  xgp(2,1)=-gploc;

  xgp(0,2)= gploc;
  xgp(1,2)= gploc;
  xgp(2,2)=-gploc;

  xgp(0,3)=-gploc;
  xgp(1,3)= gploc;
  xgp(2,3)=-gploc;

  xgp(0,4)=-gploc;
  xgp(1,4)=-gploc;
  xgp(2,4)= gploc;

  xgp(0,5)= gploc;
  xgp(1,5)=-gploc;
  xgp(2,5)= gploc;

  xgp(0,6)= gploc;
  xgp(1,6)= gploc;
  xgp(2,6)= gploc;

  xgp(0,7)=-gploc;
  xgp(1,7)= gploc;
  xgp(2,7)= gploc;

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(0,i) = x[0];
    xrefe(1,i) = x[1];
    xrefe(2,i) = x[2];

    xcurr(0,i) = xrefe(0,i) + disp[i*NODDOF_SOH8+0];
    xcurr(1,i) = xrefe(1,i) + disp[i*NODDOF_SOH8+1];
    xcurr(2,i) = xrefe(2,i) + disp[i*NODDOF_SOH8+2];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodaldisp;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodaldisp(i,0) = disp[i];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodalvel;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodalvel(i,0) = vel[i];
  }

  //vector of porosity at gp (for output only)
  std::vector<double> porosity_gp(NUMGPT_SOH8,0.0);

  //******************* FAD ************************
  /*

   // sacado data type replaces "double"
   typedef Sacado::Fad::DFad<double> FAD;  // for first derivs
   // sacado data type replaces "double" (for first+second derivs)
   typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD;

   vector<FAD> fad_disp(NUMDOF_SOH8);
   for (int i=0; i<NUMDOF_SOH8; ++i)
   {
   fad_disp[i] = disp[i];
   fad_disp[i].diff(i,NUMDOF_SOH8);   // variables to differentiate for
   }

   LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMNOD_SOH8>  fad_xrefe(false);
   LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMNOD_SOH8> fad_xcurr(false);

   for (int i=0; i<NUMNOD_SOH8; ++i)
   {
   const double* x = nodes[i]->X();

   fad_xrefe(0,i) = x[0];
   fad_xrefe(1,i) = x[1];
   fad_xrefe(2,i) = x[2];

   fad_xcurr(0,i) = fad_xrefe(0,i) + fad_disp[i*NODDOF_SOH8+0];
   fad_xcurr(1,i) = fad_xrefe(1,i) + fad_disp[i*NODDOF_SOH8+1];
   fad_xcurr(2,i) = fad_xrefe(2,i) + fad_disp[i*NODDOF_SOH8+2];
   }

   LINALG::TMatrix<FAD,NUMDOF_SOH8,1> fad_nodaldisp(false);
   LINALG::TMatrix<FAD,NUMNOD_SOH8,1> fad_epreaf(false);
   for(int i=0; i<NUMNOD_SOH8 ; i++)
   {
   fad_epreaf(i) = epreaf(i);
   // fad_epreaf(i).diff(NUMDOF_SOH8+i,NUMDOF_SOH8+NUMNOD_SOH8);
   }

   for (int i=0; i<NUMDOF_SOH8; ++i)
   fad_nodaldisp(i,0) = fad_disp[i];
   */
  //******************** FAD ***********************

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  LINALG::Matrix<6,NUMNOD_SOH8> N_XYZ2;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);
  LINALG::Matrix<6,NUMNOD_SOH8> deriv2;
  LINALG::Matrix<NUMDIM_SOH8,1> xsi;
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    //get shapefunctions and its derivatives at gausspoint
    LINALG::Matrix<NUMNOD_SOH8,1> shapefct = shapefcts[gp];
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> deriv = derivs[gp];
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invJ = invJ_[gp];

    //get coordinates of gausspoint and compute second derivatives at gausspoint
    xsi(0)=xgp(0,gp);
    xsi(1)=xgp(1,gp);
    xsi(2)=xgp(2,gp);
    DRT::UTILS::shape_function_deriv2<DRT::Element::hex8>(xsi,deriv2);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ X_,r  Y_,r  Z_,r ]^-1
     **     J^-1 = [ X_,s  Y_,s  Z_,s ]
     **            [ X_,t  Y_,t  Z_,t ]
     */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // transposed jacobian "dX/ds"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xjm0;
    xjm0.MultiplyNT(deriv,xrefe);

    // get the second derivatives of standard element at current GP w.r.t. XYZ
    DRT::UTILS::gder2<DRT::Element::hex8>(xjm0,N_XYZ,deriv2,xrefe,N_XYZ2);

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xjm;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xji;
    xjm.MultiplyNT(deriv,xcurr);
    const double det = xji.Invert(xjm);

    // determinant of deformationgradient: det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    const double J = det/detJ;

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // pressure gradient at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    // fluid velocity at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    // structure displacement and velocity at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> dispint(true);
    LINALG::Matrix<NUMDIM_SOH8,1> velint(true);

    for(int i=0; i<NUMNOD_SOH8; i++)
    for(int j=0; j<NUMDIM_SOH8; j++)
    {
      dispint(j) += nodaldisp(i*NUMDIM_SOH8+j) * shapefct(i);
      velint(j) += nodalvel(i*NUMDIM_SOH8+j) * shapefct(i);
    }

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr,N_XYZ); //  (6.17)
    //defgrd.PutScalar(0.0);
    //for (int i=0; i<3; i++)
    //  defgrd(i,i) =1;

    // set to initial state as test to receive a linear solution
    //for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

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
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
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

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    /*
     // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
     // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
     Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
     LINALG::Matrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
     glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
     glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
     glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
     glstrain(3) = cauchygreen(0,1);
     glstrain(4) = cauchygreen(1,2);
     glstrain(5) = cauchygreen(2,0);
     */

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> C_inv(false);
    C_inv.Invert(cauchygreen);

    //inverse Right Cauchy-Green tensor as vector
    LINALG::Matrix<6,1> C_inv_vec(true);
    C_inv_vec(0) = C_inv(0,0);
    C_inv_vec(1) = C_inv(1,1);
    C_inv_vec(2) = C_inv(2,2);
    C_inv_vec(3) = C_inv(0,1);
    C_inv_vec(4) = C_inv(1,2);
    C_inv_vec(5) = C_inv(2,0);

    // inverse deformation gradient F^-1
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //------------------------------------ build F^-1 as vector 9x1
    LINALG::Matrix<9,1> defgrd_inv_vec;
    defgrd_inv_vec(0)=defgrd_inv(0,0);
    defgrd_inv_vec(1)=defgrd_inv(0,1);
    defgrd_inv_vec(2)=defgrd_inv(0,2);
    defgrd_inv_vec(3)=defgrd_inv(1,0);
    defgrd_inv_vec(4)=defgrd_inv(1,1);
    defgrd_inv_vec(5)=defgrd_inv(1,2);
    defgrd_inv_vec(6)=defgrd_inv(2,0);
    defgrd_inv_vec(7)=defgrd_inv(2,1);
    defgrd_inv_vec(8)=defgrd_inv(2,2);

    //------------------------------------ build F^-T as vector 9x1
    LINALG::Matrix<9,1> defgrd_IT_vec;
    defgrd_IT_vec(0)=defgrd_inv(0,0);
    defgrd_IT_vec(1)=defgrd_inv(1,0);
    defgrd_IT_vec(2)=defgrd_inv(2,0);
    defgrd_IT_vec(3)=defgrd_inv(0,1);
    defgrd_IT_vec(4)=defgrd_inv(1,1);
    defgrd_IT_vec(5)=defgrd_inv(2,1);
    defgrd_IT_vec(6)=defgrd_inv(0,2);
    defgrd_IT_vec(7)=defgrd_inv(1,2);
    defgrd_IT_vec(8)=defgrd_inv(2,2);

    //--------------------------- build N_X operator (wrt material config)
    LINALG::Matrix<9,NUMDOF_SOH8> N_X(true); // set to zero
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      N_X(0,3*i+0) = N_XYZ(0,i);
      N_X(1,3*i+1) = N_XYZ(0,i);
      N_X(2,3*i+2) = N_XYZ(0,i);

      N_X(3,3*i+0) = N_XYZ(1,i);
      N_X(4,3*i+1) = N_XYZ(1,i);
      N_X(5,3*i+2) = N_XYZ(1,i);

      N_X(6,3*i+0) = N_XYZ(2,i);
      N_X(7,3*i+1) = N_XYZ(2,i);
      N_X(8,3*i+2) = N_XYZ(2,i);
    }

    LINALG::Matrix<NUMDIM_SOH8*NUMDIM_SOH8,NUMDIM_SOH8> F_X(true);
    for(int i=0; i<NUMDIM_SOH8; i++)
    {
      for(int n=0; n<NUMNOD_SOH8; n++)
      {
        // second derivatives w.r.t. XYZ are orderd as followed: deriv2(N,XX ; N,YY ; N,ZZ ; N,XY ; N,XZ ; N,YZ)
        F_X(i*NUMDIM_SOH8+0, 0) += N_XYZ2(0,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 0) += N_XYZ2(3,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 0) += N_XYZ2(4,n)*nodaldisp(n*NUMDIM_SOH8+i);

        F_X(i*NUMDIM_SOH8+0, 1) += N_XYZ2(3,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 1) += N_XYZ2(1,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 1) += N_XYZ2(5,n)*nodaldisp(n*NUMDIM_SOH8+i);

        F_X(i*NUMDIM_SOH8+0, 2) += N_XYZ2(4,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 2) += N_XYZ2(5,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 2) += N_XYZ2(2,n)*nodaldisp(n*NUMDIM_SOH8+i);

      }
    }

    //--------------------------- compute dJ/dX = dJ/dF : dF/dX = JF^-T : dF/dX at gausspoint
    //--------------material gradient of jacobi determinant J: GradJ = dJ/dX= dJ/dF : dF/dX = J * F^-T : dF/dX
    LINALG::Matrix<1,NUMDIM_SOH8> GradJ;
    GradJ.MultiplyTN(J, defgrd_IT_vec, F_X);

    //------linearization of jacobi determinant detF=J w.r.t. strucuture displacement   dJ/d(us) = dJ/dF : dF/dus = J * F^-T * N,X
    LINALG::Matrix<1,NUMDOF_SOH8> dJ_dus;
    dJ_dus.MultiplyTN(J,defgrd_inv_vec,N_X);

    //------linearization of material gradient of jacobi determinant GradJ  w.r.t. strucuture displacement d(GradJ)/d(us)
    //---------------------d(GradJ)/dus =  dJ/dus * F^-T . : dF/dX + J * dF^-T/dus : dF/dX + J * F^-T : N_X_X

    //dF^-T/dus : dF/dX = - (F^-1 . dN/dX . u_s . F^-1)^T : dF/dX
    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> dFinvdus_dFdX(true);
    for (int i=0; i<NUMDIM_SOH8; i++)
    for (int n =0; n<NUMNOD_SOH8; n++)
    for(int j=0; j<NUMDIM_SOH8; j++)
    {
      const int gid = NUMDIM_SOH8 * n +j;
      for (int k=0; k<NUMDIM_SOH8; k++)
      for(int l=0; l<NUMDIM_SOH8; l++)
      for(int p=0; p<NUMDIM_SOH8; p++)
      {
        dFinvdus_dFdX(p, gid) += -defgrd_inv(l,j) * N_XYZ(k,n) * defgrd_inv(k,i) * F_X(i*NUMDIM_SOH8+l,p);

      }
    }

    //dF^-T/dus
    LINALG::Matrix<NUMDIM_SOH8*NUMDIM_SOH8,NUMDOF_SOH8> dFinvdus(true);
    for (int i=0; i<NUMDIM_SOH8; i++)
    for (int n =0; n<NUMNOD_SOH8; n++)
    for(int j=0; j<NUMDIM_SOH8; j++)
    {
      const int gid = NUMDIM_SOH8 * n +j;
      for (int k=0; k<NUMDIM_SOH8; k++)
      for(int l=0; l<NUMDIM_SOH8; l++)
      dFinvdus(i*NUMDIM_SOH8+l, gid) += -defgrd_inv(l,j) * N_XYZ(k,n) * defgrd_inv(k,i);
    }

    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> dFinvdus_dFdX2;
    dFinvdus_dFdX2.MultiplyTN(F_X,dFinvdus);
    for (int i=0; i<NUMDIM_SOH8; i++)
    for (int j =0; j<NUMDOF_SOH8; j++)
    if(dFinvdus_dFdX(i,j)-dFinvdus_dFdX2(i,j)>1e-8)
    dserror("dFinvdus_dFdX falsch");

    //F^-T : N_X_X
    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> Finv_N_XYZ2(true);

    for (int n =0; n<NUMNOD_SOH8; n++)
    {
      //! second derivatives  are ordered as followed: (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      int n_dim = n*NUMDIM_SOH8;
      Finv_N_XYZ2(0, n_dim + 0) += defgrd_inv(0,0) * N_XYZ2(0,n) + defgrd_inv(1,0) * N_XYZ2(3,n)+ defgrd_inv(2,0) * N_XYZ2(4,n);
      Finv_N_XYZ2(0, n_dim + 1) += defgrd_inv(0,1) * N_XYZ2(0,n) + defgrd_inv(1,1) * N_XYZ2(3,n)+ defgrd_inv(2,1) * N_XYZ2(4,n);
      Finv_N_XYZ2(0, n_dim + 2) += defgrd_inv(0,2) * N_XYZ2(0,n) + defgrd_inv(1,2) * N_XYZ2(3,n)+ defgrd_inv(2,2) * N_XYZ2(4,n);

      Finv_N_XYZ2(1, n_dim + 0) += defgrd_inv(0,0) * N_XYZ2(3,n) + defgrd_inv(1,0) * N_XYZ2(1,n)+ defgrd_inv(2,0) * N_XYZ2(5,n);
      Finv_N_XYZ2(1, n_dim + 1) += defgrd_inv(0,1) * N_XYZ2(3,n) + defgrd_inv(1,1) * N_XYZ2(1,n)+ defgrd_inv(2,1) * N_XYZ2(5,n);
      Finv_N_XYZ2(1, n_dim + 2) += defgrd_inv(0,2) * N_XYZ2(3,n) + defgrd_inv(1,2) * N_XYZ2(1,n)+ defgrd_inv(2,2) * N_XYZ2(5,n);

      Finv_N_XYZ2(2, n_dim + 0) += defgrd_inv(0,0) * N_XYZ2(4,n) + defgrd_inv(1,0) * N_XYZ2(5,n)+ defgrd_inv(2,0) * N_XYZ2(2,n);
      Finv_N_XYZ2(2, n_dim + 1) += defgrd_inv(0,1) * N_XYZ2(4,n) + defgrd_inv(1,1) * N_XYZ2(5,n)+ defgrd_inv(2,1) * N_XYZ2(2,n);
      Finv_N_XYZ2(2, n_dim + 2) += defgrd_inv(0,2) * N_XYZ2(4,n) + defgrd_inv(1,2) * N_XYZ2(5,n)+ defgrd_inv(2,2) * N_XYZ2(2,n);
    }

    LINALG::Matrix<1,NUMDIM_SOH8> temp2;
    temp2.MultiplyTN( defgrd_IT_vec, F_X);

    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> dgradJ_dus(true);
    dgradJ_dus.MultiplyTN(temp2,dJ_dus);

    dgradJ_dus.Update(J,dFinvdus_dFdX,1.0);

    dgradJ_dus.Update(J,Finv_N_XYZ2,1.0);

    //--------------------------------------------------------------------

    const double a = ( bulkmodulus/(1-initporosity) + press - penalty/initporosity ) * J;
    const double b = -a + bulkmodulus + penalty;
    const double c = (b/a) * (b/a) + 4*penalty/a;
    double d = sqrt(c)*a;
    double sign =1.0;

    double test = 1/(2*a)*(-b+d);
    if( test >= 1.0 or test < 0.0 )
    {
      sign = -1.0;
      d = sign*d;
    }

    const double porosity = 1/(2*a)*(-b+d);

    if( porosity >= 1.0 or porosity < 0.0 )
    {
      dserror("invalid porosity!");
    }

    const double d_p = J * (-b+2*penalty)/d;
    // const double d_p_p = ( d * J + d_p * (b - 2*penalty_) ) / (d * d) * J;
    const double d_J = a/J * ( -b + 2*penalty ) / d;
    const double d_J_p = (d_p / J + ( 1-d_p*d_p/(J*J) ) / d *a);
    const double d_J_J = ( a*a/(J*J)-d_J*d_J )/ d;

    //double porosity = structmat->ComputePorosity(press, J, initporosity,gp);

    porosity_gp[gp] = porosity;

    //d(porosity) / d(p)
    const double dphi_dp= - J * porosity/a + (J+d_p)/(2*a);

    //d(porosity) / d(J)
    const double dphi_dJ= -porosity/J+ 1/(2*J) + d_J / (2*a);

    //d(porosity) / d(J)d(pressure)
    const double dphi_dJdp= -1/J*dphi_dp+ d_J_p/(2*a) - d_J*J/(2*a*a);

    //d^2(porosity) / d(J)^2
    const double dphi_dJJ= porosity/(J*J) - dphi_dJ/J - 1/(2*J*J) - d_J/(2*a*J) + d_J_J/(2*a);

    //linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1,NUMDOF_SOH8> dphi_dus;
    dphi_dus.Update( dphi_dJ , dJ_dus );

    //material porosity gradient Grad(phi) = dphi/dp * Grad(p) + dphi/dJ * Grad(J)
    LINALG::Matrix<1,NUMDIM_SOH8> grad_porosity;
    for (int idim=0; idim<NUMDIM_SOH8; ++idim)
    {
      grad_porosity(idim)=dphi_dp*Gradp(idim)+dphi_dJ*GradJ(idim);
    }

    //linearization of material porosity gradient w.r.t structure displacement
    // d ( Grad(\phi) ) / du_s = d\phi / (dJ du_s) * dJ /dX+ d\phi / dJ * dJ / (dX*du_s) + d\phi / (dp*du_s) * dp /dX
    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> dgradphi_dus;
    dgradphi_dus.MultiplyTN(dphi_dJJ, GradJ ,dJ_dus);
    dgradphi_dus.Update(dphi_dJ, dgradJ_dus, 1.0);
    dgradphi_dus.Multiply(dphi_dJdp, Gradp, dJ_dus, 1.0);

    //*****************************************************************************************
    /*
     // return gp strains (only in case of stress/strain output)
     switch (iostrain)
     {
     case INPAR::STR::strain_gl:
     {
     if (elestrain == NULL) dserror("strain data not available");
     for (int i = 0; i < 3; ++i)
     (*elestrain)(gp,i) = glstrain(i);
     for (int i = 3; i < 6; ++i)
     (*elestrain)(gp,i) = 0.5 * glstrain(i);
     }
     break;
     case INPAR::STR::strain_ea:
     {
     if (elestrain == NULL) dserror("strain data not available");
     // rewriting Green-Lagrange strains in matrix format
     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl;
     gl(0,0) = glstrain(0);
     gl(0,1) = 0.5*glstrain(3);
     gl(0,2) = 0.5*glstrain(5);
     gl(1,0) = gl(0,1);
     gl(1,1) = glstrain(1);
     gl(1,2) = 0.5*glstrain(4);
     gl(2,0) = gl(0,2);
     gl(2,1) = gl(1,2);
     gl(2,2) = glstrain(2);

     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
     temp.Multiply(gl,defgrd_inv);
     euler_almansi.MultiplyTN(defgrd_inv,temp);

     (*elestrain)(gp,0) = euler_almansi(0,0);
     (*elestrain)(gp,1) = euler_almansi(1,1);
     (*elestrain)(gp,2) = euler_almansi(2,2);
     (*elestrain)(gp,3) = euler_almansi(0,1);
     (*elestrain)(gp,4) = euler_almansi(1,2);
     (*elestrain)(gp,5) = euler_almansi(0,2);
     }
     break;
     case INPAR::STR::strain_none:
     break;
     default:
     dserror("requested strain type not available");
     }*/

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ** Here all possible material laws need to be incorporated,
     ** the stress vector, a C-matrix, and a density must be retrieved,
     ** every necessary data must be passed.
     */
    /*
     double density = 0.0;
     LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
     LINALG::Matrix<NUMSTR_SOH8,1> stress(true);
     soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
     // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc
     */

    /*
     // return gp stresses
     switch (iostress)
     {
     case INPAR::STR::stress_2pk:
     {
     if (elestress == NULL) dserror("stress data not available");
     for (int i = 0; i < NUMSTR_SOH8; ++i)
     (*elestress)(gp,i) = stress(i);
     }
     break;
     case INPAR::STR::stress_cauchy:
     {
     if (elestress == NULL) dserror("stress data not available");
     const double detF = defgrd.Determinant();

     LINALG::Matrix<3,3> pkstress;
     pkstress(0,0) = stress(0);
     pkstress(0,1) = stress(3);
     pkstress(0,2) = stress(5);
     pkstress(1,0) = pkstress(0,1);
     pkstress(1,1) = stress(1);
     pkstress(1,2) = stress(4);
     pkstress(2,0) = pkstress(0,2);
     pkstress(2,1) = pkstress(1,2);
     pkstress(2,2) = stress(2);

     LINALG::Matrix<3,3> temp;
     LINALG::Matrix<3,3> cauchystress;
     temp.Multiply(1.0/detF,defgrd,pkstress,0.0);
     cauchystress.MultiplyNT(temp,defgrd);

     (*elestress)(gp,0) = cauchystress(0,0);
     (*elestress)(gp,1) = cauchystress(1,1);
     (*elestress)(gp,2) = cauchystress(2,2);
     (*elestress)(gp,3) = cauchystress(0,1);
     (*elestress)(gp,4) = cauchystress(1,2);
     (*elestress)(gp,5) = cauchystress(0,2);
     }
     break;
     case INPAR::STR::stress_none:
     break;
     default:
     dserror("requested stress type not available");
     }
     */

    //F^-T * Grad\phi
    LINALG::Matrix<NUMDIM_SOH8,1> Finvgradphi;
    Finvgradphi.MultiplyTT(defgrd_inv, grad_porosity);

    //F^-T * d(Grad\phi)/d(u_s)
    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> Finvdgradphidus;
    Finvdgradphidus.MultiplyTN(defgrd_inv, dgradphi_dus);

    //dF^-T/du_s * Grad(\phi) = - (F^-1 . dN/dX . u_s . F^-1)^T * Grad(\phi)
    LINALG::Matrix<NUMDIM_SOH8,NUMDOF_SOH8> dFinvdus_gradphi(true);
    for (int i=0; i<NUMDIM_SOH8; i++)
    for (int n =0; n<NUMNOD_SOH8; n++)
    for(int j=0; j<NUMDIM_SOH8; j++)
    {
      const int gid = NUMDIM_SOH8 * n +j;
      for (int k=0; k<NUMDIM_SOH8; k++)
      for(int l=0; l<NUMDIM_SOH8; l++)
      dFinvdus_gradphi(i, gid) += -defgrd_inv(l,j) * N_XYZ(k,n) * defgrd_inv(k,i) * grad_porosity(l);
    }

    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> dCinv_dus (true);
    for (int n=0; n<NUMNOD_SOH8; ++n)
    for (int k=0; k<NUMDIM_SOH8; ++k)
    {
      const int gid = n*NUMDIM_SOH8+k;
      for (int i=0; i<NUMDIM_SOH8; ++i)
      {
        dCinv_dus(0,gid) += -2*C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(0,k);
        dCinv_dus(1,gid) += -2*C_inv(1,i)*N_XYZ(i,n)*defgrd_inv(1,k);
        dCinv_dus(2,gid) += -2*C_inv(2,i)*N_XYZ(i,n)*defgrd_inv(2,k);
        /* ~~~ */
        dCinv_dus(3,gid) += -C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(1,k)-defgrd_inv(0,k)*N_XYZ(i,n)*C_inv(1,i);
        dCinv_dus(4,gid) += -C_inv(1,i)*N_XYZ(i,n)*defgrd_inv(2,k)-defgrd_inv(1,k)*N_XYZ(i,n)*C_inv(2,i);
        dCinv_dus(5,gid) += -C_inv(2,i)*N_XYZ(i,n)*defgrd_inv(0,k)-defgrd_inv(2,k)*N_XYZ(i,n)*C_inv(0,i);
      }
    }

    //******************* FAD ************************
    /*
     LINALG::TMatrix<FAD,NUMNOD_SOH8,1>  fad_shapefct;
     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMNOD_SOH8> fad_N_XYZ;
     LINALG::TMatrix<FAD,6,NUMNOD_SOH8> fad_N_XYZ2;
     for (int j=0; j<NUMNOD_SOH8; ++j)
     {
     fad_shapefct(j)=shapefct(j);
     for (int i=0; i<NUMDIM_SOH8; ++i)
     fad_N_XYZ(i,j) = N_XYZ(i,j);
     for (int i=0; i<6 ; i++)
     fad_N_XYZ2(i,j) = N_XYZ2(i,j);
     }

     FAD fad_press = fad_shapefct.Dot(fad_epreaf);
     LINALG::TMatrix<FAD,NUMDIM_SOH8,1> fad_Gradp;
     fad_Gradp.Multiply(fad_N_XYZ,fad_epreaf);

     // compute F
     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMDIM_SOH8> fad_defgrd;
     fad_defgrd.MultiplyNT(fad_xcurr,fad_N_XYZ);
     FAD fad_J = Determinant3x3<FAD>(fad_defgrd);

     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMDIM_SOH8>    fad_defgrd_inv;
     fad_defgrd_inv = fad_defgrd;
     Inverse3x3(fad_defgrd_inv);

     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMDIM_SOH8> fad_cauchygreen;
     fad_cauchygreen.MultiplyTN(fad_defgrd,fad_defgrd);

     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMDIM_SOH8>    fad_C_inv;
     fad_C_inv = fad_cauchygreen;
     Inverse3x3(fad_C_inv);

     LINALG::TMatrix<FAD,6,1> fad_C_inv_vec(false);
     fad_C_inv_vec(0) = fad_C_inv(0,0);
     fad_C_inv_vec(1) = fad_C_inv(1,1);
     fad_C_inv_vec(2) = fad_C_inv(2,2);
     fad_C_inv_vec(3) = fad_C_inv(0,1);
     fad_C_inv_vec(4) = fad_C_inv(1,2);
     fad_C_inv_vec(5) = fad_C_inv(2,0);


     LINALG::TMatrix<FAD,9,1> fad_defgrd_IT_vec;
     fad_defgrd_IT_vec(0)=fad_defgrd_inv(0,0);
     fad_defgrd_IT_vec(1)=fad_defgrd_inv(1,0);
     fad_defgrd_IT_vec(2)=fad_defgrd_inv(2,0);
     fad_defgrd_IT_vec(3)=fad_defgrd_inv(0,1);
     fad_defgrd_IT_vec(4)=fad_defgrd_inv(1,1);
     fad_defgrd_IT_vec(5)=fad_defgrd_inv(2,1);
     fad_defgrd_IT_vec(6)=fad_defgrd_inv(0,2);
     fad_defgrd_IT_vec(7)=fad_defgrd_inv(1,2);
     fad_defgrd_IT_vec(8)=fad_defgrd_inv(2,2);

     LINALG::TMatrix<FAD,NUMDIM_SOH8*NUMDIM_SOH8,NUMDIM_SOH8> fad_F_X(true);
     for(int i=0; i<NUMDIM_SOH8; i++)
     {
     for(int n=0; n<NUMNOD_SOH8; n++)
     {
     fad_F_X(i*NUMDIM_SOH8+0, 0) +=   fad_N_XYZ2(0,n)*fad_nodaldisp(n*NUMDIM_SOH8+i);
     fad_F_X(i*NUMDIM_SOH8+1, 0) +=   fad_N_XYZ2(3,n)*fad_nodaldisp(n*NUMDIM_SOH8+i);
     fad_F_X(i*NUMDIM_SOH8+2, 0) +=   fad_N_XYZ2(4,n)*fad_nodaldisp(n*NUMDIM_SOH8+i);

     fad_F_X(i*NUMDIM_SOH8+0, 1) +=   fad_N_XYZ2(3,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;
     fad_F_X(i*NUMDIM_SOH8+1, 1) +=   fad_N_XYZ2(1,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;
     fad_F_X(i*NUMDIM_SOH8+2, 1) +=   fad_N_XYZ2(5,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;

     fad_F_X(i*NUMDIM_SOH8+0, 2) +=   fad_N_XYZ2(4,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;
     fad_F_X(i*NUMDIM_SOH8+1, 2) +=   fad_N_XYZ2(5,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;
     fad_F_X(i*NUMDIM_SOH8+2, 2) +=   fad_N_XYZ2(2,n)*fad_nodaldisp(n*NUMDIM_SOH8+i) ;
     }
     }

     LINALG::TMatrix<FAD,1,NUMDIM_SOH8> fad_GradJ;
     fad_GradJ.MultiplyTN(fad_J, fad_defgrd_IT_vec, fad_F_X);


     FAD fad_a     = ( bulkmodulus/(1-initporosity) + fad_press - penalty/initporosity ) * fad_J;
     FAD fad_b     = -fad_a + bulkmodulus + penalty;
     FAD fad_c   = (fad_b/fad_a) * (fad_b/fad_a) + 4*penalty/fad_a;
     FAD fad_d     = sign*sqrt(fad_c)*fad_a;

     FAD fad_porosity = 1/(2*fad_a)*(-fad_b+fad_d);

     FAD fad_d_p   =  fad_J * (-fad_b+2*penalty)/fad_d;
     FAD fad_d_J   =  fad_a/fad_J * ( -fad_b + 2*penalty ) / fad_d;

     FAD fad_dphi_dp=  - fad_J * fad_porosity/fad_a + (fad_J+fad_d_p)/(2*fad_a);
     FAD fad_dphi_dJ=  -fad_porosity/fad_J+ 1/(2*fad_J) + fad_d_J / (2*fad_a);

     LINALG::TMatrix<FAD,1,NUMDIM_SOH8>             fad_grad_porosity;
     //      fad_grad_porosity.Update(fad_dphi_dp,fad_Gradp,fad_dphi_dJ,fad_GradJ);
     for (int idim=0; idim<NUMDIM_SOH8; ++idim)
     {
     fad_grad_porosity(idim)=fad_dphi_dp*fad_Gradp(idim)+fad_dphi_dJ*fad_GradJ(idim);
     }

     for (int i=0; i<NUMDOF_SOH8; i++)
     {
     if( (dJ_dus(i)-fad_J.dx(i)) > 1e-8)
     {
     cout<<"dJdus("<<i<<"): "<<dJ_dus(i)<<endl;
     cout<<"fad_J.dx("<<i<<"): "<<fad_J.dx(i)<<endl;
     dserror("check dJdus failed!");
     }
     }
     cout<<"dJdus check done and ok"<<endl;

     for (int i=0; i<NUMDOF_SOH8; i++)
     for (int j=0; j<NUMDIM_SOH8*NUMDIM_SOH8; j++)
     {
     if( (dFinvdus(j,i)-fad_defgrd_IT_vec(j).dx(i)) > 1e-8)
     {
     cout<<"dFinvdus("<<i<<"): "<<dFinvdus(j,i)<<endl;
     cout<<"fad_defgrd_IT_vec.dx("<<i<<"): "<<fad_defgrd_IT_vec(j).dx(i)<<endl;
     dserror("check dFinvdus failed!");
     }
     }
     cout<<"dFinvdus check done and ok"<<endl;

     for (int i=0; i<NUMDOF_SOH8; i++)
     for (int j=0; j<NUMDIM_SOH8; j++)
     {
     if( (dgradJ_dus(j,i)-fad_GradJ(j).dx(i)) > 1e-8)
     {
     cout<<"dgradJ_dus("<<i<<"): "<<dgradJ_dus(j,i)<<endl;
     cout<<"fad_GradJ.dx("<<i<<"): "<<fad_GradJ(j).dx(i)<<endl;
     cout<<"GradJ:"<<endl<<GradJ<<endl;
     cout<<"fad_GradJ:"<<endl<<fad_GradJ<<endl;
     dserror("check dgradJ_dus failed!");
     }
     }
     cout<<"dgradJ_dus check done and ok"<<endl;

     for (int i=0; i<NUMDOF_SOH8; i++)
     if( (dphi_dus(i)-fad_porosity.dx(i)) > 1e-8)
     {
     cout<<"dphi_dus("<<i<<"): "<<dphi_dus(i)<<endl;
     cout<<"fad_porosity.dx("<<i<<"): "<<fad_porosity.dx(i)<<endl;
     cout<<"dphi_dus:"<<endl<<dphi_dus<<endl;
     cout<<"fad_porosity:"<<endl<<fad_porosity<<endl;
     dserror("check dgradJ_dus failed!");
     }
     cout<<"dphi_dus check done and ok"<<endl;

     for (int i=0; i<NUMDOF_SOH8; i++)
     for (int j=0; j<NUMDIM_SOH8; j++)
     {
     if( (dgradphi_dus(j,i)-fad_grad_porosity(j).dx(i)) > 1e-8)
     {
     cout<<"dgradphi_dus("<<i<<"): "<<dgradphi_dus(j,i)<<endl;
     cout<<"fad_grad_porosity.dx("<<i<<"): "<<fad_grad_porosity(j).dx(i)<<endl;
     cout<<"dgradphi_dus:"<<endl<<dgradphi_dus<<endl;
     cout<<"fad_grad_porosity:"<<endl<<fad_grad_porosity<<endl;
     dserror("check dgradphi_dus failed!");
     }
     }
     cout<<"dgradphi_dus check done and ok"<<endl;

     for (int i=0; i<NUMDOF_SOH8; i++)
     for (int j=0; j<6; j++)
     {
     if( (dCinv_dus(j,i)-fad_C_inv_vec(j).dx(i)) > 1e-8)
     {
     cout<<"dCinv_dus("<<i<<"): "<<dCinv_dus(j,i)<<endl;
     cout<<"fad_C_inv.dx("<<i<<"): "<<fad_C_inv_vec(j).dx(i)<<endl;
     cout<<"dCinv_dus:"<<endl<<dCinv_dus<<endl;
     cout<<"fad_C_inv:"<<endl<<fad_C_inv_vec<<endl;
     dserror("check dCinv_dus failed!");
     }
     }
     cout<<"dCinv_dus check done and ok"<<endl;
     */
    //******************* FAD ************************

    //B^T . C^-1
    LINALG::Matrix<NUMDOF_SOH8,1> cinvb(true);
    cinvb.MultiplyTN(bop,C_inv_vec);

    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = detJ*gpweights[gp];
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> estiff_stat(true);
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> erea_u(true);
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> erea_v(true);
    LINALG::Matrix<NUMDOF_SOH8,1> erea_force(true);
    LINALG::Matrix<NUMDOF_SOH8,1> ecoupl_force_p(true);
    LINALG::Matrix<NUMDOF_SOH8,1> ecoupl_force_v(true);
    LINALG::Matrix<NUMDOF_SOH8,NUMNOD_SOH8> erea_p(true);

    if (force != NULL or stiffmatrix != NULL or reamatrix != NULL or forcerea != NULL )
    {
      for (int k=0; k<NUMNOD_SOH8; k++)
      {
        const int fk = NUMDIM_SOH8*k;
        const double fac = detJ_w* shapefct(k);
        const double v = fac * reacoeff * porosity * porosity* J;

        for(int j=0; j<NUMDIM_SOH8; j++)
        {

          /*-------structure- fluid velocity coupling:  RHS
           "dracy-terms"
           - reacoeff * J *  phi^2 *  v^f
           */
          ecoupl_force_v(fk+j) += -v * fvelint(j);

          /* "reactive dracy-terms"
           reacoeff * J *  phi^2 *  v^s
           */
          erea_force(fk+j) += v * velint(j);

          /*-------structure- fluid pressure coupling: RHS
           *                        "porosity gradient terms"
           J *  F^-T * Grad(phi) * p
           */
          ecoupl_force_p(fk+j) += fac * J * Finvgradphi(j) * press;

          for(int i=0; i<NUMNOD_SOH8; i++)
          {
            const int fi = NUMDIM_SOH8*i;

            /* additional "reactive darcy-term"
             detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
             */
            erea_v(fk+j,fi+j) += v * shapefct(i);

            for (int l=0; l<NUMDIM_SOH8; l++)
            {
              /* additional "porosity gradient term" + "darcy-term"
               +  detJ * w(gp) * p * ( J * F^-T * d(Grad(phi))/d(us) + dJ/d(us) * F^-T * Grad(phi) + J * d(F^-T)/d(us) *Grad(phi) ) * D(us)
               - detJ * w(gp) * (  dJ/d(us) * v^f * reacoeff * phi^2 + 2* J * reacoeff * phi * d(phi)/d(us) * v^f ) * D(us)
               */
              estiff_stat(fk+j,fi+l) += fac * ( J * Finvdgradphidus(j,fi+l) * press + press * dJ_dus(fi+l) * Finvgradphi(j)
                  + press * J * dFinvdus_gradphi(j, fi+l)
                  - reacoeff * porosity * ( porosity * dJ_dus(fi+l) + 2 * J * dphi_dus(fi+l) ) * fvelint(j)
              )
              ;

              /* additional "reactive darcy-term"
               detJ * w(gp) * (  dJ/d(us) * vs * reacoeff * phi^2 + 2* J * reacoeff * phi * d(phi)/d(us) * vs ) * D(us)
               */
              erea_u(fk+j,fi+l) += fac * reacoeff * porosity * velint(j) * ( porosity * dJ_dus(fi+l) + 2 * J * dphi_dus(fi+l) );
            }
          }
        }
      }
    }

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    //   force->MultiplyTN(detJ_w, bop, stress, 1.0);

    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> tmp1;
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> tmp2;
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> tmp3;
    //LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> tmp4;

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * (1-\phi) * p^f * detJ * w(gp))
    double fac1 = -detJ_w * (1-porosity) * press;
    tmp1.Multiply(fac1,cinvb,dJ_dus);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * (1-\phi) * p^f * detJ * w(gp))
    double fac2= fac1 * J;
    tmp2.MultiplyTN(fac2,bop,dCinv_dus);

    // additional fluid stress- stiffness term (B^T .  d\phi/d(us) . C^-1  * J * p^f * detJ * w(gp))
    double fac3= detJ_w * press * J;
    tmp3.Multiply(fac3,cinvb,dphi_dus);

    if (force != NULL )
    {
      // additional fluid stress- stiffness term RHS -(B^T .  (1-phi) . C^-1  * J * p^f * detJ * w(gp))
      force->Update(fac2,cinvb,1.0);

      //stationary pressure coupling part of RHS
      // "porosity gradient terms": detJ * w(gp) * J *  F^-T * Grad(phi) * p
      force->Update(1.0,ecoupl_force_p,1.0);

      //stationary velocity coupling part of RHS
      //additional "reactive darcy-term":  - detJ * w(gp) * reacoeff * J *  phi^2 *  v^f
      force->Update(1.0,ecoupl_force_v,1.0);

      //additional "reactive term" RHS  detJ * w(gp) * ( J * reacoeff * phi^2 * v_s)
      force->Update(1.0,erea_force,1.0);

    }

    if ( reamatrix != NULL )
    {
      /* additional "reactive darcy-term"
       detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
       */
      reamatrix->Update(1.0,erea_v,1.0);
    }

    if (stiffmatrix != NULL)
    {
      // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * (1-\phi) * p^f * detJ * w(gp))
      stiffmatrix->Update(1.0,tmp1,1.0);

      // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * (1-\phi) * p^f * detJ * w(gp))
      stiffmatrix->Update(1.0,tmp2,1.0);

      // additional fluid stress- stiffness term (B^T .  d\phi/d(us) . C^-1  * J * p^f * detJ * w(gp))
      stiffmatrix->Update(1.0,tmp3,1.0);

      /* additional "porosity gradient term" + "darcy-term"
       -  detJ * w(gp) * p *( J * F^-T * d(grad(phi))/d(us) + dJ/d(us) * F^-T * grad(phi) + J * d(F^-T)/d(us) *grad(phi) ) * D(us)
       + detJ * w(gp) * ( 2 * reacoeff * phi * d(phi)/d(us) * J * v^f ) * D(us)
       */
      stiffmatrix->Update(1.0,estiff_stat,1.0);

      /* additional "reactive darcy-term"
       detJ * w(gp) * (  dJ/d(us) * vs * reacoeff * phi^2 + 2* J * reacoeff * phi * d(phi)/d(us) * vs ) * D(us)
       */
      stiffmatrix->Update(1.0,erea_u,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(C_inv_vec); // auxiliary integrated stress
      sfac.Scale(fac1); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
        + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
        + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
        + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
          bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************


      //if the reaction part is not supposed to be computed separately, we add it to the stiffness
      if ( reamatrix == NULL )
      {
        stiffmatrix->Update(1.0/dt,erea_v,1.0);
      }
    }
    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  //write porosities at GP into material (for output only)
  structmat->SetPorosityAtGP(porosity_gp);

  return;
}

    /*----------------------------------------------------------------------*
     |  evaluate only the poroelasticity fraction for the element (private)
     *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_coupling_poroelast(vector<int>& lm, // location matrix
    vector<double>& disp, // current displacements
    vector<double>& vel, // current velocities
    //  vector<double>&           residual,       // current residual displ
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> & evelnp, //current fluid velocity
    LINALG::Matrix<NUMNOD_SOH8, 1> & epreaf, //current fluid pressure
    LINALG::Matrix<NUMDOF_SOH8, (NUMDIM_SOH8 + 1) * NUMNOD_SOH8>* stiffmatrix, // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, (NUMDIM_SOH8 + 1) * NUMNOD_SOH8>* reamatrix, // element reactive matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force, // element internal force vector
    LINALG::Matrix<NUMDOF_SOH8, 1>* forcerea, // element reactive force vector
    ParameterList& params) // algorithmic parameters e.g. time
{
  /* ============================================================================*
   ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
   ** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8, 1> > shapefcts =
      soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> > derivs =
      soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  //=============================get parameters
  // get global id of the structure element
  int id = Id();
  //access fluid discretization
  RCP<DRT::Discretization> fluiddis = null;
  fluiddis = DRT::Problem::Instance()->Dis(genprob.numff, 0);
  //get corresponding fluid element
  DRT::Element* fluidele = fluiddis->gElement(id);
  if (fluidele == NULL)
    dserror("Fluid element %i not on local processor", id);

  //get fluid material
  const MAT::FluidPoro* fluidmat = static_cast<const MAT::FluidPoro*>((fluidele->Material()).get());
  if(fluidmat->MaterialType() != INPAR::MAT::m_fluidporo)
    dserror("invalid fluid material for poroelasticity");

  //get structure material
  const MAT::StructPoro* structmat = static_cast<const MAT::StructPoro*>((Material()).get());
  if(structmat->MaterialType() != INPAR::MAT::m_structporo)
    dserror("invalid structure material for poroelasticity");

  const double reacoeff = fluidmat->ComputeReactionCoeff();
  const double bulkmodulus = structmat->Bulkmodulus();
  const double penalty = structmat->Penaltyparameter();

  const double initporosity = params.get<double>("initporosity");
  double theta = params.get<double>("theta");
  //double dt   = params.get<double>("delta time");

  //=======================================================================

  // update element geometry
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xrefe; // material coord. of element
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xcurr; // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;

  // (r,s,t) gp-locations of fully integrated linear 8-node Hex
  const double gploc = 1.0/sqrt(3.0); // gp sampling point value for linear fct
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> xgp; // coord. of guass point in reference coordinates
  xgp(0,0)=-gploc;
  xgp(1,0)=-gploc;
  xgp(2,0)=-gploc;

  xgp(0,1)= gploc;
  xgp(1,1)=-gploc;
  xgp(2,1)=-gploc;

  xgp(0,2)= gploc;
  xgp(1,2)= gploc;
  xgp(2,2)=-gploc;

  xgp(0,3)=-gploc;
  xgp(1,3)= gploc;
  xgp(2,3)=-gploc;

  xgp(0,4)=-gploc;
  xgp(1,4)=-gploc;
  xgp(2,4)= gploc;

  xgp(0,5)= gploc;
  xgp(1,5)=-gploc;
  xgp(2,5)= gploc;

  xgp(0,6)= gploc;
  xgp(1,6)= gploc;
  xgp(2,6)= gploc;

  xgp(0,7)=-gploc;
  xgp(1,7)= gploc;
  xgp(2,7)= gploc;

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(0,i) = x[0];
    xrefe(1,i) = x[1];
    xrefe(2,i) = x[2];

    xcurr(0,i) = xrefe(0,i) + disp[i*NODDOF_SOH8+0];
    xcurr(1,i) = xrefe(1,i) + disp[i*NODDOF_SOH8+1];
    xcurr(2,i) = xrefe(2,i) + disp[i*NODDOF_SOH8+2];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodaldisp;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodaldisp(i,0) = disp[i];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodalvel;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodalvel(i,0) = vel[i];
  }

  LINALG::Matrix<NUMDOF_SOH8,(NUMDIM_SOH8+1)*NUMNOD_SOH8> ecoupl(true);
  LINALG::Matrix<NUMDOF_SOH8,NUMNOD_SOH8> ecoupl_p(true);
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> ecoupl_v(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  LINALG::Matrix<6,NUMNOD_SOH8> N_XYZ2;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);
  LINALG::Matrix<6,NUMNOD_SOH8> deriv2; //  second derivatives at gausspoint w.r.t. r,s,t
  LINALG::Matrix<NUMDIM_SOH8,1> xsi;
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    //get shapefunctions and its derivatives at gausspoint
    LINALG::Matrix<NUMNOD_SOH8,1> shapefct = shapefcts[gp];
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> deriv = derivs[gp];
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invJ = invJ_[gp];

    //get coordinates of gausspoint and compute second derivatives at gausspoint
    xsi(0)=xgp(0,gp);
    xsi(1)=xgp(1,gp);
    xsi(2)=xgp(2,gp);
    DRT::UTILS::shape_function_deriv2<DRT::Element::hex8>(xsi,deriv2);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ X_,r  Y_,r  Z_,r ]^-1
     **     J^-1 = [ X_,s  Y_,s  Z_,s ]
     **            [ X_,t  Y_,t  Z_,t ]
     */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // transposed jacobian "dX/ds"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xjm0;
    xjm0.MultiplyNT(deriv,xrefe);

    // get the second derivatives of standard element at current GP w.r.t. xyz
    DRT::UTILS::gder2<DRT::Element::hex8>(xjm0,N_XYZ,deriv2,xrefe,N_XYZ2);

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xjm;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> xji;
    xjm.MultiplyNT(deriv,xcurr);
    const double det = xji.Invert(xjm);

    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    const double J = det/detJ;

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr,N_XYZ); //  (6.17)
    //defgrd.PutScalar(0.0);
    //for (int i=0; i<3; i++)
    //  defgrd(i,i) =1;

    // set to initial state as test to receive a linear solution
    //for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

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
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
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

    // -----------------Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // ------------------Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::Matrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    //------------------ inverse Right Cauchy-Green tensor
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> C_inv(false);
    C_inv.Invert(cauchygreen);

    //-----------inverse Right Cauchy-Green tensor as vector
    LINALG::Matrix<6,1> C_inv_vec(true);
    C_inv_vec(0) = C_inv(0,0);
    C_inv_vec(1) = C_inv(1,1);
    C_inv_vec(2) = C_inv(2,2);
    C_inv_vec(3) = C_inv(0,1);
    C_inv_vec(4) = C_inv(1,2);
    C_inv_vec(5) = C_inv(2,0);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    //--------------------- get fluid velocity at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    //! ----------------structure displacement and velocity at integration point
    LINALG::Matrix<NUMDIM_SOH8,1> dispint(true);
    LINALG::Matrix<NUMDIM_SOH8,1> velint(true);
    for(int i=0; i<NUMNOD_SOH8; i++)
    for(int j=0; j<NUMDIM_SOH8; j++)
    {
      dispint(j) += nodaldisp(i*NUMDIM_SOH8+j) * shapefct(i);
      velint(j) += nodalvel(i*NUMDIM_SOH8+j) * shapefct(i);
    }

    // inverse deformation gradient F^-1
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //------------------------------------ build F^-1 as vector 9x1
    LINALG::Matrix<9,1> defgrd_inv_vec;
    defgrd_inv_vec(0)=defgrd_inv(0,0);
    defgrd_inv_vec(1)=defgrd_inv(0,1);
    defgrd_inv_vec(2)=defgrd_inv(0,2);
    defgrd_inv_vec(3)=defgrd_inv(1,0);
    defgrd_inv_vec(4)=defgrd_inv(1,1);
    defgrd_inv_vec(5)=defgrd_inv(1,2);
    defgrd_inv_vec(6)=defgrd_inv(2,0);
    defgrd_inv_vec(7)=defgrd_inv(2,1);
    defgrd_inv_vec(8)=defgrd_inv(2,2);

    //------------------------------------ build F^-T as vector 9x1
    LINALG::Matrix<9,1> defgrd_IT_vec;
    defgrd_IT_vec(0)=defgrd_inv(0,0);
    defgrd_IT_vec(1)=defgrd_inv(1,0);
    defgrd_IT_vec(2)=defgrd_inv(2,0);
    defgrd_IT_vec(3)=defgrd_inv(0,1);
    defgrd_IT_vec(4)=defgrd_inv(1,1);
    defgrd_IT_vec(5)=defgrd_inv(2,1);
    defgrd_IT_vec(6)=defgrd_inv(0,2);
    defgrd_IT_vec(7)=defgrd_inv(1,2);
    defgrd_IT_vec(8)=defgrd_inv(2,2);

    //--------------------------- build N_x operator (wrt material config)
    LINALG::Matrix<9,NUMDOF_SOH8> N_X(true); // set to zero
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      N_X(0,3*i+0) = N_XYZ(0,i);
      N_X(1,3*i+1) = N_XYZ(0,i);
      N_X(2,3*i+2) = N_XYZ(0,i);

      N_X(3,3*i+0) = N_XYZ(1,i);
      N_X(4,3*i+1) = N_XYZ(1,i);
      N_X(5,3*i+2) = N_XYZ(1,i);

      N_X(6,3*i+0) = N_XYZ(2,i);
      N_X(7,3*i+1) = N_XYZ(2,i);
      N_X(8,3*i+2) = N_XYZ(2,i);
    }

    LINALG::Matrix<NUMDIM_SOH8*NUMDIM_SOH8,NUMDIM_SOH8> F_X(true);
    for(int i=0; i<NUMDIM_SOH8; i++)
    {
      for(int n=0; n<NUMNOD_SOH8; n++)
      {
        F_X(i*NUMDIM_SOH8+0, 0) += N_XYZ2(0,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 0) += N_XYZ2(3,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 0) += N_XYZ2(4,n)*nodaldisp(n*NUMDIM_SOH8+i);

        F_X(i*NUMDIM_SOH8+0, 1) += N_XYZ2(3,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 1) += N_XYZ2(1,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 1) += N_XYZ2(5,n)*nodaldisp(n*NUMDIM_SOH8+i);

        F_X(i*NUMDIM_SOH8+0, 2) += N_XYZ2(4,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+1, 2) += N_XYZ2(5,n)*nodaldisp(n*NUMDIM_SOH8+i);
        F_X(i*NUMDIM_SOH8+2, 2) += N_XYZ2(2,n)*nodaldisp(n*NUMDIM_SOH8+i);
      }
    }

    //--------------material gradient of jacobi determinant J: GradJ = dJ/dX= dJ/dF : dF/dX = J * F^-T : dF/dX
    LINALG::Matrix<1,NUMDIM_SOH8> GradJ;
    GradJ.MultiplyTN(J, defgrd_IT_vec, F_X);

    //**************************************************+auxilary variables for computing the porosity and linearization
    const double a = ( bulkmodulus/(1-initporosity) + press - penalty/initporosity ) * J;
    const double b = -a + bulkmodulus + penalty;
    const double c = (b/a) * (b/a) + 4*penalty/a;
    double d = sqrt(c)*a;
    double sign =1.0;
    double test = 1/(2*a)*(-b+d);
    if( test >= 1.0 or test < 0.0 )
    {
      sign = -1.0;
      d = sign*d;
    }

    const double porosity = 1/(2*a)*(-b+d);

    if( porosity >= 1.0 or porosity < 0.0 )
    {
      dserror("invalid porosity!");
    }

    const double d_p = J * (-b+2*penalty)/d;
    const double d_p_p = ( d * J + d_p * (b - 2*penalty) ) / (d * d) * J;
    const double d_J = a/J * ( -b + 2*penalty ) / d;
    const double d_J_p = d_p / J + ( 1-d_p*d_p/(J*J) ) / d *a;

    //double porosity = structmat->ComputePorosity(press, J, initporosity,gp);

    //d(porosity) / d(p)
    const double dphi_dp= - J * porosity/a + (J+d_p)/(2*a);

    //d^2(porosity) / d(pressure)^2
    const double dphi_dpp= -J/a*dphi_dp + porosity*J*J/(a*a) - J/(2*a*a)*(J+d_p) + d_p_p/(2*a);

    //d(porosity) / d(J)
    const double dphi_dJ= -porosity/J+ 1/(2*J) + d_J / (2*a);

    //d(porosity) / d(J)d(pressure)
    double dphi_dJdp= -1/J*dphi_dp+ d_J_p/(2*a) - d_J*J/(2*a*a);

    //-----------material porosity gradient
    LINALG::Matrix<1,NUMDIM_SOH8> grad_porosity;
    for (int idim=0; idim<NUMDIM_SOH8; ++idim)
    {
      grad_porosity(idim)=dphi_dp*Gradp(idim)+dphi_dJ*GradJ(idim);
    }

    //----------------linearization of material porosity gradient w.r.t fluid pressure
    // d(Grad(phi))/dp = d^2(phi)/(dJ*dp) * GradJ * N + d^2(phi)/(dp)^2 * Gradp * N + d(phi)/dp * N,X
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> dgradphi_dp;
    dgradphi_dp.MultiplyTT(dphi_dJdp,GradJ,shapefct);
    dgradphi_dp.MultiplyNT(dphi_dpp,Gradp,shapefct,1.0);
    dgradphi_dp.Update(dphi_dp,N_XYZ,1.0);

    //******************* FAD ************************
    /*
     // sacado data type replaces "double"
     typedef Sacado::Fad::DFad<double> FAD;  // for first derivs
     // sacado data type replaces "double" (for first+second derivs)
     typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD;


     LINALG::TMatrix<FAD,NUMNOD_SOH8,1>  fad_shapefct;
     LINALG::TMatrix<FAD,NUMDIM_SOH8,NUMNOD_SOH8> fad_N_XYZ(false);
     LINALG::TMatrix<FAD,6,NUMNOD_SOH8> fad_N_XYZ2(false);
     for (int j=0; j<NUMNOD_SOH8; ++j)
     {
     fad_shapefct(j)=shapefct(j);
     for (int i=0; i<NUMDIM_SOH8; ++i)
     fad_N_XYZ(i,j) = N_XYZ(i,j);
     for (int i=0; i<6 ; i++)
     fad_N_XYZ2(i,j) = N_XYZ2(i,j);
     }

     LINALG::TMatrix<FAD,NUMDOF_SOH8,1> fad_nodaldisp(false);
     LINALG::TMatrix<FAD,NUMNOD_SOH8,1> fad_epreaf(false);
     for(int i=0; i<NUMNOD_SOH8 ; i++)
     {
     fad_epreaf(i) = epreaf(i);
     fad_epreaf(i).diff(i,NUMNOD_SOH8);
     }

     FAD fad_press = fad_shapefct.Dot(fad_epreaf);
     LINALG::TMatrix<FAD,NUMDIM_SOH8,1> fad_Gradp;
     fad_Gradp.Multiply(fad_N_XYZ,fad_epreaf);

     LINALG::TMatrix<FAD,1,NUMDIM_SOH8> fad_GradJ;
     for(int i=0; i<NUMDIM_SOH8; i++)
     fad_GradJ(i)=GradJ(i);

     FAD fad_a     = ( bulkmodulus/(1-initporosity) + fad_press - penalty/initporosity ) * J;
     FAD fad_b     = -fad_a + bulkmodulus + penalty;
     FAD fad_c   = (fad_b/fad_a) * (fad_b/fad_a) + 4*penalty/fad_a;
     FAD fad_d     = sign*sqrt(fad_c)*fad_a;

     FAD fad_porosity = 1/(2*fad_a)*(-fad_b+fad_d);

     FAD fad_d_p   =  J * (-fad_b+2*penalty)/fad_d;
     FAD fad_d_J   =  fad_a/J * ( -fad_b + 2*penalty ) / fad_d;

     FAD fad_dphi_dp=  - J * fad_porosity/fad_a + (J+fad_d_p)/(2*fad_a);
     FAD fad_dphi_dJ=  -fad_porosity/J+ 1/(2*J) + fad_d_J / (2*fad_a);

     LINALG::TMatrix<FAD,1,NUMDIM_SOH8>             fad_grad_porosity;
     for (int idim=0; idim<NUMDIM_SOH8; ++idim)
     {
     fad_grad_porosity(idim)=fad_dphi_dp*fad_Gradp(idim)+fad_dphi_dJ*fad_GradJ(idim);
     }


     for (int i=0; i<NUMNOD_SOH8; i++)
     for (int j=0; j<NUMDIM_SOH8; j++)
     {
     if( (dgradphi_dp(j,i)-fad_grad_porosity(j).dx(i)) > 1e-8)
     {
     cout<<"dgradphi_dp("<<i<<"): "<<dgradphi_dp(j,i)<<endl;
     cout<<"fad_grad_porosity.dx("<<i<<"): "<<fad_grad_porosity(j).dx(i)<<endl;
     cout<<"dgradphi_dp:"<<endl<<dgradphi_dp<<endl;
     cout<<"fad_grad_porosity:"<<endl<<fad_grad_porosity<<endl;
     dserror("check dgradphi_dus failed!");
     }
     }
     cout<<"dgradphi_dp check done and ok"<<endl;
     */
    //******************* FAD ************************

    //*****************************************************************************************

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ** Here all possible material laws need to be incorporated,
     ** the stress vector, a C-matrix, and a density must be retrieved,
     ** every necessary data must be passed.
     */
    /*
     double density = 0.0;
     LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
     LINALG::Matrix<NUMSTR_SOH8,1> stress(true);
     soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
     // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc
     */

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = detJ*gpweights[gp];

    //B^T . C^-1
    LINALG::Matrix<NUMDOF_SOH8,1> cinvb(true);
    cinvb.MultiplyTN(bop,C_inv_vec);

    //F^-T * grad\phi
    LINALG::Matrix<NUMDIM_SOH8,1> Finvgradphi;
    Finvgradphi.MultiplyTT(defgrd_inv, grad_porosity);

    //F^-T * dgrad\phi/dp
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> Finvgradphidp;
    Finvgradphidp.MultiplyTN(defgrd_inv, dgradphi_dp);

    if (force != NULL or stiffmatrix != NULL or reamatrix != NULL or forcerea != NULL )
    {
      for (int i=0; i<NUMNOD_SOH8; i++)
      {
        const int fi = NUMDIM_SOH8*i;
        const double fac = detJ_w* shapefct(i);

        for(int j=0; j<NUMDIM_SOH8; j++)
        {
          for(int k=0; k<NUMNOD_SOH8; k++)
          {
            const int fk = NUMDIM_SOH8*k;

            /*-------structure- fluid pressure coupling: "stress terms" + "porosity gradient terms"
             -B^T . ( (1-phi)*J*C^-1 - d(phi)/(dp)*p*J*C^-1 ) * Dp
             + J * F^-T * Grad(phi) * Dp + J * F^-T * d(Grad((phi))/(dp) * p * Dp
             */
            ecoupl_p(fi+j, k) += detJ_w * cinvb(fi+j) * ( -(1-porosity)
                + dphi_dp * press
            ) * J * shapefct(k)
            + fac * J * ( Finvgradphi(j) * shapefct(k) + Finvgradphidp(j,k) * press )
            ;

            /*-------structure- fluid pressure coupling:  "dracy-terms" + "reactive darcy-terms"
             - 2 * reacoeff * J * v^f * phi * d(phi)/dp  Dp
             + 2 * reacoeff * J * v^s * phi * d(phi)/dp  Dp
             */
            const double tmp = fac * reacoeff * J * 2 * porosity * dphi_dp * shapefct(k);
            ecoupl_p(fi+j, k) += -tmp * fvelint(j);

            ecoupl_p(fi+j, k) += tmp * velint(j);

            /*-------structure- fluid velocity coupling:  "dracy-terms"
             -reacoeff * J *  phi^2 *  Dv^f
             */
            ecoupl_v(fi+j, fk+j) += -fac * reacoeff * J * porosity * porosity * shapefct(k);

          }
        }
      }
    }
    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  if( forcerea != NULL )
  {
    //instationary pressure part of RHS
    //forcerea->Multiply(1.0,erea_p,epreaf,1.0);
  }

  if (force != NULL )
  {
    //all rhs terms are added in soh8_nlnstiff_poroelast
  }

  if (stiffmatrix != NULL or reamatrix != NULL)
  {
    // add structure displacement - fluid velocity part to matrix
    for (int ui=0; ui<NUMNOD_SOH8; ++ui)
    {
      const int dim_ui = NUMDIM_SOH8*ui;

      for (int jdim=0; jdim < NUMDIM_SOH8;++jdim)
      {
        const int dim_ui_jdim = dim_ui+jdim;

        for (int vi=0; vi<NUMNOD_SOH8; ++vi)
        {
          const int numdof_vi = (NUMDIM_SOH8+1)*vi;
          const int dim_vi = NUMDIM_SOH8*vi;

          for (int idim=0; idim <NUMDIM_SOH8; ++idim)
          {
            ecoupl(dim_ui_jdim , numdof_vi+idim) += ecoupl_v(dim_ui_jdim , dim_vi+idim);
          }
        }
      }
    }

  // add structure displacement - fluid pressure part to matrix
    for (int ui=0; ui<NUMNOD_SOH8; ++ui)
    {
      const int dim_ui = NUMDIM_SOH8*ui;

      for (int jdim=0; jdim < NUMDIM_SOH8;++jdim)
      {
        const int dim_ui_jdim = dim_ui+jdim;

        for (int vi=0; vi<NUMNOD_SOH8; ++vi)
        {
          ecoupl( dim_ui_jdim , (NUMDIM_SOH8+1)*vi+NUMDIM_SOH8 ) += ecoupl_p( dim_ui_jdim , vi);
        }
      }
    }
  }

  if ( reamatrix != NULL )
  {
    //reamatrix->Update(1.0,erea,1.0);
  }

  if (stiffmatrix != NULL)
  {
    // build tangent coupling matrix : effective dynamic stiffness coupling matrix
    //    K_{Teffdyn} = 1/dt C
    //                + theta K_{T}
    stiffmatrix->Update(theta,ecoupl,1.0);
  }

  return;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

