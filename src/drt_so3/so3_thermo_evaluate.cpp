/*!----------------------------------------------------------------------
\file so3_thermo_evaluate.cpp
\brief

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so3_thermo.H"
#include "so3_thermo_fwd.hpp"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include <iterator>

// headers of thermo-materials
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/robinson.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"


/*----------------------------------------------------------------------*
 | pre-evaluate the element (public)                         dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::PreEvaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la
  )
{
  // if the coupling variables are required before Evaluate() is called the 1st
  // time
  // here for Robinson's material
//  if(la.Size()>1)
//  {
//    // the temperature field has only one dof per node, disregarded by the
//    // dimension of the problem
//    const int numdofpernode_thr = NumDofPerNode(1,*(Nodes()[0]));
//    if (la[1].Size() != nen_*numdofpernode_thr)
//      dserror("Location vector length for temperatures does not match!");
//
//    if (discretization.HasState(1,"temperature"))
//    {
//      // check if you can get the temperature state
//      Teuchos::RCP<const Epetra_Vector> tempnp
//        = discretization.GetState(1,"temperature");
//      if (tempnp == Teuchos::null)
//        dserror("Cannot get state vector 'tempnp'.");
//
//      // extract local values of the global vectors
//      Teuchos::RCP<std::vector<double> >mytempnp
//        = rcp(new std::vector<double>(la[1].lm_.size()) );
//      DRT::UTILS::ExtractMyValues(*tempnp,*mytempnp,la[1].lm_);
//
//      // now set the current temperature vector in the parameter list
//      params.set<Teuchos::RCP<vector<double> > >("tempnp",mytempnp);
//    }
//  } // initial temperature dependence

  return;
}


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Thermo< so3_ele, distype>::Evaluate(
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
  // start with "none"
  typename So3_Thermo::ActionType act = So3_Thermo::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_stifftemp")   act = So3_Thermo::calc_struct_stifftemp;

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // coupling terms K_dT in stiffness matrix K^{TSI} for monolithic TSI
  case So3_Thermo::calc_struct_stifftemp:
  {
    if (la.Size()>1)
    {
      CouplToThrEvaluate(
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
  }
  break;

  //==================================================================================
  default:
  {
    // in some cases we need to write/change some data before evaluating
    // here you can pass, e.g. for Robinson's material the temperature even for
    PreEvaluate(
      params,
      discretization,
      la
      );

    // call the purely structural methods
    so3_ele::Evaluate(
      params,
      discretization,
      la[0].lm_,  // only the first column, i.e. the structural field is passed
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );

    // add the temperature-dependent terms to the structural field, i.e.
    // it's a TSI problem
    if (la.Size()>1)
    {
      CouplToThrEvaluate(
        params,
        discretization,
        la,  // coupled TSI is considered, i.e. pass the compled location array
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra,
        elevec2_epetra,
        elevec3_epetra
        );
    }

  }  // default
  break;

  }  // action

  return 0;
}  // Evaluate()


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 | here is the action for the coupling to the thermal field             |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::CouplToThrEvaluate(
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
  // start with "none"
  ActionType act = none;

  // get the required action for coupling with the thermal field
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_internalforce")  act = calc_struct_internalforce;
  else if (action=="calc_struct_nlnstiff")       act = calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")   act = calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")  act = calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stifftemp")      act = calc_struct_stifftemp;
  else if (action=="calc_struct_stress")         act = calc_struct_stress;
  else if (action=="calc_struct_reset_istep")    act = calc_struct_reset_istep;
  else if (action=="calc_struct_update_istep")   act = calc_struct_update_istep;
  else dserror("Unknown type of action for So3_Thermo: %s",action.c_str());

  // what should the element do
  switch(act)
  {
  //============================================================================
  // internal force vector for TSI only
  case calc_struct_internalforce:
  {
    // internal force vector
    LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);
    // elemat1+2, elevec2+3 are not used anyway

    // need current displacement and residual/incremental displacements
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0,"residual displacement");

    if (disp==Teuchos::null || res==Teuchos::null)
      dserror("Cannot get state vectors 'displacement' and/or residual");

    // build the location vector only for the structure field
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);
    // create a dummy element matrix to apply linearised EAS-stuff onto
    LINALG::Matrix<numdofperelement_,numdofperelement_> myemat(true);

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage the
    // coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector/nsd_
    std::vector<double> mytempnp( ( (la[0].lm_).size() )/nsd_, 0.0 );

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
        = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_thr = NumDofPerNode(1,*(Nodes()[0]));
      if (la[1].Size() != nen_*numdofpernode_thr)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
    }

    // default: geometrically non-linear analysis with Total Lagrangean approach
    if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_nonlinear)
    {
      dserror("will follow soon!");
    }  // kintype_==nonlinear
    // geometric linear
    else if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_linear)
    {
      // calculate the THERMOmechanical term for fint
      so3_thermo_fint_lin(
        la,
        mydisp,
        myres,
        mytempnp,
        &elevec1,
        NULL,
        NULL,
        params,
        INPAR::STR::stress_none,
        INPAR::STR::strain_none
        );
    }  // kintype_==linear
  }  // calc_struct_internalforce
  break;

  //============================================================================
  // (non)linear stiffness for TSI
  case calc_struct_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<numdofperelement_,numdofperelement_> elemat1(elemat1_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);
    // elemat2, elevec2+3 are not used anyway

    // need current displacement and residual/incremental displacements
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0,"residual displacement");

    if (disp==Teuchos::null || res==Teuchos::null)
      dserror("Cannot get state vectors 'displacement' and/or residual");

    // build the location vector only for the structure field
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage the
    // coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector/nsd_
    std::vector<double> mytempnp( ( (la[0].lm_).size() )/nsd_, 0.0 );

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
        = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_thr = NumDofPerNode(1,*(Nodes()[0]));
      if (la[1].Size() != nen_*numdofpernode_thr)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // default: geometrically non-linear analysis with Total Lagrangean approach
      if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_nonlinear)
      {
        dserror("will follow soon!");
      }  // kintype_==nonlinear
      // geometric linear
      else if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_linear)
      {
        // calculate the THERMOmechanical term for fint
        so3_thermo_fint_lin(
          la,
          mydisp,
          myres,
          mytempnp,
          &elevec1,
          NULL,
          NULL,
          params,
          INPAR::STR::stress_none,
          INPAR::STR::strain_none
          );
      }  // kintype_==linear
    }
  }  // calc_struct_nlnstiff
  break;

  //============================================================================
  // (non)linear stiffness, mass matrix and internal force vector for TSI
  case calc_struct_nlnstiffmass:
  case calc_struct_nlnstifflmass:
  {
    // stiffness
    LINALG::Matrix<numdofperelement_,numdofperelement_> elemat1(elemat1_epetra.A(),true);
    // massmatrix
    LINALG::Matrix<numdofperelement_,numdofperelement_> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);
    // elevec2+3 are not used anyway

    // need current displacement and residual/incremental displacements
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0,"residual displacement");

    if (disp==Teuchos::null || res==Teuchos::null)
      dserror("Cannot get state vectors 'displacement' and/or residual");

    // build the location vector only for the structure field
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage the
    // coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector/nsd_
    std::vector<double> mytempnp( ( (la[0].lm_).size() )/nsd_, 0.0 );

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
        = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_thr = NumDofPerNode(1,*(Nodes()[0]));
      if (la[1].Size() != nen_*numdofpernode_thr)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // default: geometrically non-linear analysis with Total Lagrangean approach
      if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_nonlinear)
      {
        dserror("will follow soon!");
      }  // kintype_==nonlinear
      // geometric linear
      else if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_linear)
      {
        // build the current temperature vector
        LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[1]),true);  // view only!
        // calculate the THERMOmechanical term for fint
        so3_thermo_fint_lin(
          la,
          mydisp,
          myres,
          mytempnp,
          &elevec1,
          NULL,
          NULL,
          params,
          INPAR::STR::stress_none,
          INPAR::STR::strain_none
          );
      }  // kintype_==linear
    }

    if (act==calc_struct_nlnstifflmass)
      so3_thermo_lumpmass(&elemat2);

  }  // calc_struct_nlnstiff(l)mass
  break;

  //==================================================================================
  // evaluate stresses and strains at gauss points
  case calc_struct_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway
    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==so3_ele::Owner())
    {
      Teuchos::RCP<const Epetra_Vector> disp
        = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res
        = discretization.GetState(0,"residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null)
        dserror("Cannot get state vectors 'displacement'");

      vector<double> mydisp((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
      vector<double> myres((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

      Teuchos::RCP<vector<char> > stressdata
        = params.get<Teuchos::RCP<vector<char> > >("stress", Teuchos::null);
      if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
//      Teuchos::RCP<vector<char> > straindata
//        = params.get<Teuchos::RCP<vector<char> > >("strain", Teuchos::null);
//      Teuchos::RCP<vector<char> > plstraindata
//        = params.get<Teuchos::RCP<vector<char> > >("plstrain", Teuchos::null);
//      if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
//      if (plstraindata==Teuchos::null) dserror("Cannot get 'plastic strain' data");

      LINALG::Matrix<numgpt_1,nstr_> stress;
      // strains are required for robinson's material
//      LINALG::Matrix<numgpt_1,nstr_> strain;
//      LINALG::Matrix<numgpt_1,nstr_> plstrain;

      INPAR::STR::StressType iostress
        = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress",
            INPAR::STR::stress_none);
//      INPAR::STR::StrainType iostrain
//        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain",
//            INPAR::STR::strain_none);

      // initialise the vectors
      // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage the
      // coupling field is not yet known. Pass coupling vectors filled with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp( ( (la[0].lm_).size() )/nsd_, 0.0 );

      // call the purely structural method
      // If a visco-plastic Robinson's material is used evaluate the element using
      // the current temperature
      // that is NOT a beautiful implementation, but it works
//TODO 2012-10-26      if (mat->MaterialType() != INPAR::MAT::m_vp_robinson)
//      {
//        soh8_linstiffmass(la[0].lm_,mydisp,myres,NULL,NULL,NULL,NULL,&stress,
//          &strain,&plstrain,params,iostress,iostrain,ioplstrain);
//      }

      // need current temperature state,
      // call the temperature discretization: thermo equates 2nd dofset
      // disassemble temperature
      if (discretization.HasState(1,"temperature"))
      {
        // check if you can get the temperature state
        Teuchos::RCP<const Epetra_Vector> tempnp
          = discretization.GetState(1,"temperature");
        if (tempnp==Teuchos::null)
          dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = NumDofPerNode(1,*(Nodes()[0]));
        if (la[1].Size() != nen_*numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");

        // extract the current temperatures
        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

        // purely structural method, this is the coupled routine, i.e., a 2nd
        // discretisation exists, i.e., --> we always have a temperature state
        // If a visco-plastic Robinson's material is used evaluate the element using
        // the current temperature
        // that is NOT a beautiful implementation, but it works
        // get the temperature dependent stress
        LINALG::Matrix<numgpt_1,nstr_> stresstemp(true);
//TODO 2012-10-26        if (mat->MaterialType() == INPAR::MAT::m_vp_robinson)
//        {
//          soh8_linstiffmass(la[0].lm_,mydisp,myres,&mytempnp,NULL,NULL,NULL,
//            &stress,&strain,&plstrain,params,iostress,iostrain,ioplstrain);
//        }
//        else
//        {

        // default: geometrically non-linear analysis with Total Lagrangean approach
        if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_nonlinear)
        {
          dserror("will follow soon!");
        }  // kintype_==nonlinear
        // geometric linear
        else if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_linear)
        {
          // calculate the THERMOmechanical term for fint: temperature stresses
          so3_thermo_fint_lin(
            la,
            mydisp,
            myres,
            mytempnp,
            NULL,
            &stresstemp,
            NULL,
            params,
            iostress,
            INPAR::STR::strain_none);
        }  // kintype_==linear
//        }

        // total stress
        // add stresstemp to the mechanical stress
        // stress = stress_d + stress_T
        stress.Update(1.0,stresstemp,1.0);
      }

      {
        DRT::PackBuffer data;
        so3_ele::AddtoPack(data, stress);
        data.StartPacking();
        so3_ele::AddtoPack(data, stress);
        std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
      }

//      {
//        DRT::PackBuffer data;
//        so3_ele::AddtoPack(data, strain);
//        data.StartPacking();
//        so3_ele::AddtoPack(data, strain);
//        std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
//      }
//
//      {
//        DRT::PackBuffer data;
//        so3_ele::AddtoPack(data, plstrain);
//        data.StartPacking();
//        so3_ele::AddtoPack(data, plstrain);
//        std::copy(data().begin(),data().end(),std::back_inserter(*plstraindata));
//      }
    }  // end proc Owner
  }  // calc_struct_stress
  break;

  //============================================================================
  // required for predictor TangDis --> can be helpful in compressible case!
  case calc_struct_reset_istep:
  {
    // do nothing;
  }
  break;

  //============================================================================
  case calc_struct_update_istep:
  {
    // TODO 2012-10-31 check if Update has to be called here again
    Teuchos::RCP<MAT::Material> mat = so3_ele::Material();

    // Update of history for thermo-(visco-)plastic material if they exist
    if (mat->MaterialType() == INPAR::MAT::m_thermopllinelast)
    {
      MAT::ThermoPlasticLinElast* thrpllinelast = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      thrpllinelast->Update();
    }
    // incremental update of internal variables/history
    else if (mat->MaterialType() == INPAR::MAT::m_vp_robinson)
    {
      MAT::Robinson* robinson = static_cast <MAT::Robinson*>(mat.get());
      bool imrlike = false;
      robinson->Update(imrlike, 0.0);
    }
  }  // calc_struct_update_istep
  break;

  //============================================================================
  // coupling term k_dT of stiffness matrix for monolithic TSI
  case calc_struct_stifftemp:
  {
    // mechanical-thermal system matrix
    LINALG::Matrix<numdofperelement_,nen_> stiffmatrixcoupl(elemat1_epetra.A(),true);
    // elemat2,elevec1-3 are not used anyway

    // default: geometrically non-linear analysis with Total Lagrangean approach
    if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_nonlinear)
    {
      dserror("will follow soon!");
    }  // kintype_==linear
    // geometric linear
    else if (kintype_ == DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_linear)
    {
      // need current displacement and residual/incremental displacements
      Teuchos::RCP<const Epetra_Vector> disp
        = discretization.GetState(0,"displacement");
      if (disp==Teuchos::null)
        dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp((la[0].lm_).size());
      // build the location vector only for the structure field
      DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);

      // calculate the mechanical-thermal system matrix
      so3_thermo_kdt_lin(la,mydisp,&stiffmatrixcoupl);

    }  // kintype_==linear
  }  // calc_struct_stifftemp
  break;

  //============================================================================
  default:
    dserror("Unknown type of action for So3_Thermo");
  } // action

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 05/10 |
 | contribution to r_d (private)                                        |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_fint_lin(
  DRT::Element::LocationArray& la,  // location array
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  vector<double>& temp, // current temperature
  LINALG::Matrix<numdofperelement_,1>* force,  // element internal force vector
  LINALG::Matrix<numgpt_1,nstr_>* elestress,  // stresses at GP
  LINALG::Matrix<numgpt_1,nstr_>* elestrain,  // strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{
  // update element geometry hex8, 3D: (8x3)
  LINALG::Matrix<nen_,nsd_> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // x, current  coord. of element
  // vector of the current element temperatures
  LINALG::Matrix<nen_,1> etemp(true);

  DRT::Node** nodes = Nodes();
  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];

    etemp(i,0) = temp[i+0];
  }

  // we need the (residual) displacement -- current increment of displacement
  LINALG::Matrix<numdofperelement_,1> res_d;
  for (int i = 0; i<numdofperelement_; ++i)
  {
    res_d(i) = residual[i];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  LINALG::Matrix<nsd_,nen_> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<nsd_,nsd_> defgrd(true);
  // shape functions and their first derivatives
  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<nsd_,nen_> deriv;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // geometrically linear, i.e. reference == current state, i.e. F == I
    // set to initial state (defgrd == identity)
    for (int i=0; i<3; ++i)
      defgrd(i,i) = 1.0;

    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<nstr_,numdofperelement_> boplin;
    for (int i=0; i<nen_; ++i)
    {
      boplin(0,numdofpernode_*i+0) = N_XYZ(0,i);
      boplin(0,numdofpernode_*i+1) = 0.0;
      boplin(0,numdofpernode_*i+2) = 0.0;
      boplin(1,numdofpernode_*i+0) = 0.0;
      boplin(1,numdofpernode_*i+1) = N_XYZ(1,i);
      boplin(1,numdofpernode_*i+2) = 0.0;
      boplin(2,numdofpernode_*i+0) = 0.0;
      boplin(2,numdofpernode_*i+1) = 0.0;
      boplin(2,numdofpernode_*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,numdofpernode_*i+0) = N_XYZ(1,i);
      boplin(3,numdofpernode_*i+1) = N_XYZ(0,i);
      boplin(3,numdofpernode_*i+2) = 0.0;
      boplin(4,numdofpernode_*i+0) = 0.0;
      boplin(4,numdofpernode_*i+1) = N_XYZ(2,i);
      boplin(4,numdofpernode_*i+2) = N_XYZ(1,i);
      boplin(5,numdofpernode_*i+0) = N_XYZ(2,i);
      boplin(5,numdofpernode_*i+1) = 0.0;
      boplin(5,numdofpernode_*i+2) = N_XYZ(0,i);
    }
    // copy structural shape functions needed for the thermo field
    // identical shapefunctions for the displacements and the temperatures

    // product of shapefunctions and element temperatures for stresstemp
    // N_T . T
    LINALG::Matrix<1,1> Ntemp(true);
    Ntemp.MultiplyTN(shapefunct,etemp);
    const double scalartemp  = shapefunct.Dot(etemp);

    // build iterative strains
    LINALG::Matrix<nstr_,1> straininc(true);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    // calculate the stress part dependent on the temperature in the material
    LINALG::Matrix<nstr_,1> ctemp(true);
    LINALG::Matrix<nstr_,1> stresstemp(true);
    LINALG::Matrix<nstr_,nstr_> cmat(true);
    LINALG::Matrix<nstr_,1> glstrain(true);
    LINALG::Matrix<nstr_,1> plglstrain(true);
    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) stresstemp = C . Delta T
    // do not call the material for Robinson's material
//    // TODO 2012-10-30
//    if ( !(Material()->MaterialType() == INPAR::MAT::m_vp_robinson) )
    so3_thermo_materialize(&stresstemp,&ctemp,&Ntemp,&cmat,&defgrd,&glstrain,&plglstrain,
      straininc,scalartemp,&density,gp,params);

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

      // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress==NULL) dserror("stress data not available");
      for (int i=0; i<nstr_; ++i)
        (*elestress)(gp,i) = stresstemp(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress==NULL) dserror("stress data not available");

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

      LINALG::Matrix<nsd_,nsd_> temp;
      LINALG::Matrix<nsd_,nsd_> cauchystresstemp;
      temp.Multiply(1.0,defgrd,pkstresstemp);
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

    // integrate internal force vector r_d
    // f = f + (B^T . sigma_temp) * detJ * w(gp)
    if (force != NULL)
    {
      // old implementation hex8_thermo double detJ_w = detJ*gpweights[gp];
      double detJ_w = detJ*intpoints_.Weight(gp);//gpweights[gp];
      force->MultiplyTN(detJ_w, boplin, stresstemp, 1.0);
    }  // if (force != NULL)

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So3_Thermo::so3_thermo_fint_lin


/*----------------------------------------------------------------------*
 | evaluate only the mechanical-thermal stiffness term       dano 03/11 |
 | for monolithic TSI, contribution to k_dT (private)                   |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_kdt_lin(
  DRT::Element::LocationArray& la,
  vector<double>& disp,  // current displacement
  LINALG::Matrix<numdofperelement_,nen_>* stiffmatrixcoupl  // (nsd_*nen_ x nen_)
  )
{
  // update element geometry (8x3)
  LINALG::Matrix<nen_,nsd_> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // x, current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  LINALG::Matrix<nsd_,nen_> N_XYZ;
  // shape functions and their first derivatives
  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<nsd_,nen_> deriv;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // linear B-operator B_L = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    //
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<nstr_,numdofperelement_> boplin;
    for (int i=0; i<nen_; ++i)
    {
      boplin(0,numdofpernode_*i+0) = N_XYZ(0,i);
      boplin(0,numdofpernode_*i+1) = 0.0;
      boplin(0,numdofpernode_*i+2) = 0.0;
      boplin(1,numdofpernode_*i+0) = 0.0;
      boplin(1,numdofpernode_*i+1) = N_XYZ(1,i);
      boplin(1,numdofpernode_*i+2) = 0.0;
      boplin(2,numdofpernode_*i+0) = 0.0;
      boplin(2,numdofpernode_*i+1) = 0.0;
      boplin(2,numdofpernode_*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,numdofpernode_*i+0) = N_XYZ(1,i);
      boplin(3,numdofpernode_*i+1) = N_XYZ(0,i);
      boplin(3,numdofpernode_*i+2) = 0.0;
      boplin(4,numdofpernode_*i+0) = 0.0;
      boplin(4,numdofpernode_*i+1) = N_XYZ(2,i);
      boplin(4,numdofpernode_*i+2) = N_XYZ(1,i);
      boplin(5,numdofpernode_*i+0) = N_XYZ(2,i);
      boplin(5,numdofpernode_*i+1) = 0.0;
      boplin(5,numdofpernode_*i+2) = N_XYZ(0,i);
    }

    // copy structural shape functions needed for the thermo field
    // identical shapefunctions for the displacements and the temperatures

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated
    */
    // get the thermal material tangent
    LINALG::Matrix<nstr_,1> ctemp(true);
    Ctemp(&ctemp);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double detJ_w = detJ*intpoints_.Weight(gp);
    // update linear coupling matrix K_dT
    if (stiffmatrixcoupl != NULL)
    {
      // C_temp . N_temp
      LINALG::Matrix<nstr_,nen_> cn(true);
      cn.MultiplyNT(ctemp,shapefunct); // (6x8)=(6x1)(1x8)
      // integrate stiffness term
      // k_dT = k_dT + (B^T . C_temp . N_temp) * detJ * w(gp)
      stiffmatrixcoupl->MultiplyTN(detJ_w, boplin, cn, 1.0);
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
}  // so3_thermo_kdt_lin()


/*----------------------------------------------------------------------*
 | material law with temperature part for So3_thermo         dano 05/10 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_materialize(
  LINALG::Matrix<nstr_,1>* stresstemp,
  LINALG::Matrix<nstr_,1>* ctemp,
  LINALG::Matrix<1,1>* Ntemp,  // temperature of element
  LINALG::Matrix<nstr_,nstr_>* cmat,
  LINALG::Matrix<nsd_,nsd_>* defgrd, //
  LINALG::Matrix<nstr_,1>* glstrain,
  LINALG::Matrix<nstr_,1>* plglstrain,
  LINALG::Matrix<nstr_,1>& straininc,
  const double& scalartemp,
  double* density,
  const int gp,
  Teuchos::ParameterList& params
  )
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stresstemp) dserror("No stress vector supplied");
  if (!ctemp) dserror("No material tangent matrix supplied");
//  if (!Ntemp) dserror("No temperature supplied");
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
    case INPAR::MAT::m_vp_robinson: /*-- visco-plastic Robinson's material */
    {
      MAT::Robinson* robinson = static_cast <MAT::Robinson*>(mat.get());
      robinson->Evaluate(*glstrain,*plglstrain,straininc,scalartemp,gp,params,*cmat,*stresstemp);
      *density = robinson->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of temperature dependent material");
    break;
  } // switch (mat->MaterialType())

  return;
}  // so3_thermo_materialize()


/*----------------------------------------------------------------------*
 | get the constant temperature fraction for stresstemp      dano 05/10 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Ctemp(
  LINALG::Matrix<nstr_,1>* ctemp
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
    case INPAR::MAT::m_vp_robinson: /*-- visco-plastic Robinson's material */
    {
      // so far: do nothing, because the displacement-dependent coupling term
      // is neglected
      return;
      break;
    }
    default:
      dserror("Cannot ask material for the temperature-dependent material tangent");
      break;
  } // switch (mat->MaterialType())

}  // Ctemp()


/*----------------------------------------------------------------------*
 | lump mass matrix (private)                               bborn 07/08 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::so3_thermo_lumpmass(
  LINALG::Matrix<numdofperelement_,numdofperelement_>* emass
  )
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
}  // so3_thermo_lumpmass


/*----------------------------------------------------------------------*
 | initialise Jacobian                                       dano 08/12 |
 | in called once in Initialize() in so3_thermo_eletypes.cpp            |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::InitJacobianMapping()
{
  // get the material coordinates
  LINALG::Matrix<nen_,nsd_> xrefe;
  for (int i=0; i<nen_; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  // initialise the derivatives of the shape functions
  LINALG::Matrix<nsd_,nen_> deriv;

  // coordinates of the current integration point (xsi_)
  for (int gp=0; gp<numgpt_; ++gp)
  {
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim=0;idim<nsd_;idim++)
    {
      xsi_[gp](idim) = gpcoord[idim];
    }
    // first derivatives of shape functions (deriv)
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    // compute Jacobian matrix and determinant
    // actually compute its transpose....
    /*
      +-            -+ T      +-            -+
      | dx   dx   dx |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dr   dr   dr |
      |              |        |              |
      | dy   dy   dy |        | dx   dy   dz |
      | --   --   -- |   =    | --   --   -- |
      | dr   ds   dt |        | ds   ds   ds |
      |              |        |              |
      | dz   dz   dz |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dt   dt   dt |
      +-            -+        +-            -+
     */
    // derivatives of coordinates w.r.t material coordinates xjm_ = dx/ds
    invJ_[gp].Multiply(deriv,xrefe);
    // xij_ = ds/dx
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] < 1.0E-16)
      dserror("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",detJ_[gp]);
  }

  return;
}  // InitJacobianMapping


/*----------------------------------------------------------------------*/
