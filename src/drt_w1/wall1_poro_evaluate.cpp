/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_evaluate.cpp

 \brief evaluate methods for 2D wall element for structure part of porous medium

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

\level 2
 *----------------------------------------------------------------------*/


#include "wall1_poro.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include <iterator>

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/matlist_reactions.H"
#include "../drt_mat/fluidporo_multiphase.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"
//#include "Sacado.hpp"

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                  vuong 12/12      |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::PreEvaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization&      discretization,
                                        DRT::Element::LocationArray& la)
{
  if(scatracoupling_)
  {
    if(la.Size()>2)
    {
      if (discretization.HasState(2,"scalar"))
      {
        // check if you can get the scalar state
        Teuchos::RCP<const Epetra_Vector> scalarnp
          = discretization.GetState(2,"scalar");

        // extract local values of the global vectors
        std::vector<double> myscalar (la[2].lm_.size());
        DRT::UTILS::ExtractMyValues(*scalarnp,myscalar,la[2].lm_);

        if(NumMaterial()<3)
          dserror("no third material defined for Wall poro element!");
        Teuchos::RCP<MAT::Material> scatramat = Material(2);

        int numscal=1;
        if( scatramat->MaterialType() == INPAR::MAT::m_matlist or
            scatramat->MaterialType() == INPAR::MAT::m_matlist_reactions
            )
        {
          Teuchos::RCP<MAT::MatList> matlist = Teuchos::rcp_dynamic_cast<MAT::MatList>(scatramat);
          numscal = matlist->NumMat();
        }

        Teuchos::RCP<std::vector<double> > scalar = Teuchos::rcp( new std::vector<double>(numscal,0.0) );
        if((int)myscalar.size() != numscal*numnod_)
          dserror("sizes do not match!");

        for(int i=0; i<numnod_; i++)
          for(int j=0; j<numscal; j++)
            scalar->at(j) += myscalar[numscal*i+j]/numnod_;

        params.set("scalar",scalar);
      }
    }
    else
    {
//      const double time = params.get("total time",0.0);
//    // find out whether we will use a time curve and get the factor
//      int num = 0; // TO BE READ FROM INPUTFILE AT EACH ELEMENT!!!
//      std::vector<double> xrefe; xrefe.resize(3);
//      DRT::Node** nodes = Nodes();
//      // get displacements of this element
//    //  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
//     for (int i=0; i<numnod_; ++i){
//        const double* x = nodes[i]->X();
//        xrefe [0] +=  x[0]/numnod_;
//        xrefe [1] +=  x[1]/numnod_;
//        xrefe [2] +=  x[2]/numnod_;
//
//      }
//      const double* coordgpref = &xrefe[0];
//      double functfac = DRT::Problem::Instance()->Funct(num).Evaluate(0,coordgpref,time);
//      params.set<double>("scalar",functfac);

//      Teuchos::RCP<std::vector<double> > scalar = Teuchos::rcp( new std::vector<double>(4,0.0) );
//      for(int i=0; i<numnod_; i++)
//      {
//          scalar->at(0) = 0.0;
//          scalar->at(1) = 17.6;
//          scalar->at(2) = 0.0;
//          scalar->at(3) = 0.0;
//      }
//
//      params.set("scalar",scalar);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                         vuong 12/12  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_Poro<distype>::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  if(not init_)
    dserror("internal element data not initialized!");

  this->SetParamsInterfacePtr(params);
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (this->IsParamsInterface())
  {
    act = this->ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action=="struct_poro_calc_fluidcoupling")       act = ELEMENTS::struct_poro_calc_fluidcoupling;
    else if (action=="struct_poro_calc_scatracoupling")      act = ELEMENTS::struct_poro_calc_scatracoupling;
    else if (action=="struct_poro_calc_prescoupling")  act = ELEMENTS::struct_poro_calc_prescoupling;
  }

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // off diagonal terms in stiffness matrix for monolithic coupling
  case ELEMENTS::struct_poro_calc_fluidcoupling:
  case ELEMENTS::struct_poro_calc_scatracoupling:
  case ELEMENTS::struct_poro_calc_prescoupling:
  {
    //in some cases we need to write/change some data before evaluating
    PreEvaluate(params,
                discretization,
                la);

    MyEvaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;
  //==================================================================================
  default:
  {
    //in some cases we need to write/change some data before evaluating
    PreEvaluate(params,
                discretization,
                la);

    //evaluate parent solid element
    DRT::ELEMENTS::Wall1::Evaluate(params,
                      discretization,
                      la[0].lm_,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);

    //add volume coupling specific terms
   MyEvaluate(params,
              discretization,
              la,
              elemat1_epetra,
              elemat2_epetra,
              elevec1_epetra,
              elevec2_epetra,
              elevec3_epetra);
  }
  break;
  } // action

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        vuong 12/12    |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_Poro<distype>::MyEvaluate(
                                    Teuchos::ParameterList&      params,
                                    DRT::Discretization&         discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix&    elemat1_epetra,
                                    Epetra_SerialDenseMatrix&    elemat2_epetra,
                                    Epetra_SerialDenseVector&    elevec1_epetra,
                                    Epetra_SerialDenseVector&    elevec2_epetra,
                                    Epetra_SerialDenseVector&    elevec3_epetra
                                    )
{
  this->SetParamsInterfacePtr(params);
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (this->IsParamsInterface())
  {
    act = this->ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action=="calc_struct_internalforce")            act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_nlnstiff")                 act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_nlnstiffmass")             act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="struct_poro_calc_fluidcoupling")       act = ELEMENTS::struct_poro_calc_fluidcoupling;
    else if (action=="calc_struct_stress")                   act = ELEMENTS::struct_calc_stress;
    else if (action=="struct_poro_calc_prescoupling")  act = ELEMENTS::struct_poro_calc_prescoupling;
    //else if (action=="postprocess_stress")                act = postprocess_stress;
    //else dserror("Unknown type of action for Wall1_Poro: %s",action.c_str());
  }


  // --------------------------------------------------
  // Now do the nurbs specific stuff
  if (Shape() == DRT::Element::nurbs4 || Shape() == DRT::Element::nurbs9)
  {
    myknots_.resize(2);

    switch (act)
    {
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    case ELEMENTS::struct_poro_calc_prescoupling:
    case ELEMENTS::struct_calc_stress:
    {
      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

      bool zero_sized = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots_,Id());

      // skip zero sized elements in knot span --- they correspond to interpolated nodes
      if (zero_sized)
      {
        return (0);
      }

      break;
    }
    default:
      myknots_.clear();
      break;
    }
  }

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // nonlinear stiffness, damping and internal force vector for poroelasticity
  case ELEMENTS::struct_calc_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<numdof_,numdof_> elemat1(elemat1_epetra.A(),true);
    //damping
    LINALG::Matrix<numdof_,numdof_> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);
    // elemat2,elevec2+3 are not used anyway

    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdim_,numnod_> mydisp(true);
    ExtractValuesFromGlobalVector(discretization,0,lm, &mydisp, NULL,"displacement");

    LINALG::Matrix<numdof_,numdof_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;

    enum INPAR::STR::DampKind damping = params.get<enum INPAR::STR::DampKind>("damping",INPAR::STR::damp_none);
    LINALG::Matrix<numdof_,numdof_>* matptr2 = NULL;
    if (elemat2.IsInitialized() and (damping==INPAR::STR::damp_material) ) matptr2 = &elemat2;

    if(la.Size()>1)
    {
      if (discretization.HasState(1,"fluidvel"))
      {
        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures
        LINALG::Matrix<numdim_,numnod_> myvel(true);
        LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
        LINALG::Matrix<numnod_,1> myepreaf(true);

        if (discretization.HasState(0,"velocity"))
          ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

        // extract local values of the global vectors
        ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

        //calculate tangent stiffness matrix
        nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,matptr,matptr2,&elevec1,//NULL,
            //NULL,NULL,
            params);
      }
      else if (la.Size()>2)
        if(discretization.HasState(2,"solid_pressure"))
        {
          // extract local values of the global vectors
          LINALG::Matrix<numnod_,1> myepreaf(true);
          ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myepreaf,"solid_pressure");

          //calculate tangent stiffness matrix
          nlnstiff_poroelast_presbased(lm,mydisp,myepreaf,matptr,&elevec1,params);
        }

    }
  }
  break;

  //==================================================================================
  // nonlinear stiffness, mass matrix and internal force vector for poroelasticity
  case ELEMENTS::struct_calc_nlnstiffmass:
  {

    // stiffness
    LINALG::Matrix<numdof_,numdof_> elemat1(elemat1_epetra.A(),true);
    // mass
    LINALG::Matrix<numdof_,numdof_> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);
    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdim_,numnod_> mydisp(true);
    ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");

    LINALG::Matrix<numdof_,numdof_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;

    // we skip this evaluation if the coupling is not setup yet, i.e.
    // if the secondary dofset or the secondary material was not set
    // this can happen during setup of the time integrator or restart
    // TODO: there might be a better way. For instance do not evaluate
    //       before the setup of the multiphysics problem is completed.
    if(la.Size()>1 and NumMaterial() > 1)
    {
      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures

      LINALG::Matrix<numdim_,numnod_> myvel(true);
      LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
      LINALG::Matrix<numnod_,1> myepreaf(true);

      if (discretization.HasState(0,"velocity"))
        ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      //TODO: this is kind of a hack. Find a better way! (e.g. move the pressure based variant into own element)
      if (discretization.HasState(1,"fluidvel"))
      {
        // extract local values of the global vectors
        ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

        nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,matptr,NULL,&elevec1,//NULL,
            //NULL,NULL,
            params);
      }
      else if (la.Size()>2)
        if(discretization.HasState(2,"solid_pressure"))
        {
          // extract local values of the global vectors
          ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myepreaf,"solid_pressure");

          nlnstiff_poroelast_presbased(lm,mydisp,myepreaf,matptr,&elevec1,params);
        }
    }

  }
  break;

  //==================================================================================
  // coupling terms in force-vector and stiffness matrix for poroelasticity
  case ELEMENTS::struct_poro_calc_fluidcoupling:
  {
    // stiffness
    LINALG::Matrix<numdof_,(numdim_+1)*numnod_> elemat1(elemat1_epetra.A(),true);
    //LINALG::Matrix<numdof_,(numdim_+1)*numnod_> elemat2(elemat2_epetra.A(),true);

    // internal force vector
    //LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    //LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);

    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;

    // need current fluid state,
    // call the fluid discretization: fluid equates 2nd dofset
    // disassemble velocities and pressures
    if (discretization.HasState(1,"fluidvel"))
    {
      LINALG::Matrix<numdim_,numnod_> myvel(true);
      LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
      LINALG::Matrix<numnod_,1> myepreaf(true);

      LINALG::Matrix<numdim_,numnod_> mydisp(true);
      ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");

      if (discretization.HasState(0,"velocity"))
        ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      // extract local values of the global vectors
      ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

        coupling_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,matptr,//NULL,
          NULL,NULL,params);
    }
    else if (la.Size()>2)
    {
      if(discretization.HasState(2,"solid_pressure") && discretization.HasState(1,"porofluid"))
      {
        //TODO: move this to function, once split is performed
        std::vector<double> myephi(la[1].Size());
        Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"porofluid");
        DRT::UTILS::ExtractMyValues(*matrix_state,myephi,la[1].lm_);
        GetMaterials_presbased();
        const int numphases = fluidmultimat_->NumMat();
        Epetra_SerialDenseMatrix ephi(numphases, numnod_,true);

        for (int i = 0; i < numnod_; i++)
        {
          for (int j = 0; j < numphases; j++)
          {
            ephi(j,i) = myephi[i*numphases+j];
          }
        }
        LINALG::Matrix<numdim_,numnod_> mydisp(true);
        ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");

        coupling_poroelast_presbased(lm,mydisp,ephi,elemat1_epetra,params);
      }
      else
        dserror("cannot find global states displacement or solidpressure");
    }

  }
  break;

  //==================================================================================
  // nonlinear stiffness and internal force vector for poroelasticity
  case ELEMENTS::struct_calc_internalforce:
  {
    // stiffness
    LINALG::Matrix<numdof_,numdof_> elemat1(elemat1_epetra.A(),true);
    LINALG::Matrix<numdof_,numdof_> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);
    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdim_,numnod_> mydisp(true);
    ExtractValuesFromGlobalVector(discretization,0,lm, &mydisp, NULL,"displacement");

    // need current fluid state,
    // call the fluid discretization: fluid equates 2nd dofset
    // disassemble velocities and pressures
    if (discretization.HasState(1,"fluidvel"))
    {
      // extract local values of the global vectors
      LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
      LINALG::Matrix<numnod_,1> myepreaf(true);
      ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

      LINALG::Matrix<numdim_,numnod_> myvel(true);
      ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      //calculate tangent stiffness matrix
      nlnstiff_poroelast(lm,mydisp,myvel,myfluidvel,myepreaf,NULL,NULL,&elevec1,//NULL,
          //NULL,NULL,
          params);
    }
    else if (la.Size()>2)
      if(discretization.HasState(2,"solid_pressure"))
      {
        // extract local values of the global vectors
        LINALG::Matrix<numnod_,1> myepreaf(true);
        ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myepreaf,"solid_pressure");

        //calculate tangent stiffness matrix
        nlnstiff_poroelast_presbased(lm,mydisp,myepreaf,NULL,&elevec1,params);
      }
  }
  break;
  //==================================================================================
  // evaluate stresses and strains at gauss points
  case ELEMENTS::struct_calc_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==Owner())
    {
      // get the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      LINALG::Matrix<numdim_,numnod_> mydisp(true);
      ExtractValuesFromGlobalVector(discretization,0,lm, &mydisp, NULL,"displacement");

      Teuchos::RCP<std::vector<char> > couplstressdata
        = params.get<Teuchos::RCP<std::vector<char> > >("couplstress", Teuchos::null);

      if (couplstressdata==Teuchos::null) dserror("Cannot get 'couplstress' data");

      Epetra_SerialDenseMatrix couplstress(numgpt_,Wall1::numstr_);

      INPAR::STR::StressType iocouplstress
        = DRT::INPUT::get<INPAR::STR::StressType>(params, "iocouplstress",
            INPAR::STR::stress_none);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1,"fluidvel"))
      {
        // extract local values of the global vectors
        LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
        LINALG::Matrix<numnod_,1> myepreaf(true);
        ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

        couplstress_poroelast(mydisp,
                              myfluidvel,
                              myepreaf,
                              &couplstress,
                              NULL,
                              params,
                              iocouplstress);
      }
      else if (la.Size()>2)
        if(discretization.HasState(2,"solid_pressure"))
        {
          // extract local values of the global vectors
          LINALG::Matrix<numdim_,numnod_> myfluidvel(true);
          LINALG::Matrix<numnod_,1> myepreaf(true);
          ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myepreaf,"solid_pressure");

          couplstress_poroelast(mydisp,
                                myfluidvel,
                                myepreaf,
                                &couplstress,
                                NULL,
                                params,
                                iocouplstress);
        }

      // pack the data for postprocessing
      {
        DRT::PackBuffer data;
        // get the size of stress
        Wall1::AddtoPack(data, couplstress);
        data.StartPacking();
        // pack the stresses
        Wall1::AddtoPack(data, couplstress);
        std::copy(data().begin(),data().end(),std::back_inserter(*couplstressdata));
      }
    }  // end proc Owner
  }  // ELEMENTS::struct_calc_stress
  break;

  //==================================================================================
  default:
    //do nothing
    break;
  } // action
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (private) vuong 12/12 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::nlnstiff_poroelast(
    std::vector<int>&                  lm,           // location matrix
    LINALG::Matrix<numdim_, numnod_>&  disp,         // current displacements
    LINALG::Matrix<numdim_, numnod_>&  vel,          // current velocities
    LINALG::Matrix<numdim_, numnod_> & evelnp,       // current fluid velocities
    LINALG::Matrix<numnod_, 1> &       epreaf,       // current fluid pressure
    LINALG::Matrix<numdof_, numdof_>*  stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_, numdof_>*  reamatrix,    // element reactive matrix
    LINALG::Matrix<numdof_, 1>*        force,        // element internal force vector
    Teuchos::ParameterList&            params        // algorithmic parameters e.g. time
    )
{
  GetMaterials();

  // update element geometry
  LINALG::Matrix<numdim_,numnod_> xrefe; // material coord. of element
  LINALG::Matrix<numdim_,numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  LINALG::Matrix<numdof_,numdof_> erea_v(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoop(  params,
                        xrefe,
                        xcurr,
                        disp,
                        vel,
                        evelnp,
                        epreaf,
                        NULL,
                        erea_v,
                        stiffmatrix,
                        reamatrix,
                        force);

  if ( reamatrix != NULL )
  {
    /* additional "reactive darcy-term"
     detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
     */
    reamatrix->Update(1.0,erea_v,1.0);
  }

  return;
}  // nlnstiff_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (private) vuong 12/12 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::nlnstiff_poroelast_presbased(
    std::vector<int>&                  lm,           // location matrix
    LINALG::Matrix<numdim_, numnod_>&  disp,         // current displacements
    LINALG::Matrix<numnod_, 1> &       epreaf,       // current fluid pressure
    LINALG::Matrix<numdof_, numdof_>*  stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_, 1>*        force,        // element internal force vector
    Teuchos::ParameterList&            params        // algorithmic parameters e.g. time
    )
{
  GetMaterials_presbased();

  // update element geometry
  LINALG::Matrix<numdim_,numnod_> xrefe; // material coord. of element
  LINALG::Matrix<numdim_,numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoop_presbased(
                        params,
                        xrefe,
                        xcurr,
                        disp,
                        epreaf,
                        stiffmatrix,
                        force);

  return;
}  // nlnstiff_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::GaussPointLoop(
                                    Teuchos::ParameterList& params,
                                    const LINALG::Matrix<numdim_,numnod_>& xrefe,
                                    const LINALG::Matrix<numdim_,numnod_>& xcurr,
                                    const LINALG::Matrix<numdim_,numnod_>& nodaldisp,
                                    const LINALG::Matrix<numdim_,numnod_>& nodalvel,
                                    const LINALG::Matrix<numdim_,numnod_> & evelnp,
                                    const LINALG::Matrix<numnod_,1> & epreaf,
                                    const LINALG::Matrix<numnod_, 1>*  porosity_dof,
                                    LINALG::Matrix<numdof_,numdof_>& erea_v,
                                    LINALG::Matrix<numdof_, numdof_>*  stiffmatrix,
                                    LINALG::Matrix<numdof_,numdof_>* reamatrix,
                                    LINALG::Matrix<numdof_,1>*              force
                                        )
{

  /*--------------------------------- get node weights for nurbs elements */
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnod_; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights_(inode) = cp->W();
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  // first derivatives N_XYZ at gp w.r.t. material coordinates
  LINALG::Matrix<numdim_,numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<numdim_,numdim_> defgrd(true);
  // shape function at gp w.r.t. reference coordinates
  LINALG::Matrix<numnod_,1> shapefct;
  // first derivatives at gp w.r.t. reference coordinates
  LINALG::Matrix<numdim_,numnod_> deriv ;

  LINALG::Matrix<numstr_,1> fstress(true);

  for (int gp=0; gp<numgpt_; ++gp)
  {

    //evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    //compute deformation gradient
    ComputeDefGradient(defgrd,N_XYZ,xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<numdim_,numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) = dJ/dF : dF/dus = J * F^-T * N,X
    static LINALG::Matrix<1,numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static LINALG::Matrix<1,numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J,
        volchange,
        dJ_dus,
        dvolchange_dus,
        defgrd,
        defgrd_inv,
        N_XYZ,
        nodaldisp);

    // non-linear B-operator
    LINALG::Matrix<numstr_,numdof_> bop;
    ComputeBOperator(bop,defgrd,N_XYZ);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // pressure gradient at integration point
    LINALG::Matrix<numdim_,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    // fluid velocity at integration point
    LINALG::Matrix<numdim_,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<numdim_,numdim_>              fvelder;
    fvelder.MultiplyNT(evelnp,N_XYZ);

    // structure displacement and velocity at integration point
    LINALG::Matrix<numdim_,1> velint(true);

    for(int i=0; i<numnod_; i++)
      for(int j=0; j<numdim_; j++)
        velint(j) += nodalvel(j,i) * shapefct(i);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<numdim_,numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<numdim_,numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // compute some auxiliary matrixes for computation of linearization
    //dF^-T/dus
    LINALG::Matrix<numdim_*numdim_,numdof_> dFinvTdus(true);
    //F^-T * Grad p
    LINALG::Matrix<numdim_,1> Finvgradp;
    //dF^-T/dus * Grad p
    LINALG::Matrix<numdim_,numdof_> dFinvdus_gradp(true);
    //dC^-1/dus * Grad p
    LINALG::Matrix<numstr_,numdof_> dCinv_dus (true);

    ComputeAuxiliaryValues(N_XYZ,defgrd_inv,C_inv,Gradp,dFinvTdus,Finvgradp,dFinvdus_gradp,dCinv_dus);

    //--------------------------------------------------------------------

    //linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1,numdof_> dphi_dus;
    double porosity=0.0;

    ComputePorosityAndLinearization(params,press,volchange,gp,shapefct,porosity_dof,dvolchange_dus,porosity,dphi_dus);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    if(fluidmat_->Type() == MAT::PAR::darcy_brinkman )
    {
      FillMatrixAndVectorsBrinkman(
                                    gp,
                                    J,
                                    porosity,
                                    fvelder,
                                    defgrd_inv,
                                    bop,
                                    C_inv,
                                    dphi_dus,
                                    dJ_dus,
                                    dCinv_dus,
                                    dFinvTdus,
                                    stiffmatrix,
                                    force,
                                    fstress);
    }

    FillMatrixAndVectors(   gp,
                            shapefct,
                            N_XYZ,
                            J,
                            press,
                            porosity,
                            velint,
                            fvelint,
                            fvelder,
                            defgrd_inv,
                            bop,
                            C_inv,
                            Finvgradp,
                            dphi_dus,
                            dJ_dus,
                            dCinv_dus,
                            dFinvdus_gradp,
                            dFinvTdus,
                            erea_v,
                            stiffmatrix,
                            force,
                            fstress);
  }//end of gaussloop
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::GaussPointLoop_presbased(
                                    Teuchos::ParameterList& params,
                                    const LINALG::Matrix<numdim_,numnod_>&  xrefe,
                                    const LINALG::Matrix<numdim_,numnod_>&  xcurr,
                                    const LINALG::Matrix<numdim_,numnod_>&  nodaldisp,
                                    const LINALG::Matrix<numnod_,1> &       epreaf,
                                    LINALG::Matrix<numdof_, numdof_>*       stiffmatrix,
                                    LINALG::Matrix<numdof_,1>*              force
                                        )
{

  /*--------------------------------- get node weights for nurbs elements */
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnod_; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights_(inode) = cp->W();
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  // first derivatives N_XYZ at gp w.r.t. material coordinates
  LINALG::Matrix<numdim_,numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<numdim_,numdim_> defgrd(true);
  // shape function at gp w.r.t. reference coordinates
  LINALG::Matrix<numnod_,1> shapefct;
  // first derivatives at gp w.r.t. reference coordinates
  LINALG::Matrix<numdim_,numnod_> deriv ;

  LINALG::Matrix<numstr_,1> fstress(true);

  for (int gp=0; gp<numgpt_; ++gp)
  {

    //evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    //compute deformation gradient
    ComputeDefGradient(defgrd,N_XYZ,xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<numdim_,numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) = dJ/dF : dF/dus = J * F^-T * N,X
    static LINALG::Matrix<1,numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static LINALG::Matrix<1,numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J,
        volchange,
        dJ_dus,
        dvolchange_dus,
        defgrd,
        defgrd_inv,
        N_XYZ,
        nodaldisp);

    // non-linear B-operator
    LINALG::Matrix<numstr_,numdof_> bop;
    ComputeBOperator(bop,defgrd,N_XYZ);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<numdim_,numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<numdim_,numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // compute some auxiliary matrixes for computation of linearization
    //dC^-1/dus
    LINALG::Matrix<numstr_,numdof_> dCinv_dus (true);
    for (int n=0; n<numnod_; ++n)
      for (int k=0; k<numdim_; ++k)
      {
        const int gid = n*numdim_+k;
        for (int i=0; i<numdim_; ++i)
        {
          dCinv_dus(0,gid) += -2*C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(0,k);
          dCinv_dus(1,gid) += -2*C_inv(1,i)*N_XYZ(i,n)*defgrd_inv(1,k);
          /* ~~~ */
          dCinv_dus(2,gid) += -C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(1,k)-defgrd_inv(0,k)*N_XYZ(i,n)*C_inv(1,i);
        }
      }

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    FillMatrixAndVectors_presbased(
                            gp,
                            shapefct,
                            N_XYZ,
                            J,
                            press,
                            bop,
                            C_inv,
                            dJ_dus,
                            dCinv_dus,
                            stiffmatrix,
                            force);
  }//end of gaussloop
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectors(
    const int &                                     gp,
    const LINALG::Matrix<numnod_,1>&                shapefct,
    const LINALG::Matrix<numdim_,numnod_>&          N_XYZ,
    const double&                                   J,
    const double&                                   press,
    const double&                                   porosity,
    const LINALG::Matrix<numdim_,1>&                velint,
    const LINALG::Matrix<numdim_,1>&                fvelint,
    const LINALG::Matrix<numdim_,numdim_>&          fvelder,
    const LINALG::Matrix<numdim_,numdim_>&          defgrd_inv,
    const LINALG::Matrix<numstr_,numdof_>&          bop,
    const LINALG::Matrix<numdim_,numdim_>&          C_inv,
    const LINALG::Matrix<numdim_,1>&                Finvgradp,
    const LINALG::Matrix<1,numdof_>&                dphi_dus,
    const LINALG::Matrix<1,numdof_>&                dJ_dus,
    const LINALG::Matrix<numstr_,numdof_>&          dCinv_dus,
    const LINALG::Matrix<numdim_,numdof_>&          dFinvdus_gradp,
    const LINALG::Matrix<numdim_*numdim_,numdof_>&  dFinvTdus,
    LINALG::Matrix<numdof_,numdof_>&                erea_v,
    LINALG::Matrix<numdof_, numdof_>*               stiffmatrix,
    LINALG::Matrix<numdof_,1>*                      force,
    LINALG::Matrix<numstr_,1>&                      fstress)
{
  //const double reacoeff = fluidmat_->ComputeReactionCoeff();

  static LINALG::Matrix<numdim_,numdim_> matreatensor(true);
  static LINALG::Matrix<numdim_,numdim_> reatensor(true);
  static LINALG::Matrix<numdim_,numdim_> linreac_dphi(true);
  static LINALG::Matrix<numdim_,numdim_> linreac_dJ(true);
  static LINALG::Matrix<numdim_,1> reafvel(true);
  static LINALG::Matrix<numdim_,1> reavel(true);
  {
    static LINALG::Matrix<numdim_,numdim_> temp(false);
    fluidmat_->ComputeReactionTensor(matreatensor,J,porosity);
    fluidmat_->ComputeLinMatReactionTensor(linreac_dphi,linreac_dJ,J,porosity);
    temp.Multiply(1.0,matreatensor,defgrd_inv);
    reatensor.MultiplyTN(defgrd_inv,temp);
    reavel.Multiply(reatensor,velint);
    reafvel.Multiply(reatensor,fvelint);
  }

  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;

  {
    for (int k=0; k<numnod_; k++)
    {
      const int fk = numdim_*k;
      const double fac = detJ_w* shapefct(k);
      const double v = fac * porosity * porosity* J * J;

      for(int j=0; j<numdim_; j++)
      {
        /*-------structure- velocity coupling:  RHS
         "darcy-terms"
         - reacoeff * J^2 *  phi^2 *  v^f
         */
        (*force)(fk+j) += -v * reafvel(j);

        /* "reactive darcy-terms"
         reacoeff * J^2 *  phi^2 *  v^s
         */
        (*force)(fk+j) += v * reavel(j);

        /*-------structure- fluid pressure coupling: RHS
         *                        "pressure gradient terms"
         - J *  F^-T * Grad(p) * phi
         */
        (*force)(fk+j) += fac * J * Finvgradp(j) * ( - porosity);

        for(int i=0; i<numnod_; i++)
        {
          const int fi = numdim_*i;

          for (int l=0; l<numdim_; l++)
          {
            /* additional "reactive darcy-term"
             detJ * w(gp) * ( J^2 * reacoeff * phi^2  ) * D(v_s)
             */
            erea_v(fk+j,fi+l) += v * reatensor(j,l) * shapefct(i);

            /* additional "pressure gradient term"
             -  detJ * w(gp) * phi *  ( dJ/d(us) * F^-T * Grad(p) - J * d(F^-T)/d(us) *Grad(p) ) * D(us)
             - detJ * w(gp) * d(phi)/d(us) * J * F^-T * Grad(p) * D(us)
             */
            (*stiffmatrix)(fk+j,fi+l) += fac * (
                                              - porosity * dJ_dus(fi+l) * Finvgradp(j)
                                              - porosity * J * dFinvdus_gradp(j, fi+l)
                                              - dphi_dus(fi+l) * J * Finvgradp(j)
                                            )
            ;

            /* additional "reactive darcy-term"
               detJ * w(gp) * 2 * ( dJ/d(us) * vs * reacoeff * phi^2 + J * reacoeff * phi * d(phi)/d(us) * vs ) * D(us)
             - detJ * w(gp) *  2 * ( J * dJ/d(us) * v^f * reacoeff * phi^2 + J * reacoeff * phi * d(phi)/d(us) * v^f ) * D(us)
             */
            (*stiffmatrix)(fk+j,fi+l) += fac * J * porosity *  2.0 * ( reavel(j) - reafvel(j) ) *
                                            ( porosity * dJ_dus(fi+l) + J * dphi_dus(fi+l) );

              for (int m=0; m<numdim_; ++m)
                for (int n=0; n<numdim_; ++n)
                  for (int p=0; p<numdim_; ++p)
                    (*stiffmatrix)(fk+j,fi+l) += v * ( velint(p) - fvelint(p) ) * (
                                                    dFinvTdus(j*numdim_+m,fi+l) * matreatensor(m,n) * defgrd_inv(n,p)
                                                  + defgrd_inv(m,j) * matreatensor(m,n) * dFinvTdus(p*numdim_+n,fi+l)
                                                  );
            //check if derivatives of reaction tensor are zero --> significant speed up
            if (fluidmat_->PermeabilityFunction() != MAT::PAR::const_)
            {
              for (int m=0; m<numdim_; ++m)
                for (int n=0; n<numdim_; ++n)
                  for (int p=0; p<numdim_; ++p)
                    (*stiffmatrix)(fk+j,fi+l) += v * ( velint(p) - fvelint(p) ) * (
                                                  + defgrd_inv(m,j) * (
                                                  linreac_dphi(m,n) * dphi_dus(fi+l) + linreac_dJ(m,n) * dJ_dus(fi+l)
                                                  ) * defgrd_inv(n,p)
                                                  );
            }
          }
        }
      }
    }
  }

  //inverse Right Cauchy-Green tensor as vector
  static LINALG::Matrix<numstr_,1> C_inv_vec;
  C_inv_vec(0) = C_inv(0,0);
  C_inv_vec(1) = C_inv(1,1);
  C_inv_vec(2) = C_inv(0,1);

  //B^T . C^-1
  static LINALG::Matrix<numdof_,1> cinvb(true);
  cinvb.MultiplyTN(bop,C_inv_vec);

  const double fac1 = -detJ_w * press;
  const double fac2= fac1 * J;

  // update internal force vector
  if (force != NULL )
  {
    // additional fluid stress- stiffness term RHS -(B^T .  C^-1  * J * p^f * detJ * w(gp))
    force->Update(fac2,cinvb,1.0);
  }  // if (force != NULL )

  // update stiffness matrix
  if (stiffmatrix != NULL)
  {
    static LINALG::Matrix<numdof_,numdof_> tmp;

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * p^f * detJ * w(gp))
    tmp.Multiply(fac1,cinvb,dJ_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * p^f * detJ * w(gp))
    tmp.MultiplyTN(fac2,bop,dCinv_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    LINALG::Matrix<numstr_,1> sfac(C_inv_vec); // auxiliary integrated stress

    //scale and add viscous stress
    sfac.Update(detJ_w,fstress,fac2); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]

    std::vector<double> SmB_L(2); // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<numnod_; ++inod)
    {
      SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(2) * N_XYZ(1, inod);
      SmB_L[1] = sfac(2) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod);
      for (int jnod=0; jnod<numnod_; ++jnod)
      {
        double bopstrbop = 0.0; // intermediate value
        for (int idim=0; idim<numdim_; ++idim)
          bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
        (*stiffmatrix)(numdim_*inod+0,numdim_*jnod+0) += bopstrbop;
        (*stiffmatrix)(numdim_*inod+1,numdim_*jnod+1) += bopstrbop;
      }
    } // end of integrate `geometric' stiffness******************************
  }

}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectors_presbased(
    const int &                                     gp,
    const LINALG::Matrix<numnod_,1>&                shapefct,
    const LINALG::Matrix<numdim_,numnod_>&          N_XYZ,
    const double&                                   J,
    const double&                                   press,
    const LINALG::Matrix<numstr_,numdof_>&          bop,
    const LINALG::Matrix<numdim_,numdim_>&          C_inv,
    const LINALG::Matrix<1,numdof_>&                dJ_dus,
    const LINALG::Matrix<numstr_,numdof_>&          dCinv_dus,
    LINALG::Matrix<numdof_, numdof_>*               stiffmatrix,
    LINALG::Matrix<numdof_,1>*                      force)
{

  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;

  //inverse Right Cauchy-Green tensor as vector
  static LINALG::Matrix<numstr_,1> C_inv_vec;
  C_inv_vec(0) = C_inv(0,0);
  C_inv_vec(1) = C_inv(1,1);
  C_inv_vec(2) = C_inv(0,1);

  //B^T . C^-1
  static LINALG::Matrix<numdof_,1> cinvb(true);
  cinvb.MultiplyTN(bop,C_inv_vec);

  const double fac1 = -detJ_w * press;
  const double fac2= fac1 * J;

  // update internal force vector
  if (force != NULL )
  {
    // additional fluid stress- stiffness term RHS -(B^T .  C^-1  * J * p^f * detJ * w(gp))
    force->Update(fac2,cinvb,1.0);
  }  // if (force != NULL )

  // update stiffness matrix
  if (stiffmatrix != NULL)
  {
    static LINALG::Matrix<numdof_,numdof_> tmp;

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * p^f * detJ * w(gp))
    tmp.Multiply(fac1,cinvb,dJ_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * p^f * detJ * w(gp))
    tmp.MultiplyTN(fac2,bop,dCinv_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    LINALG::Matrix<numstr_,1> sfac(C_inv_vec); // auxiliary integrated stress
    sfac.Scale(fac2);

    std::vector<double> SmB_L(2); // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<numnod_; ++inod)
    {
      SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(2) * N_XYZ(1, inod);
      SmB_L[1] = sfac(2) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod);
      for (int jnod=0; jnod<numnod_; ++jnod)
      {
        double bopstrbop = 0.0; // intermediate value
        for (int idim=0; idim<numdim_; ++idim)
          bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
        (*stiffmatrix)(numdim_*inod+0,numdim_*jnod+0) += bopstrbop;
        (*stiffmatrix)(numdim_*inod+1,numdim_*jnod+1) += bopstrbop;
      }
    } // end of integrate `geometric' stiffness******************************
  }

}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectorsBrinkman(
    const int &                                     gp,
    const double&                                   J,
    const double&                                   porosity,
    const LINALG::Matrix<numdim_,numdim_>&          fvelder,
    const LINALG::Matrix<numdim_,numdim_>&          defgrd_inv,
    const LINALG::Matrix<numstr_,numdof_>&          bop,
    const LINALG::Matrix<numdim_,numdim_>&          C_inv,
    const LINALG::Matrix<1,numdof_>&                dphi_dus,
    const LINALG::Matrix<1,numdof_>&                dJ_dus,
    const LINALG::Matrix<numstr_,numdof_>&          dCinv_dus,
    const LINALG::Matrix<numdim_*numdim_,numdof_>&  dFinvTdus,
    LINALG::Matrix<numdof_, numdof_>*               stiffmatrix,
    LINALG::Matrix<numdof_,1>*                      force,
    LINALG::Matrix<numstr_,1>&                      fstress)
{
  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;

  const double visc = fluidmat_->Viscosity();
  LINALG::Matrix<numdim_,numdim_> CinvFvel;
  LINALG::Matrix<numdim_,numdim_> tmp;
  CinvFvel.Multiply(C_inv,fvelder);
  tmp.MultiplyNT(CinvFvel,defgrd_inv);
  LINALG::Matrix<numdim_,numdim_> tmp2(tmp);
  tmp.UpdateT(1.0,tmp2,1.0);

  fstress(0) = tmp(0,0);
  fstress(1) = tmp(1,1);
  fstress(2) = tmp(0,1);

  fstress.Scale(detJ_w * visc * J * porosity);

  //B^T . C^-1
  LINALG::Matrix<numdof_,1> fstressb(true);
  fstressb.MultiplyTN(bop,fstress);

  if (force != NULL )
    force->Update(1.0,fstressb,1.0);

  //evaluate viscous terms (for darcy-brinkman flow only)
  if (stiffmatrix != NULL)
  {
    LINALG::Matrix<numdim_,numdim_> tmp4;
    tmp4.MultiplyNT(fvelder,defgrd_inv);

    double fac = detJ_w * visc;

    LINALG::Matrix<numstr_,numdof_> fstress_dus (true);
    for (int n=0; n<numnod_; ++n)
      for (int k=0; k<numdim_; ++k)
      {
        const int gid = n*numdim_+k;

        fstress_dus(0,gid) += 2*( dCinv_dus(0,gid)*tmp4(0,0) + dCinv_dus(2,gid)*tmp4(1,0) );
        fstress_dus(1,gid) += 2*( dCinv_dus(2,gid)*tmp4(0,1) + dCinv_dus(1,gid)*tmp4(1,1) );
        /* ~~~ */
        fstress_dus(2,gid) += + dCinv_dus(0,gid)*tmp4(0,1) + dCinv_dus(2,gid)*tmp4(1,1)
                              + dCinv_dus(2,gid)*tmp4(0,0) + dCinv_dus(1,gid)*tmp4(1,0);

        for(int j=0; j<numdim_; j++)
        {
          fstress_dus(0,gid) += 2*CinvFvel(0,j) * dFinvTdus(j*numdim_  ,gid);
          fstress_dus(1,gid) += 2*CinvFvel(1,j) * dFinvTdus(j*numdim_+1,gid);
          /* ~~~ */
          fstress_dus(2,gid) += + CinvFvel(0,j) * dFinvTdus(j*numdim_+1,gid)
                                + CinvFvel(1,j) * dFinvTdus(j*numdim_  ,gid);
        }
      }

    LINALG::Matrix<numdof_,numdof_> tmp;

    // additional viscous fluid stress- stiffness term (B^T . fstress . dJ/d(us) * porosity * detJ * w(gp))
    tmp.Multiply(fac*porosity,fstressb,dJ_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // additional fluid stress- stiffness term (B^T .  d\phi/d(us) . fstress  * J * w(gp))
    tmp.Multiply(fac*J,fstressb,dphi_dus);
    stiffmatrix->Update(1.0,tmp,1.0);

    // additional fluid stress- stiffness term (B^T .  phi . dfstress/d(us)  * J * w(gp))
    tmp.MultiplyTN(detJ_w * visc * J * porosity,bop,fstress_dus);
    stiffmatrix->Update(1.0,tmp,1.0);
  }
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (private) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::coupling_poroelast(
    std::vector<int>&                                 lm,            // location matrix
    LINALG::Matrix<numdim_, numnod_>&                 disp,          // current displacements
    LINALG::Matrix<numdim_, numnod_>&                 vel,           // current velocities
    LINALG::Matrix<numdim_, numnod_> &                evelnp,        //current fluid velocity
    LINALG::Matrix<numnod_, 1> &                      epreaf,        //current fluid pressure
    LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix,   // element stiffness matrix
    LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* reamatrix,     // element reactive matrix
    LINALG::Matrix<numdof_, 1>*                       force,         // element internal force vector
    Teuchos::ParameterList&                           params)        // algorithmic parameters e.g. time
{

  GetMaterials();

  //=======================================================================

  // update element geometry
  LINALG::Matrix<numdim_,numnod_> xrefe; // material coord. of element
  LINALG::Matrix<numdim_,numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  if(stiffmatrix!=NULL)
    GaussPointLoopOD(  params,
                       xrefe,
                       xcurr,
                       disp,
                       vel,
                       evelnp,
                       epreaf,
                       NULL,
                       *stiffmatrix);

  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (private) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::coupling_poroelast_presbased(
    std::vector<int>&                                 lm,            // location matrix
    LINALG::Matrix<numdim_, numnod_>&                 disp,          // current displacements
    Epetra_SerialDenseMatrix&                         phi,           // current primary variable for poro-multiphase flow
    Epetra_SerialDenseMatrix&                         couplmat,   // element stiffness matrix
    Teuchos::ParameterList&                           params)        // algorithmic parameters e.g. time
{

  GetMaterials_presbased();

  //=======================================================================

  // update element geometry
  LINALG::Matrix<numdim_,numnod_> xrefe; // material coord. of element
  LINALG::Matrix<numdim_,numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/

    GaussPointLoopOD_presbased(  params,
                       xrefe,
                       xcurr,
                       disp,
                       phi,
                       couplmat);

  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::GaussPointLoopOD(
                                    Teuchos::ParameterList&                 params,
                                    const LINALG::Matrix<numdim_,numnod_>&  xrefe,
                                    const LINALG::Matrix<numdim_,numnod_>&  xcurr,
                                    const LINALG::Matrix<numdim_,numnod_>&  nodaldisp,
                                    const LINALG::Matrix<numdim_,numnod_>&  nodalvel,
                                    const LINALG::Matrix<numdim_,numnod_>&  evelnp,
                                    const LINALG::Matrix<numnod_,1> &       epreaf,
                                    const LINALG::Matrix<numnod_, 1>*       porosity_dof,
                                    LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl
                                        )
{

  /*--------------------------------- get node weights for nurbs elements */
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnod_; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights_(inode) = cp->W();
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<numdim_,numnod_> N_XYZ;       //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<numdim_,numdim_> defgrd(true); //  deformation gradiant evaluated at gauss point
  LINALG::Matrix<numnod_,1> shapefct;           //  shape functions evalulated at gauss point
  LINALG::Matrix<numdim_,numnod_> deriv(true);  //  first derivatives at gausspoint w.r.t. r,s,t
  //LINALG::Matrix<numdim_,1> xsi;
  for (int gp=0; gp<numgpt_; ++gp)
  {
    //evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd,N_XYZ,xcurr);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChange(
        J,
        volchange,
        defgrd,
        N_XYZ,
        nodaldisp);

    // non-linear B-operator (may so be called, meaning
    LINALG::Matrix<numstr_,numdof_> bop;
    ComputeBOperator(bop,defgrd,N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<numdim_,numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<numdim_,numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // inverse deformation gradient F^-1
    LINALG::Matrix<numdim_,numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    LINALG::Matrix<numdim_,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    //--------------------- get fluid velocity at integration point
    LINALG::Matrix<numdim_,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<numdim_,numdim_>              fvelder;
    fvelder.MultiplyNT(evelnp,N_XYZ);

    //! ----------------structure displacement and velocity at integration point
    LINALG::Matrix<numdim_,1> velint(true);
    for(int i=0; i<numnod_; i++)
      for(int j=0; j<numdim_; j++)
        velint(j) += nodalvel(j,i) * shapefct(i);

    //**************************************************+auxilary variables for computing the porosity and linearization
    double dphi_dp=0.0;
    double porosity=0.0;

    ComputePorosityAndLinearizationOD(params,
                                      press,
                                      volchange,
                                      gp,
                                      shapefct,
                                      porosity_dof,
                                      porosity,
                                      dphi_dp);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    FillMatrixAndVectorsOD(
                              gp,
                              shapefct,
                              N_XYZ,
                              J,
                              porosity,
                              dphi_dp,
                              velint,
                              fvelint,
                              defgrd_inv,
                              Gradp,
                              bop,
                              C_inv,
                              ecoupl);

    if(fluidmat_->Type() == MAT::PAR::darcy_brinkman )
    {
      FillMatrixAndVectorsBrinkmanOD(
                                      gp,
                                      shapefct,
                                      N_XYZ,
                                      J,
                                      porosity,
                                      dphi_dp,
                                      fvelder,
                                      defgrd_inv,
                                      bop,
                                      C_inv,
                                      ecoupl);
    }//darcy-brinkman
    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::GaussPointLoopOD_presbased(
                                    Teuchos::ParameterList&                 params,
                                    const LINALG::Matrix<numdim_,numnod_>&  xrefe,
                                    const LINALG::Matrix<numdim_,numnod_>&  xcurr,
                                    const LINALG::Matrix<numdim_,numnod_>&  nodaldisp,
                                    Epetra_SerialDenseMatrix&               phi,
                                    Epetra_SerialDenseMatrix&               couplmat
                                        )
{

  /*--------------------------------- get node weights for nurbs elements */
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnod_; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights_(inode) = cp->W();
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<numdim_,numnod_> N_XYZ;       //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<numdim_,numdim_> defgrd(true); //  deformation gradiant evaluated at gauss point
  LINALG::Matrix<numnod_,1> shapefct;           //  shape functions evalulated at gauss point
  LINALG::Matrix<numdim_,numnod_> deriv(true);  //  first derivatives at gausspoint w.r.t. r,s,t

  const int numphases = phi.RowDim();
  Epetra_SerialDenseMatrix solpressderiv(numnod_,numphases);
  // compute derivative of solid pressure w.r.t primary variable phi at node
  ComputeSolPressureDeriv(phi,numphases,solpressderiv);

  for (int gp=0; gp<numgpt_; ++gp)
  {
    //evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd,N_XYZ,xcurr);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChange(
        J,
        volchange,
        defgrd,
        N_XYZ,
        nodaldisp);

    // non-linear B-operator (may so be called, meaning
    LINALG::Matrix<numstr_,numdof_> bop;
    ComputeBOperator(bop,defgrd,N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<numdim_,numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<numdim_,numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    FillMatrixAndVectorsOD_presbased(
                              gp,
                              shapefct,
                              N_XYZ,
                              J,
                              bop,
                              C_inv,
                              solpressderiv,
                              couplmat);

    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}


/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectorsOD(
    const int &                             gp,
    const LINALG::Matrix<numnod_,1>&        shapefct,
    const LINALG::Matrix<numdim_,numnod_>&  N_XYZ,
    const double&                           J,
    const double&                           porosity,
    const double&                           dphi_dp,
    const LINALG::Matrix<numdim_,1>&        velint,
    const LINALG::Matrix<numdim_,1>&        fvelint,
    const LINALG::Matrix<numdim_,numdim_>&  defgrd_inv,
    const LINALG::Matrix<numdim_,1>&        Gradp,
    const LINALG::Matrix<numstr_,numdof_>&  bop,
    const LINALG::Matrix<numdim_,numdim_>&  C_inv,
    LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl)
{
  LINALG::Matrix<numdim_,numdim_> matreatensor(true);
  LINALG::Matrix<numdim_,numdim_> reatensor(true);
  LINALG::Matrix<numdim_,numdim_> linreac_dphi(true);
  LINALG::Matrix<numdim_,numdim_> linreac_dJ(true);
  LINALG::Matrix<numdim_,1> reafvel(true);
  LINALG::Matrix<numdim_,1> reavel(true);
  {
    LINALG::Matrix<numdim_,numdim_> temp(true);
    fluidmat_->ComputeReactionTensor(matreatensor,J,porosity);
    fluidmat_->ComputeLinMatReactionTensor(linreac_dphi,linreac_dJ,J,porosity);
    temp.Multiply(1.0,matreatensor,defgrd_inv);
    reatensor.MultiplyTN(defgrd_inv,temp);
    reavel.Multiply(reatensor,velint);
    reafvel.Multiply(reatensor,fvelint);
  }

  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<numstr_,1> C_inv_vec;
  C_inv_vec(0) = C_inv(0,0);
  C_inv_vec(1) = C_inv(1,1);
  C_inv_vec(2) = C_inv(0,1);

  //B^T . C^-1
  LINALG::Matrix<numdof_,1> cinvb(true);
  cinvb.MultiplyTN(bop,C_inv_vec);

  //F^-T * Grad p
  LINALG::Matrix<numdim_,1> Finvgradp;
  Finvgradp.MultiplyTN(defgrd_inv, Gradp);

  //F^-T * N_XYZ
  LINALG::Matrix<numdim_,numnod_> FinvNXYZ;
  FinvNXYZ.MultiplyTN(defgrd_inv, N_XYZ);

  {
    for (int i=0; i<numnod_; i++)
    {
      const int fi = numdim_*i;
      const double fac = detJ_w* shapefct(i);

      for(int j=0; j<numdim_; j++)
      {
        for(int k=0; k<numnod_; k++)
        {
          const int fk = (numdim_+1)*k;
          const int fk_press = fk+numdim_;

          /*-------structure- fluid pressure coupling: "stress terms" + "pressure gradient terms"
           -B^T . ( -1*J*C^-1 ) * Dp
           - J * F^-T * Grad(p) * dphi/dp * Dp - J * F^-T * d(Grad((p))/(dp) * phi * Dp
           */
          ecoupl(fi+j, fk_press) +=   detJ_w * cinvb(fi+j) * ( -1.0) * J * shapefct(k)
                                    - fac * J * (   dphi_dp  * Finvgradp(j) * shapefct(k)
                                                  + porosity * FinvNXYZ(j,k)
                                                )
                                            ;

          /*-------structure- fluid pressure coupling:  "darcy-terms" + "reactive darcy-terms"
           - 2 * reacoeff * J * v^f * phi * d(phi)/dp  Dp
           + 2 * reacoeff * J * v^s * phi * d(phi)/dp  Dp
           + J * J * phi * phi * defgrd_^-T * d(mat_reacoeff)/d(phi) * defgrd_^-1 * (v^s-v^f) * d(phi)/dp Dp
           */
          const double tmp = fac * J * J * 2 * porosity * dphi_dp * shapefct(k);
          ecoupl(fi+j, fk_press ) += -tmp * reafvel(j);

          ecoupl(fi+j, fk_press ) += tmp * reavel(j);

          //check if derivatives of reaction tensor are zero --> significant speed up
          if (fluidmat_->PermeabilityFunction() != MAT::PAR::const_)
          {
            const double tmp2 = 0.5 * tmp * porosity;
            for (int m=0; m<numdim_; ++m)
            {
              for (int n=0; n<numdim_; ++n)
              {
                for (int p=0; p<numdim_; ++p)
                {
                  ecoupl(fi+j, fk_press ) +=  tmp2 * defgrd_inv(m,j) * linreac_dphi(m,n) * defgrd_inv(n,p) *
                                              ( velint(p) - fvelint(p) );
                }
              }
            }
          }

          /*-------structure- fluid velocity coupling:  "darcy-terms"
           -reacoeff * J * J *  phi^2 *  Dv^f
           */
          const double v = fac * J * J * porosity * porosity;
          for(int l=0; l<numdim_; l++)
            ecoupl(fi+j, fk+l) += -v * reatensor(j,l) * shapefct(k);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectorsOD_presbased(
    const int &                             gp,
    const LINALG::Matrix<numnod_,1>&        shapefct,
    const LINALG::Matrix<numdim_,numnod_>&  N_XYZ,
    const double&                           J,
    const LINALG::Matrix<numstr_,numdof_>&  bop,
    const LINALG::Matrix<numdim_,numdim_>&  C_inv,
    const Epetra_SerialDenseMatrix&         solpressderiv,
    Epetra_SerialDenseMatrix&               couplmat)
{

  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<numstr_,1> C_inv_vec;
  C_inv_vec(0) = C_inv(0,0);
  C_inv_vec(1) = C_inv(1,1);
  C_inv_vec(2) = C_inv(0,1);

  //B^T . C^-1
  LINALG::Matrix<numdof_,1> cinvb(true);
  cinvb.MultiplyTN(bop,C_inv_vec);

  const int numphases = solpressderiv.ColDim();

  {
    for (int i=0; i<numnod_; i++)
    {
      const int fi = numdim_*i;

      for(int j=0; j<numdim_; j++)
      {
        for(int k=0; k<numnod_; k++)
        {
          for(int iphase = 0; iphase < numphases; iphase++)
          {
            int fk_press = k*numphases+iphase;

            /*-------structure- fluid pressure coupling: "stress term"
             -B^T . ( -1*J*C^-1 ) * Dp
             */
            couplmat(fi+j, fk_press) +=   detJ_w * cinvb(fi+j) * ( -1.0) * J * shapefct(k) * solpressderiv(k,iphase)
                                                ;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::FillMatrixAndVectorsBrinkmanOD(
    const int &                             gp,
    const LINALG::Matrix<numnod_,1>&        shapefct,
    const LINALG::Matrix<numdim_,numnod_>&  N_XYZ,
    const double&                           J,
    const double&                           porosity,
    const double&                           dphi_dp,
    const LINALG::Matrix<numdim_,numdim_>&  fvelder,
    const LINALG::Matrix<numdim_,numdim_>&  defgrd_inv,
    const LINALG::Matrix<numstr_,numdof_>&  bop,
    const LINALG::Matrix<numdim_,numdim_>&  C_inv,
    LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl)
{
  const double detJ_w = detJ_[gp]*intpoints_.Weight(gp)*thickness_;
  const double visc = fluidmat_->Viscosity();

  LINALG::Matrix<numstr_,1> fstress;

  LINALG::Matrix<numdim_,numdim_> CinvFvel;
  LINALG::Matrix<numdim_,numdim_> tmp;
  CinvFvel.Multiply(C_inv,fvelder);
  tmp.MultiplyNT(CinvFvel,defgrd_inv);
  LINALG::Matrix<numdim_,numdim_> tmp2(tmp);
  tmp.UpdateT(1.0,tmp2,1.0);

  fstress(0) = tmp(0,0);
  fstress(1) = tmp(1,1);
  fstress(2) = tmp(0,1);

  //B^T . \sigma
  LINALG::Matrix<numdof_,1> fstressb;
  fstressb.MultiplyTN(bop,fstress);
  LINALG::Matrix<numdim_,numnod_> N_XYZ_Finv;
  N_XYZ_Finv.Multiply(defgrd_inv,N_XYZ);

  //dfstress/dv^f
  LINALG::Matrix<numstr_,numdof_> dfstressb_dv;
  for (int i=0; i<numnod_; i++)
  {
    const int fi = numdim_*i;
    for(int j=0; j<numdim_; j++)
    {
      int k = fi+j;
      dfstressb_dv(0,k) = 2 * N_XYZ_Finv(0,i) * C_inv(0,j);
      dfstressb_dv(1,k) = 2 * N_XYZ_Finv(1,i) * C_inv(1,j);
      //**********************************
      dfstressb_dv(2,k) = N_XYZ_Finv(0,i) * C_inv(1,j) + N_XYZ_Finv(1,i) * C_inv(0,j);
    }
  }

  //B^T . dfstress/dv^f
  LINALG::Matrix<numdof_,numdof_> dfstressb_dv_bop(true);
  dfstressb_dv_bop.MultiplyTN(bop,dfstressb_dv);

  for (int i=0; i<numnod_; i++)
  {
    const int fi = numdim_*i;

    for(int j=0; j<numdim_; j++)
    {
      for(int k=0; k<numnod_; k++)
      {
        const int fk_sub = numdim_*k;
        const int fk = (numdim_+1)*k;
        const int fk_press = fk+numdim_;

        /*-------structure- fluid pressure coupling: "darcy-brinkman stress terms"
         B^T . ( \mu*J - d(phi)/(dp) * fstress ) * Dp
         */
        ecoupl(fi+j, fk_press) += detJ_w * fstressb(fi+j) * dphi_dp * visc * J * shapefct(k);
        for(int l=0; l<numdim_; l++)
        {
          /*-------structure- fluid velocity coupling: "darcy-brinkman stress terms"
           B^T . ( \mu*J - phi * dfstress/dv^f ) * Dp
           */
          ecoupl(fi+j, fk+l) += detJ_w * visc * J * porosity * dfstressb_dv_bop(fi+j, fk_sub+l);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::couplstress_poroelast(
    LINALG::Matrix<numdim_, numnod_>&  disp,         // current displacements
    LINALG::Matrix<numdim_, numnod_> & evelnp,       // current fluid velocities
    LINALG::Matrix<numnod_, 1> &       epreaf,       // current fluid pressure
    Epetra_SerialDenseMatrix*          elestress,    // stresses at GP
    Epetra_SerialDenseMatrix*          elestrain,    // strains at GP
    Teuchos::ParameterList&            params,       // algorithmic parameters e.g. time
    const INPAR::STR::StressType       iostress      // stress output option
    )
{
  // update element geometry
  LINALG::Matrix<numdim_,numnod_> xrefe; // material coord. of element
  LINALG::Matrix<numdim_,numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }
  LINALG::Matrix<numnod_,1> shapefct;
  LINALG::Matrix<numdim_,numdim_> defgrd(true);
  LINALG::Matrix<numdim_,numnod_> N_XYZ;
  LINALG::Matrix<numdim_,numnod_> deriv ;

  for (int gp=0; gp<numgpt_; ++gp)
  {

    //evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChange(
        J,
        volchange,
        defgrd,
        N_XYZ,
        disp);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd,N_XYZ,xcurr);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<numdim_,numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<numdim_,numdim_> C_inv;
    C_inv.Invert(cauchygreen);

    //inverse Right Cauchy-Green tensor as vector
    LINALG::Matrix<numstr_,1> C_inv_vec;
    for(int i =0, k=0;i<numdim_; i++)
      for(int j =0;j<numdim_-i; j++,k++)
        C_inv_vec(k)=C_inv(i+j,j);

    LINALG::Matrix<Wall1::numstr_,1> couplstress(true);
    couplstress(0) = -1.0*J*press*C_inv_vec(0);
    couplstress(1) = -1.0*J*press*C_inv_vec(1);
    couplstress(2) = 0.0; // this is needed to be compatible with the implementation of the wall element
    couplstress(3) = -1.0*J*press*C_inv_vec(2);

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress==NULL) dserror("stress data not available");
      for (int i=0; i<numstr_; ++i)
        (*elestress)(gp,i) = couplstress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress==NULL) dserror("stress data not available");

      // push forward of material stress to the spatial configuration
      LINALG::Matrix<numdim_,numdim_> cauchycouplstress;
      PK2toCauchy(couplstress,defgrd,cauchycouplstress);

      (*elestress)(gp,0) = cauchycouplstress(0,0);
      (*elestress)(gp,1) = cauchycouplstress(1,1);
      (*elestress)(gp,2) = 0.0;
      (*elestress)(gp,3) = cauchycouplstress(0,1);
    }
    break;
    case INPAR::STR::stress_none:
      break;

    default:
      dserror("requested stress type not available");
      break;
    }
  }
}//couplstress_poroelast

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::InitElement()
{
  LINALG::Matrix<numdim_,numnod_> deriv ;
  LINALG::Matrix<numnod_,numdim_> xrefe;
  for (int i=0; i<numnod_; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    for(int j=0; j<numdim_; ++j)
      xrefe(i,j) = Nodes()[i]->X()[j];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  for (int gp=0; gp<numgpt_; ++gp)
  {
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim=0;idim<numdim_;idim++)
    {
       xsi_[gp](idim) = gpcoord[idim];
    }
  }

  if (distype != DRT::Element::nurbs4 and distype != DRT::Element::nurbs9)
  {
    for (int gp=0; gp<numgpt_; ++gp)
    {
      DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

      invJ_[gp].Multiply(deriv,xrefe);
      detJ_[gp] = invJ_[gp].Invert();
      if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0",detJ_[gp]);
    }
  }

  scatracoupling_=false;

  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  if(probtype == prb_poroscatra or probtype == prb_immersed_cell)
    scatracoupling_=true;

  init_=true;

  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeJacobianDeterminantVolumeChangeAndLinearizations(
    double & J,
    double & volchange,
    LINALG::Matrix<1,numdof_>& dJ_dus,
    LINALG::Matrix<1,numdof_>& dvolchange_dus,
    const LINALG::Matrix<numdim_,numdim_>& defgrd,
    const LINALG::Matrix<numdim_,numdim_>& defgrd_inv,
    const LINALG::Matrix<numdim_,numnod_>& N_XYZ,
    const LINALG::Matrix<numdim_,numnod_>& nodaldisp
    )
{
  //compute J
  J=defgrd.Determinant();
  //compute linearization of J
  ComputeLinearizationOfJacobian(dJ_dus,J,N_XYZ,defgrd_inv);

  if(kintype_==INPAR::STR::kinem_nonlinearTotLag) //total lagrange (nonlinear)
  {
    //for nonlinear kinematics the Jacobian of the deformation gradient is the volume change
    volchange=J;
    dvolchange_dus=dJ_dus;
  }
  else if(kintype_==INPAR::STR::kinem_linear) //linear kinematics
  {
    //for linear kinematics the volume change is the trace of the linearized strains

    //gradient of displacements
    static LINALG::Matrix<numdim_,numdim_> dispgrad;
    dispgrad.Clear();
    //gradient of displacements
    dispgrad.MultiplyNT(nodaldisp,N_XYZ);

    volchange=1.0;
    //volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
    for(int i=0; i<numdim_;++i)
      volchange += dispgrad(i,i);

    for(int i=0; i<numdim_;++i)
      for(int j=0; j<numnod_;++j)
        dvolchange_dus(numdim_*j+i)=N_XYZ(i,j);
  }
  else
    dserror("invalid kinematic type!");

  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeJacobianDeterminantVolumeChange(
    double & J,
    double & volchange,
    const LINALG::Matrix<numdim_,numdim_>& defgrd,
    const LINALG::Matrix<numdim_,numnod_>& N_XYZ,
    const LINALG::Matrix<numdim_,numnod_>& nodaldisp
    )
{
  //compute J
  J=defgrd.Determinant();

  if(kintype_==INPAR::STR::kinem_nonlinearTotLag) //total lagrange (nonlinear)
  {
    //for nonlinear kinematics the Jacobian of the deformation gradient is the volume change
    volchange=J;
  }
  else if(kintype_==INPAR::STR::kinem_linear) //linear kinematics
  {
    //for linear kinematics the volume change is the trace of the linearized strains

    //gradient of displacements
    static LINALG::Matrix<numdim_,numdim_> dispgrad;
    dispgrad.Clear();
    //gradient of displacements
    dispgrad.MultiplyNT(nodaldisp,N_XYZ);

    volchange=1.0;
    //volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
    for(int i=0; i<numdim_;++i)
      volchange += dispgrad(i,i);
  }
  else
    dserror("invalid kinematic type!");
}

/*-----------------------------------------------------------------------------*
 *  transform of 2. Piola Kirchhoff stresses to cauchy stresses     vuong 12/12|
 *----------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::PK2toCauchy(
  LINALG::Matrix<Wall1::numstr_,1>& stress,
  LINALG::Matrix<numdim_,numdim_>& defgrd,
  LINALG::Matrix<numdim_,numdim_>& cauchystress
  )
{
  // calculate the Jacobi-deterinant
  const double detF = (defgrd).Determinant();

  // sigma = 1/J . F . S . F^T
  LINALG::Matrix<numdim_,numdim_> pkstress;
  pkstress(0,0) = (stress)(0);
  pkstress(0,1) = (stress)(2);
  pkstress(1,0) = pkstress(0,1);
  pkstress(1,1) = (stress)(1);

  LINALG::Matrix<numdim_,numdim_> temp;
  temp.Multiply((1.0/detF),(defgrd),pkstress);
  (cauchystress).MultiplyNT(temp,(defgrd));

}  // PK2toCauchy()

/*-----------------------------------------------------------------------------*
 * compute deformation gradient                                     vuong 03/15|
 *----------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeDefGradient(
    LINALG::Matrix<numdim_,numdim_>&       defgrd,   ///<<    (i) deformation gradient at gausspoint
    const LINALG::Matrix<numdim_,numnod_>& N_XYZ,    ///<<    (i) derivatives of shape functions w.r.t. reference coordinates
    const LINALG::Matrix<numdim_,numnod_>& xcurr     ///<<    (i) current position of gausspoint
  )
{
  if(kintype_==INPAR::STR::kinem_nonlinearTotLag) //total lagrange (nonlinear)
  {
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr,N_XYZ); //  (6.17)
  }
  else if(kintype_==INPAR::STR::kinem_linear) //linear kinmatics
  {
    defgrd.Clear();
    for(int i=0;i<numdim_;i++)
      defgrd(i,i) = 1.0;
  }
  else
    dserror("invalid kinematic type!");

  return;

}  // ComputeDefGradient


/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
inline void
DRT::ELEMENTS::Wall1_Poro<distype>::ComputeBOperator(
                                                        LINALG::Matrix<numstr_,numdof_>& bop,
                                                        const LINALG::Matrix<numdim_,numdim_>& defgrd,
                                                        const LINALG::Matrix<numdim_,numnod_>& N_XYZ)
{
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
  for (int i=0; i<numnod_; ++i)
  {
    bop(0,noddof_*i+0) = defgrd(0,0)*N_XYZ(0,i);
    bop(0,noddof_*i+1) = defgrd(1,0)*N_XYZ(0,i);
    bop(1,noddof_*i+0) = defgrd(0,1)*N_XYZ(1,i);
    bop(1,noddof_*i+1) = defgrd(1,1)*N_XYZ(1,i);
    /* ~~~ */
    bop(2,noddof_*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
    bop(2,noddof_*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
  }

}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeShapeFunctionsAndDerivatives(
    const int & gp,
    LINALG::Matrix<numnod_,1>& shapefct,
    LINALG::Matrix<numdim_,numnod_>& deriv ,
    LINALG::Matrix<numdim_,numnod_>& N_XYZ)
{

  // get values of shape functions and derivatives in the gausspoint
  if(distype != DRT::Element::nurbs4
     and
     distype != DRT::Element::nurbs9)
  {
  // shape functions and their derivatives for polynomials
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);
  }
  else
  {
    // nurbs version
    DRT::NURBS::UTILS::nurbs_get_funct_deriv
    (shapefct  ,
      deriv  ,
      xsi_[gp],
      myknots_,
      weights_,
      distype );

    LINALG::Matrix<numnod_,numdim_> xrefe;
    for (int i=0; i<numnod_; ++i)
    {
      Node** nodes=Nodes();
      if(!nodes) dserror("Nodes() returned null pointer");
      xrefe(i,0) = Nodes()[i]->X()[0];
      xrefe(i,1) = Nodes()[i]->X()[1];
    }
    invJ_[gp].Multiply(deriv,xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0",detJ_[gp]);
  }

  /* get the inverse of the Jacobian matrix which looks like:
   **            [ X_,r  Y_,r  Z_,r ]^-1
   **     J^-1 = [ X_,s  Y_,s  Z_,s ]
   **            [ X_,t  Y_,t  Z_,t ]
   */

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::Wall1_Poro<distype>::ComputeJacobianDeterminant(
    const int & gp,
    const LINALG::Matrix<numdim_,numnod_>&   xcurr,
    const   LINALG::Matrix<numdim_,numnod_>& deriv)
{
  // get Jacobian matrix and determinant w.r.t. spatial configuration
  //! transposed jacobian "dx/ds"
  LINALG::Matrix<numdim_,numdim_> xjm;
  //! inverse of transposed jacobian "ds/dx"
  LINALG::Matrix<numdim_,numdim_> xji;
  xjm.MultiplyNT(deriv,xcurr);
  const double det = xji.Invert(xjm);

  // determinant of deformationgradient: det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  const double J = det/detJ_[gp];

  return J;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
inline   void
DRT::ELEMENTS::Wall1_Poro<distype>::ComputeLinearizationOfJacobian(
    LINALG::Matrix<1,numdof_>& dJ_dus,
    const double& J,
    const LINALG::Matrix<numdim_,numnod_>& N_XYZ,
    const LINALG::Matrix<numdim_,numdim_>& defgrd_inv)
{
  //--------------------------- build N_X operator (wrt material config)
  LINALG::Matrix<numdim_*numdim_,numdof_> N_X(true); // set to zero
  for (int i=0; i<numnod_; ++i)
  {
    N_X(0,numdim_*i+0) = N_XYZ(0,i);
    N_X(1,numdim_*i+1) = N_XYZ(0,i);

    N_X(2,numdim_*i+0) = N_XYZ(1,i);
    N_X(3,numdim_*i+1) = N_XYZ(1,i);
  }

  //------------------------------------ build F^-1 as vector 4x1
  LINALG::Matrix<numdim_*numdim_,1> defgrd_inv_vec;
  defgrd_inv_vec(0)=defgrd_inv(0,0);
  defgrd_inv_vec(1)=defgrd_inv(0,1);
  defgrd_inv_vec(2)=defgrd_inv(1,0);
  defgrd_inv_vec(3)=defgrd_inv(1,1);

  //------linearization of jacobi determinant detF=J w.r.t. strucuture displacement   dJ/d(us) = dJ/dF : dF/dus = J * F^-T * N,X
  dJ_dus.MultiplyTN(J,defgrd_inv_vec,N_X);
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeAuxiliaryValues(const LINALG::Matrix<numdim_,numnod_>& N_XYZ,
    const LINALG::Matrix<numdim_,numdim_>& defgrd_inv,
    const LINALG::Matrix<numdim_,numdim_>& C_inv,
    const LINALG::Matrix<numdim_,1>&  Gradp,
    LINALG::Matrix<numdim_*numdim_,numdof_>& dFinvTdus,
    LINALG::Matrix<numdim_,1>& Finvgradp,
    LINALG::Matrix<numdim_,numdof_>& dFinvdus_gradp,
    LINALG::Matrix<numstr_,numdof_>& dCinv_dus)
{
  //F^-T * Grad p
  Finvgradp.MultiplyTN(defgrd_inv, Gradp);

  if(kintype_!=INPAR::STR::kinem_linear)
  {
    //dF^-T/dus
    for (int i=0; i<numdim_; i++)
      for (int n =0; n<numnod_; n++)
        for(int j=0; j<numdim_; j++)
        {
          const int gid = numdim_ * n +j;
          for (int k=0; k<numdim_; k++)
            for(int l=0; l<numdim_; l++)
              dFinvTdus(i*numdim_+l, gid) += -defgrd_inv(l,j) * N_XYZ(k,n) * defgrd_inv(k,i);
        }

    //dF^-T/dus * Grad p
    for (int i=0; i<numdim_; i++)
      for (int n =0; n<numnod_; n++)
        for(int j=0; j<numdim_; j++)
        {
          const int gid = numdim_ * n +j;
          for(int l=0; l<numdim_; l++)
            dFinvdus_gradp(i, gid) += dFinvTdus(i*numdim_+l, gid)  * Gradp(l);
        }
  }

  for (int n=0; n<numnod_; ++n)
    for (int k=0; k<numdim_; ++k)
    {
      const int gid = n*numdim_+k;
      for (int i=0; i<numdim_; ++i)
      {
        dCinv_dus(0,gid) += -2*C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(0,k);
        dCinv_dus(1,gid) += -2*C_inv(1,i)*N_XYZ(i,n)*defgrd_inv(1,k);
        /* ~~~ */
        dCinv_dus(2,gid) += -C_inv(0,i)*N_XYZ(i,n)*defgrd_inv(1,k)-defgrd_inv(0,k)*N_XYZ(i,n)*C_inv(1,i);
      }
    }
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputePorosityAndLinearization
(   Teuchos::ParameterList&             params,
    const double&                       press,
    const double&                       J,
    const int&                          gp,
    const LINALG::Matrix<numnod_,1>&    shapfct,
    const LINALG::Matrix<numnod_,1>*    myporosity,
    const LINALG::Matrix<1,numdof_>&    dJ_dus,
    double &                            porosity,
    LINALG::Matrix<1,numdof_>&          dphi_dus)
{
  double dphi_dJ=0.0;

  structmat_->ComputePorosity( params,
                              press,
                              J,
                              gp,
                              porosity,
                              NULL,
                              &dphi_dJ,
                              NULL,
                              NULL,
                              NULL        //dphi_dpp not needed
                              );

  //linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
  dphi_dus.Update( dphi_dJ , dJ_dus );

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputePorosityAndLinearizationOD
(   Teuchos::ParameterList&          params,
    const double&                    press,
    const double&                    J,
    const int&                       gp,
    const LINALG::Matrix<numnod_,1>&       shapfct,
    const LINALG::Matrix<numnod_,1>*       myporosity,
    double &                         porosity,
    double &                         dphi_dp)
{
  structmat_->ComputePorosity( params,
                              press,
                              J,
                              gp,
                              porosity,
                              &dphi_dp,
                              NULL,       //dphi_dJ not needed
                              NULL,       //dphi_dJdp not needed
                              NULL,       //dphi_dJJ not needed
                              NULL        //dphi_dpp not needed
                              );

  return;
}

/*----------------------------------------------------------------------*
 * derivative of sol. pres. at node for multiphase flow kremheller  03/17|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ComputeSolPressureDeriv
(   const Epetra_SerialDenseMatrix&  phi,
    const int&                       numphases,
    Epetra_SerialDenseMatrix&        solpressderiv
)
{
  //TODO: do we really have to do this inside the solid element
  //      maybe give the solpressurederive directly to solid element
  solpressderiv.Scale(0.0);

  for (int inode=0; inode < numnod_;inode++)
  {
    std::vector<double> phiVec(numphases);
    std::vector<double> genpress(numphases);

    std::vector<double> press(numphases);
    std::vector<double> sat(numphases);
    Epetra_SerialDenseMatrix helpderiv(numphases,numphases,true);
    Epetra_SerialDenseMatrix satderiv(numphases,numphases,true);
    Epetra_SerialDenseMatrix pressderiv(numphases,numphases,true);

    // compute phi at node
    for (int j=0; j < numphases;j++)
    {
      phiVec[j]=(phi(j,inode));
    }

    // evaluate the pressures
    fluidmultimat_->EvaluateGenPressure(genpress,phiVec);

    //! transform generalized pressures to true pressure values
    fluidmultimat_->TransformGenPresToTruePres(genpress,press);

    // explicit evaluation of saturation
    fluidmultimat_->EvaluateSaturation(sat,phiVec,press);

    // calculate the derivative of the pressure (actually first its inverse)
    fluidmultimat_->EvaluateDerivOfDofWrtPressure(pressderiv,phiVec);

    // now invert the derivatives of the dofs w.r.t. pressure to get the derivatives
    // of the pressure w.r.t. the dofs
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(pressderiv);
      int err = inverse.Invert();
      if (err != 0)
        dserror("Inversion of matrix for pressure derivative failed with error code %d.",err);
    }

    // calculate derivatives of saturation w.r.t. pressure
    fluidmultimat_->EvaluateDerivOfSaturationWrtPressure(helpderiv,phiVec);

    // chain rule: the derivative of saturation w.r.t. dof =
    // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
    satderiv.Multiply('N','N',1.0,helpderiv,pressderiv,0.0);

    // compute derivative of solid pressure w.r.t. dofs with product rule
    for(int iphase=0; iphase<numphases; iphase++)
      for(int jphase=0; jphase<numphases; jphase++)
        solpressderiv(inode,iphase) +=   pressderiv(jphase,iphase)*sat[jphase]
                                        + satderiv(jphase,iphase)*press[jphase];
  }

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::ExtractValuesFromGlobalVector(
                                    const DRT::Discretization&          discretization, ///< discretization
                                    const int&                          dofset,
                                    const std::vector<int>&             lm,             ///<
                                    LINALG::Matrix<numdim_,numnod_> *   matrixtofill,   ///< vector field
                                    LINALG::Matrix<numnod_,1> *         vectortofill,   ///< scalar field
                                    const std::string                   state          ///< state of the global vector
)
{
  //todo put on higher level
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(dofset,state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  const int numdofpernode = discretization.NumDof(dofset,Nodes()[0]);

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  if(numdofpernode==numdim_+1)
  {
    for (int inode=0; inode<numnod_; ++inode)  // number of nodes
    {
      // fill a vector field via a pointer
      if (matrixtofill != NULL)
      {
        for(int idim=0; idim<numdim_; ++idim) // number of dimensions
        {
          (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
        }  // end for(idim)
      }
      // fill a scalar field via a pointer
      if (vectortofill != NULL)
        (*vectortofill)(inode,0) = mymatrix[numdim_+(inode*numdofpernode)];
    }
  }
  else if (numdofpernode==numdim_)
  {
    for (int inode=0; inode<numnod_; ++inode)  // number of nodes
    {
      // fill a vector field via a pointer
      if (matrixtofill != NULL)
      {
        for(int idim=0; idim<numdim_; ++idim) // number of dimensions
        {
          (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
        }  // end for(idim)
      }
    }
  }
  else if (numdofpernode==1)
    for (int inode=0; inode<numnod_; ++inode)  // number of nodes
    {
      if (vectortofill != NULL)
        (*vectortofill)(inode,0) = mymatrix[inode*numdofpernode];
    }
  else
  {
    for (int inode=0; inode<numnod_; ++inode)  // number of nodes
    {
      if (vectortofill != NULL)
        (*vectortofill)(inode,0) = mymatrix[inode*numdofpernode];
    }
  }
   // dserror("wrong number of dofs for extract");
}

/*----------------------------------------------------------------------*
 *                                                            vuong 12/12|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad9>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::nurbs9>;
