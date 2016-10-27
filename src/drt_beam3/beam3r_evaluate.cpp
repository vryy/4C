/*----------------------------------------------------------------------------*/
/*!
\file beam3r_evaluate.cpp

\brief evaluation methods for 3D nonlinear Reissner beam element

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

#include "beam3r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../headers/FAD_utils.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_structure_new/str_model_evaluator_data.H"

#include <iostream>
#include <iomanip>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1, //nonlinear stiffness matrix
    Epetra_SerialDenseMatrix& elemat2, //nonlinear mass matrix
    Epetra_SerialDenseVector& elevec1, //nonlinear internal (elastic) forces
    Epetra_SerialDenseVector& elevec2, //nonlinear inertia forces
    Epetra_SerialDenseVector& elevec3)
{
  // Set structure params interface pointer
  SetParamsInterfacePtr(params);
  // Set statmech params interface pointer
  if (IsParamsInterface())
    SetStatMechParamsInterfacePtr();

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if (action == "calc_none")                    dserror("No action supplied");
    else if (action=="calc_struct_linstiff")                               act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")                               act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce")                          act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")                           act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")                           act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass")                          act = ELEMENTS::struct_calc_nlnstifflmass; //with lumped mass matrix
    else if (action=="calc_struct_stress")                                 act = ELEMENTS::struct_calc_stress;
    else if (action=="calc_struct_eleload")                                act = ELEMENTS::struct_calc_eleload;
    else if (action=="calc_struct_fsiload")                                act = ELEMENTS::struct_calc_fsiload;
    else if (action=="calc_struct_update_istep")                           act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")                            act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_ptcstiff")                               act = ELEMENTS::struct_calc_ptcstiff;
    else if (action=="calc_struct_energy")                                 act = ELEMENTS::struct_calc_energy;
    else dserror("Unknown type of action for Beam3r");
  }

  // nnodetriad: number of nodes used for interpolation of triad field
  const int nnodetriad = NumNode();

  switch(act)
  {
    case ELEMENTS::struct_calc_ptcstiff:
    {
      switch(nnodetriad)
      {
        case 2:EvaluatePTC<2>(params, elemat1); break;
        case 3:EvaluatePTC<3>(params, elemat1); break;
        case 4:EvaluatePTC<4>(params, elemat1); break;
        case 5:EvaluatePTC<5>(params, elemat1); break;
        default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
      }
      break;
    }

    case ELEMENTS::struct_calc_linstiff:
    {
      // only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
      break;
    }

    case ELEMENTS::struct_calc_energy:
    {
      if(elevec1 != Teuchos::null)
      {
        if(elevec1.M()!=1)
          dserror("energy vector of invalid size %i, expected row dimension 1 (total elastic energy of element)!", elevec1.M());
        elevec1(0)=Eint_;
      }
      break;
    }

    // nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    case ELEMENTS::struct_calc_internalinertiaforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom

      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      if (act == ELEMENTS::struct_calc_nlnstiffmass)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,1>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            }
            else
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,2>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            }
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<3,3,1>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<3,2,2>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<4,4,1>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<4,2,2>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<5,5,1>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<5,2,2>(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
        }
      }

      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        // TODO there is a method 'Beam3r::lumpmass'; check generality and functionality and enable action here
        dserror("Lumped mass matrix not implemented for beam3r elements so far!");
      }

      else if (act == ELEMENTS::struct_calc_nlnstiff)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,1>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            }
            else
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,2>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            }
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<3,3,1>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<3,2,2>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<4,4,1>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<4,2,2>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<5,5,1>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<5,2,2>(params,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }

      }
      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,1>(params,mydisp,NULL,NULL,&elevec1,NULL);
            }
            else
            {
              CalcInternalAndInertiaForcesAndStiff<2,2,2>(params,mydisp,NULL,NULL,&elevec1,NULL);
            }
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<3,3,1>(params,mydisp,NULL,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<3,2,2>(params,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<4,4,1>(params,mydisp,NULL,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<4,2,2>(params,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<5,5,1>(params,mydisp,NULL,NULL,&elevec1,NULL);
            else
              CalcInternalAndInertiaForcesAndStiff<5,2,2>(params,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }

      else if (act == ELEMENTS::struct_calc_internalinertiaforce)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<2,2,1>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<2,2,2>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<3,3,1>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<3,2,2>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<4,4,1>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<4,2,2>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcInternalAndInertiaForcesAndStiff<5,5,1>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            else
              CalcInternalAndInertiaForcesAndStiff<5,2,2>(params,mydisp,NULL,NULL,&elevec1,&elevec2);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }

      break;
    }

    case ELEMENTS::struct_calc_update_istep:
    {
      /* the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam at the end of the time step has to be stored*/
      Qconvnode_ = Qnewnode_;
      QconvGPmass_ = QnewGPmass_;
      wconvGPmass_ = wnewGPmass_;
      aconvGPmass_ = anewGPmass_;
      amodconvGPmass_ = amodnewGPmass_;
      rttconvGPmass_ = rttnewGPmass_;
      rttmodconvGPmass_ = rttmodnewGPmass_;
      rtconvGPmass_ = rtnewGPmass_;
      rconvGPmass_ = rnewGPmass_;
      dispthetaconvnode_ = dispthetanewnode_;
      QconvGPdampstoch_ = QnewGPdampstoch_;
      break;
    }

    case ELEMENTS::struct_calc_reset_istep:
    {
      /* the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      Qnewnode_ = Qconvnode_;
      dispthetanewnode_ = dispthetaconvnode_;
      QnewGPmass_=QconvGPmass_;
      wnewGPmass_=wconvGPmass_;
      anewGPmass_=aconvGPmass_;
      amodnewGPmass_=amodconvGPmass_;
      rttnewGPmass_=rttconvGPmass_;
      rttmodnewGPmass_=rttmodconvGPmass_;
      rtnewGPmass_=rtconvGPmass_;
      rnewGPmass_=rconvGPmass_;
      QnewGPdampstoch_=QconvGPdampstoch_;
      break;
    }

    case ELEMENTS::struct_calc_brownianforce:
    case ELEMENTS::struct_calc_brownianstiff:
    {
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get element velocity
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      if (act == ELEMENTS::struct_calc_brownianforce)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<2,2,1>(params,myvel,mydisp,NULL,&elevec1);
            else
              CalcBrownianForcesAndStiff<2,2,2>(params,myvel,mydisp,NULL,&elevec1);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<3,3,1>(params,myvel,mydisp,NULL,&elevec1);
            else
              CalcBrownianForcesAndStiff<3,2,2>(params,myvel,mydisp,NULL,&elevec1);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<4,4,1>(params,myvel,mydisp,NULL,&elevec1);
            else
              CalcBrownianForcesAndStiff<4,2,2>(params,myvel,mydisp,NULL,&elevec1);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<5,5,1>(params,myvel,mydisp,NULL,&elevec1);
            else
              CalcBrownianForcesAndStiff<5,2,2>(params,myvel,mydisp,NULL,&elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }

      }
      else if (act == ELEMENTS::struct_calc_brownianstiff)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<2,2,1>(params,myvel,mydisp,&elemat1,&elevec1);
            else
              CalcBrownianForcesAndStiff<2,2,2>(params,myvel,mydisp,&elemat1,&elevec1);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<3,3,1>(params,myvel,mydisp,&elemat1,&elevec1);
            else
              CalcBrownianForcesAndStiff<3,2,2>(params,myvel,mydisp,&elemat1,&elevec1);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<4,4,1>(params,myvel,mydisp,&elemat1,&elevec1);
            else
              CalcBrownianForcesAndStiff<4,2,2>(params,myvel,mydisp,&elemat1,&elevec1);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              CalcBrownianForcesAndStiff<5,5,1>(params,myvel,mydisp,&elemat1,&elevec1);
            else
              CalcBrownianForcesAndStiff<5,2,2>(params,myvel,mydisp,&elemat1,&elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }
      else
        dserror("You shouldn't be here.");

      break;
    }

    case ELEMENTS::struct_calc_stress:
    {
      dserror("No stress output implemented for beam3r elements");
      break;
    }

    case ELEMENTS::struct_calc_recover:
    {
      // do nothing here
      break;
    }

    default:
      std::cout << "\ncalled element with action type " << ActionType2String(act);
      dserror("This action type is not implemented for Beam3eb");
    break;
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::EvaluateNeumann(Teuchos::ParameterList& params,
                                           DRT::Discretization& discretization,
                                           DRT::Condition& condition,
                                           std::vector<int>& lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  // find out whether we will use a time curve
  bool usetime = true;
  double time = -1.0;

  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get<double>("total time",-1.0);

  if (time<0.0) usetime = false;

  // nnodetriad: number of nodes used for interpolation of triad field
  const int nnodetriad = NumNode();
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  int nnodecl = nnodetriad;
  if (centerline_hermite_) nnodecl=2;

  // vpernode: number of interpolated values per node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))
  int vpernode=1;
  if (centerline_hermite_)
    vpernode=2;

  // number of DOFs per node depending on type of node
  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  const DiscretizationType distype = this->Shape();

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  // amplitude of load curve at current time called, 6 components (3 forces, 3 moments)
  std::vector<double> curvefac(6,1.0);

  for (int i=0; i<6; ++i)
  {
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[i];

    if (curvenum>=0 && usetime)
      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(neumann_lineload));

  // declaration of variables in order to store shape functions
  // used for interpolation of triad field
  Epetra_SerialDenseVector I_i(nnodetriad);
  // used for interpolation of centerline
  Epetra_SerialDenseVector H_i(vpernode*nnodecl);

  // get values and switches from the condition

  // onoff is related to the first numdf flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the numdf "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");
  // funct is related to the numdf "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

  // integration points in parameter space and weights
  double xi=0.0;
  double wgt=0.0;

  // integration loops
  for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  {
    xi = intpoints.qxg[numgp][0];
    wgt = intpoints.qwgt[numgp];

    // evaluation of shape functions at Gauss points
    DRT::UTILS::shape_function_1D(I_i,xi,distype);
    if (centerline_hermite_)
      DRT::UTILS::shape_function_hermite_1D(H_i,xi,reflength_,line2);
    else
      DRT::UTILS::shape_function_1D(H_i,xi,distype);

    // position vector at the gauss point at reference configuration needed for function evaluation
    std::vector<double> X_ref(3,0.0);

    // calculate coordinates of corresponding Gauss point in reference configuration
    for (int node=0;node<nnodecl;node++)
    {
      for (int dim=0;dim<3;dim++)
      {
        X_ref[dim]+=H_i[vpernode*node]*Nodes()[node]->X()[dim];

        if (centerline_hermite_)
          X_ref[dim]+=H_i[vpernode*node+1]*(Trefnode_[node])(dim);
      }
    }

    double fac=0;
    fac = wgt * jacobiGPneumannline_[numgp];

    // load vector ar
    double ar[6];

    // loop the relevant dofs of a node
    for (int dof=0; dof<6; ++dof)
      ar[dof] = fac * (*onoff)[dof] * (*val)[dof] * curvefac[dof];
    double functionfac = 1.0;
    int functnum = -1;

    //sum up load components
    for (int dof=0; dof<6; ++dof)
    {
      if (functions)
        functnum = (*functions)[dof];
      else
        functnum = -1;

      // evaluate function at the position of the current GP
      if (functnum>0)
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof, &X_ref[0], time, NULL);   // TODO X_ref[0] is only x-coordinate. is this done on purpose?
      else
        functionfac = 1.0;

      for (int node=0; node<nnodecl; ++node)
      {
        if (dof<3)
        {
          elevec1[dofpercombinode*node+dof] += H_i[vpernode*node] *ar[dof] *functionfac;

          if (centerline_hermite_)
            elevec1[dofpercombinode*node+6+dof] += H_i[vpernode*node+1] *ar[dof] *functionfac;
        }
        else // dof<6
          elevec1[dofpercombinode*node+dof] += I_i[node] *ar[dof] *functionfac;
      }

      for (int node=nnodecl; node<nnodetriad; ++node)
        if (dof>2 && dof<6)
          elevec1[dofperclnode*nnodecl+dofpertriadnode*node+dof-3] += I_i[node] *ar[dof] *functionfac;
    }
  } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

  return 0;
}

/*----------------------------------------------------------------------------------------------------------------------*
 | get constitutive matrices from material law                                                               grill 03/16|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::GetConstitutiveMatrices(LINALG::TMatrix<T,3,3>& CN,
                                                           LINALG::TMatrix<T,3,3>& CM) const
{
  // first of all we get the material law
   Teuchos::RCP<const MAT::Material> currmat = Material();
   double ym = 0.0;
   double sm = 0.0;

   // assignment of material parameters; only St.Venant material is accepted for this beam
   switch(currmat->MaterialType())
   {
     case INPAR::MAT::m_stvenant:// only linear elastic material supported
     {
       const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
       ym = actmat->Youngs();
       sm = ym / (2*(1 + actmat->PoissonRatio()));
     }
     break;
     default:
       dserror("unknown or improper type of material law");
     break;
   }

   // defining material constitutive matrix CN between Gamma and N according to Jelenic 1999, section 2.4
   CN.PutScalar(0);
   CN(0,0) = ym*crosssec_;
   CN(1,1) = sm*crosssecshear_;
   CN(2,2) = sm*crosssecshear_;

   // defining material constitutive matrix CM between curvature and moment according to Jelenic 1999, section 2.4
   CM.PutScalar(0);
   CM(0,0) = sm*Irr_;
   CM(1,1) = ym*Iyy_;
   CM(2,2) = ym*Izz_;

  return;
}

/*----------------------------------------------------------------------------------------------------------------------*
 |push forward material stress vector and constitutive matrix to their spatial counterparts by rotation matrix Lambda   |
 |according to Romero 2004, eq. (3.10)                                                                       cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::pushforward(const LINALG::TMatrix<T,3,3>& Lambda,
                                                const LINALG::TMatrix<T,3,1>& stress_mat,
                                                const LINALG::TMatrix<T,3,3>& C_mat,
                                                LINALG::TMatrix<T,3,1>& stress_spatial,
                                                LINALG::TMatrix<T,3,3>& c_spatial) const
{
  // introduce auxiliary variable for pushforward of rotational matrices
  LINALG::TMatrix<T,3,3> temp;

  // push forward stress vector
  stress_spatial.Multiply(Lambda,stress_mat);

  // push forward constitutive matrix according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,C_mat);
  c_spatial.MultiplyNT(temp,Lambda);

   return;
}

/*------------------------------------------------------------------------------------------------------------*
 | calculate internal and inertia forces and their contributions to stiffmatrix                    cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::CalcInternalAndInertiaForcesAndStiff(Teuchos::ParameterList&        params,
                                                                      std::vector<double>&      disp,
                                                                      Epetra_SerialDenseMatrix* stiffmatrix,
                                                                      Epetra_SerialDenseMatrix* massmatrix,
                                                                      Epetra_SerialDenseVector* force,
                                                                      Epetra_SerialDenseVector* inertia_force)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  /********************************** Initialize/resize variables **************************************
   *****************************************************************************************************/

  //********************************* statmech periodic boundary conditions ****************************

  // unshift node positions, i.e. manipulate element displacement vector
  // as if there where no periodic boundary conditions
  if(StatMechParamsInterfacePtr() != Teuchos::null)
    UnShiftNodePosition(disp,nnodecl);

  //********************************** quantities valid for entire element *****************************
  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  // internal force vector
  LINALG::TMatrix<FADordouble, dofperclnode*nnodecl+dofpertriadnode*nnodetriad, 1> f_int(true);

  // reference triad Lambda_r and corresponding quaternion Q_r
  LINALG::TMatrix<FADordouble,3,3> Lambda_r(true);
  LINALG::TMatrix<FADordouble,4,1> Q_r(true);

  // angle of relative rotation between node I and J according to (3.10), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Phi_IJ(true);

  // clear internal (elastic) energy
  Eint_=0.0;

  //**************************************** nodal quantities *******************************************

  // current nodal DOFs relevant for centerline interpolation in total Lagrangian style, i.e. initial values + displacements
  LINALG::TMatrix<FADordouble,3*vpernode*nnodecl,1> disp_totlag_centerline(true);

  // quaternions of all nodal triads
  std::vector<LINALG::TMatrix<FADordouble,4,1> > Q_i(nnodetriad);

  // rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
  std::vector<LINALG::TMatrix<FADordouble,3,1> > Psi_li(nnodetriad);

  //*************************** physical quantities evaluated at a certain GP ***************************

  // derivation of beam centerline with respect to arc-length parameter: r'(x) from (2.12), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> r_s;
  // spin matrix related to vector r_s
  LINALG::TMatrix<FADordouble,3,3> r_s_hat;
  // interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Psi_l;
  /* derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter
   * at a certain Gauss point according to (3.11), Jelenic 1999*/
  LINALG::TMatrix<FADordouble,3,1> Psi_l_s;
  // triad at GP
  LINALG::TMatrix<FADordouble,3,3> Lambda;

  // 3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> K;
  // 3D vector of material axial and shear strains from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Gamma;

  // convected stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressN;
  LINALG::TMatrix<FADordouble,3,1> stressM;
  LINALG::TMatrix<FADordouble,3,3> CN;
  LINALG::TMatrix<FADordouble,3,3> CM;

  // spatial stresses n and m according to (3.10), Romero 2004 and spatial constitutive matrices c_n and c_m according to page 148, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressn;
  LINALG::TMatrix<FADordouble,3,1> stressm;
  LINALG::TMatrix<FADordouble,3,3> cn;
  LINALG::TMatrix<FADordouble,3,3> cm;

  //********************************** (generalized) shape functions ************************************
  /* Note: index i refers to the i-th shape function (i = 0 ... nnode*vpernode-1)
   * the vectors store individual shape functions, NOT an assembled matrix of shape functions)*/

  /* vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial shape functions evaluated at the numgp-th Gauss point
   * these shape functions are used for the interpolation of the triad field*/
  std::vector<LINALG::Matrix<1,nnodetriad> > I_i;
  // same for the derivatives
  std::vector<LINALG::Matrix<1,nnodetriad> > I_i_xi;

  /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite) shape functions evaluated at the numgp-th GP
   * these shape functions are used for the interpolation of the beam centerline*/
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i;
  // same for the derivatives
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i_xi;

  // vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at a certain Gauss point according to (3.18), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itilde(nnodetriad);

  // vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^nnode at a certain Gauss point according to (3.19), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itildeprime(nnodetriad);

  /*************************** update/compute quantities valid for entire element **********************
   *****************************************************************************************************/

  // update disp_totlag
  UpdateDispTotLagAndNodalTriads<nnodetriad,nnodecl,vpernode>(disp,disp_totlag_centerline,Q_i);

  // compute reference triad Lambda_r according to (3.9), Jelenic 1999
  CalcRefQuaternion<FADordouble>(Q_i[nodeI_],Q_i[nodeJ_],Q_r,Phi_IJ);
  LARGEROTATIONS::quaterniontotriad(Q_r,Lambda_r);

  // setup constitutive matrices
  GetConstitutiveMatrices<FADordouble>(CN,CM);

  /* compute nodal local rotations according to (3.8), Jelenic 1999
   * this is done individually for each node in order to avoid function argument std::vector<LINALG::TMatrix<...> >
   * a function with this argument type cannot be called with an argument of type std::vector<LINALG::Matrix<...> > (needed e.g. in SetUpReferenceGeometry) */
  for (unsigned int node=0; node<nnodetriad; ++node)
    CalcPsi_li<FADordouble>(Q_i[node],Q_r,Psi_li[node]);

  /******************************* elasticity: compute fint and stiffmatrix ****************************
   *****************************************************************************************************/

  //************************* residual and stiffmatrix contributions from forces ***********************
  // for these contributions, reduced integration is applied to avoid locking

  // get integration points for elasticity
  DRT::UTILS::IntegrationPoints1D gausspoints_elast_force(MyGaussRule(res_elastic_force));

  // reuse variables for individual shape functions and resize to new numgp
  I_i.resize(gausspoints_elast_force.nquad);
  H_i_xi.resize(gausspoints_elast_force.nquad);

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints_elast_force,I_i,this->Shape());

  EvaluateShapeFunctionDerivsAllGPs<nnodecl,vpernode>(gausspoints_elast_force,H_i_xi,this->Shape());


  // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints_elast_force.nquad; numgp++)
  {
    // weight of GP in parameter space
    const double wgt = gausspoints_elast_force.qwgt[numgp];

    Calc_r_s<nnodecl,vpernode,FADordouble>(disp_totlag_centerline, H_i_xi[numgp], jacobiGPelastf_[numgp], r_s);

    Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[numgp], Psi_l);
    Calc_Lambda<FADordouble>(Psi_l,Q_r,Lambda);

    // compute spin matrix related to vector rprime for later use
    LARGEROTATIONS::computespin<FADordouble>(r_s_hat,r_s);

    // compute material strains Gamma and K
    computeGamma<FADordouble>(r_s,Lambda,GammarefGP_[numgp],Gamma);

    // compute material stresses by multiplying strains with constitutive matrix
    stressN.Multiply(CN,Gamma);

    /* compute spatial stresses and constitutive matrices from convected ones according to Jelenic 1999, page 148, paragraph
     * between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward<FADordouble>(Lambda,stressN,CN,stressn,cn);

    /* computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     * to single blocks of matrix in eq. (4.3)*/
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      /* upper left block
       * note: jacobi factor cancels out because it is defined by ds=(ds/dxi)*dxi
       *       and I^{i'} in Jelenic1999 is derivative with respect to arc-length parameter in reference configuration s
       *       which can be computed from I_i_xi by multiplication with the inverse determinant: I^{i'}=I_i_s=I_i_xi*(dxi/ds) */
      for (int k=0; k<3; ++k)
      {
        f_int(dofpercombinode*node+k) += H_i_xi[numgp](vpernode*node)*stressn(k)*wgt;
        if (centerline_hermite_)
          f_int(dofpercombinode*node+6+k) += H_i_xi[numgp](vpernode*node+1)*stressn(k)*wgt;
      }

      // lower left block
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          f_int(dofpercombinode*node+3+i) -= r_s_hat(i,j)*stressn(j)*I_i[numgp](node)*wgt*jacobiGPelastf_[numgp];
    }
    for (unsigned int node=nnodecl; node<nnodetriad; ++node)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      // lower left block
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          f_int(dofperclnode*nnodecl+dofpertriadnode*node+i) -= r_s_hat(i,j)*stressn(j)*I_i[numgp](node)*wgt*jacobiGPelastf_[numgp];
    }

#ifndef BEAM3RAUTOMATICDIFF

    if (stiffmatrix != NULL)
    {
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[numgp]);

      /* computation of stiffness matrix according to Jelenic 1999, eq. (4.7); computation split up with respect
       * to single blocks of matrix in eq. (4.7).
       * note: again, jacobi factor cancels out in terms whith I^{i'}=I_i_s=I_i_xi*(dxi/ds) (see comment above)
       *       but be careful: Itildeprime and rprime are indeed derivatives with respect to arc-length parameter in reference configuration s */

      // auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;
      LINALG::Matrix<3,3> auxmatrix3;

      for (unsigned int nodei=0; nodei<nnodecl; nodei++)
      {
        for (unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // upper left block
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+j) += H_i_xi[numgp](vpernode*nodei)*H_i_xi[numgp](vpernode*nodej)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
              if (centerline_hermite_)
              {
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+j) += H_i_xi[numgp](vpernode*nodei+1)*H_i_xi[numgp](vpernode*nodej)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
                (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+6+j) += H_i_xi[numgp](vpernode*nodei)*H_i_xi[numgp](vpernode*nodej+1)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+6+j) += H_i_xi[numgp](vpernode*nodei+1)*H_i_xi[numgp](vpernode*nodej+1)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
              }
            }

          // lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(r_s_hat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+6+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej+1)*wgt;
            }

          // upper right block
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix3.Multiply(auxmatrix2,Itilde[nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+3+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+3+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;

        }
        for (unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // upper right block
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix3.Multiply(auxmatrix2,Itilde[nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;
        }
      }
      for (unsigned int nodei=nnodecl; nodei<nnodetriad; nodei++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        for (unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(r_s_hat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+6+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;

        }
        for (unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)
        {
          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;
        }
      }

    } // if (stiffmatrix != NULL)
#endif


    // add elastic energy from forces at this GP
    for(int dim=0; dim<3; dim++)
    {
      Eint_ += 0.5 * FADUTILS::CastToDouble(Gamma(dim)) * FADUTILS::CastToDouble(stressN(dim)) *jacobiGPelastf_[numgp]*wgt;
    }

  }//for(int numgp=0; numgp < gausspoints_elast_force.nquad; numgp++)


  //************************* residual and stiffmatrix contributions from moments ***********************

  // get integration points for elasticity
  DRT::UTILS::IntegrationPoints1D gausspoints_elast_moment(MyGaussRule(res_elastic_moment));

  // reuse variables for individual shape functions and resize to new numgp
  I_i.resize(gausspoints_elast_moment.nquad);
  I_i_xi.resize(gausspoints_elast_moment.nquad);

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAndDerivsAllGPs<nnodetriad,1>(gausspoints_elast_moment,I_i,I_i_xi,this->Shape());

  // reset norm of maximal bending curvature
  Kmax_=0.0;

  // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp<gausspoints_elast_moment.nquad; numgp++)
  {
    // weight of GP in parameter space
    const double wgt = gausspoints_elast_moment.qwgt[numgp];

    Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[numgp], Psi_l);
    Calc_Psi_l_s<nnodetriad,FADordouble>(Psi_li, I_i_xi[numgp], jacobiGPelastm_[numgp], Psi_l_s);
    Calc_Lambda<FADordouble>(Psi_l,Q_r,Lambda);

    // compute material curvature K
    computeK<FADordouble>(Psi_l,Psi_l_s,KrefGP_[numgp],K);

    // determine norm of maximal bending curvature at this GP and store in class variable if needed
    double Kmax = std::sqrt(FADUTILS::CastToDouble(K(1))*FADUTILS::CastToDouble(K(1)) + FADUTILS::CastToDouble(K(2))*FADUTILS::CastToDouble(K(2)) );
    if(Kmax > Kmax_)
      Kmax_=Kmax;

    // compute material stresses by multiplying curvature with constitutive matrix
    stressM.Multiply(CM,K);

    /* compute spatial stresses and constitutive matrix from material ones according to Jelenic 1999, page 148, paragraph
     * between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward<FADordouble>(Lambda,stressM,CM,stressm,cm);

    /* computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     * to single blocks of matrix in eq. (4.3)*/
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      // lower right block
      for (int i=0; i<3; ++i)
        f_int(dofpercombinode*node+3+i) += I_i_xi[numgp](node)*stressm(i)*wgt;
    }
    for (unsigned int node=nnodecl; node<nnodetriad; ++node)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      // lower right block
      for (int i=0; i<3; ++i)
        f_int(dofperclnode*nnodecl+dofpertriadnode*node+i) += I_i_xi[numgp](node)*stressm(i)*wgt;
    }

#ifndef BEAM3RAUTOMATICDIFF

    if (stiffmatrix != NULL)
    {
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[numgp]);
      computeItildeprime<nnodetriad,double>(Psi_l,Psi_l_s,Itildeprime,Phi_IJ,Lambda_r,Psi_li,I_i[numgp],I_i_xi[numgp],jacobiGPelastm_[numgp]);

      /* computation of stiffness matrix according to Jelenic 1999, eq. (4.7)*/

      // auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;

      for(unsigned int nodei=0; nodei<nnodecl; nodei++)
      {
        for(unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) -= auxmatrix1(i,j)*wgt;
        }
        for(unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix1(i,j)*wgt;
        }
      }

      for(unsigned int nodei=nnodecl; nodei<nnodetriad; nodei++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        for(unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) -= auxmatrix1(i,j)*wgt;
        }
        for(unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix1(i,j)*wgt;
        }
      }

    }//if (stiffmatrix != NULL)
#endif

    // add elastic energy from moments at this GP
    for(int dim=0; dim<3; dim++)
    {
      Eint_ += 0.5 * FADUTILS::CastToDouble(K(dim)) * FADUTILS::CastToDouble(stressM(dim)) *jacobiGPelastm_[numgp]*wgt;
    }

  }//for(int numgp=0; numgp < gausspoints_elast_moment.nquad; numgp++)

  if(force!=NULL)
  {
    for (unsigned int i=0; i<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; i++)
    {
      (*force)(i)= FADUTILS::CastToDouble(f_int(i));
    }
  }

#ifdef BEAM3RAUTOMATICDIFF

  if (stiffmatrix != NULL)
  {
    // compute stiffness matrix with FAD
    for (unsigned int i=0; i<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; i++)
    {
      for (unsigned int j=0; j<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }

    /* we need to transform the stiffmatrix because its entries are derivatives with respect to additive rotational increments
     * we want a stiffmatrix containing derivatives with respect to multiplicative rotational increments
     * therefore apply a trafo matrix to all those 3x3 blocks in stiffmatrix which correspond to derivation with respect to rotational DOFs
     * the trafo matrix is simply the T-Matrix (see Jelenic1999, (2.4)): \Delta_{mult} \vec \theta_{inode} = \mat T(\vec \theta_{inode} * \Delta_{addit} \vec \theta_{inode}*/

    LINALG::TMatrix<FAD,3,3> tempmat(true);
    LINALG::TMatrix<FAD,3,3> newstiffmat(true);
    LINALG::TMatrix<FAD,3,3> Tmat(true);
    LINALG::TMatrix<FAD,3,1> theta_totlag_j(true);

    for (unsigned int jnode=0; jnode<nnodecl; jnode++)
    {
      // compute physical total angle theta_totlag
      LARGEROTATIONS::quaterniontoangle(Q_i[jnode],theta_totlag_j);

      // compute Tmatrix of theta_totlag_i
      Tmat = LARGEROTATIONS::Tmatrix(theta_totlag_j);

      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        // block1: derivative of nodal positions with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

        // block3: derivative of nodal tangents with respect to theta (rotational DOFs)
        if(centerline_hermite_)
        {
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+6+i,dofpercombinode*jnode+3+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*inode+6+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();
        }

      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

      }
    }

    for (unsigned int jnode=nnodecl; jnode<nnodetriad; jnode++)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      // compute physical total angle theta_totlag
      LARGEROTATIONS::quaterniontoangle(Q_i[jnode],theta_totlag_j);

      // compute Tmatrix of theta_totlag_i
      Tmat = LARGEROTATIONS::Tmatrix(theta_totlag_j);

      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        // block1: derivative of nodal positions with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

        // block3: derivative of nodal tangents with respect to theta (rotational DOFs)
        if(centerline_hermite_)
        {
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+6+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*inode+6+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();
        }

      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)
      {
        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

      }
    }

  }
#endif


  /******************************* inertia: compute fint and massmatrix ********************************
   *****************************************************************************************************/

  // calculation of inertia forces/moments and massmatrix; neglect inertia in case of Brownian dynamics (Statmech) simulation
  if ( (massmatrix != NULL or inertia_force != NULL) and (!statmechprob_) )
  {

    /* calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
     * finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
     * needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
     * steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
     * consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
     * a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
     * our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
     * according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)*/

    /* Update: we now use a multiplicative update of rotational DOFs on time integrator level. Moreover, a new Lie group GenAlpha has been implemented
     *         that consistently updates the discrete TRANSLATIONAL velocity and acceleration vectors according to this element-internal scheme. This would allow us to
     *         use the global vel and acc vector at least for translational inertia contributions. Nevertheless, we stick to this completely element-internal
     *         temporal discretization of spatially continuous variables (angular velocity and acceleration) because the reverse order of discretization (spatial -> temporal)
     *         is much more intricate basically because of the triad interpolation. See also the discussion in Christoph Meier's Dissertation on this topic. (Maximilian Grill, 08/16)*/

    double dt = 1000.0;
    double beta = -1.0;
    double gamma = -1.0;
    double alpha_f = -1.0;
    double alpha_m = -1.0;

    if (this->IsParamsInterface())
    {
      dt = ParamsInterface().GetDeltaTime();
      beta = ParamsInterface().GetBeamParamsInterfacePtr()->GetBeta();
      gamma = ParamsInterface().GetBeamParamsInterfacePtr()->GetGamma();
      alpha_f = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlphaf();
      alpha_m = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlpham();
    }
    else
    {
      beta = params.get<double>("rot_beta",1000);
      gamma = params.get<double>("rot_gamma",1000);
      alpha_f = params.get<double>("rot_alphaf",1000);
      alpha_m = params.get<double>("rot_alpham",1000);
      dt = params.get<double>("delta time",1000);
    }

    const bool materialintegration=true;        // TODO unused?
    const double diff_factor_vel = gamma/(beta*dt);
    const double diff_factor_acc = (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f));

    LINALG::Matrix<3,3> Lambdanewmass(true);
    LINALG::Matrix<3,3> Lambdaconvmass(true);

    // get the material law
    Teuchos::RCP<const MAT::Material> currmat = Material();
    double rho = 0.0;

    // assignment of material parameters; only St.Venant material is accepted for this beam
    switch(currmat->MaterialType())
    {
      case INPAR::MAT::m_stvenant:// only linear elastic material supported
      {
        const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
        rho = actmat->Density();
      }
      break;
      default:
        dserror("unknown or improper type of material law");
      break;
    }

    /* tensor of mass moments of inertia and cross-section value. These values are used in order to artificially scale
     * the the translational and rotational inertia terms with given input parameters if necessary*/
    LINALG::Matrix<3,3> Jp(true);
    Jp(0,0)=inertscalerot1_*(Iyy_+Izz_);
    Jp(1,1)=inertscalerot2_*Iyy_;
    Jp(2,2)=inertscalerot2_*Izz_;
    Jp.Scale(rho);

    const double scaledcrosssec=inertscaletrans_*crosssec_;

    // get integration scheme for inertia forces and mass matrix
    DRT::UTILS::IntegrationPoints1D gausspoints_mass(MyGaussRule(res_inertia));
    // reuse variables for individual shape functions and resize to new numgp
    I_i.resize(gausspoints_mass.nquad);
    H_i.resize(gausspoints_mass.nquad);

    // evaluate all shape functions at all specified Gauss points
    EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints_mass,I_i,this->Shape());
    EvaluateShapeFunctionsAllGPs<nnodecl,vpernode>(gausspoints_mass,H_i,this->Shape());

    // Calculate current centerline position at gauss points (needed for element intern time integration)
    for (int gp=0; gp<gausspoints_mass.nquad; gp++)//loop through Gauss points
      Calc_r<nnodecl,vpernode,double>(FADUTILS::CastToDouble<FADordouble,3*vpernode*nnodecl,1>(disp_totlag_centerline),H_i[gp],rnewGPmass_[gp]);

    Ekin_=0.0;
    L_=0.0;
    P_=0.0;

    for (int gp=0; gp<gausspoints_mass.nquad; gp++)//loop through Gauss points
    {
      // weight of GP in parameter space
      const double wgtmass = gausspoints_mass.qwgt[gp];

      LINALG::Matrix<3,3> Jp_bar(Jp);
      Jp_bar.Scale(diff_factor_acc);

      LINALG::Matrix<3,1> dL(true);

      // update quaternions at GPs for exact Gauss quadrature
      Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[gp], Psi_l);
      Calc_Qgauss<double>(FADUTILS::CastToDouble<FADordouble,3,1>(Psi_l),FADUTILS::CastToDouble<FADordouble,4,1>(Q_r),QnewGPmass_[gp]);
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[gp]);

      Lambdanewmass.Clear();
      Lambdaconvmass.Clear();
      // compute current and old triad at Gauss point
      LARGEROTATIONS::quaterniontotriad<double>(QnewGPmass_[gp],Lambdanewmass);
      LARGEROTATIONS::quaterniontotriad<double>(QconvGPmass_[gp],Lambdaconvmass);

      // rotation between last converged position and current position expressed as a quaternion
      LINALG::Matrix<4,1>  deltaQ(true);
      LARGEROTATIONS::quaternionproduct<double>(LARGEROTATIONS::inversequaternion<double>(QconvGPmass_[gp]),QnewGPmass_[gp],deltaQ);

      // spatial rotation between last converged position and current position expressed as a three element rotation vector
      LINALG::Matrix<3,1> deltatheta(true);
      LARGEROTATIONS::quaterniontoangle<double>(deltaQ,deltatheta);

      // compute material counterparts of spatial vectors
      LINALG::Matrix<3,1> deltaTHETA(true);
      LINALG::Matrix<3,1> Wconvmass(true);
      LINALG::Matrix<3,1> Wnewmass(true);
      LINALG::Matrix<3,1> Aconvmass(true);
      LINALG::Matrix<3,1> Anewmass(true);
      LINALG::Matrix<3,1> Amodconvmass(true);
      LINALG::Matrix<3,1> Amodnewmass(true);
      deltaTHETA.MultiplyTN(Lambdanewmass,deltatheta);
      Wconvmass.MultiplyTN(Lambdaconvmass,wconvGPmass_[gp]);
      Aconvmass.MultiplyTN(Lambdaconvmass,aconvGPmass_[gp]);
      Amodconvmass.MultiplyTN(Lambdaconvmass,amodconvGPmass_[gp]);

      /* update angular velocities and accelerations according to Newmark time integration scheme either in
       * material description (see Jelenic, 1999, p. 146, equations (2.8) and (2.9)) or in spatial description
       * (for testing purposes, not recommended by Jelenic). The corresponding equations are adapted according to
       * the gen-alpha Lie group time integration scheme proposed in [Arnold, Brüls (2007)], [Brüls, Cardona, 2010]
       * and [Brüls, Cardona, Arnold (2012)]. In the predictor step of the time integration the following
       * formulas automatically deliver a constant displacement (deltatheta=0), consistent velocity and consistent acceleration
       * predictor. This fact has to be reflected in a consistent manner by the choice of the predictor in the input file*/
      if (materialintegration)
      {
        for (int i=0;i<3;i++)
        {
          Anewmass(i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltaTHETA(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*Wconvmass(i)      \
                        -alpha_f/(1.0-alpha_f)*Aconvmass(i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*Amodconvmass(i);

          Wnewmass(i)=gamma/(beta*dt)*deltaTHETA(i)+(1-gamma/beta)*Wconvmass(i)+dt*(1-gamma/(2*beta))*Amodconvmass(i);

          Amodnewmass(i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*Anewmass(i) + alpha_f*Aconvmass(i) - alpha_m*Amodconvmass(i));
        }
        wnewGPmass_[gp].Multiply(Lambdanewmass,Wnewmass);
        anewGPmass_[gp].Multiply(Lambdanewmass,Anewmass);
        amodnewGPmass_[gp].Multiply(Lambdanewmass,Amodnewmass);
      }
      else
      {
        for (int i=0;i<3;i++)
        {
          wnewGPmass_[gp](i)=gamma/(beta*dt)*deltatheta(i)+(1-gamma/beta)*wconvGPmass_[gp](i)+dt*(1-gamma/(2*beta))*amodconvGPmass_[gp](i);

          anewGPmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltatheta(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*wconvGPmass_[gp](i)      \
                        -alpha_f/(1.0-alpha_f)*aconvGPmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*amodconvGPmass_[gp](i);

          amodnewGPmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*anewGPmass_[gp](i) + alpha_f*aconvGPmass_[gp](i) - alpha_m*amodconvGPmass_[gp](i) );
        }
        Wnewmass.MultiplyTN(Lambdanewmass,wnewGPmass_[gp]);
        Anewmass.MultiplyTN(Lambdanewmass,anewGPmass_[gp]);
        Amodnewmass.MultiplyTN(Lambdanewmass,amodnewGPmass_[gp]);
      }

      LINALG::Matrix<3,1> deltar(true);
      for (int i=0;i<3;i++)
      {
        deltar(i)=rnewGPmass_[gp](i)-rconvGPmass_[gp](i);
      }
      for (int i=0;i<3;i++)
      {
        rttnewGPmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltar(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*rtconvGPmass_[gp](i)      \
                      -alpha_f/(1.0-alpha_f)*rttconvGPmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*rttmodconvGPmass_[gp](i);

        rtnewGPmass_[gp](i)=gamma/(beta*dt)*deltar(i)+(1-gamma/beta)*rtconvGPmass_[gp](i)+dt*(1-gamma/(2*beta))*rttmodconvGPmass_[gp](i);

        rttmodnewGPmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*rttnewGPmass_[gp](i) + alpha_f*rttconvGPmass_[gp](i) - alpha_m*rttmodconvGPmass_[gp](i) );
      }

      // spin matrix of the material angular velocity, i.e. S(W)
      LINALG::Matrix<3,3> SWnewmass(true);
      LARGEROTATIONS::computespin<double>(SWnewmass,Wnewmass);
      LINALG::Matrix<3,1> Jp_Wnewmass(true);
      LINALG::Matrix<3,1> auxvector1(true);
      LINALG::Matrix<3,1> Pi_t(true);
      Jp_Wnewmass.Multiply(Jp,Wnewmass);
      for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
          auxvector1(i)+=SWnewmass(i,j)*Jp_Wnewmass(j)+Jp(i,j)*Anewmass(j);

      Pi_t.Multiply(Lambdanewmass,auxvector1);
      LINALG::Matrix<3,1> r_tt(true);
      LINALG::Matrix<3,1> r_t(true);
      LINALG::Matrix<3,1> r(true);

      r_tt = rttnewGPmass_[gp];
      r_t = rtnewGPmass_[gp];
      r = rnewGPmass_[gp];

      LINALG::Matrix<3,3> S_r(true);
      LARGEROTATIONS::computespin<double>(S_r,r);
      dL.Multiply(S_r,r_t);
      dL.Scale(rho*scaledcrosssec);
      LINALG::Matrix<3,1> Lambdanewmass_Jp_Wnewmass(true);
      Lambdanewmass_Jp_Wnewmass.Multiply(Lambdanewmass,Jp_Wnewmass);
      dL.Update(1.0,Lambdanewmass_Jp_Wnewmass,1.0);
      for (int i=0;i<3;i++)
      {
        L_(i)+=wgtmass*jacobiGPmass_[gp]*dL(i);
        P_(i)+=wgtmass*jacobiGPmass_[gp]*rho*scaledcrosssec*r_t(i);
      }

      LINALG::Matrix<3,3> S_Pit(true);
      LARGEROTATIONS::computespin<double>(S_Pit,Pi_t);
      LINALG::Matrix<3,3> SJpWnewmass(true);
      LARGEROTATIONS::computespin<double>(SJpWnewmass,Jp_Wnewmass);
      LINALG::Matrix<3,3> SWnewmass_Jp(true);
      SWnewmass_Jp.Multiply(SWnewmass,Jp);
      Jp_bar.Update(diff_factor_vel,SWnewmass_Jp,1.0);
      Jp_bar.Update(-diff_factor_vel,SJpWnewmass,1.0);

      LINALG::Matrix<3,3> Tmatrix(true);
      Tmatrix=LARGEROTATIONS::Tmatrix(deltatheta);

      LINALG::Matrix<3,3> Lambdanewmass_Jpbar(true);
      Lambdanewmass_Jpbar.Multiply(Lambdanewmass, Jp_bar);
      LINALG::Matrix<3,3> LambdaconvmassT_Tmatrix(true);
      LambdaconvmassT_Tmatrix.MultiplyTN(Lambdaconvmass, Tmatrix);
      LINALG::Matrix<3,3> Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix(true);
      Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix.Multiply(Lambdanewmass_Jpbar, LambdaconvmassT_Tmatrix);
      LINALG::Matrix<3,3> auxmatrix1(true);
      auxmatrix1.Update(-1.0,S_Pit,1.0);
      auxmatrix1.Update(1.0,Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix,1.0);

      if (inertia_force != NULL)
      {
        // inertia forces
        for (int i=0; i<3; i++)
        {
          for (unsigned int node=0; node<nnodecl; node++)
          {
            // translational contribution
            (*inertia_force)(dofpercombinode*node+i) += jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*node)*r_tt(i);
            if (centerline_hermite_)
              (*inertia_force)(dofpercombinode*node+6+i) += jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*node+1)*r_tt(i);
            //rotational contribution
            (*inertia_force)(dofpercombinode*node+3+i) += jacobiGPmass_[gp]*wgtmass*I_i[gp](node)*Pi_t(i);
          }
          for (unsigned int node=nnodecl; node<nnodetriad; node++)    // this loop is only entered in case of nnodetriad>nnodecl
          {
            // rotational contribution
            (*inertia_force)(dofperclnode*nnodecl+dofpertriadnode*node+i) += jacobiGPmass_[gp]*wgtmass*I_i[gp](node)*Pi_t(i);
          }
        }
      }

      if (massmatrix != NULL)
      {
        // linearization of inertia forces: massmatrix
        for (unsigned int jnode=0; jnode<nnodecl; jnode++)
        {
          // translational contribution
          for (unsigned int inode=0; inode<nnodecl; inode++)
            for (int k=0;k<3;k++)
            {
              (*massmatrix)(dofpercombinode*inode+k,dofpercombinode*jnode+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode);
              if (centerline_hermite_)
              {
                (*massmatrix)(dofpercombinode*inode+6+k,dofpercombinode*jnode+6+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode+1);
                (*massmatrix)(dofpercombinode*inode+k,dofpercombinode*jnode+6+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode+1);
                (*massmatrix)(dofpercombinode*inode+6+k,dofpercombinode*jnode+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode);
              }
            }

          // rotational contribution
          LINALG::Matrix<3,3> auxmatrix2(true);
          auxmatrix2.Multiply(auxmatrix1,Itilde[jnode]);
          for (unsigned int inode=0; inode<nnodecl; inode++)
          {
            for (int i=0; i<3; i++)
              for (int j=0; j<3; j++)
                (*massmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
          }
          for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
          {
            for (int i=0; i<3; i++)
              for (int j=0; j<3; j++)
                (*massmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
          }
        }
        for (unsigned int jnode=nnodecl; jnode<nnodetriad; ++jnode)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // rotational contribution
          LINALG::Matrix<3,3> auxmatrix2(true);
          auxmatrix2.Multiply(auxmatrix1,Itilde[jnode]);
          for (unsigned int inode=0; inode<nnodecl; inode++)
          {
            for (int i=0; i<3; i++)
              for (int j=0; j<3; j++)
                (*massmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
          }
          for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
          {
            for (int i=0; i<3; i++)
              for (int j=0; j<3; j++)
                (*massmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
          }
        }
      }

      // Calculation of kinetic energy
      LINALG::Matrix<1,1> ekinrot(true);
      LINALG::Matrix<1,1> ekintrans(true);
      ekinrot.MultiplyTN(Wnewmass,Jp_Wnewmass);
      ekintrans.MultiplyTN(r_t,r_t);
      Ekin_+=0.5*(ekinrot.Norm2() + rho*scaledcrosssec*ekintrans.Norm2())*jacobiGPmass_[gp]*wgtmass;
      Ekintorsion_+=0.5* Wnewmass(0)*Jp_Wnewmass(0)*jacobiGPmass_[gp]*wgtmass;
      Ekinbending_+=0.5* Wnewmass(1)*Jp_Wnewmass(1)*jacobiGPmass_[gp]*wgtmass;
      Ekinbending_+=0.5* Wnewmass(2)*Jp_Wnewmass(2)*jacobiGPmass_[gp]*wgtmass;
      Ekintrans_+=0.5*rho*scaledcrosssec*ekintrans.Norm2()*jacobiGPmass_[gp]*wgtmass;

      Jp_Wnewmass.Multiply(Jp,Wnewmass);
    }//for (int gp=0; gp<gausspoints_mass.nquad; gp++)

    // In Lie group GenAlpha algorithm, the mass matrix is multiplied with factor (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)) later.
    // so we apply inverse factor here because the correct prefactors for displacement/velocity/acceleration dependent terms have been applied individually above
    if (massmatrix!=NULL)
      massmatrix->Scale(beta*dt*dt*(1.0-alpha_f)/(1.0-alpha_m));
  }//if (massmatrix != NULL or inertia_force != NULL)

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | calculation of thermal (i.e. stochastic) and damping forces according to Brownian dynamics      grill 06/16|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::CalcBrownianForcesAndStiff(Teuchos::ParameterList&   params,
                                                       std::vector<double>&      vel,
                                                       std::vector<double>&      disp,
                                                       Epetra_SerialDenseMatrix* stiffmatrix,
                                                       Epetra_SerialDenseVector* force)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  //********************************* statmech periodic boundary conditions ****************************

  // unshift node positions, i.e. manipulate element displacement vector
  // as if there where no periodic boundary conditions
  if(StatMechParamsInterfacePtr() != Teuchos::null)
    UnShiftNodePosition(disp,nnodecl);

  /****** update/compute key variables describing displacement and velocity state of this element *****/

  // current nodal DOFs relevant for centerline interpolation in total Lagrangian style, i.e. initial values + displacements
  LINALG::Matrix<3*vpernode*nnodecl,1> disp_totlag_centerline(true);

  // discrete centerline (i.e. translational) velocity vector
  LINALG::Matrix<3*vpernode*nnodecl,1> vel_centerline(true);

  // quaternions of all nodal triads
  std::vector<LINALG::Matrix<4,1> > Q_i(nnodetriad);

  // update disp_totlag_centerline and nodal triads
  UpdateDispTotLagAndNodalTriads<nnodetriad,nnodecl,vpernode>(disp,disp_totlag_centerline,Q_i);

  // update current values of centerline (i.e. translational) velocity
  ExtractCenterlineDofValues<nnodecl,vpernode,double>(vel,vel_centerline);

  /****** compute and assemble force and stiffness contributions from viscous damping and stochastic forces *****/

  // add stiffness and forces due to translational damping effects
  EvaluateTranslationalDamping<nnodecl,vpernode,3>(params,vel_centerline,disp_totlag_centerline,stiffmatrix,force);

  // add stiffness and forces (i.e. moments) due to rotational damping effects
  EvaluateRotationalDamping<nnodetriad,nnodecl,vpernode,3>(params,Q_i,stiffmatrix,force);

  // add stochastic forces and (if required) resulting stiffness
  EvaluateStochasticForces<nnodecl,vpernode,3,3>(params,disp_totlag_centerline,stiffmatrix,force);

}

/*------------------------------------------------------------------------------------------------------------*
 | update (total) displacement vector and set nodal triads (as quaternions)                        grill 03/16|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads(const std::vector<double>&            disp,
                                                           LINALG::Matrix<3*vpernode*nnodecl,1>& disp_totlag_centerline,
                                                           std::vector<LINALG::Matrix<4,1> >&    Q_i)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  // get current values of translational nodal DOFs in total Lagrangean manner (initial value + disp)
  // rotational DOFs need different handling, depending on whether FAD is used or not (see comment below)
  ExtractCenterlineDofValues<nnodecl,vpernode,double>(disp,disp_totlag_centerline);
  AddRefValuesDispCenterline<nnodecl,vpernode,double>(disp_totlag_centerline);

  // get current displacement values of rotational DOFs (i.e. relative rotation with respect to reference config)
  for(int dim=0; dim<3; dim++)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
      dispthetanewnode_[node](dim) = disp[dofpercombinode*node+3+dim];

    for (unsigned int node=nnodecl; node<nnodetriad; ++node)
      dispthetanewnode_[node](dim) = disp[dofperclnode*nnodecl+dofpertriadnode*node+dim];
  }


  // rotational displacement at a certain node in quaternion form
  LINALG::Matrix<4,1> deltaQ;
  // initial nodal rotation vector in quaternion form
  LINALG::Matrix<4,1> Q0;

  // Compute current nodal triads
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    // get initial nodal rotation vectors and transform to quaternions
    LARGEROTATIONS::angletoquaternion(theta0node_[node],Q0);

    // rotate initial triads by relative rotation vector from displacement vector (via quaternion product)
    LARGEROTATIONS::angletoquaternion(dispthetanewnode_[node],deltaQ);
    LARGEROTATIONS::quaternionproduct(Q0,deltaQ,Qnewnode_[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
    Qnewnode_[node].Scale(1/Qnewnode_[node].Norm2());

    // copy quaternions of nodal triads to TMatrix FADordouble
    for (unsigned int i=0; i<4; ++i)
      Q_i[node](i) = Qnewnode_[node](i);
  }

}

/*------------------------------------------------------------------------------------------------------------*
 | update (total) displacement vector and set nodal triads (as quaternions)                        grill 03/16|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads(const std::vector<double>&                         disp,
                                                           LINALG::TMatrix<FADordouble,3*vpernode*nnodecl,1>& disp_totlag_centerline,
                                                           std::vector<LINALG::TMatrix<FADordouble,4,1> >&    Q_i)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  // get current values of translational nodal DOFs in total Lagrangean manner (initial value + disp)
  // rotational DOFs need different handling, depending on whether FAD is used or not (see comment below)
  ExtractCenterlineDofValues<nnodecl,vpernode,FADordouble>(disp,disp_totlag_centerline);
  AddRefValuesDispCenterline<nnodecl,vpernode,FADordouble>(disp_totlag_centerline);

  // get current displacement values of rotational DOFs (i.e. relative rotation with respect to reference config)
  for(int dim=0; dim<3; dim++)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
      dispthetanewnode_[node](dim) = disp[dofpercombinode*node+3+dim];

    for (unsigned int node=nnodecl; node<nnodetriad; ++node)
      dispthetanewnode_[node](dim) = disp[dofperclnode*nnodecl+dofpertriadnode*node+dim];
  }


  // rotational displacement at a certain node in quaternion form
  LINALG::Matrix<4,1> deltaQ;
  // initial nodal rotation vector in quaternion form
  LINALG::Matrix<4,1> Q0;

  // Compute current nodal triads
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    // get initial nodal rotation vectors and transform to quaternions
    LARGEROTATIONS::angletoquaternion(theta0node_[node],Q0);

    // rotate initial triads by relative rotation vector from displacement vector (via quaternion product)
    LARGEROTATIONS::angletoquaternion(dispthetanewnode_[node],deltaQ);
    LARGEROTATIONS::quaternionproduct(Q0,deltaQ,Qnewnode_[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
    Qnewnode_[node].Scale(1/Qnewnode_[node].Norm2());

    // copy quaternions of nodal triads to TMatrix FADordouble
    for (unsigned int i=0; i<4; ++i)
      Q_i[node](i) = Qnewnode_[node](i);
  }

#ifdef BEAM3RAUTOMATICDIFF
  // set differentiation variables for FAD: translational DOFs
  for (int dim=0; dim<3; ++dim)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      disp_totlag_centerline(dofperclnode*node+dim).diff(dofpercombinode*node+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);

      // have Hermite interpolation? then set tangent DOFs as well
      if(vpernode==2)
        disp_totlag_centerline(dofperclnode*node+3+dim).diff(dofpercombinode*node+6+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);
    }
  }

  // rotation vector theta at a specific node in a total Lagrangean manner (with respect to global reference coordinate system)
  std::vector<LINALG::TMatrix<FAD,3,1> > theta_totlag_i(nnodetriad);

  // compute nodal quaternions based on multiplicative increments of rotational DOFs
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    // compute physical total angle theta_totlag
    LARGEROTATIONS::quaterniontoangle(Q_i[node],theta_totlag_i[node]);
  }

  // set differentiation variables for FAD: rotational DOFs
  for (unsigned int dim=0; dim<3; ++dim)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
      theta_totlag_i[node](dim).diff(dofpercombinode*node+3+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);

    for (unsigned int node=nnodecl; node<nnodetriad; ++node)
      theta_totlag_i[node](dim).diff(dofperclnode*nnodecl+dofpertriadnode*node+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);
  }

  /* Attention: although the nodal quaternions Q_i have already been computed correctly, we need the following step
   *            in order to track the dependency of subsequently calculated quantities via FAD */
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    Q_i[node].PutScalar(0.0);
    LARGEROTATIONS::angletoquaternion(theta_totlag_i[node],Q_i[node]);
  }

#endif //#ifndef BEAM3RAUTOMATICDIFF

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                   cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::lumpmass(Epetra_SerialDenseMatrix* massmatrix)
{
  // lump mass matrix
  if (massmatrix != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*massmatrix).N(); ++c) // parse columns
    {
      double d = 0.0;
      for (int r=0; r<(*massmatrix).M(); ++r) // parse rows
      {
        d += (*massmatrix)(r,c); // accumulate row entries
        (*massmatrix)(r,c) = 0.0;
      }

      (*massmatrix)(c,c) = d; // apply sum of row entries on diagonal
    }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::EvaluatePTC(Teuchos::ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  //apply PTC rotation damping term using a Lobatto integration rule; implemented for 2 nodes only
  if(nnode>2 or centerline_hermite_)
    dserror("PTC was originally implemented for 2-noded Reissner beam element only. Check functionality for "
        "numnodes>2 and/or Hermite interpolation and extend if needed!");

  for (unsigned int node=0; node<nnode; node++)
  {
    //computing angle increment from current position in comparison with last converged position for damping
    LINALG::Matrix<4,1> deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconvnode_[node]),Qnewnode_[node],deltaQ);
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    //scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale( params.get<double>("crotptc",0.0) );

    //each node gets a block diagonal damping term; the Lobatto integration weight is 0.5 for 2-noded elements
    //jacobi determinant is constant and equals 0.5*refelelength for 2-noded elements
    for(int k=0; k<3; k++)
      for (int l=0; l<3; l++)
        elemat1(node*6+3+k,node*6+3+l) += artstiff(k,l)*0.5*0.5*reflength_;

    //PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      elemat1(node*6+k,node*6+k) += params.get<double>("ctransptc",0.0)*0.5*0.5*reflength_;
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::HowManyRandomNumbersINeed() const
{
  // get Gauss rule for evaluation of stochastic force contributions
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  /* at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   * random numbers for the translational degrees of freedom */
#ifndef BEAM3RCONSTSTOCHFORCE
  return (3*gausspoints.nquad);
#else
  return (3);
#endif
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public)                                    cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3r::EvaluateRotationalDamping(Teuchos::ParameterList&          params,  //!<parameter list
                                              const std::vector<LINALG::Matrix<4,1> >& Qnode,
                                              Epetra_SerialDenseMatrix*                stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                force)  //!< element internal force vector
{
  const unsigned int dofperclnode = 3*vpernode;
  const unsigned int dofpertriadnode = 3;
  const unsigned int dofpercombinode = dofperclnode + dofpertriadnode;

  // get time step size
  double dt = 1000;
  if (IsParamsInterface())
    dt = ParamsInterface().GetDeltaTime();
  else
    dt = params.get<double>("delta time",1000);

  // get damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  GetDampingCoefficients(gamma);

  // get Gauss points and weights for evaluation of viscous damping contributions
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  // reference triad Lambda_r and corresponding quaternion Q_r
  LINALG::Matrix<3,3> Lambda_r(true);
  LINALG::Matrix<4,1> Q_r(true);

  // angle of relative rotation between node I and J according to (3.10), Jelenic 1999
  LINALG::Matrix<3,1> Phi_IJ(true);

  // rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,1> > Psi_li(nnodetriad);

  //*************************** physical quantities evaluated at a certain GP ***************************

  // interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
  LINALG::Matrix<3,1> Psi_l;

  // material triad and corresponding quaternion at a certain Gauss point
  LINALG::Matrix<3,3> LambdaGP;
  LINALG::Matrix<4,1> QnewGP;

  //********************************** (generalized) shape functions ************************************
  /* Note: index i refers to the i-th shape function (i = 0 ... nnodetriad-1)
   * the vectors store individual shape functions, NOT an assembled matrix of shape functions)*/

  /* vector whose numgp-th element is a 1xnnodetriad-matrix with all Lagrange polynomial shape functions evaluated at the numgp-th Gauss point
   * these shape functions are used for the interpolation of the triad field */
  std::vector<LINALG::Matrix<1,nnodetriad> > I_i(gausspoints.nquad);

  // evaluate all shape functions at all specified Gauss points
  EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints,I_i,this->Shape());

  /* vector with nnodetriad elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode according to (3.19), Jelenic 1999*/
  std::vector<LINALG::TMatrix<double,3,3> > Itilde(nnodetriad);


  // compute reference triad Lambda_r according to (3.9), Jelenic 1999
  CalcRefQuaternion<double>(Qnode[nodeI_],Qnode[nodeJ_],Q_r,Phi_IJ);
  LARGEROTATIONS::quaterniontotriad(Q_r,Lambda_r);

  // compute nodal local rotations according to (3.8), Jelenic 1999
  for (unsigned int node=0; node<nnodetriad; ++node)
    CalcPsi_li(Qnode[node],Q_r,Psi_li[node]);


  for (int gp=0; gp<gausspoints.nquad; gp++)//loop through Gauss points         // ToDo cleanup of auxiliary variables, comments, double-check linearization !!!
  {
    // update quaternions at GPs for exact Gauss quadrature
    Calc_Psi_l<nnodetriad>(Psi_li, I_i[gp], Psi_l);
    Calc_Qgauss<double>(Psi_l,Q_r,QnewGP);

    // store in class variable in order to get QconvGPmass_ in subsequent time step
    QnewGPdampstoch_[gp]=QnewGP;

    // compute triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(QnewGP,LambdaGP);

    // rotation between last converged position and current position expressed as a quaternion
    LINALG::Matrix<4,1> deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(QconvGPdampstoch_[gp]),QnewGP,deltaQ);

    // extract rotation vector from quaternion
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    // angular velocity at this Gauss point according to backward Euler scheme
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);

    // compute matrix Lambda*[gamma(2) 0 0 \\ 0 0 0 \\ 0 0 0]*Lambda^t = gamma(2) * g_1 \otimes g_1
    // where g_1 is first base vector, i.e. first column of Lambda
    LINALG::Matrix<3,3> g1g1gamma;
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        g1g1gamma(k,j) = (LambdaGP)(k,0)*(LambdaGP)(j,0)*gamma(2);

    // compute vector gamma(2) * g_1 \otimes g_1 * \omega
    LINALG::Matrix<3,1> g1g1gammaomega;
    g1g1gammaomega.Multiply(g1g1gamma,omega);

    if (force != NULL)
    {
      //loop over all nodes
      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        //loop over three dimensions in line direction
        for (unsigned int idim=0; idim<ndim; idim++)
          (*force)(dofpercombinode*inode+3+idim) += g1g1gammaomega(idim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)
      {
        //loop over three dimensions in line direction
        for (unsigned int idim=0; idim<ndim; idim++)
          (*force)(dofperclnode*nnodecl+dofpertriadnode*inode+idim) += g1g1gammaomega(idim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
      }
    }

    if (stiffmatrix != NULL)
    {
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[gp]);

      // compute matrix gamma(2) * g_1 \otimes g_1 * \omega * Tmat
      LINALG::Matrix<3,3> g1g1gammaTmat;
      g1g1gammaTmat.Multiply(g1g1gamma,LARGEROTATIONS::Tmatrix(deltatheta));      // ToDo: check this term: why do we need this Tmatrix? should be multiplicative increment anyway?!

      // compute spin matrix S(\omega)
      LINALG::Matrix<3,3> Sofomega;
      LARGEROTATIONS::computespin(Sofomega,omega);

      // compute matrix gamma(2) * g_1 \otimes g_1 *S(\omega)
      LINALG::Matrix<3,3> g1g1gammaSofomega;
      g1g1gammaSofomega.Multiply(g1g1gamma,Sofomega);

      // compute spin matrix S(gamma(2) * g_1 \otimes g_1 *\omega)
      LINALG::Matrix<3,3> Sofg1g1gammaomega;
      LARGEROTATIONS::computespin(Sofg1g1gammaomega,g1g1gammaomega);

      // auxiliary matrices
      LINALG::Matrix<3,3> sum(true);
      LINALG::Matrix<3,3> auxmatrix(true);

      sum += g1g1gammaTmat;
      sum.Scale(1/dt);
      sum += g1g1gammaSofomega;
      sum -= Sofg1g1gammaomega;

      // loop over first nnodecl line nodes
      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        // loop over first nnodecl column nodes
        for (unsigned int jnode=0; jnode<nnodecl; jnode++)
        {
          auxmatrix.Multiply(sum,Itilde[jnode]);

          // loop over three dimensions in line and column direction
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<3; jdim++)
            {
              (*stiffmatrix)(dofpercombinode*inode+3+idim,dofpercombinode*jnode+3+jdim) += auxmatrix(idim,jdim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
            }

        }
        for (unsigned int jnode=nnodecl; jnode<nnodetriad; jnode++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          auxmatrix.Multiply(sum,Itilde[jnode]);

          // loop over three dimensions in line and column direction
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<3; jdim++)
            {
              (*stiffmatrix)(dofpercombinode*inode+3+idim,dofperclnode*nnodecl+dofpertriadnode*jnode+jdim) += auxmatrix(idim,jdim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
            }
        }
      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        // loop over all column nodes
        for (unsigned int jnode=0; jnode<nnodecl; jnode++)
        {
          auxmatrix.Multiply(sum,Itilde[jnode]);

          // loop over three dimensions in line and column direction
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<3; jdim++)
            {
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+idim,dofpercombinode*jnode+3+jdim) += auxmatrix(idim,jdim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
            }

        }
        for (unsigned int jnode=nnodecl; jnode<nnodetriad; jnode++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          auxmatrix.Multiply(sum,Itilde[jnode]);

          // loop over three dimensions in line and column direction
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<3; jdim++)
            {
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+idim,dofperclnode*nnodecl+dofpertriadnode*jnode+jdim) += auxmatrix(idim,jdim)*(I_i[gp])(inode)*gausspoints.qwgt[gp]*jacobiGPdampstoch_[gp];
            }
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3r::EvaluateTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const LINALG::Matrix<ndim*vpernode*nnodecl,1>& vel_centerline,
                                                  const LINALG::Matrix<ndim*vpernode*nnodecl,1>& disp_totlag_centerline,
                                                  Epetra_SerialDenseMatrix*        stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector*        force)//!< element internal force vector
{
  /* only nodes for centerline interpolation are considered here (= first nnodecl nodes of this element);
     each of these nodes holds 3*vpernode translational DoFs AND 3 rotational DoFs */
  const unsigned int dofpernode = 3*vpernode+3;

  // get time step size
  double dt = 1000;
  if (IsParamsInterface())
    dt = ParamsInterface().GetDeltaTime();
  else
    dt = params.get<double>("delta time",1000);

  // velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;

  // position of beam centerline point corresponding to a certain Gauss point
  LINALG::Matrix<ndim,1> r(true);
  // tangent vector (derivative of beam centerline curve r with respect to arc-length parameter s)
  LINALG::Matrix<ndim,1> r_s(true);
  // velocity of beam centerline point relative to background fluid velocity
  LINALG::Matrix<ndim,1> vel_rel(true);

  // damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  GetDampingCoefficients(gamma);

  // viscous force vector per unit length at current GP
  LINALG::Matrix<ndim,1> f_visc(true);
  // damping matrix
  LINALG::Matrix<ndim,ndim> damp_mat(true);

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite) shape functions evaluated at the numgp-th GP
   * these shape functions are used for the interpolation of the beam centerline*/
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i(gausspoints.nquad);
  // same for the derivatives
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i_xi(gausspoints.nquad);

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAndDerivsAllGPs<nnodecl,vpernode>(gausspoints,H_i,H_i_xi,this->Shape());

  for(int gp=0; gp<gausspoints.nquad; gp++)
  {
    // compute position vector r of point in physical space corresponding to Gauss point
    Calc_r<nnodecl,vpernode,double>(disp_totlag_centerline, H_i[gp], r);

    // compute tangent vector t_{\par}=r' at current Gauss point
    Calc_r_s<nnodecl,vpernode,double>(disp_totlag_centerline, H_i_xi[gp], jacobiGPdampstoch_[gp], r_s);

    // compute velocity and gradient of background flow field at point r
    GetBackgroundVelocity<ndim,double>(params,r,velbackground,velbackgroundgrad);

    // compute velocity vector at this Gauss point via same interpolation as for centerline position vector
    CalcInterpolation<nnodecl,vpernode,3,double>(vel_centerline, H_i[gp], vel_rel);
    vel_rel -= velbackground;

    // loop over lines and columns of damping matrix
    for (unsigned int idim=0; idim<ndim; idim++)
      for (unsigned int jdim=0; jdim<ndim; jdim++)
        damp_mat(idim,jdim) = (idim==jdim)*gamma(1) + (gamma(0) - gamma(1))*r_s(idim)*r_s(jdim);

    // compute viscous force vector per unit length at current GP
    f_visc.Multiply(damp_mat,vel_rel);

    if (force != NULL)
    {
      // loop over all nodes used for centerline interpolation
      for (unsigned int inode=0; inode<nnodecl; inode++)
        // loop over dimensions
        for (unsigned int idim=0; idim<ndim; idim++)
        {
          (*force)(inode*dofpernode+idim) += H_i[gp](vpernode*inode)*jacobiGPdampstoch_[gp]*gausspoints.qwgt[gp]*f_visc(idim);
          if (centerline_hermite_)
            (*force)(inode*dofpernode+6+idim) += H_i[gp](vpernode*inode+1)*jacobiGPdampstoch_[gp]*gausspoints.qwgt[gp]*f_visc(idim);
        }
    }

    if (stiffmatrix != NULL)
    {
      // compute matrix product of damping matrix and gradient of background velocity
      LINALG::Matrix<ndim,ndim> dampmatvelbackgroundgrad(true);
      dampmatvelbackgroundgrad.Multiply(damp_mat,velbackgroundgrad);


      // loop over all nodes used for centerline interpolation
      for (unsigned int inode=0; inode<nnodecl; inode++)
        //loop over all column nodes
        for (unsigned int jnode=0; jnode<nnodecl; jnode++)
        {
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<ndim; jdim++)
            {
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode)*jacobiGPdampstoch_[gp]* damp_mat(idim,jdim) / dt;
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+jdim) -= gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode)*jacobiGPdampstoch_[gp]* dampmatvelbackgroundgrad(idim,jdim);
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+idim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode)* (gamma(0) - gamma(1))*r_s(jdim) * vel_rel(jdim);
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode)* (gamma(0) - gamma(1))*r_s(idim) * vel_rel(jdim);

              if (centerline_hermite_)
              {
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode)*jacobiGPdampstoch_[gp]* damp_mat(idim,jdim) / dt;
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+jdim) -= gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode)*jacobiGPdampstoch_[gp]* dampmatvelbackgroundgrad(idim,jdim);
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+idim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode)* (gamma(0) - gamma(1))*r_s(jdim) * vel_rel(jdim);
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode)* (gamma(0) - gamma(1))*r_s(idim) * vel_rel(jdim);

                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode+1)*jacobiGPdampstoch_[gp]* damp_mat(idim,jdim) / dt;
                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+jdim) -= gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode+1)*jacobiGPdampstoch_[gp]* dampmatvelbackgroundgrad(idim,jdim);
                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+idim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode+1)* (gamma(0) - gamma(1))*r_s(jdim) * vel_rel(jdim);
                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode+1)* (gamma(0) - gamma(1))*r_s(idim) * vel_rel(jdim);

                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode+1)*jacobiGPdampstoch_[gp]* damp_mat(idim,jdim) / dt;
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+jdim) -= gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode+1)*jacobiGPdampstoch_[gp]* dampmatvelbackgroundgrad(idim,jdim);
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+idim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)* (gamma(0) - gamma(1))*r_s(jdim) * vel_rel(jdim);
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+jdim) += gausspoints.qwgt[gp]*H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)* (gamma(0) - gamma(1))*r_s(idim) * vel_rel(jdim);
              }
            }
        }
    }

  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodecl,unsigned int vpernode, unsigned int ndim, unsigned int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
void DRT::ELEMENTS::Beam3r::EvaluateStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const LINALG::Matrix<ndim*vpernode*nnodecl,1>& disp_totlag_centerline,
                                              Epetra_SerialDenseMatrix*        stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*        force)//!< element internal force vector
{
  /* only nodes for centerline interpolation are considered here (= first nnodecl nodes of this element);
     each of these nodes holds 3*vpernode translational DoFs AND 3 rotational DoFs */
  const unsigned int dofpernode = 3*vpernode+3;

  // damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  GetDampingCoefficients(gamma);

  /* get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5 */
  Teuchos::RCP<Epetra_MultiVector> randomforces = StatMechParamsInterface().GetRandomForces();

  // my random number vector at current GP
  LINALG::Matrix<ndim,1> randnumvec(true);

  // tangent vector (derivative of beam centerline curve r with respect to arc-length parameter s)
  LINALG::Matrix<ndim,1> r_s(true);

  // stochastic force vector per unit length at current GP
  LINALG::Matrix<ndim,1> f_stoch(true);

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite) shape functions evaluated at the numgp-th GP
   * these shape functions are used for the interpolation of the beam centerline*/
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i(gausspoints.nquad);
  // same for the derivatives
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i_xi(gausspoints.nquad);

  // evaluate all shape function derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAndDerivsAllGPs<nnodecl,vpernode>(gausspoints,H_i,H_i_xi,this->Shape());


  for (int gp=0; gp < gausspoints.nquad; gp++)
  {
    // compute tangent vector t_{\par}=r' at current Gauss point
    Calc_r_s<nnodecl,vpernode,FADordouble>(disp_totlag_centerline, H_i_xi[gp], jacobiGPdampstoch_[gp], r_s);

    // extract random numbers from global vector
    for (unsigned int idim=0; idim<ndim; idim++)
    {
#ifndef BEAM3RCONSTSTOCHFORCE
      randnumvec(idim) = (*randomforces)[gp*randompergauss+idim][LID()];
#else
      randnumvec(idim) = (*randomforces)[idim][LID()];
#endif
    }

    // compute stochastic force vector per unit length at current GP
    f_stoch.Clear();
    for (unsigned int idim=0; idim<ndim; idim++)
      for (unsigned int jdim=0; jdim<ndim; jdim++)
        f_stoch(idim) += (std::sqrt(gamma(1))*(idim==jdim) + (std::sqrt(gamma(0)) - std::sqrt(gamma(1)))*r_s(idim)*r_s(jdim))*randnumvec(jdim);


    if (force != NULL)
    {
      // loop over all nodes
      for (unsigned int inode=0; inode<nnodecl; inode++)
        //loop dimensions with respect to lines
        for (unsigned int idim=0; idim<ndim; idim++)
        {
          (*force)(inode*dofpernode+idim) -= H_i[gp](vpernode*inode)*f_stoch(idim)*std::sqrt(jacobiGPdampstoch_[gp]*gausspoints.qwgt[gp]);
          if (centerline_hermite_)
            (*force)(inode*dofpernode+6+idim) -= H_i[gp](vpernode*inode+1)*f_stoch(idim)*std::sqrt(jacobiGPdampstoch_[gp]*gausspoints.qwgt[gp]);
        }
    }


    if (stiffmatrix != NULL)
    {
      // loop over all nodes used for centerline interpolation
      for (unsigned int inode=0; inode<nnodecl; inode++)
        //loop over all column nodes
        for (unsigned int jnode=0; jnode<nnodecl; jnode++)
        {
          for (unsigned int idim=0; idim<ndim; idim++)
            for (unsigned int jdim=0; jdim<ndim; jdim++)
            {
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+idim) -= H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode)*r_s(jdim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));
              (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+jdim) -= H_i[gp](vpernode*inode)*H_i_xi[gp](vpernode*jnode)*r_s(idim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));

              if (centerline_hermite_)
              {
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+idim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode)*r_s(jdim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+jdim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode)*r_s(idim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));

                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+idim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)*r_s(jdim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));
                (*stiffmatrix)(inode*dofpernode+idim,jnode*dofpernode+6+jdim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)*r_s(idim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));

                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+idim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)*r_s(jdim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));
                (*stiffmatrix)(inode*dofpernode+6+idim,jnode*dofpernode+6+jdim) -= H_i[gp](vpernode*inode+1)*H_i_xi[gp](vpernode*jnode+1)*r_s(idim)*randnumvec(jdim)*std::sqrt(gausspoints.qwgt[gp]/ jacobiGPdampstoch_[gp])*(std::sqrt(gamma(0)) - std::sqrt(gamma(1)));
              }
            }
        }
    }

  } // end: loop GPs

}
