/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.cpp

\brief optimizer of the topology optimization

\maintainer Martin Winklmaier

\level 3
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "opti_GCMMA.H"
#include "opti_resulttest.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/optimization_density.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
TOPOPT::Optimizer::Optimizer(Teuchos::RCP<DRT::Discretization> optidis,
    Teuchos::RCP<DRT::Discretization> fluiddis, const Teuchos::ParameterList& params,
    Teuchos::RCP<IO::DiscretizationWriter>& output)
    : optidis_(optidis),
      fluiddis_(fluiddis),
      params_(params),
      gradienttype_(
          DRT::INPUT::IntegralValue<INPAR::TOPOPT::GradientType>(params_, "GRADIENT_TYPE"))
{
  const Teuchos::ParameterList& optimizer_params = params_.sublist("TOPOLOGY OPTIMIZER");

  // fluid fields for optimization
  fluidvel_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_Vector>>);
  adjointvel_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_Vector>>);

  // topology density fields
  dens_ = Teuchos::rcp(new Epetra_Vector(*ColMap(), false));


  // value of the objective function
  obj_ = 0.0;
  // gradient of the objective function
  obj_der_ = Teuchos::rcp(new Epetra_Vector(*RowMap()));

  /// number of constraints
  switch (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(optimizer_params, "TESTCASE"))
  {
    case INPAR::TOPOPT::optitest_channel:
    case INPAR::TOPOPT::optitest_channel_with_step:
    case INPAR::TOPOPT::optitest_cornerflow:
    case INPAR::TOPOPT::optitest_lin_poro:
    case INPAR::TOPOPT::optitest_quad_poro:
    case INPAR::TOPOPT::optitest_cub_poro:
    case INPAR::TOPOPT::optitest_no:
    case INPAR::TOPOPT::optitest_workflow_without_fluiddata:
    case INPAR::TOPOPT::optitest_snake_one_constr:
    {
      num_constr_ = 1;
      break;
    }
    case INPAR::TOPOPT::optitest_snake_multiple_constr:
    {
      num_constr_ = 5;
      break;
    }
    default:
    {
      dserror("unknown optimization case");
      break;
    }
  }

  // value of the constraint(s);
  constr_ = Teuchos::rcp(new Epetra_SerialDenseVector(num_constr_));
  // gradient of the constraint(s)
  constr_der_ = Teuchos::rcp(new Epetra_MultiVector(*RowMap(), num_constr_, false));

  int numFDPoints = 0;  // number of additional points for finite differences
  switch (gradienttype_)
  {
    case INPAR::TOPOPT::gradientByAdjoints:
    {
      numFDPoints = 1;
      break;
    }  // dummy
    case INPAR::TOPOPT::gradientByFD1:
    {
      numFDPoints = 1;
      break;
    }
    case INPAR::TOPOPT::gradientByFD2:
    {
      numFDPoints = 2;
      break;
    }
    default:
    {
      dserror("unknown type of gradient computation");
      break;
    }
  }
  objective_FD_ = Teuchos::rcp(new Epetra_SerialDenseVector(numFDPoints));
  constraints_FD_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(num_constr_, numFDPoints));

  // set initial density field if present
  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialDensityField>(
                             optimizer_params, "INITIALFIELD"),
      optimizer_params.get<int>("INITFUNCNO"));

  // topology density fields
  Teuchos::RCP<Epetra_Vector> dens = Teuchos::rcp(new Epetra_Vector(*RowMap(), false));
  LINALG::Export(*dens_, *dens);

  if (1 == 1)  // use this when different optimizers are present
  {
    optimizer_ = Teuchos::rcp(new OPTI::GCMMA(
        optidis_, params_, dens, num_constr_, Teuchos::null, Teuchos::null, output));
  }

  // write output using the derived optimizer -> optimizer must be present here
  Output();

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeValues(double& objective, Epetra_SerialDenseVector& constraints)
{
  // initialize local values processor-wise
  double loc_obj = 0.0;
  Epetra_SerialDenseVector loc_constr = Epetra_SerialDenseVector(num_constr_);

  const bool doAdjoint = false;

  // check if all data is present
  CheckData(doAdjoint);

  // Transform fluid field to colmap
  TransformFlowFields(doAdjoint, true);

  Teuchos::ParameterList params;

  params.set("action", "compute_values");

  params.set("objective_value", &loc_obj);
  params.set("constraint_values", &loc_constr);

  params.set("fluidvel", fluidvel_);
  params.set("fluiddis", fluiddis_);

  params.set("topopt_density", Teuchos::rcp_const_cast<const Epetra_Vector>(dens_));

  optidis_->ClearState();
  optidis_->Evaluate(params, Teuchos::null, Teuchos::null);

  optidis_->ClearState();

  // communicate local values over processors
  optidis_->Comm().SumAll(&loc_obj, &objective, 1);
  optidis_->Comm().SumAll(loc_constr.Values(), constraints.Values(), num_constr_);

  // Re-Transform fluid field to rowmap
  TransformFlowFields(doAdjoint, false);  // TODO is this required (back-mapping + 2nd argument)
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradients(
    Teuchos::RCP<Epetra_Vector> obj_der, Teuchos::RCP<Epetra_MultiVector> constr_der)
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective gradient");
  obj_der->PutScalar(0.0);
  constr_der->PutScalar(0.0);

  bool doAdjoint = true;
  if (gradienttype_ == INPAR::TOPOPT::gradientByAdjoints)
    doAdjoint = true;
  else if (gradienttype_ == INPAR::TOPOPT::gradientByFD1)
    doAdjoint = false;
  else
    dserror("unknown type of gradient computation");

  // check if all data is present
  CheckData(doAdjoint);

  // Transform fluid field to colmap
  TransformFlowFields(doAdjoint, true);

  Teuchos::ParameterList params;

  params.set("action", "compute_gradients");

  params.set("constraints_derivations", constr_der);
  params.set("objective_derivations", obj_der);

  params.set("fluidvel", fluidvel_);
  params.set("adjointvel", adjointvel_);
  params.set("fluiddis", fluiddis_);

  //  for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=fluidvel_->begin();
  //      i!=fluidvel_->end();i++)
  //    std::cout << "in optimizer at gradient computation fluidvel of step " << i->first << " is "
  //    << *i->second << std::endl;
  //  for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=adjointvel_->begin();
  //      i!=adjointvel_->end();i++)
  //    std::cout << "in optimizer at gradient computation adjointvel of step " << i->first << " is
  //    " << *i->second << std::endl;

  params.set("topopt_density", Teuchos::rcp_const_cast<const Epetra_Vector>(dens_));

  optidis_->ClearState();
  optidis_->Evaluate(params, Teuchos::null, Teuchos::null);
  optidis_->ClearState();

  // Re-Transform fluid field to rowmap
  TransformFlowFields(doAdjoint, false);  // TODO is this required (back-mapping + 2nd argument)
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradientDirectionForFD(
    const double fac, const int GID, const int index, const int numFDPoints)
{
  double& objective = (*objective_FD_)(index);
  Epetra_SerialDenseVector constraints(View, (*constraints_FD_)[index], num_constr_);

  ComputeValues(objective, constraints);

  if ((index == numFDPoints - 1) and (dens_->Map().MyGID(GID)))
  {
    switch (gradienttype_)
    {
      case INPAR::TOPOPT::gradientByFD1:
      {
        obj_der_->ReplaceGlobalValue(GID, 0, ((*objective_FD_)(0) - obj_) / fac);

        for (int i = 0; i < constr_->Length(); i++)
          constr_der_->ReplaceGlobalValue(GID, i, ((*constraints_FD_)(i, 0) - (*constr_)(i)) / fac);

        break;
      }
      case INPAR::TOPOPT::gradientByFD2:
      {
        obj_der_->ReplaceGlobalValue(GID, 0,
            ((*objective_FD_)(0) - (*objective_FD_)(1)) /
                (-2.0 * fac));  // second, here used direction is -c, so factor -1 to handle this

        for (int i = 0; i < constr_->Length(); i++)
          constr_der_->ReplaceGlobalValue(
              GID, i, ((*constraints_FD_)(i, 0) - (*constraints_FD_)(i, 1)) / (-2.0 * fac));

        break;
      }
      case INPAR::TOPOPT::gradientByAdjoints:
      default:
      {
        dserror("wrong type of gradient computation");
        break;
      }
    }
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::Output() { optimizer_->Output(); }



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::SetInitialDensityField(
    const INPAR::TOPOPT::InitialDensityField initfield, const int startfuncno)
{
  if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
      INPAR::TOPOPT::dens_node_based)
  {
    // loop all nodes on the processor
    for (int lnodeid = 0; lnodeid < optidis_->NumMyColNodes(); lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = optidis_->lColNode(lnodeid);
      const double* coords = lnode->X();

      switch (initfield)
      {
        case INPAR::TOPOPT::initdensfield_zero_field:
        {
          dens_->PutScalar(0.0);
          break;
        }
        case INPAR::TOPOPT::initdensfield_field_by_function:
        {
          // evaluate component k of spatial function
          double initialval = DRT::Problem::Instance()
                                  ->Funct(startfuncno - 1)
                                  .Evaluate(0, lnode->X(), 0.0);       // scalar
          int err = dens_->ReplaceMyValues(1, &initialval, &lnodeid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        case INPAR::TOPOPT::initdensfield_channelflow0:
        case INPAR::TOPOPT::initdensfield_channelflow05:
        case INPAR::TOPOPT::initdensfield_channelflow1:
        {
          // case with border value 0
          double initialval = 0.0;
          if (abs(coords[1]) < 0.2 - 1.0e-12)
            initialval = 1.0;
          else
            initialval = 0.0;

          // case with border value 0.5
          if ((initfield == INPAR::TOPOPT::initdensfield_channelflow05) and
              (abs(abs(coords[1]) - 0.2) < 1.0e-12))
            initialval = 0.5;

          // case with border value 1.0
          if ((initfield == INPAR::TOPOPT::initdensfield_channelflow1) and
              (abs(abs(coords[1]) - 0.2) < 1.0e-12))
            initialval = 1.0;

          // evaluate component k of spatial function
          int err = dens_->ReplaceMyValues(1, &initialval, &lnodeid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        case INPAR::TOPOPT::initdensfield_channelstepflow0:
        case INPAR::TOPOPT::initdensfield_channelstepflow05:
        case INPAR::TOPOPT::initdensfield_channelstepflow1:
        {
          // case with border value 0
          double initialval = 1.0;  // default case with density one
          if ((coords[1] < 0.4 + 1.0e-12) and (coords[0] > 1.5 - 1.0e-12) and
              (coords[0] < 1.9 + 1.0e-12))
            initialval = 0.0;  // edge area with val 0

          if (((coords[1] > 0.4 - 1.0e-12) and (coords[1] < 0.4 + 1.0e-12) and
                  (coords[0] > 1.5 - 1.0e-12) and
                  (coords[0] < 1.9 + 1.0e-12)) or  // upper boundary line of step
              ((coords[0] > 1.5 - 1.0e-12) and (coords[0] < 1.5 + 1.0e-12) and
                  (coords[1] < 0.4 + 1.0e-12)) or  // left boundary line of step
              ((coords[0] > 1.9 - 1.0e-12) and (coords[0] < 1.9 + 1.0e-12) and
                  (coords[1] < 0.4 + 1.0e-12)))  // right boundary line of step
          {
            if (initfield == INPAR::TOPOPT::initdensfield_channelstepflow05)
              initialval = 0.5;  // case with border value 0.5
            if (initfield == INPAR::TOPOPT::initdensfield_channelstepflow1)
              initialval = 1.0;  // case with border value 1.0
          }

          // evaluate component k of spatial function
          int err = dens_->ReplaceMyValues(1, &initialval, &lnodeid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        default:
        {
          dserror("unknown initial field");
          break;
        }
      }
    }
  }
  else if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
           INPAR::TOPOPT::dens_ele_based)
  {
    const int nsd = DRT::Problem::Instance()->NDim();

    // loop all nodes on the processor
    for (int leleid = 0; leleid < optidis_->NumMyColElements(); leleid++)
    {
      const DRT::Element* ele = optidis_->lColElement(leleid);

      // get "approximation of element center"
      Epetra_SerialDenseMatrix coords(nsd, 1);
      const DRT::Node* const* nodes = ele->Nodes();
      for (int i = 0; i < ele->NumNode(); i++)
      {
        Epetra_SerialDenseMatrix nodecoords(nsd, 1);
        for (int idim = 0; idim < nsd; idim++) nodecoords(idim, 0) = nodes[i]->X()[idim];

        coords += nodecoords;
      }
      coords.Scale(1.0 / ele->NumNode());

      switch (initfield)
      {
        case INPAR::TOPOPT::initdensfield_zero_field:
        {
          dens_->PutScalar(0.0);
          break;
        }
        case INPAR::TOPOPT::initdensfield_field_by_function:
        {
          // evaluate component k of spatial function
          double initialval = DRT::Problem::Instance()
                                  ->Funct(startfuncno - 1)
                                  .Evaluate(0, coords.A(), 0.0);      // scalar
          int err = dens_->ReplaceMyValues(1, &initialval, &leleid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        case INPAR::TOPOPT::initdensfield_channelflow0:
        case INPAR::TOPOPT::initdensfield_channelflow05:
        case INPAR::TOPOPT::initdensfield_channelflow1:
        {
          // case with border value 0
          double initialval = 0.0;
          if (abs(coords(1, 0)) < 0.2 - 1.0e-12)
            initialval = 1.0;
          else
            initialval = 0.0;

          // case with border value 0.5
          if ((initfield == INPAR::TOPOPT::initdensfield_channelflow05) and
              (abs(abs(coords(1, 0)) - 0.2) < 1.0e-12))
            initialval = 0.5;

          // case with border value 1.0
          if ((initfield == INPAR::TOPOPT::initdensfield_channelflow1) and
              (abs(abs(coords(1, 0)) - 0.2) < 1.0e-12))
            initialval = 1.0;

          // evaluate component k of spatial function
          int err = dens_->ReplaceMyValues(1, &initialval, &leleid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        case INPAR::TOPOPT::initdensfield_channelstepflow0:
        case INPAR::TOPOPT::initdensfield_channelstepflow05:
        case INPAR::TOPOPT::initdensfield_channelstepflow1:
        {
          // case with border value 0
          double initialval = 1.0;  // default case with density one
          if ((coords(1, 0) < 0.4 + 1.0e-12) and (coords(0, 0) > 1.5 - 1.0e-12) and
              (coords(0, 0) < 1.9 + 1.0e-12))
            initialval = 0.0;  // edge area with val 0

          if (((coords(1, 0) > 0.4 - 1.0e-12) and (coords(1, 0) < 0.4 + 1.0e-12) and
                  (coords(0, 0) > 1.5 - 1.0e-12) and
                  (coords(0, 0) < 1.9 + 1.0e-12)) or  // upper boundary line of step
              ((coords(0, 0) > 1.5 - 1.0e-12) and (coords(0, 0) < 1.5 + 1.0e-12) and
                  (coords(1, 0) < 0.4 + 1.0e-12)) or  // left boundary line of step
              ((coords(0, 0) > 1.9 - 1.0e-12) and (coords(0, 0) < 1.9 + 1.0e-12) and
                  (coords(1, 0) < 0.4 + 1.0e-12)))  // right boundary line of step
          {
            if (initfield == INPAR::TOPOPT::initdensfield_channelstepflow05)
              initialval = 0.5;  // case with border value 0.5
            if (initfield == INPAR::TOPOPT::initdensfield_channelstepflow1)
              initialval = 1.0;  // case with border value 1.0
          }

          // evaluate component k of spatial function
          int err = dens_->ReplaceMyValues(1, &initialval, &leleid);  // lnodeid = ldofid
          if (err != 0) dserror("dof not on proc");
          break;
        }
        default:
        {
          dserror("unknown initial field");
          break;
        }
      }
    }
  }
  else
    dserror("not implemented type of density function");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportFlowParams(Teuchos::RCP<Teuchos::ParameterList>& fluidParams)
{
  // save the fluid parameter
  fluidParams_ = fluidParams;

  // set the general parameter on element level
  Teuchos::ParameterList opti_ele_params;
  opti_ele_params.set("action", "set_general_optimization_parameter");

  // set material parameters of fluid and optimization material
  {
    // get the material
    const int nummat = DRT::Problem::Instance()->Materials()->Num();
    for (int id = 1; id - 1 < nummat; ++id)
    {
      Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

      if (imat == Teuchos::null)
        dserror("Could not find material Id %d", id);
      else
      {
        switch (imat->Type())
        {
          case INPAR::MAT::m_fluid:
          {
            const MAT::PAR::Parameter* matparam = imat->Parameter();
            const MAT::PAR::NewtonianFluid* mat =
                static_cast<const MAT::PAR::NewtonianFluid*>(matparam);

            opti_ele_params.set("density", mat->density_);
            opti_ele_params.set("viscosity", mat->viscosity_);
            break;
          }
          case INPAR::MAT::m_opti_dens:
          {
            const MAT::PAR::Parameter* matparam = imat->Parameter();
            const MAT::PAR::TopOptDens* mat = static_cast<const MAT::PAR::TopOptDens*>(matparam);

            opti_ele_params.set("MIN_PORO", mat->PoroBdDown());
            opti_ele_params.set("MAX_PORO", mat->PoroBdUp());

            const INPAR::TOPOPT::OptiCase testcase =
                (INPAR::TOPOPT::OptiCase)(DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(
                    params_.sublist("TOPOLOGY OPTIMIZER"), "TESTCASE"));
            switch (testcase)
            {
              case INPAR::TOPOPT::optitest_channel:
              case INPAR::TOPOPT::optitest_channel_with_step:
              case INPAR::TOPOPT::optitest_cornerflow:
              case INPAR::TOPOPT::optitest_lin_poro:
              case INPAR::TOPOPT::optitest_quad_poro:
              case INPAR::TOPOPT::optitest_cub_poro:
              {
                opti_ele_params.set("SMEAR_FAC", (double)(-(int)testcase));
                break;
              }
              default:
              {
                opti_ele_params.set("SMEAR_FAC", mat->SmearFac());
                break;
              }
            }

            break;
          }
          default:
          {
            dserror("unknown material %s", imat->Name().c_str());
            break;
          }
        }
      }
    }
  }

  // flow parameter
  opti_ele_params.set("timealgo", fluidParams_->get<int>("time int algo"));
  // parameter for stabilization
  opti_ele_params.sublist("RESIDUAL-BASED STABILIZATION") =
      fluidParams_->sublist("RESIDUAL-BASED STABILIZATION");
  opti_ele_params.set("dt", fluidParams_->get<double>("time step size"));
  opti_ele_params.set("maxtimesteps", fluidParams_->get<int>("max number timesteps"));
  opti_ele_params.set("theta", fluidParams_->get<double>("theta"));

  // objective parameter
  opti_ele_params.set("dissipation",
      fluidParams_->get<INPAR::TOPOPT::ObjectiveDissipation>("OBJECTIVE_DISSIPATION"));
  opti_ele_params.set("pres_drop", fluidParams_->get<bool>("OBJECTIVE_PRESSURE_DROP"));

  opti_ele_params.set("dissipation_fac", fluidParams_->get<double>("DISSIPATION_FAC"));
  opti_ele_params.set("pres_drop_fac", fluidParams_->get<double>("PRESSURE_DROP_FAC"));

  opti_ele_params.set("theta_obj", params_.sublist("TOPOLOGY OPTIMIZER").get<double>("THETA"));

  opti_ele_params.set(
      "vol_bd", params_.sublist("TOPOLOGY OPTIMIZER").get<double>("VOLUME_BOUNDARY"));

  opti_ele_params.set<INPAR::TOPOPT::DensityField>(
      "dens_type", DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE"));

  opti_ele_params.set<INPAR::TOPOPT::OptiCase>(
      "opti_case", DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(
                       params_.sublist("TOPOLOGY OPTIMIZER"), "TESTCASE"));

  optidis_->Evaluate(opti_ele_params, Teuchos::null, Teuchos::null);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportFluidData(Teuchos::RCP<Epetra_Vector> vel, int step)
{
  if (fluidvel_->find(step) != fluidvel_->end())
    *(fluidvel_->find(step)->second) = *vel;
  else
  {
    // a copy is required here, otherwise the values are changed within the fluid time integration
    Teuchos::RCP<Epetra_Vector> new_vel = Teuchos::rcp(new Epetra_Vector(*vel));  // copy
    fluidvel_->insert(std::pair<int, Teuchos::RCP<Epetra_Vector>>(step, new_vel));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportAdjointFluidData(Teuchos::RCP<Epetra_Vector> vel, int step)
{
  if (adjointvel_->find(step) != adjointvel_->end())
    *(adjointvel_->find(step)->second) = *vel;
  else
  {
    // a copy is required here, otherwise the values are changed within the adjoint time integration
    Teuchos::RCP<Epetra_Vector> new_vel = Teuchos::rcp(new Epetra_Vector(*vel));
    adjointvel_->insert(std::pair<int, Teuchos::RCP<Epetra_Vector>>(step, new_vel));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::Iterate(bool& doGradient)
{
  Teuchos::RCP<Epetra_Vector> dens = Teuchos::null;

  if (doGradient)
  {
    dens = optimizer_->Iterate(obj_, obj_der_, constr_, constr_der_);
  }
  else
  {
    dens = optimizer_->Iterate(obj_, Teuchos::null, constr_, Teuchos::null);
  }

  LINALG::Export(*dens, *dens_);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::FinishIteration(bool& doGradient)
{
  return optimizer_->FinishIteration(obj_, constr_, doGradient);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::Converged(bool& doGradient)
{
  bool converged = false;

  // dont use gradients if not needed (safety)
  if (doGradient)
  {
    converged = optimizer_->Converged(obj_, obj_der_, constr_, constr_der_);
  }
  else
  {
    converged = optimizer_->Converged(obj_, Teuchos::null, constr_, Teuchos::null);
  }

  return converged;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::CheckData(const bool doAdjoint)
{
  size_t num_sols = fluidParams_->get<int>("max number timesteps");

  if (fluidParams_->get<int>("time int algo") == INPAR::FLUID::timeint_stationary)
  {
    if (num_sols != fluidvel_->size())
      dserror("fluid field and time step numbers do not fit: n_f = %i, n_t = %i", fluidvel_->size(),
          num_sols);

    if (doAdjoint == true)
    {
      if (num_sols != adjointvel_->size())
        dserror("adjoint field and time step numbers do not fit: n_a = %i, n_t = %i",
            adjointvel_->size(), num_sols);
    }
  }
  else
  {
    // n timesteps
    // -> solutions at time t^0,t^1,...,t^n
    // -> n+1 solutions
    size_t num_sols = fluidParams_->get<int>("max number timesteps") + 1;

    if (num_sols != fluidvel_->size())
      dserror("fluid field and time step numbers do not fit: n_f = %i, n_t = %i", fluidvel_->size(),
          num_sols);

    if (doAdjoint == true)
    {
      if (num_sols - 1 != adjointvel_->size())  // no adjoint solution of time step 0 required
        dserror("adjoint field and time step numbers do not fit: n_a = %i, n_t = %i",
            adjointvel_->size(), num_sols);
    }
  }

  return true;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::TransformFlowFields(const bool doAdjoint, const bool rowToCol)
{
  const Epetra_Map* targetmap = NULL;
  if (rowToCol)
    targetmap = fluiddis_->DofColMap();
  else
    targetmap = fluiddis_->DofRowMap();

  for (std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator i = fluidvel_->begin();
       i != fluidvel_->end(); i++)
  {
    // if maps have been mapped before, dont change them
    if (i->second->Map().PointSameAs(*targetmap)) continue;

    Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*targetmap, false));

    // export vector to target map
    LINALG::Export(*i->second, *vec);
    i->second = vec;
  }

  // if maps have been mapped before, dont change them
  if (doAdjoint == true)
  {
    for (std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator i = adjointvel_->begin();
         i != adjointvel_->end(); i++)
    {
      if (i->second->Map().PointSameAs(*targetmap)) continue;

      Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*targetmap, false));

      // export vector to target map
      LINALG::Export(*i->second, *vec);
      i->second = vec;
    }
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::RowMap()
{
  if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
      INPAR::TOPOPT::dens_node_based)
    return optidis_->NodeRowMap();
  else if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
           INPAR::TOPOPT::dens_ele_based)
    return optidis_->ElementRowMap();
  else
    dserror("not implemented type of density field");

  return NULL;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::ColMap()
{
  if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
      INPAR::TOPOPT::dens_node_based)
    return optidis_->NodeColMap();
  else if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params_, "DENS_TYPE") ==
           INPAR::TOPOPT::dens_ele_based)
    return optidis_->ElementColMap();
  else
    dserror("not implemented type of density field");

  return NULL;
}


void TOPOPT::Optimizer::UpdateOptimizationParameter()
{
  if (DRT::INPUT::IntegralValue<bool>(params_.sublist("TOPOLOGY OPTIMIZER"), "update_smooth") ==
      true)
  {
    // get the material
    const int nummat = DRT::Problem::Instance()->Materials()->Num();
    for (int id = 1; id - 1 < nummat; ++id)
    {
      Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

      if (imat == Teuchos::null)
        dserror("Could not find material Id %d", id);
      else
      {
        switch (imat->Type())
        {
          case INPAR::MAT::m_opti_dens:
          {
            MAT::PAR::Parameter* matparam = imat->Parameter();
            MAT::PAR::TopOptDens* mat = static_cast<MAT::PAR::TopOptDens*>(matparam);

            const INPAR::TOPOPT::OptiCase testcase =
                (INPAR::TOPOPT::OptiCase)(DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(
                    params_.sublist("TOPOLOGY OPTIMIZER"), "TESTCASE"));
            switch (testcase)
            {
              case INPAR::TOPOPT::optitest_channel:
              case INPAR::TOPOPT::optitest_channel_with_step:
              case INPAR::TOPOPT::optitest_cornerflow:
              case INPAR::TOPOPT::optitest_lin_poro:
              case INPAR::TOPOPT::optitest_quad_poro:
              case INPAR::TOPOPT::optitest_cub_poro:
                break;
              default:
              {
                const int stepnum =
                    params_.sublist("TOPOLOGY OPTIMIZER").get<int>("update_smooth_every_iter");
                const double fac =
                    params_.sublist("TOPOLOGY OPTIMIZER").get<double>("update_smooth_fac");
                if ((optimizer_->OuterIter() % stepnum == 0) and (optimizer_->Iter() > 0) and
                    (optimizer_->InnerIter() == 1))
                {
                  // Update material
                  mat->UpdateSmearFac(fac * mat->SmearFac());

                  // update the according parameter on element level
                  Teuchos::ParameterList opti_ele_params;
                  opti_ele_params.set("action", "update_general_optimization_parameter");
                  opti_ele_params.set("SMEAR_FAC", mat->SmearFac());
                  optidis_->Evaluate(opti_ele_params, Teuchos::null, Teuchos::null);
                }

                if (mat->SmearFac() > 2000000001)
                  dserror("optimization should already be converged");
                break;
              }
            }
            break;
          }
          case INPAR::MAT::m_fluid:
            break;
          default:
          {
            dserror("unknown material %s", imat->Name().c_str());
            break;
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> TOPOPT::Optimizer::CreateFieldTest()
{
  return Teuchos::rcp(new OPTI::OptiResultTest(*optimizer_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::AdoptDensityForFD(const double value, const int GID)
{
  if (dens_->Map().MyGID(GID)) dens_->SumIntoGlobalValue(GID, 0, value);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int TOPOPT::Optimizer::Iter() const { return optimizer_->Iter(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ReadRestart(const int step)
{
  // required here for objective derivation
  Teuchos::RCP<IO::DiscretizationReader> reader = optimizer_->ReadRestart(step);

  Teuchos::RCP<Epetra_Vector> dens = optimizer_->X();
  LINALG::Export(*dens, *dens_);

  *obj_der_ = *optimizer_->ObjDeriv();
  *constr_der_ = *optimizer_->ConstrDeriv();

  return;
}
