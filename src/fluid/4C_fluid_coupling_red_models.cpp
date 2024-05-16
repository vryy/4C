/*-----------------------------------------------------------*/
/*! \file

\brief evaluation of 3D/redD coupled vascular bc


\level 3
*/
/*-----------------------------------------------------------*/

// #ifdef D_COUPLED_ARTNET

#include "4C_fluid_coupling_red_models.hpp"

#include "4C_fluid_ele_action.hpp"

#include <math.h>
#include <stdio.h>

FOUR_C_NAMESPACE_OPEN
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 11/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidCouplingWrapperBase::FluidCouplingWrapperBase(
    Teuchos::RCP<DRT::Discretization> dis_3D, Teuchos::RCP<DRT::Discretization> dis_redD,
    //                                                         Teuchos::RCP<red_D_time_int>
    //                                                         RedD_Time_integ,
    IO::DiscretizationWriter& output, double dt_3D, double dt_redD)
    :  // call constructor for "nontrivial" objects
      discret3_d_(dis_3D),
      discret_red_d_(dis_redD),
      //  reduced_D_time_integ_(RedD_Time_integ),
      output_(output)
{
  // ---------------------------------------------------------------------
  // Read in the time step
  // ---------------------------------------------------------------------
  dt_f3_ = dt_3D;
  dt_rm_ = dt_redD;

  // check whether the "reduceD time step" is a deviser  of "3D time step"
  //    This is essential since it is not advisable to change the time
  //    step of the reduced-D problem during the simulation
  int quotient = int(floor(dt_f3_ / dt_rm_));
  double remainder = dt_f3_ - double(quotient) * dt_rm_;
  if (remainder != 0.0)
  {
    FOUR_C_THROW("\"Fluid 3D\" must have a time step multiple of that of \"reduced-D\" problem");
  }

  // ---------------------------------------------------------------------
  // Read in all conditions
  // ---------------------------------------------------------------------
  std::vector<CORE::Conditions::Condition*> couplingcond;
  discret3_d_->GetCondition("Art_3D_redD_CouplingCond", couplingcond);

  std::vector<CORE::Conditions::Condition*> couplingcond2;
  discret_red_d_->GetCondition("Art_redD_3D_CouplingCond", couplingcond2);

  // the number of lines of coupling boundary conditions found in the input
  // note that some of these lines could belong to the same physical condition
  // which is then marked by the same 'ConditionID'
  // The corresponding resorting of result values has to be done later
  unsigned int numcondlines = couplingcond.size();

  // ---------------------------------------------------------------------
  // Check whether the 2nd conditions are similar in number as the first
  // ---------------------------------------------------------------------
  if (numcondlines != couplingcond2.size())
  {
    FOUR_C_THROW(
        "coupled problem beween reduced-D and 3D must have equal number of condition on both sides "
        "of the discretization boundaries");
  }

  if (numcondlines > 0)  // if there is at least one coupling bc
  {
    map3_dnp_ = Teuchos::rcp(new std::map<std::string, double>);
    map3_dn_ = Teuchos::rcp(new std::map<std::string, double>);
    map_red_dnp_ = Teuchos::rcp(new std::map<std::string, double>);
    map_red_dn_ = Teuchos::rcp(new std::map<std::string, double>);
    // -------------------------------------------------------------------
    // get the maximum allowable number of iterations at the boundary
    // which should be the same!
    // -------------------------------------------------------------------

    int N_iter = couplingcond[0]->parameters().Get<int>("MaximumIterations");

    // -------------------------------------------------------------------
    // make sure that each coupling has two conditions of same ID
    // -------------------------------------------------------------------
    for (unsigned int i = 0; i < numcondlines; i++)
    {
      bool CondIsFine = false;
      int condid = couplingcond[i]->parameters().Get<int>("ConditionID");
      for (unsigned int j = 0; j < numcondlines; j++)
      {
        int condid2 = couplingcond2[j]->parameters().Get<int>("ConditionID");
        if (condid2 == condid)
        {
          CondIsFine = true;
          break;
        }
      }
      if (!CondIsFine)
      {
        FOUR_C_THROW(
            "[3D/Reduced-D COUPLING] condition [%d] is defined only on the 3D side", condid);
      }
    }

    // -------------------------------------------------------------------
    // now care for the fact that there could be more than one input line
    // belonging to the same coupling boundary condition
    // -------------------------------------------------------------------
    for (unsigned int i = 0; i < numcondlines; i++)
    {
      int condid = couplingcond[i]->parameters().Get<int>("ConditionID");

      // -----------------------------------------------------------------
      // Find the second condition number which has the same id as the
      // first
      // -----------------------------------------------------------------
      unsigned int j = 0;
      for (j = 0; j < numcondlines; j++)
      {
        if (condid == couplingcond2[j]->parameters().Get<int>("ConditionID"))
        {
          break;
        }
      }
      int thisN_iter = couplingcond[i]->parameters().Get<int>("MaximumIterations");
      if (thisN_iter != N_iter)
        FOUR_C_THROW(
            "all maximum number of iterations on the coupling boundary between 3-D and reduced-D "
            "boundary should be the same!!!");


      // ------------------------------------------------------------------
      // allocate the coupling bc class members for every case
      // ------------------------------------------------------------------
      Teuchos::RCP<FluidCouplingBc> couplingbc = Teuchos::rcp(
          new FluidCouplingBc(discret3_d_, discret_red_d_, output_, dt_f3_, dt_rm_, condid, i, j));

      // -----------------------------------------------------------------
      // sort coupling bc's in map and test, if one condition ID appears
      // more than once. Currently this case is forbidden.
      // -----------------------------------------------------------------
      bool inserted = coup_map3_d_.insert(std::make_pair(condid, couplingbc)).second;
      if (!inserted)
        FOUR_C_THROW(
            "There are more than one 3D-to-OneD coupling condition lines with the same ID. This "
            "can not yet be handled.");
    }  // end loop over condition lines from input

    // -------------------------------------------------------------------
    // Fill the coupled variables boundary condition
    //    +-----------------------------+--------------------+
    //    |  BoundaryVariable           |  BoundaryValue     |
    //    +-----------------------------+--------------------+
    //    |         pressure1           |       300.25       |
    //    |         .......             .                    .
    //    |         flow10              .        20.43       |
    //    +-----------------------------+--------------------+
    // -------------------------------------------------------------------

    // -------------------------------------------------------------------
    // Fill the 3D coupling variable
    // -------------------------------------------------------------------
    for (unsigned int i = 0; i < numcondlines; i++)
    {
      // Get condition ID
      int id = couplingcond[i]->parameters().Get<int>("ConditionID");

      // Get returned coupling variable
      std::string variable = ((couplingcond[i])->parameters().Get<std::string>("ReturnedVariable"));

      // Build a new std::string from [coupling Variable name][Condition Id]
      std::stringstream VariableWithId;
      VariableWithId << variable << "_" << id;
      double value = 0.0;

      // Build the map

      map3_dnp_->insert(std::make_pair(VariableWithId.str(), value));
      map3_dn_->insert(std::make_pair(VariableWithId.str(), value));
    }

    // ------------------------------------------------------------------
    // Fill the reduced-D coupling variable
    // ------------------------------------------------------------------
    for (unsigned int i = 0; i < numcondlines; i++)
    {
      // Get condition ID
      int id = couplingcond2[i]->parameters().Get<int>("ConditionID");

      // Get returned coupling variable
      std::string variable =
          ((couplingcond2[i])->parameters().Get<std::string>("ReturnedVariable"));

      // Build a new std::string from [coupling Variable name][Condition Id]
      std::stringstream VariableWithId;
      VariableWithId << variable << "_" << id;
      double value = 0.0;

      // Build the map
      map_red_dnp_->insert(std::make_pair(VariableWithId.str(), value));
      map_red_dn_->insert(std::make_pair(VariableWithId.str(), value));
    }

  }  // end if there were conditions

  return;
}  // end FluidCouplingWrapperBase


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap flow rate calculation                             ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase ::FlowRateCalculation(double time, double dta)
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    mapiter->second->FluidCouplingBc ::FlowRateCalculation(time, dta, mapiter->first);
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap flow rate calculation                             ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::PressureCalculation(double time, double dta)
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    mapiter->second->FluidCouplingBc::PressureCalculation(time, dta, mapiter->first);
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap outflow boundary pressure application             ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::ApplyBoundaryConditions(
    double time, double dta, double theta)
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  // ---------------------------------------------------------------------
  // Read in all conditions
  // ---------------------------------------------------------------------

  // Read in the 3D coupling conditions
  std::vector<CORE::Conditions::Condition*> conds3D;
  discret3_d_->GetCondition("Art_3D_redD_CouplingCond", conds3D);

  // Read in the reduced-D coupling conditions
  std::vector<CORE::Conditions::Condition*> conds_redD;
  discret_red_d_->GetCondition("Art_redD_3D_CouplingCond", conds_redD);

  int condID;
  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    // Get condition ID
    condID = mapiter->first;

    //----------------------------------------------------------------
    // find the parameters that must be calculated and returned by
    // the 3D problem
    //----------------------------------------------------------------
    for (auto& i : conds3D)
    {
      if (i->parameters().Get<int>("ConditionID") == condID)
      {
        // get returned value name from 3D boundary
        std::string variable_str = (i->parameters().Get<std::string>("ReturnedVariable"));

        // concatenate the variable name with the variable id
        std::stringstream CouplingVariable;
        CouplingVariable << variable_str << "_" << condID;

        if (variable_str == "flow")
        {
          auto itr = map3_dnp_->find(CouplingVariable.str());
          if (itr == map3_dnp_->end())
          {
            FOUR_C_THROW("[3D/Reduced-D COUPLING] 3D map has no variable %s for condition [%d]",
                variable_str.c_str(), condID);
          }
          (*map3_dnp_)[CouplingVariable.str()] =
              mapiter->second->FluidCouplingBc::FlowRateCalculation(time, dta, condID);
        }
        else if (variable_str == "pressure")
        {
          std::map<std::string, double>::iterator itr = map3_dnp_->find(CouplingVariable.str());
          if (itr == map3_dnp_->end())
          {
            FOUR_C_THROW("[3D/Reduced-D COUPLING] 3D map has no variable %s for condition [%d]",
                variable_str.c_str(), condID);
          }
          double density = 0.0;
          double viscosity = 0.0;
          double area = 0.0;
          area = mapiter->second->FluidCouplingBc::Area(density, viscosity, condID);
          (*map3_dnp_)[CouplingVariable.str()] =
              mapiter->second->FluidCouplingBc::PressureCalculation(time, dta, condID);
          (*map3_dnp_)[CouplingVariable.str()] /= area;
        }
        else
        {
          FOUR_C_THROW("(%s): No such coupling variable on the 3D side is defined yet",
              variable_str.c_str());
        }
        if (discret3_d_->Comm().MyPID() == 0)
        {
          std::cout << "3D condition "
                    << " [" << condID << "] returns " << variable_str << " "
                    << (*map3_dnp_)[CouplingVariable.str()] << " at time " << time << std::endl;
        }
        break;
      }
    }
  }

  // -------------------------------------------------------------------
  // Solve the reduced-D problem
  // -------------------------------------------------------------------

  int NumOfSteps = int(dt_f3_ / dt_rm_);

  for (int N = 0; N < NumOfSteps; N++)
  {
    for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
    {
      // Get condition ID
      condID = mapiter->first;
      //----------------------------------------------------------------
      // find the parameters that must be calculated and returned by
      // the reducedD problem
      //----------------------------------------------------------------
      for (unsigned int i = 0; i < conds3D.size(); i++)
      {
        if (conds_redD[i]->parameters().Get<int>("ConditionID") == condID)
        {
          // get returned value name from 3D boundary
          std::string variable_str =
              (conds_redD[i]->parameters().Get<std::string>("ReturnedVariable"));

          // concatenate the variable name with the variable id
          std::stringstream CouplingVariable;
          CouplingVariable << variable_str << "_" << condID;

          if (variable_str == "pressure")
          {
            auto itr = map_red_dnp_->find(CouplingVariable.str());
            if (itr == map_red_dnp_->end())
            {
              FOUR_C_THROW(
                  "[3D/Reduced-D COUPLING] reduced-D map has no variable %s for condition [%d]",
                  variable_str.c_str(), condID);
            }
            (*map_red_dnp_)[CouplingVariable.str()] = 0.0;
          }
          else if (variable_str == "flow")
          {
            std::map<std::string, double>::iterator itr =
                map_red_dnp_->find(CouplingVariable.str());
            if (itr == map_red_dnp_->end())
            {
              FOUR_C_THROW(
                  "[3D/Reduced-D COUPLING] reduced-D map has no variable %s for condition [%d]",
                  variable_str.c_str(), condID);
            }
            (*map_red_dnp_)[CouplingVariable.str()] = 0.0;
          }

          else
          {
            FOUR_C_THROW("(%s): No such coupling variable on the 3D side is defined yet",
                variable_str.c_str());
          }
          break;
        }
      }
    }
    // -----------------------------------------------------------------
    // Define a map that will have the interpolated values at the
    // reduced-D time subscale
    // -----------------------------------------------------------------
    Teuchos::RCP<std::map<std::string, double>> map3D_inter_to_Red =
        Teuchos::rcp(new std::map<std::string, double>);
    double dstep = 1.0 / double(NumOfSteps);

    // -----------------------------------------------------------------
    // Calculate the variables with in the reduced-D time subscale
    //
    //
    //    ^
    // V  |
    //    |                                      +
    //    |                                   .  .
    //    |                               o
    //    |                           .   .      .
    //    |                       //
    //    |                   .   .       .      .
    //    |               o
    //    |           .   .       .       .      .
    //    |       +
    //    |       .       .       .       .      .
    //    |
    //    |       .       .       .       .      .
    //    +-------+-------+-------//------+------+--------->
    //    |       0       1       i      N-2     N-1     Step
    //
    //
    //  ds = 1/N
    //                /  i  \                                           .
    //  V|   =  V|  + |-----|*
    //    i       1   \ N-1 /
    //
    //
    //
    // -----------------------------------------------------------------

    std::map<std::string, double>::iterator itr_sub;
    for (itr_sub = map3_dnp_->begin(); itr_sub != map3_dnp_->end(); itr_sub++)
    {
      std::string var_str = itr_sub->first;
      double var = (*map3_dn_)[var_str];
      double dvar = itr_sub->second - (*map3_dn_)[var_str];
      var = var + double(N + 1) * dstep * (dvar);

      (*map3D_inter_to_Red)[var_str] = var;
    }

    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    //    params->set("3D map of values", map3_Dnp_);
    params->set("3D map of values", map3D_inter_to_Red);
    params->set("reducedD map of values", map_red_dnp_);
    double subscale_time = time - (dt_rm_ * double(NumOfSteps - N - 1));
    params->set("time", subscale_time);
    // #endif

    //    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp( new Teuchos::ParameterList);
    //    params->set("3D map of values", map3_Dnp_);
    //    params->set("reducedD map of values", mapRed_Dnp_);


    this->Integrate(true, params);
  }

  // -------------------------------------------------------------------
  // Get the reduced-D results from all of the processors
  // if a processor doesn't have any reduced-D results then it must
  // return values equivelant to zero
  // -------------------------------------------------------------------

  for (unsigned int i = 0; i < conds3D.size(); i++)
  {
    // Get condition ID
    int ID = conds_redD[i]->parameters().Get<int>("ConditionID");

    // Concatenate the returned value with the condition ID
    std::string ReturnedVariable =
        (conds_redD[i]->parameters().Get<std::string>("ReturnedVariable"));

    std::stringstream VariableWithId;
    VariableWithId << ReturnedVariable << "_" << ID;

    // Get the variable with is returned
    double var = (*map_red_dnp_)[VariableWithId.str()];

    // update the variable on all of the processors
    double par_var = 0.0;
    discret3_d_->Comm().SumAll(&var, &par_var, 1);
    (*map_red_dnp_)[VariableWithId.str()] = par_var;

    // Apply the boundary condition to the outlets
    if (ReturnedVariable == "pressure")
    {
      coup_map3_d_[ID]->OutflowBoundary(par_var, time, dta, theta, ID);
    }
    else if (ReturnedVariable == "flow")
    {
      coup_map3_d_[ID]->InflowBoundary(par_var, time, dta, theta, ID);
    }
    else
    {
      FOUR_C_THROW(
          "Reduced-dimensional problem, returned a value of type [%s] at the condition (%d)",
          ReturnedVariable.c_str(), ID);
      exit(0);
    }
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap update of residual                                ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::UpdateResidual(Teuchos::RCP<Epetra_Vector> residual)
{
  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  (*map_red_dn_) = (*map_red_dnp_);
  (*map3_dn_) = (*map3_dnp_);

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    mapiter->second->FluidCouplingBc::UpdateResidual(residual);
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap update of residual                                ismail 05/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::EvaluateDirichlet(
    Teuchos::RCP<Epetra_Vector> velnp, const Epetra_Map& condmap, double time)
{
  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  (*map_red_dn_) = (*map_red_dnp_);
  (*map3_dn_) = (*map3_dnp_);

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    mapiter->second->FluidCouplingBc::EvaluateDirichlet(velnp, condmap, time);
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap restart writing                                   ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::WriteRestart(IO::DiscretizationWriter& output)
{
  std::map<std::string, double>::iterator it;
  //! map of coupling variables returned by the 3-D model at time step n+1
  for (it = map3_dnp_->begin(); it != map3_dnp_->end(); it++)
  {
    std::stringstream stream;
    stream << it->first << "_3D_np";
    output.WriteDouble(stream.str(), it->second);
  }
  //! map of coupling variables returned by the 3-D model at time step n
  for (it = map3_dn_->begin(); it != map3_dn_->end(); it++)
  {
    std::stringstream stream;
    stream << it->first << "_3D_n";
    output.WriteDouble(stream.str(), it->second);
  }
  //! map of coupling variables returned by the reduced-D model at time step n+1
  for (it = map3_dnp_->begin(); it != map3_dnp_->end(); it++)
  {
    std::stringstream stream;
    stream << it->first << "_Red_np";
    output.WriteDouble(stream.str(), it->second);
  }
  //! map of coupling variables returned by the reduced-D model at time step n
  for (it = map_red_dn_->begin(); it != map_red_dn_->end(); it++)
  {
    std::stringstream stream;
    stream << it->first << "_Red_n";
    output.WriteDouble(stream.str(), it->second);
  }

  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
  {
    mapiter->second->FluidCouplingBc::WriteRestart(output, mapiter->first);
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap restart reading                                   ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingWrapperBase::ReadRestart(IO::DiscretizationReader& reader)
{
  std::map<std::string, double>::iterator it;
  //! map of coupling variables returned by the 3-D model at time step n+1
  for (it = map3_dnp_->begin(); it != map3_dnp_->end(); it++)
  {
    std::stringstream stream;
    double val = 0.0;
    stream << it->first << "_3D_np";
    val = reader.ReadDouble(stream.str());
    it->second = val;
  }
  //! map of coupling variables returned by the 3-D model at time step n
  for (it = map3_dn_->begin(); it != map3_dn_->end(); it++)
  {
    std::stringstream stream;
    double val = 0.0;
    stream << it->first << "_3D_n";
    val = reader.ReadDouble(stream.str());
    it->second = val;
  }
  //! map of coupling variables returned by the reduced-D model at time step n+1
  for (it = map3_dnp_->begin(); it != map3_dnp_->end(); it++)
  {
    std::stringstream stream;
    double val = 0.0;
    stream << it->first << "_Red_np";
    val = reader.ReadDouble(stream.str());
    it->second = val;
  }
  //! map of coupling variables returned by the reduced-D model at time step n
  for (it = map_red_dn_->begin(); it != map_red_dn_->end(); it++)
  {
    std::stringstream stream;
    double val = 0.0;
    stream << it->first << "_Red_n";
    val = reader.ReadDouble(stream.str());
    it->second = val;
  }

  std::map<const int, Teuchos::RCP<class FluidCouplingBc>>::iterator mapiter;

  for (mapiter = coup_map3_d_.begin(); mapiter != coup_map3_d_.end(); mapiter++)
    mapiter->second->FluidCouplingBc::ReadRestart(reader, mapiter->first);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 12/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidCouplingBc::FluidCouplingBc(Teuchos::RCP<DRT::Discretization> dis_3D,
    Teuchos::RCP<DRT::Discretization> dis_redD, IO::DiscretizationWriter& output, double dt_3d,
    double dt_rm, int condid, int numcond, int numcond2)
    :  // call constructor for "nontrivial" objects
      condid_(condid),
      discret_3_d_(dis_3D),
      discret_red_d_(dis_redD),
      output_(output)
{
  // ---------------------------------------------------------------------
  // read in all 3D to reducedD boundary conditions
  // ---------------------------------------------------------------------
  std::vector<CORE::Conditions::Condition*> couplingcond;
  dis_redD->GetCondition("Art_redD_3D_CouplingCond", couplingcond);


  // ---------------------------------------------------------------------
  // get time steps size
  // ---------------------------------------------------------------------
  dt_f3_ = dt_3d;
  dt_rm_ = dt_rm;

  // ---------------------------------------------------------------------
  // get the processor ID from the communicator
  // ---------------------------------------------------------------------
  myrank_ = discret_3_d_->Comm().MyPID();

  // ---------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // ---------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_3_d_->DofRowMap();
  couplingbc_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  flowrate_ = 0.0;

  // ---------------------------------------------------------------------
  // Calculate alfa, the velocity correction factor:
  //
  //
  //  -------------+                     -------------+<------
  //               +<------                           +<------
  //       3D      +<------                   3D      +<------
  //    GEOMETRY   +<------                GEOMETRY   +<------
  //               +<------                           +<------
  //               +<------                           +<------
  //  -------------+                     -------------+<------
  //                   ^                                  ^
  //                   |                                  |
  //           Actual applied                      Actual calculated
  //              velocity                             velocity
  //
  //              Q|
  //               |calcutated
  //      alfa = --------
  //              Q|
  //               |applied
  // ---------------------------------------------------------------------

  if ((couplingcond[numcond2]->parameters().Get<std::string>("ReturnedVariable")) == "flow")
  {
    double density = 0.0;
    double viscosity = 0.0;
    double time = 0.0;
    double flowrate = this->FlowRateCalculation(time, dt_3d, condid);
    double area = this->Area(density, viscosity, condid);
    if (flowrate == 0.0)
    {
      FOUR_C_THROW(
          "3D SURF condition (%d) expects a flowrate from 1D problem,\nthus it must have a "
          "Dirichlet BC of 1 in the direction of flow",
          condid);
    }
    else
    {
      alfa_ = fabs(area / flowrate);
    }
    if (myrank_ == 0)
    {
      std::cout << "Velocity correction factor cond(" << condid << ") is: " << alfa_ << std::endl;
    }
  }
  else
  {
    alfa_ = 1.0;
  }

  velocity_ = 0.0;
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Restart writing                                        ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingBc::WriteRestart(IO::DiscretizationWriter& output, int condnum)
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!
  std::stringstream stream1, stream2, stream3;

  // also write vector couplingbc_
  stream3 << "couplingbc" << condnum;
  output.WriteVector(stream3.str(), couplingbc_);


  // write time steps size
  output.WriteDouble("dta_3D", dt_f3_);
  output.WriteDouble("reduced_D_dta", dt_rm_);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Restart reading                                        ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidCouplingBc::ReadRestart(IO::DiscretizationReader& reader, int condnum)
{
  std::stringstream stream1, stream2, stream3;

  // also read vector couplingbc_
  stream3 << "couplingbc" << condnum;
  reader.ReadVector(couplingbc_, stream3.str());

  // read time steps size
  dt_f3_ = reader.ReadDouble("dta_3D");
  dt_rm_ = reader.ReadDouble("reduced_D_dta");


  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Area calculation                                        ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!

*/
double FLD::UTILS::FluidCouplingBc::Area(double& density, double& viscosity, int condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", FLD::calc_area);
  eleparams.set<double>("area", 0.0);
  eleparams.set<double>("viscosity", 0.0);
  eleparams.set<double>("density", 0.0);

  const std::string condstring("Art_3D_redD_CouplingCond");

  discret_3_d_->EvaluateCondition(eleparams, condstring, condid);

  double actarea = eleparams.get<double>("area");
  density = eleparams.get<double>("density");
  viscosity = eleparams.get<double>("viscosity");

  // find the lowest proc number that knows the material data
  int numproc = discret_3_d_->Comm().NumProc();
  int theproc = -1;  // the lowest proc that has the desired information
  std::vector<double> alldens(numproc);

  discret_3_d_->Comm().GatherAll(&density, alldens.data(), 1);
  for (int i = 0; i < numproc; i++)
    if (alldens[i] > 0.0)
    {
      theproc = i;
      break;
    }
  if (theproc < 0) FOUR_C_THROW("Something parallel went terribly wrong!");

  // do the actual communication of density ...
  discret_3_d_->Comm().Broadcast(&density, 1, theproc);
  // ... and viscosity
  discret_3_d_->Comm().Broadcast(&viscosity, 1, theproc);

  // get total area in parallel case
  double pararea = 0.0;
  discret_3_d_->Comm().SumAll(&actarea, &pararea, 1);

  if (myrank_ == 0)
  {
    std::cout << "3D/Reduced-D coupling condition Id: " << condid << " area = " << pararea
              << std::endl;
  }
  return pararea;
}  // FluidImplicitTimeInt::Area



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Flow rate calculation                                  ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
  modified by chfoe 04/08

  Calculate the flow rate across an impedance boundary surface

  Flow rates are
  (1) calculated for single element surfaces
  (2) added up over the elements of the single procs
  (3) communicated and added over the procs
  (4) and finally stored within the vector 'flowrates_'

  The vector of the flowrates holds the flow rate history of the
  very last cycle!

*/
double FLD::UTILS::FluidCouplingBc::FlowRateCalculation(double time, double dta, int condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", FLD::calc_flowrate);
  eleparams.set("total time", time);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_3_d_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> flowrates = CORE::LINALG::CreateVector(*dofrowmap, true);
  const std::string condstring("Art_3D_redD_CouplingCond");
  discret_3_d_->EvaluateCondition(eleparams, flowrates, condstring, condid);

  double local_flowrate = 0.0;
  for (int i = 0; i < dofrowmap->NumMyElements(); i++)
  {
    local_flowrate += ((*flowrates)[i]);
  }

  double flowrate = 0.0;
  dofrowmap->Comm().SumAll(&local_flowrate, &flowrate, 1);

  return flowrate;
}  // FluidImplicitTimeInt::FlowRateCalculation


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Pressure calculation                                   ismail 04/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
  Calculate the pressure across a coupling boundary surface

  Flow rates are
  (1) calculated for single element surfaces
  (2) added up over the elements of the single procs
  (3) communicated and added over the procs
  (4) devide the integrated pressure over the cross-sectional area
  (4) and finally stored within the vector 'pressures_'

  The vector of the flowrates holds the flow rate history of the
  very last cycle!

*/
double FLD::UTILS::FluidCouplingBc::PressureCalculation(double time, double dta, int condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", FLD::calc_pressure_bou_int);
  eleparams.set<double>("pressure boundary integral", 0.0);
  eleparams.set("total time", time);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_3_d_->DofRowMap();

  // get elemental flowrates ...
  Teuchos::RCP<Epetra_Vector> myStoredPressures = Teuchos::rcp(new Epetra_Vector(*dofrowmap, 100));
  const std::string condstring("Art_3D_redD_CouplingCond");
  discret_3_d_->EvaluateCondition(eleparams, myStoredPressures, condstring, condid);

  // ... as well as actual total flowrate on this proc
  double actpressure = eleparams.get<double>("pressure boundary integral");

  // get total flowrate in parallel case
  double parpressure = 0.0;
  discret_3_d_->Comm().SumAll(&actpressure, &parpressure, 1);

  return parpressure;
}  // FluidImplicitTimeInt::PressureCalculation


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply outflow boundary to the coupled surface          ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*

  This routine contains major parts of the following paper:

  Olufsen et al.: "Numerical Simulation and Experimental Validation of
  Blood Flow in Arteries with Structured-Tree Outflow Conditions",
  Annals of Biomedical Eingineering, Vol. 28, pp. 1281--1299, 2000.

  Basic Idea:
  (1) Evaluate convolution integral over one cycle to obtain outflow
      pressure
  (2) Apply this pressure as a Neumann-load type at the outflow boundary

*/
void FLD::UTILS::FluidCouplingBc::OutflowBoundary(
    double pressure, double time, double dta, double theta, int condid)
{
  // call the element to apply the pressure
  Teuchos::ParameterList eleparams;
  // action for elements
  // the reason we have Outlet impedance as action is because we don't
  // want to rewrite the implimented code
  eleparams.set<int>("action", FLD::Outletimpedance);

  eleparams.set("total time", time);
  eleparams.set("delta time", dta);
  eleparams.set("thsl", theta * dta);
  eleparams.set("WindkesselPressure", pressure);

  if (myrank_ == 0)
    printf(
        "3D/reduced-D coupling condition Id: %d Pressure %f at time %f\n", condid, pressure, time);


  couplingbc_->PutScalar(0.0);
  const std::string condstring("Art_3D_redD_CouplingCond");
  discret_3_d_->EvaluateCondition(eleparams, couplingbc_, condstring, condid);

  //  discret_3D_->ClearState();

  return;
}  // FluidImplicitTimeInt::OutflowBoundary

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply inflow boundary to the coupled surface           ismail 05/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*
  (2) Apply this flowrate as a Dirichlet boundary at the inflow boundary

*/
void FLD::UTILS::FluidCouplingBc::InflowBoundary(
    double flowrate, double time, double dta, double theta, int condid)
{
  // call the element to apply the pressure
  Teuchos::ParameterList eleparams;
  // action for elements
  // the reason we have Outlet impedance as action is because we don't
  // want to rewrite the implimented code


  if (myrank_ == 0)
    printf("3D/reduced-D coupling condition Id: %d flowrate = %f\n", condid, flowrate);


  std::vector<CORE::Conditions::Condition*> cond3D;
  discret_3_d_->GetCondition("Art_3D_redD_CouplingCond", cond3D);

  double area = 0.0;
  double density = 0.0;
  double viscosity = 0.0;

  area = Area(density, viscosity, condid_);
  velocity_ = flowrate / area;

  return;
}  // FluidCouplingBc::InflowBoundary



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update of residual vector                              ismail 12/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
 */
void FLD::UTILS::FluidCouplingBc::UpdateResidual(Teuchos::RCP<Epetra_Vector> residual)
{
  residual->Update(1.0, *couplingbc_, 1.0);
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update dirichlet values                                ismail 05/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
 */
void FLD::UTILS::FluidCouplingBc::EvaluateDirichlet(
    Teuchos::RCP<Epetra_Vector> velnp, const Epetra_Map& condmap, double time)
{
  return;
  std::cout << "Evaluating Dirich!" << std::endl;
  // FOUR_C_THROW("Dirichlet coupling is not fixed yet, if you see the message then something is
  // wrong!");
  //  std::cout<<"3D discretization:"<<std::endl<<*discret_3D_<<std::endl;
  std::vector<CORE::Conditions::Condition*> conds_red;
  discret_red_d_->GetCondition("Art_redD_3D_CouplingCond", conds_red);

  CORE::Conditions::Condition* cond_red = nullptr;
  for (unsigned int i = 0; i != conds_red.size(); i++)
  {
    if (conds_red[i]->parameters().Get<int>("ConditionID") == condid_)
    {
      cond_red = conds_red[i];
      break;
    }
  }

  if ((cond_red->parameters().Get<std::string>("ReturnedVariable")) != "flow")
  {
    return;
  }

  //  residual->Update(1.0,*couplingbc_,1.0);
  std::vector<CORE::Conditions::Condition*> conds;
  discret_3_d_->GetCondition("Art_3D_redD_CouplingCond", conds);

  CORE::Conditions::Condition* cond;
  for (unsigned int i = 0; i != conds.size(); i++)
  {
    if (conds[i]->parameters().Get<int>("ConditionID") == condid_)
    {
      cond = conds[i];
      break;
    }
  }

  double area = 0.0;
  double density = 0.0;
  double viscosity = 0.0;

  area = this->Area(density, viscosity, condid_);

  double Dflowrate = this->FlowRateCalculation(time, dt_f3_, condid_);
  alfa_ = fabs(area / Dflowrate);

  alfa_ = 1.0;
  velocity_ *= alfa_;

  std::cout << "velocity: " << std::endl;
  std::cout << "Dflowrate: " << Dflowrate << std::endl;
  std::cout << "area: " << area << std::endl;
  const std::vector<int>* nodes = cond->GetNodes();

  for (unsigned int i = 0; i < nodes->size(); i++)
  {
    int gid = (*nodes)[i];
    std::cout << "Node(" << gid << "): ";

    if (discret_3_d_->HaveGlobalNode(gid))
    {
      DRT::Node* node = discret_3_d_->gNode(gid);
      unsigned int numDof = discret_3_d_->NumDof(node);
      std::cout << "(" << numDof << ") dof --> ";
      for (unsigned int dof = 0; dof < numDof - 1; dof++)
      {
        int dof_gid = discret_3_d_->Dof(node, dof);
        //        std::cout<<"("<<dof<<")+>["<<dof_gid<<"]\t";
        if (condmap.MyGID(dof_gid))
        {
          int lid = discret_3_d_->DofRowMap()->LID(dof_gid);

          double val = (*velnp)[lid] * velocity_;
          //        std::cout<<"Vel["<<gid<<"]: "<<(*velnp) [lid]<<std::endl;
          if ((*velnp)[lid] > 1.0)
          {
            FOUR_C_THROW("coupled 3D/Reduced-D must have Dirichlet BC = 1");
            exit(1);
          }
          std::cout << "[" << dof_gid << "]\t|" << val << "\t<-<" << (*velnp)[lid] << "|\t";
          velnp->ReplaceGlobalValues(1, &val, &dof_gid);
        }
      }
    }
    std::cout << std::endl;
  }
  //  exit(1);
}

FOUR_C_NAMESPACE_CLOSE
