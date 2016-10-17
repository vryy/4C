/*!
\file xfem_coupling_base.cpp

\brief is the base for the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/



#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_base.H"

#include "../drt_lib/drt_condition_utils.H"


INPAR::XFEM::EleCouplingCondType XFEM::CondType_stringToEnum(const std::string& condname)
{
  if     (condname == "XFEMSurfFSIPart")            return INPAR::XFEM::CouplingCond_SURF_FSI_PART;
  else if(condname == "XFEMSurfFSIMono")            return INPAR::XFEM::CouplingCond_SURF_FSI_MONO;
  else if(condname == "XFEMSurfCrackFSIPart")       return INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART;
  else if(condname == "XFEMSurfFluidFluid")         return INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID;
  else if(condname == "XFEMLevelsetWeakDirichlet")  return INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET;
  else if(condname == "XFEMLevelsetNeumann")        return INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN;
  else if(condname == "XFEMLevelsetNavierSlip")     return INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP;
  else if(condname == "XFEMLevelsetTwophase")       return INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE;
  else if(condname == "XFEMLevelsetCombustion")     return INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION;
  else if(condname == "XFEMSurfWeakDirichlet")      return INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET;
  else if(condname == "XFEMSurfNeumann")            return INPAR::XFEM::CouplingCond_SURF_NEUMANN;
  else if(condname == "XFEMSurfNavierSlip")         return INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP;
  //else dserror("condition type not supported: %s", condname.c_str());

  return INPAR::XFEM::CouplingCond_NONE;
}

/*--------------------------------------------------------------------------*
* constructor
*--------------------------------------------------------------------------*/
XFEM::CouplingBase::CouplingBase(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step       ///< time step
) :
bg_dis_(bg_dis),
cond_name_(cond_name),
cond_dis_(cond_dis),
coupling_id_(coupling_id),
cutter_dis_(Teuchos::null),
coupl_dis_(Teuchos::null),
averaging_strategy_(INPAR::XFEM::invalid),
myrank_(bg_dis_->Comm().MyPID()),
dt_(-1.0),
time_(time),
step_(step)
{
  // initialize element level configuration map (no evaluation)
  InitConfigurationMap();
}

/*--------------------------------------------------------------------------*
 * Initialize Configuration Map --> No Terms are evaluated at the interface
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::InitConfigurationMap()
{
  //Configuration of Consistency Terms
  //all components:
  configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_Col] = std::pair<bool,double>(false,0.0);
  //normal terms:
  configuration_map_[INPAR::XFEM::F_Con_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Con_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_n_Col] = std::pair<bool,double>(false,0.0);
  //tangential terms:
  configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Con_t_Col] = std::pair<bool,double>(false,0.0);

  //Configuration of Adjount Consistency Terms
  //all components:
  configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Adj_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Adj_Col] = std::pair<bool,double>(false,0.0);
  //normal terms:
  configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Adj_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Adj_n_Col] = std::pair<bool,double>(false,0.0);
  //tangential terms:
  configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Adj_t_Col] = std::pair<bool,double>(false,0.0);

  //Configuration of Penalty Terms
  //all components:
  configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Pen_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Pen_Col] = std::pair<bool,double>(false,0.0);
  //normal terms:
  configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Pen_n_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Pen_n_Col] = std::pair<bool,double>(false,0.0);
  //tangential terms:
  configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::FStr_Pen_t_Col] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::XStr_Pen_t_Col] = std::pair<bool,double>(false,0.0);

  //Starting from here are some special Terms
  configuration_map_[INPAR::XFEM::F_LB_Rhs] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_LB_Rhs] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::F_TJ_Rhs] = std::pair<bool,double>(false,0.0);
  configuration_map_[INPAR::XFEM::X_TJ_Rhs] = std::pair<bool,double>(false,0.0);
  return;
}


void XFEM::CouplingBase::SetElementConditions()
{
  // number of column cutter boundary elements
  int nummycolele = cutter_dis_->NumMyColElements();

  cutterele_conds_.clear();
  cutterele_conds_.reserve(nummycolele);

  // initialize the vector invalid coupling-condition type "NONE"
  EleCoupCond init_pair = EleCoupCond(INPAR::XFEM::CouplingCond_NONE,NULL);
  for(int lid=0; lid<nummycolele; lid++) cutterele_conds_.push_back(init_pair);

  //-----------------------------------------------------------------------------------
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    DRT::Element* cutele = cutter_dis_->lColElement(lid);

    // loop all possible XFEM-coupling conditions
    for(size_t cond=0; cond < conditions_to_copy_.size(); cond++)
    {
      INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(conditions_to_copy_[cond]);

      // non-coupling condition found (e.g. FSI coupling)
      if(cond_type == INPAR::XFEM::CouplingCond_NONE) continue;

      // get all conditions with given condition name
      std::vector<DRT::Condition*> mycond;
      DRT::UTILS::FindElementConditions(cutele, conditions_to_copy_[cond], mycond);

      std::vector<DRT::Condition*> mynewcond;
      GetConditionByCouplingId(mycond, coupling_id_, mynewcond);

      DRT::Condition* cond_unique = NULL;

      // safety checks
      if(mynewcond.size() == 0)
      {
        continue; // try the next condition type
      }
      else if(mynewcond.size() == 1) // unique condition found
      {
        cond_unique = mynewcond[0];
      }
      else if(mynewcond.size()>1)
      {
        // get the right condition
        dserror("%i conditions of the same name with coupling id %i, for element %i! %s coupling-condition not unique!", mynewcond.size(), coupling_id_, cutele->Id(), conditions_to_copy_[cond].c_str());
      }

      // non-unique conditions for one cutter element
      if( cutterele_conds_[lid].first != INPAR::XFEM::CouplingCond_NONE )
      {
        dserror("There are two different condition types for the same cutter dis element with id %i: 1st %i, 2nd %i. Make the XFEM coupling conditions unique!",
            cutele->Id(), cutterele_conds_[lid].first, cond_type);
      }

      // store the unique condition pointer to the cutting element
      cutterele_conds_[lid] = EleCoupCond(cond_type, cond_unique);
    }
  }

  //-----------------------------------------------------------------------------------
  // check if all column cutter elements have a valid condition type
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    if(cutterele_conds_[lid].first == INPAR::XFEM::CouplingCond_NONE)
      dserror("cutter element with local id %i has no valid coupling-condition", lid);
  }

}

void XFEM::CouplingBase::GetConditionByCouplingId(
    const std::vector<DRT::Condition*> & mycond,
    const int coupling_id,
    std::vector<DRT::Condition*> & mynewcond
)
{
  mynewcond.clear();

  // select the conditions with specified "couplingID"
  for(size_t i=0; i< mycond.size(); ++i)
  {
    DRT::Condition* cond = mycond[i];
    const int id = cond->GetInt("label");

    if(id == coupling_id)
      mynewcond.push_back(cond);
  }
}

void XFEM::CouplingBase::Status(
    const int coupling_idx,
    const int side_start_gid
)
{
  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    printf("   +----------+-----------+-----------------------------+---------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   | %8i | %9i | %27s | %7i | %27s | %27s | %27s | %27s |\n",
        coupling_idx ,
        side_start_gid,
        TypeToStringForPrint(CondType_stringToEnum(cond_name_)).c_str(),
        coupling_id_,
        DisNameToString(cutter_dis_).c_str(),
        DisNameToString(cond_dis_).c_str(),
        DisNameToString(coupl_dis_).c_str(),
        AveragingToStringForPrint(averaging_strategy_).c_str());
  }
}



void XFEM::CouplingBase::SetAveragingStrategy()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    averaging_strategy_ = INPAR::XFEM::Xfluid_Sided;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // ask the first cutter element
    const int lid=0;
    const int val = cutterele_conds_[lid].second->GetInt("COUPSTRATEGY");
    averaging_strategy_ = static_cast<INPAR::XFEM::AveragingStrategy>(val);
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    averaging_strategy_ = INPAR::XFEM::Harmonic;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP:
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
  {
    averaging_strategy_ = INPAR::XFEM::Xfluid_Sided;
    break;
  }
  default: dserror("which is the averaging strategy for this type of coupling %i?", cond_type); break;
  }
}


void XFEM::CouplingBase::SetCouplingDiscretization()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    coupl_dis_ = cutter_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // depending on the weighting strategy
    if(averaging_strategy_==INPAR::XFEM::Xfluid_Sided)
    {
      coupl_dis_ = cutter_dis_;
    }
    else if(averaging_strategy_==INPAR::XFEM::Embedded_Sided or
        averaging_strategy_==INPAR::XFEM::Mean )
    {
      coupl_dis_ = cond_dis_;
    }
    else dserror("invalid coupling strategy for fluid-fluid application");
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    coupl_dis_ = bg_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET: // set this to Teuchos::null when the values are read from the function instead of the ivelnp vector
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP:
  {
    coupl_dis_ = cutter_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
  {
    coupl_dis_ = Teuchos::null;
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
  }
}

void XFEM::CouplingBase::EvaluateDirichletFunction(
    LINALG::Matrix<3,1>& ivel,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond,
    double time
)
{
  std::vector<double> final_values(3,0.0);

  EvaluateFunction(final_values, x.A(), cond, time);

  ivel(0,0) = final_values[0];
  ivel(1,0) = final_values[1];
  ivel(2,0) = final_values[2];
}

void XFEM::CouplingBase::EvaluateNeumannFunction(
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond,
    double time
)
{
  std::vector<double> final_values(3,0.0);

  //---------------------------------------
  const std::string* condtype = cond->Get<std::string>("type");

  // get usual body force
  if (!(*condtype == "neum_dead" or *condtype == "neum_live"))
    dserror("Unknown Neumann condition");
  //---------------------------------------

  EvaluateFunction(final_values, x.A(), cond, time);

  itraction(0,0) = final_values[0];
  itraction(1,0) = final_values[1];
  itraction(2,0) = final_values[2];
}

void XFEM::CouplingBase::EvaluateFunction(
    std::vector<double>& final_values,
    const double* x,
    const DRT::Condition* cond,
    const double time
)
{
  if(cond == NULL) dserror("invalid condition");

  const int numdof = cond->GetInt("numdof");

  if(numdof != (int)final_values.size())
    dserror("you specified NUMDOF %i in the input file, however, only %i dofs allowed!", numdof, (int)final_values.size());

  //---------------------------------------
  // get values and switches from the condition
  const std::vector<int>*    curve     = cond->Get<std::vector<int>    >("curve");
  const std::vector<int>*    onoff     = cond->Get<std::vector<int>    >("onoff");
  const std::vector<double>* val       = cond->Get<std::vector<double> >("val"  );
  const std::vector<int>*    functions = cond->Get<std::vector<int>    >("funct");

  // uniformly distributed random noise

  DRT::Condition& secondary = const_cast<DRT::Condition&>(*cond);
  const std::vector<double>* percentage = secondary.GetMutable<std::vector<double> >("randnoise");

  if(time < -1e-14) dserror("Negative time in curve/function evaluation: time = %f", time);

  //---------------------------------------
  // set this condition
  //---------------------------------------
  for(int dof=0;dof<numdof;++dof)
  {
    // get factor given by spatial function
    int functnum = -1;
    if (functions) functnum = (*functions)[dof];

    // check for potential time curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[dof];

    // initialization of time-curve factor and function factor
    double functionfac = 1.0;
    double curvefac = 1.0;

    // compute potential time curve or set time-curve factor to one
    if (curvenum >= 0)
    {
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    }

    double num = (*onoff)[dof]*(*val)[dof]*curvefac;

    if (functnum>0)
    {
      functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof%numdof,x,time,NULL);
    }

    // uniformly distributed noise
    double noise = 0.0;
    if( percentage != NULL )
    {
      if(percentage->size() != 1) dserror("expect vector of length one!");
      const double perc = percentage->at(0);

      if(fabs(perc)> 1e-14 )
      {
        const double randomnumber = DRT::Problem::Instance()->Random()->Uni(); // uniformly distributed between -1.0, 1.0
        noise = perc * randomnumber;
      }
    }

    final_values[dof] = num*(functionfac+noise);
  } // loop dofs
}

void XFEM::CouplingBase::EvaluateScalarFunction(
    double & final_values,
    const double* x,
    const double& val,
    const DRT::Condition* cond,
    const double time
)
{
  if(cond == NULL) dserror("invalid condition");

  const int numdof = 1;

  //---------------------------------------
  // get values and switches from the condition
  const std::vector<int>*    curve     = cond->Get<std::vector<int>    >("curve");
  const std::vector<int>*    functions = cond->Get<std::vector<int>    >("funct");

  // uniformly distributed random noise

  if((*functions).size()!=1 or (*curve).size()!=1)
    dserror("Do not call EvaluateScalarFunction with more than one function/value/curve provided");

  DRT::Condition& secondary = const_cast<DRT::Condition&>(*cond);
  const std::vector<double>* percentage = secondary.GetMutable<std::vector<double> >("randnoise");

  if(time < -1e-14) dserror("Negative time in curve/function evaluation: time = %f", time);

  //---------------------------------------
  // set this condition
  //---------------------------------------
  for(int dof=0;dof<numdof;++dof)
  {
    // get factor given by spatial function
    int functnum = -1;
    if (functions) functnum = (*functions)[dof];

    // check for potential time curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[dof];

    // initialization of time-curve factor and function factor
    double functionfac = 1.0;
    double curvefac = 1.0;

    // compute potential time curve or set time-curve factor to one
    if (curvenum >= 0)
    {
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    }

    double num = val*curvefac;

    if (functnum>0)
    {
      functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof%numdof,x,time,NULL);
    }

    // uniformly distributed noise
    double noise = 0.0;
    if( percentage != NULL )
    {
      if(percentage->size() != 1) dserror("expect vector of length one!");
      const double perc = percentage->at(0);

      if(fabs(perc)> 1e-14 )
      {
        const double randomnumber = DRT::Problem::Instance()->Random()->Uni(); // uniformly distributed between -1.0, 1.0
        noise = perc * randomnumber;
      }
    }

    final_values = num*(functionfac+noise);
  } // loop dofs
}
