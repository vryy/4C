/*!
\file xfem_coupling_base.cpp

\brief is the base for the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

<pre>
Maintainer: Benedikt Schott
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
  else if(condname == "XFEMLevelsetTwophase")       return INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE;
  else if(condname == "XFEMLevelsetCombustion")     return INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION;
  else if(condname == "XFEMSurfWeakDirichlet")      return INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET;
  else if(condname == "XFEMSurfNeumann")            return INPAR::XFEM::CouplingCond_SURF_NEUMANN;
  //else dserror("condition type not supported: %s", condname.c_str());

  return INPAR::XFEM::CouplingCond_NONE;
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
      // get all conditions with given condition name
      std::vector<DRT::Condition*> mycond;
      DRT::UTILS::FindElementConditions(cutele, conditions_to_copy_[cond], mycond);

      // safety checks
      if (mycond.size()>1)
      {
//        dserror("%i conditions of the same name for element %i! %s coupling-condition not unique!", mycond.size(), cutele->Id(), conditions_to_copy_[cond].c_str());
      }
      else if(mycond.size() == 0)
      {
        continue; // try the next condition type
      }
      else
      {
//        std::cout << "unique condition found!" << std::endl;
      }

      INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(conditions_to_copy_[cond]);

      // non-coupling condition found (e.g. FSI coupling)
      if(cond_type == INPAR::XFEM::CouplingCond_NONE) continue;

      // non-unique conditions for one cutter element
      if( cutterele_conds_[lid].first != INPAR::XFEM::CouplingCond_NONE )
      {
        dserror("There are two different condition types for the same cutter dis element with id %i: 1st %i, 2nd %i. Make the XFEM coupling conditions unique!",
            cutele->Id(), cutterele_conds_[lid].first, cond_type);
      }

      // store the unique condition pointer to the cutting element
      cutterele_conds_[lid] = EleCoupCond(cond_type, mycond[0]);
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
    printf("   +----------+-----------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   | %8i | %9i | %27s | %27s | %27s | %27s | %27s |\n",
        coupling_idx ,
        side_start_gid,
        TypeToStringForPrint(CondType_stringToEnum(cond_name_)).c_str(),
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
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    averaging_strategy_ = INPAR::XFEM::Xfluid_Sided;
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
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
  {
    coupl_dis_ = cutter_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
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

    final_values[dof] = num*functionfac;
  } // loop dofs
}
