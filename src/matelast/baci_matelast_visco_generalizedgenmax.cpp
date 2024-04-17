/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the viscous contribution a generalized maxwell model

\level 2
*/
/*----------------------------------------------------------------------*/
#include "baci_matelast_visco_generalizedgenmax.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::ELASTIC::PAR::GeneralizedGenMax::GeneralizedGenMax(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      numbranch_(*matdata->Get<int>("NUMBRANCH")),
      matids_(matdata->Get<std::vector<int>>("MATIDS")),
      solve_(*matdata->Get<std::string>("SOLVE"))

{
}

MAT::ELASTIC::GeneralizedGenMax::GeneralizedGenMax(MAT::ELASTIC::PAR::GeneralizedGenMax* params)
    : params_(params), branchespotsum_(0), internalpotsum_(0)
{
  // loop over materials of GeneralizedGenMax (branches)
  std::vector<int>::const_iterator m;
  for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
  {
    // make sure the summands of the current branch is empty
    internalpotsum_.clear();
    // get parameters of each branch
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> ViscoBranch = MAT::ELASTIC::Summand::Factory(matid);

    double nummat = -1.0;
    const std::vector<int>* branchmatids = nullptr;

    ViscoBranch->ReadMaterialParameters(nummat, branchmatids);

    // loop over materials of ViscoBranch (components of the viscoelastic branch)
    for (int i = 0; i < nummat; ++i)
    {
      // get parameters of each component
      int curmatid = branchmatids->at(i);
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(curmatid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      // write summand in the vector of summands of each branch
      internalpotsum_.push_back(sum);
    }

    // write into vector of summands of the GeneralizedGenMax material
    branchespotsum_.push_back(internalpotsum_);

  }  // end for-loop over branches
}

void MAT::ELASTIC::GeneralizedGenMax::ReadMaterialParameters(
    int& numbranch, const std::vector<int>*& matids, std::string& solve)
{
  numbranch = params_->numbranch_;
  matids = params_->matids_;
  solve = params_->solve_;
}

// Viscobranch
MAT::ELASTIC::PAR::ViscoBranch::ViscoBranch(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      nummat_(*matdata->Get<int>("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS"))
{
}

MAT::ELASTIC::ViscoBranch::ViscoBranch(MAT::ELASTIC::PAR::ViscoBranch* params) : params_(params) {}

void MAT::ELASTIC::ViscoBranch::ReadMaterialParameters(
    double& nummat, const std::vector<int>*& matids)
{
  nummat = params_->nummat_;
  matids = params_->matids_;
}

// Viscopart
MAT::ELASTIC::PAR::ViscoPart::ViscoPart(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), tau_(*matdata->Get<double>("TAU"))
{
}

MAT::ELASTIC::ViscoPart::ViscoPart(MAT::ELASTIC::PAR::ViscoPart* params) : params_(params) {}

void MAT::ELASTIC::ViscoPart::ReadMaterialParameters(double& tau) { tau = params_->tau_; }

FOUR_C_NAMESPACE_CLOSE
