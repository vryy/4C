/*----------------------------------------------------------------------*/
/*! \file

\brief generalized maxwell model
\level 2
\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "visco_generalizedgenmax.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | Constructor Parameter                                                        |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::GeneralizedGenMax::GeneralizedGenMax(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numbranch_(matdata->GetInt("NUMBRANCH")),
      matids_(matdata->Get<std::vector<int>>("MATIDS")),
      solve_(*matdata->Get<std::string>("SOLVE"))

{
}


/*----------------------------------------------------------------------*
 |  Constructor Material
 *----------------------------------------------------------------------*/
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
    ;
    const std::vector<int>* branchmatids = NULL;

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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::GeneralizedGenMax::ReadMaterialParameters(
    int& numbranch, const std::vector<int>*& matids, std::string& solve)
{
  numbranch = params_->numbranch_;
  matids = params_->matids_;
  solve = params_->solve_;

  return;
}


/*----------------------------------------------------------------------*/
/*-----------------VISCOBRANCH------------------------------------------*/

MAT::ELASTIC::PAR::ViscoBranch::ViscoBranch(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor
 *----------------------------------------------------------------------*/
MAT::ELASTIC::ViscoBranch::ViscoBranch(MAT::ELASTIC::PAR::ViscoBranch* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::ViscoBranch::ReadMaterialParameters(
    double& nummat, const std::vector<int>*& matids)
{
  nummat = params_->nummat_;
  matids = params_->matids_;

  return;
}

/*----------------------------------------------------------------------*/
/*-----------------VISCOPART--------------------------------------------*/
MAT::ELASTIC::PAR::ViscoPart::ViscoPart(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), tau_(matdata->GetDouble("TAU"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor
 *----------------------------------------------------------------------*/
MAT::ELASTIC::ViscoPart::ViscoPart(MAT::ELASTIC::PAR::ViscoPart* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::ViscoPart::ReadMaterialParameters(double& tau)
{
  tau = params_->tau_;
  return;
}
