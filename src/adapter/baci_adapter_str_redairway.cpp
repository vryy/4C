/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for coupling with reduced-D airway trees

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_adapter_str_redairway.hpp"

#include "baci_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
ADAPTER::StructureRedAirway::StructureRedAirway(Teuchos::RCP<Structure> stru)
    : StructureWrapper(stru)
{
  //----------------------------------------------------------------------
  // make sure
  if (structure_ == Teuchos::null)
    FOUR_C_THROW("Failed to create the underlying structural adapter");

  // first get all Neumann conditions on structure

  std::vector<DRT::Condition*> surfneumcond;
  Discretization()->GetCondition("SurfaceNeumann", surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) FOUR_C_THROW("no Neumann conditions on structure");

  // now filter those Neumann conditions that are due to the coupling
  std::vector<int> tmp;
  for (unsigned int i = 0; i < numneumcond; ++i)
  {
    DRT::Condition* actcond = surfneumcond[i];
    if (actcond->Type() == DRT::Condition::RedAirwayTissue)
    {
      int condID = *actcond->Get<int>("coupling id");
      coupcond_[condID] = actcond;
      tmp.push_back(condID);
      vn_[condID] = 0.0;
      vnp_[condID] = 0.0;
    }
  }
  unsigned int numcond = tmp.size();
  if (numcond == 0) FOUR_C_THROW("no coupling conditions found");
  coupmap_ =
      Teuchos::rcp(new Epetra_Map(tmp.size(), tmp.size(), tmp.data(), 0, Discretization()->Comm()));
}


/*======================================================================*/
void ADAPTER::StructureRedAirway::SetPressure(Teuchos::RCP<Epetra_Vector> couppres)
{
  const Epetra_BlockMap& condmap = couppres->Map();

  for (int i = 0; i < condmap.NumMyElements(); ++i)
  {
    int condID = condmap.GID(i);
    DRT::Condition* cond = coupcond_[condID];
    std::vector<double> newval(6, 0.0);
    newval[0] = (*couppres)[i];
    cond->Add("val", newval);
  }
}


/*======================================================================*/
void ADAPTER::StructureRedAirway::CalcVol(std::map<int, double>& V)
{
  if (!(Discretization()->Filled())) FOUR_C_THROW("FillComplete() was not called");

  Teuchos::ParameterList params;
  params.set("action", "calc_struct_constrvol");

  // set displacements
  Discretization()->ClearState();
  Discretization()->SetState("displacement", Dispnp());

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (int i = 0; i < coupmap_->NumMyElements(); ++i)
  {
    int condID = coupmap_->GID(i);
    DRT::Condition& cond = *(coupcond_[condID]);
    double tmp = 0.;
    params.set("ConditionID", condID);
    params.set<Teuchos::RCP<DRT::Condition>>("condition", Teuchos::rcp(&cond, false));

    // define element matrices and vectors
    CORE::LINALG::SerialDenseMatrix elematrix1;
    CORE::LINALG::SerialDenseMatrix elematrix2;
    CORE::LINALG::SerialDenseVector elevector1;
    CORE::LINALG::SerialDenseVector elevector2;
    CORE::LINALG::SerialDenseVector elevector3;

    std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*Discretization(), lm, lmowner, lmstride);

      // reshape element matrices and vectors and init to zero
      elevector3.size(1);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params, *Discretization(), lm, elematrix1, elematrix2,
          elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("error while evaluating elements");

      // "assembly
      tmp += elevector3[0];
    }
    double vol = 0.;
    Discretization()->Comm().SumAll(&tmp, &vol, 1);
    if (vol < 0.0) vol *= -1.0;
    V[condID] = vol;
  }
}


/*======================================================================*/
void ADAPTER::StructureRedAirway::InitVol()
{
  CalcVol(vn_);

  if (Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "------------------------ Initial tissue volumes ----------------------"
              << std::endl;
    for (unsigned int i = 0; i < vn_.size(); ++i)
      std::cout << "ID:  " << coupmap_->GID(i) << "     V:  " << vn_[coupmap_->GID(i)] << std::endl;
    std::cout << "----------------------------------------------------------------------"
              << std::endl;
  }
}


/*======================================================================*/
void ADAPTER::StructureRedAirway::CalcFlux(
    Teuchos::RCP<Epetra_Vector> coupflux, Teuchos::RCP<Epetra_Vector> coupvol, double dt)
{
  CalcVol(vnp_);

  for (int i = 0; i < coupmap_->NumMyElements(); ++i)
  {
    int coupID = coupmap_->GID(i);
    int lid = coupflux->Map().LID(coupID);
    (*coupflux)[lid] = (vnp_[coupID] - vn_[coupID]) / dt;
    (*coupvol)[lid] = vnp_[coupID];
  }
}


/*======================================================================*/
void ADAPTER::StructureRedAirway::Update()
{
  StructureWrapper::Update();

  for (int i = 0; i < coupmap_->NumMyElements(); ++i)
  {
    int coupID = coupmap_->GID(i);
    vn_[coupID] = vnp_[coupID];
  }
}

FOUR_C_NAMESPACE_CLOSE
