/*!-----------------------------------------------------------------------------------------------------------
\file beam3_base.cpp

\brief base class for all beam elements

\level 2

\maintainer Christoph Meier
 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3_base.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(int id, int owner) :
DRT::Element(id,owner)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(const DRT::ELEMENTS::Beam3Base& old) :
 DRT::Element(old)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ =
        Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>
        (p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> >("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Beam3Base::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetAdditiveDofGIDs(
    const DRT::Discretization& discret,
    const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->PositionDofIndices(dofindices,node);
  this->TangentDofIndices(dofindices,node);
  this->Rotation1DDofIndices(dofindices,node);
  this->TangentLengthDofIndices(dofindices,node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i=0; i<dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0,&node,dofindices[i]));

  return dofgids;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetRotVecDofGIDs(
    const DRT::Discretization& discret,
    const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->RotationVecDofIndices(dofindices,node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i=0; i<dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0,&node,dofindices[i]));

  return dofgids;
}
