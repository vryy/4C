/*!-----------------------------------------------------------------------------------------------------------
\file beam3_base.cpp

\brief base class for all beam elements

\level 2

\maintainer Christoph Meier
 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3_base.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(int id, int owner) :
DRT::Element(id,owner),
interface_ptr_(Teuchos::null),
sm_interface_ptr_(Teuchos::null)

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
void DRT::ELEMENTS::Beam3Base::SetStatMechParamsInterfacePtr()
{
  sm_interface_ptr_ = interface_ptr_->GetStatMechParamInterface();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Beam3Base::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<STATMECH::ParamsInterface> DRT::ELEMENTS::Beam3Base::StatMechParamsInterfacePtr() const
{
  return sm_interface_ptr_;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetRefPosAtXi(LINALG::Matrix<3,1>& refpos,
                                             const double& xi) const
{
  const int numclnodes = this->NumCenterlineNodes();
  int numnodalvalues = 0;
  this->HermiteCenterlineInterpolation() ? numnodalvalues = 2 : numnodalvalues=1;
  std::vector<double> refdofval;
  refdofval.reserve(3*numnodalvalues*numclnodes);
  const DRT::Node* const* nodes = this->Nodes();

  if (nodes[0]==NULL) dserror("Cannot get nodes of this element");

  // centerline nodes are always the first numclnodes nodes of the element
  for (int n=0; n<numclnodes; ++n)
  {
    std::vector<int> dofindices;

    this->PositionDofIndices(dofindices, *(nodes[n]));

    for (std::vector<int>::const_iterator it=dofindices.begin(); it!=dofindices.end(); ++it)
      refdofval.push_back((nodes[n]->X())[*it]);

    if (this->HermiteCenterlineInterpolation())
    {
      LINALG::Matrix<3,1> Tref_node;
      this->GetRefTangentAtNode(Tref_node,n);
      for (int d=0; d<3; ++d)
        refdofval.push_back(Tref_node(d));
    }
  }

  this->GetPosAtXi(refpos,xi,refdofval);
}

/*-----------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of          |
 | periodic boundary conditions; if two nodes within one element are se-       |
 | parated by a periodic boundary, one of them is shifted such that the final  |
 | distance in R^3 is the same as the initial distance in the periodic         |
 | space; the shift affects computation on element level within that very      |
 | iteration step, only (no change in global variables performed)              |
 *-----------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(
    std::vector<double>& disp, unsigned int nnode)
{
  /* note: crosslinker are set only if both bindingspots of the filaments that are
   * crosslinked lie inside the periodic bounding box. Therefore the crosslinker
   * elements need the shifted configuration on element level (in contrary to the
   * filament elements). As they were set in the bounding box they get not un-
   * shifted here, which is correct. They only get shifted if during crosslinker
   * diffusion one node exits the periodic bounding box.
   *
   * author Jonas Eichinger                                       08/16 */

  // get period length
  Teuchos::RCP<std::vector<double> > periodlength =
      StatMechParamsInterface().GetPeriodLength();

  // get dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // do nothing in case periodic boundary conditions are inactive
  if(periodlength->at(0) <= 0.0)
    return;

  /* get number of degrees of freedom per node; note:
   * the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  // loop through all nodes except for the first node which remains
  // fixed as reference node
  for(unsigned int i=1;i<nnode;i++)
    for(int dof= ndim - 1; dof > -1; dof--)
    {
      /* if the distance in some coordinate direction between some node and the
       * first node becomes smaller by adding or subtracting the period length,
       * the respective node has obviously been shifted due to periodic boundary
       * conditions and should be shifted back for evaluation of element
       * matrices and vectors; this way of detecting shifted nodes works as long
       * as the element length is smaller than half the periodic length*/
      if(std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) +
                   periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof])) <
         std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                   (Nodes()[0]->X()[dof] + disp[numdof*0+dof]))
      )
      {
        disp[numdof*i+dof] += periodlength->at(dof);
      }
      else if(std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                        periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof])) <
              std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                        (Nodes()[0]->X()[dof] + disp[numdof*0+dof]))
      )
      {
        disp[numdof*i+dof] -= periodlength->at(dof);
      }
    }

  // that's it
  return;

} // NodeShift()

