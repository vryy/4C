/*!-----------------------------------------------------------------------------------------------------------
\file beam3_base.cpp

\brief base class for all beam elements

\level 2

\maintainer Christoph Meier
 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3_base.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_inpar/inpar_statmech.H"

#include <Sacado.hpp>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(int id, int owner) :
DRT::Element(id,owner),
interface_ptr_(Teuchos::null),
sm_interface_ptr_(Teuchos::null)

{
  // todo: this is a temporary hack, should of course be set from outside
  bspotposxi_.push_back(-1.0);
//  bspotposxi_.push_back(1.0);
  bspotstatus_[0] = -1;
//  bspotstatus_[1] = -1;

  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(const DRT::ELEMENTS::Beam3Base& old) :
 DRT::Element(old)
{
  // todo: this is a temporary hack, should of course be set from outside
  bspotposxi_.push_back(-1.0);
//  bspotposxi_.push_back(1.0);
  bspotstatus_[0] = -1;
//  bspotstatus_[1] = -1;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  // bspotposxi_
  AddtoPack(data,bspotposxi_);
  // bspotstatus_
  AddtoPack(data,bspotstatus_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // bspotposxi_
  ExtractfromPack(position,data,bspotposxi_);
  // bspotstatus_
  ExtractfromPack(position,data,bspotstatus_);

  return;
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
  const int numnodalvalues = this->HermiteCenterlineInterpolation() ? 2 : 1;

  std::vector<double> zerovec;
  zerovec.resize(3*numnodalvalues*numclnodes);

  this->GetPosAtXi(refpos,xi,zerovec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetDampingCoefficients(LINALG::Matrix<3,1>& gamma) const
{
  /* These are coefficients for a straight cylindrical rod taken from
   * Howard, p. 107, table 6.2. The order is as follows:
   * (0) damping of translation parallel to axis,
   * (1) damping of translation orthogonal to axis,
   * (2) damping of rotation around its own axis */

  gamma(0) = 2*PI*StatMechParamsInterface().GetEta();
  gamma(1) = 4*PI*StatMechParamsInterface().GetEta();
  gamma(2) = 4*PI*StatMechParamsInterface().GetEta() * std::sqrt(4*this->Iyy()/PI);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<unsigned int ndim, typename T> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity(Teuchos::ParameterList&      params,  //!<parameter list
                                                const LINALG::TMatrix<T,ndim,1>&  evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::TMatrix<T,ndim,1>&        velbackground,  //!< velocity of background fluid
                                                LINALG::TMatrix<T,ndim,ndim>&     velbackgroundgrad) const//!<gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in direction params.get<int>("DBCDISPDIR",0) and orthogonal to z-axis.
   * In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  //velocity at upper boundary of domain
  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  double time = -1.0;
  double dt = 1000;
  if (IsParamsInterface())
  {
    time = ParamsInterface().GetTotalTime();
    dt = ParamsInterface().GetDeltaTime();
  }
  else
  {
    time = params.get<double>("total time",-1.0);
    dt = params.get<double>("delta time");
  }
  double starttime = StatMechParamsInterface().GetStartTimeAction();
  double shearamplitude = StatMechParamsInterface().GetShearAmplitude();
  int curvenumber = StatMechParamsInterface().GetCurveNumber() - 1;
  int dbcdispdir = StatMechParamsInterface().GetDbcDispDir() - 1;

  Teuchos::RCP<std::vector<double> > periodlength = StatMechParamsInterface().GetPeriodLength();
  INPAR::STATMECH::DBCType dbctype = StatMechParamsInterface().GetDbcType();
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
     dbctype==INPAR::STATMECH::dbctype_affineshear||
     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
    shearflow = true;

  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
  if(periodlength->at(0) > 0.0)
    if(shearflow && time>starttime && fabs(time-starttime)>dt/1e4 && curvenumber >=  0 && dbcdispdir >= 0 )
    {
      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1))[1];

      //compute background velocity
      velbackground(dbcdispdir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;

      //compute gradient of background velocity
      velbackgroundgrad(dbcdispdir,ndim-1) = uppervel / periodlength->at(ndim-1);
    }
}

/*-----------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of          |
 | periodic boundary conditions; if two nodes within one element are se-       |
 | parated by a periodic boundary, one of them is shifted such that the final  |
 | distance in R^3 is the same as the initial distance in the periodic         |
 | space; the shift affects computation on element level within that very      |
 | iteration step, only (no change in global variables performed)              |
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(std::vector<double>& disp,
    const std::vector<double>& periodlength, unsigned int nnode) const
{
  /* note: crosslinker are set only if both bindingspots of the filaments that are
   * crosslinked lie inside the periodic bounding box. Therefore the crosslinker
   * elements need the shifted configuration on element level (in contrary to the
   * filament elements). As they were set in the bounding box they get not un-
   * shifted here, which is correct. They only get shifted if during crosslinker
   * diffusion one node exits the periodic bounding box.
   *
   * author Jonas Eichinger                                       08/16 */

  // get dimension
  const int ndim = 3;

  // do nothing in case periodic boundary conditions are inactive
  if(periodlength[0] <= 0.0)
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
                   periodlength[dof] - (Nodes()[0]->X()[dof]+disp[numdof*0+dof])) <
         std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                   (Nodes()[0]->X()[dof] + disp[numdof*0+dof]))
      )
      {
        disp[numdof*i+dof] += periodlength[dof];
      }
      else if(std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                        periodlength[dof] - (Nodes()[0]->X()[dof]+disp[numdof*0+dof])) <
              std::fabs((Nodes()[i]->X()[dof]+disp[numdof*i+dof]) -
                        (Nodes()[0]->X()[dof] + disp[numdof*0+dof]))
      )
      {
        disp[numdof*i+dof] -= periodlength[dof];
      }
    }

  // that's it
  return;

} // NodeShift()

//! shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(std::vector<double>& disp,
                                        const unsigned int nnode) const
{
  this->UnShiftNodePosition(disp,*(StatMechParamsInterface().GetPeriodLength()),nnode);
}
/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetPosOfBindingSpot(LINALG::Matrix<3,1>&       pos,
                                                   std::vector<double>&       disp,
                                                   const int&                 bspotlocn,
                                                   const std::vector<double>& periodlength) const
{
  // unshift node position to get correct postion at xi
  UnShiftNodePosition(disp,periodlength,NumCenterlineNodes());

  const double xi = bspotposxi_[bspotlocn];
  // get position
  GetPosAtXi(pos,xi,disp);

  // check if xi lies outside the periodic box, if it does, shift it back in
  for(int dim=0; dim<3; ++dim)
  {
    if(periodlength[dim]>0)
    {
      if(pos(dim) < 0)
        pos(dim) += periodlength[dim];
      else if(pos(dim) > periodlength[dim])
        pos(dim) -= periodlength[dim];
    }
  }

  // that is it
  return;

}

// explicit template instantiations
template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3,double>(Teuchos::ParameterList&,
                                                                        const LINALG::TMatrix<double,3,1>&,
                                                                        LINALG::TMatrix<double,3,1>&,
                                                                        LINALG::TMatrix<double,3,3>&) const;
template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3,Sacado::Fad::DFad<double> >(Teuchos::ParameterList&,
                                                                        const LINALG::TMatrix<Sacado::Fad::DFad<double>,3,1>&,
                                                                        LINALG::TMatrix<Sacado::Fad::DFad<double>,3,1>&,
                                                                        LINALG::TMatrix<Sacado::Fad::DFad<double>,3,3>&) const;

