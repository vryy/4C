/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a truss element used as mechanical link
       between two beam elements

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_beaminteraction_link_truss.hpp"

#include "4C_beaminteraction_link.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_truss3.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


BEAMINTERACTION::BeamLinkTrussType BEAMINTERACTION::BeamLinkTrussType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::COMM::ParObject* BEAMINTERACTION::BeamLinkTrussType::Create(const std::vector<char>& data)
{
  BEAMINTERACTION::BeamLinkTruss* my_truss_linker = new BEAMINTERACTION::BeamLinkTruss();
  my_truss_linker->Unpack(data);
  return my_truss_linker;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkTruss::BeamLinkTruss()
    : BeamLinkPinJointed(),
      linkele_(Teuchos::null),
      bspotforces_(2, CORE::LINALG::SerialDenseVector(true))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkTruss::BeamLinkTruss(const BEAMINTERACTION::BeamLinkTruss& old)
    : BEAMINTERACTION::BeamLinkPinJointed(old),
      bspotforces_(2, CORE::LINALG::SerialDenseVector(true))
{
  if (linkele_ != Teuchos::null)
    linkele_ =
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Truss3>(Teuchos::rcp(old.linkele_->Clone(), true));
  else
    linkele_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLink> BEAMINTERACTION::BeamLinkTruss::Clone() const
{
  Teuchos::RCP<BEAMINTERACTION::BeamLinkTruss> newlinker =
      Teuchos::rcp(new BEAMINTERACTION::BeamLinkTruss(*this));
  return newlinker;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::Init(int id, const std::vector<std::pair<int, int>>& eleids,
    const std::vector<CORE::LINALG::Matrix<3, 1>>& initpos,
    const std::vector<CORE::LINALG::Matrix<3, 3>>& inittriad,
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset)
{
  issetup_ = false;

  BeamLinkPinJointed::Init(id, eleids, initpos, inittriad, linkertype, timelinkwasset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::Setup(const int matnum)
{
  CheckInit();

  // call setup of base class first
  BeamLinkPinJointed::Setup(matnum);

  /* the idea is to use a truss element as auxiliary object that provides us with a
   * response force depending on the position of the two material points on the
   *  parent elements (i.e. binding spots) it is connected to;
   *
   * note: the element instance created in this way can only be used in a limited way
   *       because it is not embedded in a discretization. For example,
   *       Nodes() and other methods are not functional because the
   *       pointers to nodes are not set. Same for reference position of nodes via X() ...
   *
   *       We really only use it as a calculation routine for a sophisticated
   *       (displacement-reaction force) relation here! */
  linkele_ = Teuchos::rcp(new DRT::ELEMENTS::Truss3(-1, 0));

  // set material
  linkele_->SetMaterial(matnum);

  // set cross-section area Fixme hard-coded dummy value for now
  linkele_->SetCrossSec(1.0);

  // set dummy node Ids, in order to make NumNodes() method of element return the correct number of
  // nodes
  int nodeids[] = {-1, -1};
  linkele_->SetNodeIds(2, nodeids);

  std::vector<double> refpos(6, 0.0);

  for (unsigned int i = 0; i < 3; ++i)
  {
    refpos[i] = GetBindSpotPos1()(i);
    refpos[3 + i] = GetBindSpotPos2()(i);
  }

  linkele_->SetUpReferenceGeometry(refpos);

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::Pack(CORE::COMM::PackBuffer& data) const
{
  CheckInitSetup();

  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class
  BeamLinkPinJointed::Pack(data);

  // pack linker element
  if (linkele_ != Teuchos::null) linkele_->Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  BeamLinkPinJointed::Unpack(basedata);

  // Unpack data of sub material (these lines are copied from element.cpp)
  std::vector<char> dataele;
  ExtractfromPack(position, data, dataele);
  if (dataele.size() > 0)
  {
    CORE::COMM::ParObject* object = CORE::COMM::Factory(dataele);  // Unpack is done here
    DRT::ELEMENTS::Truss3* linkele = dynamic_cast<DRT::ELEMENTS::Truss3*>(object);
    if (linkele == nullptr) FOUR_C_THROW("failed to unpack Truss3 object within BeamLinkTruss");
    linkele_ = Teuchos::rcp(linkele);
  }
  else
    linkele_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::EvaluateForce(
    CORE::LINALG::SerialDenseVector& forcevec1, CORE::LINALG::SerialDenseVector& forcevec2)
{
  CheckInitSetup();

  std::map<std::string, std::vector<double>> ele_state;
  GetDispForElementEvaluation(ele_state);

  CORE::LINALG::SerialDenseVector force(6, true);
  CORE::LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(ele_state, force, stiffmat);

  std::copy(&force(0), &force(0) + 3, &forcevec1(0));
  std::copy(&force(0) + 3, &force(0) + 6, &forcevec2(0));

  bspotforces_[0] = forcevec1;
  bspotforces_[1] = forcevec2;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::EvaluateStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
    CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
    CORE::LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  std::map<std::string, std::vector<double>> ele_state;
  GetDispForElementEvaluation(ele_state);

  CORE::LINALG::SerialDenseVector force(6, true);
  CORE::LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(ele_state, force, stiffmat);

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      stiffmat11(i, j) = stiffmat(i, j);
      stiffmat12(i, j) = stiffmat(i, 3 + j);
      stiffmat21(i, j) = stiffmat(3 + i, j);
      stiffmat22(i, j) = stiffmat(3 + i, 3 + j);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::EvaluateForceStiff(CORE::LINALG::SerialDenseVector& forcevec1,
    CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
    CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
    CORE::LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  std::map<std::string, std::vector<double>> ele_state;
  GetDispForElementEvaluation(ele_state);

  CORE::LINALG::SerialDenseVector force(6, true);
  CORE::LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(ele_state, force, stiffmat);

  std::copy(&force(0), &force(0) + 3, &forcevec1(0));
  std::copy(&force(0) + 3, &force(0) + 6, &forcevec2(0));

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      stiffmat11(i, j) = stiffmat(i, j);
      stiffmat12(i, j) = stiffmat(i, 3 + j);
      stiffmat21(i, j) = stiffmat(3 + i, j);
      stiffmat22(i, j) = stiffmat(3 + i, 3 + j);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::ResetState(std::vector<CORE::LINALG::Matrix<3, 1>>& bspotpos,
    std::vector<CORE::LINALG::Matrix<3, 3>>& bspottriad)
{
  CheckInitSetup();

  BeamLinkPinJointed::ResetState(bspotpos, bspottriad);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::FillStateVariablesForElementEvaluation(
    CORE::LINALG::Matrix<6, 1, double>& absolute_nodal_positions) const
{
  for (unsigned int i = 0; i < 3; ++i)
  {
    absolute_nodal_positions(i) = GetBindSpotPos1()(i);
    absolute_nodal_positions(3 + i) = GetBindSpotPos2()(i);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::GetDispForElementEvaluation(
    std::map<std::string, std::vector<double>>& ele_state) const
{
  const auto ref_position = linkele_->X();
  std::vector<double> disp(6, 0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    disp[i] = GetBindSpotPos1()(i) - ref_position(i, 0);
    disp[3 + i] = GetBindSpotPos2()(i) - ref_position(i + 3, 0);
  }

  ele_state.emplace(std::make_pair("disp", disp));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::ScaleLinkerReferenceLength(double scalefac)
{
  linkele_->ScaleReferenceLength(scalefac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::GetBindingSpotForce(
    int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const
{
  bspotforce = bspotforces_[bspotid];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetCurrentLinkerLength() const
{
  CORE::LINALG::Matrix<6, 1> xcurr;

  FillStateVariablesForElementEvaluation(xcurr);

  CORE::LINALG::Matrix<6, 1> curr_nodal_coords;
  curr_nodal_coords(0) = xcurr(0) - xcurr(3);
  curr_nodal_coords(1) = xcurr(1) - xcurr(4);
  curr_nodal_coords(2) = xcurr(2) - xcurr(5);
  curr_nodal_coords(3) = curr_nodal_coords(0);
  curr_nodal_coords(4) = curr_nodal_coords(1);
  curr_nodal_coords(5) = curr_nodal_coords(2);

  return linkele_->CurrLength(curr_nodal_coords);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetInternalEnergy() const
{
  return linkele_->GetInternalEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetKineticEnergy() const { return 0.0; }

FOUR_C_NAMESPACE_CLOSE
