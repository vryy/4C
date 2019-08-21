/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a truss element used as mechanical link
       between two beam elements

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#include "beam_link_truss.H"

#include "../drt_truss3/truss3.H"


#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_factory.H"

#include <Teuchos_RCP.hpp>
#include "beam_link.H"


BEAMINTERACTION::BeamLinkTrussType BEAMINTERACTION::BeamLinkTrussType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::ParObject* BEAMINTERACTION::BeamLinkTrussType::Create(const std::vector<char>& data)
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
      bspotforces_(2, LINALG::SerialDenseVector(true))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkTruss::BeamLinkTruss(const BEAMINTERACTION::BeamLinkTruss& old)
    : BEAMINTERACTION::BeamLinkPinJointed(old), bspotforces_(2, LINALG::SerialDenseVector(true))
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
    const std::vector<LINALG::Matrix<3, 1>>& initpos,
    const std::vector<LINALG::Matrix<3, 3>>& inittriad,
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
  int nodeids[2];
  for (unsigned int i = 0; i < 2; ++i) nodeids[i] = -1;
  linkele_->SetNodeIds(2, &nodeids[0]);


  std::vector<double> refpos(6, 0.0);
  std::vector<double> refrotvec_dummy(6, 0.0);

  for (unsigned int i = 0; i < 3; ++i)
  {
    refpos[i] = GetBindSpotPos1()(i);
    refpos[3 + i] = GetBindSpotPos2()(i);
  }

  linkele_->SetUpReferenceGeometry(refpos, refrotvec_dummy);

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::Pack(DRT::PackBuffer& data) const
{
  CheckInitSetup();

  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  BeamLinkPinJointed::Unpack(basedata);

  // Unpack data of sub material (these lines are copied from drt_element.cpp)
  std::vector<char> dataele;
  ExtractfromPack(position, data, dataele);
  if (dataele.size() > 0)
  {
    DRT::ParObject* object = DRT::UTILS::Factory(dataele);  // Unpack is done here
    DRT::ELEMENTS::Truss3* linkele = dynamic_cast<DRT::ELEMENTS::Truss3*>(object);
    if (linkele == NULL) dserror("failed to unpack Truss3 object within BeamLinkTruss");
    linkele_ = Teuchos::rcp(linkele);
  }
  else
    linkele_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::EvaluateForce(
    LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseVector& forcevec2)
{
  CheckInitSetup();

  LINALG::Matrix<6, 1> absolute_nodal_positions;

  FillStateVariablesForElementEvaluation(absolute_nodal_positions);

  LINALG::SerialDenseVector force(6, true);
  LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(absolute_nodal_positions, force, stiffmat);

  std::copy(&force(0), &force(0) + 3, &forcevec1(0));
  std::copy(&force(0) + 3, &force(0) + 6, &forcevec2(0));

  bspotforces_[0] = forcevec1;
  bspotforces_[1] = forcevec2;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::EvaluateStiff(LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  LINALG::Matrix<6, 1> absolute_nodal_positions;

  FillStateVariablesForElementEvaluation(absolute_nodal_positions);

  LINALG::SerialDenseVector force(6, true);
  LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(absolute_nodal_positions, force, stiffmat);

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
bool BEAMINTERACTION::BeamLinkTruss::EvaluateForceStiff(LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2, LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  LINALG::Matrix<6, 1> absolute_nodal_positions;

  FillStateVariablesForElementEvaluation(absolute_nodal_positions);

  LINALG::SerialDenseVector force(6, true);
  LINALG::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->CalcInternalForceStiffTotLag(absolute_nodal_positions, force, stiffmat);

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
void BEAMINTERACTION::BeamLinkTruss::ResetState(
    std::vector<LINALG::Matrix<3, 1>>& bspotpos, std::vector<LINALG::Matrix<3, 3>>& bspottriad)
{
  CheckInitSetup();

  BeamLinkPinJointed::ResetState(bspotpos, bspottriad);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::FillStateVariablesForElementEvaluation(
    LINALG::TMatrix<double, 6, 1>& absolute_nodal_positions) const
{
  for (unsigned int i = 0; i < 3; ++i)
  {
    absolute_nodal_positions(i) = GetBindSpotPos1()(i);
    absolute_nodal_positions(3 + i) = GetBindSpotPos2()(i);
  }
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
    int bspotid, LINALG::SerialDenseVector& bspotforce) const
{
  bspotforce = bspotforces_[bspotid];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetCurrentLinkerLength() const { return linkele_->Lcurr(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetInternalEnergy() const
{
  return linkele_->GetInternalEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::GetKineticEnergy() const { return 0.0; }
