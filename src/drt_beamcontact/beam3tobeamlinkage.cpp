/*----------------------------------------------------------------------*/
/*!
\file beam3tobeamlinkage.cpp

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_RCP.hpp>

#include "beam3tobeamlinkage.H"

#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3_base.H"

#include "../drt_fem_general/largerotations.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_factory.H"

BEAMINTERACTION::BeamToBeamLinkageType BEAMINTERACTION::BeamToBeamLinkageType::instance_;
BEAMINTERACTION::Beam3rLin2LinkageType BEAMINTERACTION::Beam3rLin2LinkageType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::ParObject* BEAMINTERACTION::Beam3rLin2LinkageType::Create( const std::vector<char> & data )
{
  BEAMINTERACTION::Beam3rLin2Linkage * my_beam3rlin2 = new BEAMINTERACTION::Beam3rLin2Linkage();
  my_beam3rlin2->Unpack(data);
  return my_beam3rlin2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamLinkage::BeamToBeamLinkage() :
    ParObject(),
    isinit_(false),
    issetup_(false),
    //element1_(Teuchos::null),
    //element2_(Teuchos::null),
    //xi1_(0.0),
    //xi2_(0.0),
    id_(-1),
    bspotpos1_(true),
    bspotpos2_(true),
    bspottriad1_(true),
    bspottriad2_(true),
    Lambdarel1_(true),
    Lambdarel2_(true)
{
  // empty constructor
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamLinkage::Init(
//    DRT::ELEMENTS::Beam3Base*  element1,
//    DRT::ELEMENTS::Beam3Base*  element2,
//    const double&  xi1,
//    const double&  xi2,
    const int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix<3,1> >& initpos,
    const std::vector<LINALG::Matrix<3,3> >& inittriad)
{
  issetup_ = false;

//  xi1_ = xi1;
//  xi2_ = xi2;

  id_ = id;
  eleids_ = eleids;

  bspotpos1_ = initpos[0];
  bspotpos2_ = initpos[1];


  /* the initial triads of the connecting element are chosen such that the first base
   * vector points in the direction of the distance vector of the two connection sites;
   * second and third base vector are arbitrarily constructed with help of "smallest rotation"
   * operation */
  LINALG::Matrix<3,3> linkeletriad(true);

  LINALG::Matrix<3,1> distvec(true);
  distvec.Update(1.0, GetBindSpotPos2(), -1.0, GetBindSpotPos1());

  LARGEROTATIONS::CalculateSRTriads<double>(distvec,inittriad[0],linkeletriad);

  LARGEROTATIONS::triadtoquaternion(linkeletriad,bspottriad1_);
  bspottriad2_=bspottriad1_;

  /* store relative rotation matrix between triads of connecting element and
   * the material triads of the "parent elements"; these remain constant over
   * the entire life of this connection */
  Lambdarel1_.MultiplyTN(inittriad[0],linkeletriad);
  Lambdarel2_.MultiplyTN(inittriad[1],linkeletriad);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamLinkage::Setup()
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamLinkage::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // isinit_
  AddtoPack(data,isinit_);
  // issetup
  AddtoPack(data,issetup_);
  // add id
  int id = Id();
  AddtoPack(data,id);

  // add eleids_
  AddtoPack(data,eleids_);
  // add involved procs
  AddtoPack(data,involvedprocs_);
  // bspotpos1_
  AddtoPack(data,bspotpos1_);
  // bspotpos2_
  AddtoPack(data,bspotpos2_);
  // bspottriad1_
  AddtoPack(data,bspottriad1_);
  // bspottriad2_
  AddtoPack(data,bspottriad2_);
  // Lambdarel1_
  AddtoPack(data,Lambdarel1_);
  // Lambdarel2_
  AddtoPack(data,Lambdarel2_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamLinkage::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // isinit_
  isinit_ = DRT::ParObject::ExtractInt(position,data);
  // issetup
  issetup_ = DRT::ParObject::ExtractInt(position,data);
  // id_
  ExtractfromPack(position,data,id_);

  // eleids_
  ExtractfromPack(position,data,eleids_);
  // involved procs
  ExtractfromPack(position,data,involvedprocs_);
  // bspotpos1_
  ExtractfromPack(position,data,bspotpos1_);
  // bspotpos2_
  ExtractfromPack(position,data,bspotpos2_);
  // bspottriad1_
  ExtractfromPack(position,data,bspottriad1_);
  // bspottriad2_
  ExtractfromPack(position,data,bspottriad2_);
  // Lambdarel1_
  ExtractfromPack(position,data,Lambdarel1_);
  // Lambdarel2_
  ExtractfromPack(position,data,Lambdarel2_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamLinkage::ResetState(
    std::vector<LINALG::Matrix<3,1> >& bspotpos,
    std::vector<LINALG::Matrix<3,3> >& bspottriad)
{
  CheckInitSetup();

  bspotpos1_=bspotpos[0];
  bspotpos2_=bspotpos[1];

  LINALG::TMatrix<double,3,3> currenttriad(true);
  currenttriad.Multiply(bspottriad[0],Lambdarel1_);
  LARGEROTATIONS::triadtoquaternion<double>(currenttriad,bspottriad1_);

  currenttriad.Clear();
  currenttriad.Multiply(bspottriad[1],Lambdarel2_);
  LARGEROTATIONS::triadtoquaternion<double>(currenttriad,bspottriad2_);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>
BEAMINTERACTION::BeamToBeamLinkage::Create()
{
  // for now, we always use a 2-noded linear Reissner element
  return Teuchos::rcp(new BEAMINTERACTION::Beam3rLin2Linkage());

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::Beam3rLin2Linkage::Beam3rLin2Linkage() :
    BeamToBeamLinkage(),
    linkele_(Teuchos::null)
{
  // empty constructor
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::Beam3rLin2Linkage::Setup()
{
  CheckInit();

  // call setup of base class first
  BeamToBeamLinkage::Setup();

  /* the idea is to use a beam element as auxiliary object that provides us with a
   * response force (and moment) depending on the position and orientation of the
   * two material cross-sections (binding spots) it is connected to;
   *
   * note: the element instance created in this way can only be used in a limited way
   *       because it is not embedded in a discretization. For example,
   *       Nodes() and other methods are not functional because the
   *       pointers to nodes are not set. Same for reference position of nodes via X() ...
   *
   *       We really only use it as a calculation routine for a sophisticated
   *       (displacement-reaction force) relation here! */
  linkele_ = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(-1,0));

  linkele_->SetCrossSec(4.751658e-06);
  linkele_->SetCrossSecShear(4.751658e-06);
  linkele_->SetMaterial(1);
  linkele_->SetIyy(4.4918e-11);
  linkele_->SetIzz(4.4918e-11);
  linkele_->SetIrr(8.9836e-11);
  linkele_->SetCenterlineHermite(false);

  // set dummy node Ids, in order to make NumNodes() method of element return the correct number of nodes
  int nodeids[2];
  for (unsigned int i=0; i<2; ++i) nodeids[i]=-1;
  linkele_->SetNodeIds(2,&nodeids[0]);


  // the triads at the two connection sites are chosen identical initially, so we only use the first one
  LINALG::Matrix<3,1> linkelerotvec(true);
  LARGEROTATIONS::quaterniontoangle(GetBindSpotQuaternion1(),linkelerotvec);

  std::vector<double> refpos(6,0.0);
  std::vector<double> refrotvec(6,0.0);

  for (unsigned int i=0; i<3; ++i)
  {
    refpos[i] = GetBindSpotPos1()(i);
    refpos[3+i] = GetBindSpotPos2()(i);

    refrotvec[i] = linkelerotvec(i);
    refrotvec[3+i] = linkelerotvec(i);
  }

  linkele_->SetUpReferenceGeometry<2,2,1>(refpos,refrotvec);

  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Beam3rLin2Linkage::Pack(DRT::PackBuffer& data) const
{
  CheckInitSetup();

  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Node
  BeamToBeamLinkage::Pack(data);

  // pack linker element
  if (linkele_!=Teuchos::null)
    linkele_->Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Beam3rLin2Linkage::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  BeamToBeamLinkage::Unpack(basedata);

  // Unpack data of sub material (these lines are copied from drt_element.cpp)
  std::vector<char> dataele;
  ExtractfromPack(position,data,dataele);
  if (dataele.size()>0)
  {
    DRT::ParObject* object = DRT::UTILS::Factory(dataele);  // Unpack is done here
    DRT::ELEMENTS::Beam3r* linkele = dynamic_cast<DRT::ELEMENTS::Beam3r*>(object);
    if (linkele==NULL)
      dserror("failed to unpack elastic material");
    linkele_ = Teuchos::rcp(linkele);
  }
  else linkele_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::Beam3rLin2Linkage::EvaluateForce(
    LINALG::TMatrix<double,6,1>& forcevec1,
    LINALG::TMatrix<double,6,1>& forcevec2)
{
  CheckInitSetup();

  LINALG::TMatrix<double,6,1> disp_totlag_centerline;

  for (unsigned int i=0; i<3; ++i)
  {
    disp_totlag_centerline(i) = GetBindSpotPos1()(i);
    disp_totlag_centerline(3+i) = GetBindSpotPos2()(i);
  }

  std::vector<LINALG::TMatrix<double,4,1> > Qnode;
  Qnode.push_back(GetBindSpotQuaternion1());
  Qnode.push_back(GetBindSpotQuaternion2());


  Epetra_SerialDenseVector force(12);

  linkele_->CalcInternalAndInertiaForcesAndStiff<2,2,1,double>(
      disp_totlag_centerline,
      Qnode,
      NULL,
      NULL,
      &force,
      NULL);

  for (unsigned int i=0; i<6; ++i)
  {
    forcevec1(i) = force(i);
    forcevec2(i) = force(6+i);
  }

  // *************************** DEBUG *******************************************
//  std::cout <<"\nBeam3rLin2Linkage::EvaluateForce: f_ele=\n";
//  force.Print(std::cout);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::Beam3rLin2Linkage::EvaluateStiff(
    LINALG::TMatrix<double,6,6>& stiffmat11,
    LINALG::TMatrix<double,6,6>& stiffmat12,
    LINALG::TMatrix<double,6,6>& stiffmat21,
    LINALG::TMatrix<double,6,6>& stiffmat22)
{
  CheckInitSetup();

  LINALG::TMatrix<double,6,1> disp_totlag_centerline;

  for (unsigned int i=0; i<3; ++i)
  {
    disp_totlag_centerline(i) = GetBindSpotPos1()(i);
    disp_totlag_centerline(3+i) = GetBindSpotPos2()(i);
  }

  std::vector<LINALG::TMatrix<double,4,1> > Qnode;
  Qnode.push_back(GetBindSpotQuaternion1());
  Qnode.push_back(GetBindSpotQuaternion2());


  Epetra_SerialDenseMatrix stiffmat(12,12);

  linkele_->CalcInternalAndInertiaForcesAndStiff<2,2,1,double>(
      disp_totlag_centerline,
      Qnode,
      &stiffmat,
      NULL,
      NULL,
      NULL);

  for (unsigned int i=0; i<6; ++i)
    for (unsigned int j=0; j<6; ++j)
    {
      stiffmat11(i,j) = stiffmat(i,j);
      stiffmat12(i,j) = stiffmat(i,6+j);
      stiffmat21(i,j) = stiffmat(6+i,j);
      stiffmat22(i,j) = stiffmat(6+i,6+j);
    }


  // *************************** DEBUG *******************************************
//  std::cout <<"\nBeam3rLin2Linkage::EvaluateStiff: stiff_ele=\n";
//  stiffmat.Print(std::cout);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::Beam3rLin2Linkage::EvaluateForceStiff(
    LINALG::TMatrix<double,6,1>& forcevec1,
    LINALG::TMatrix<double,6,1>& forcevec2,
    LINALG::TMatrix<double,6,6>& stiffmat11,
    LINALG::TMatrix<double,6,6>& stiffmat12,
    LINALG::TMatrix<double,6,6>& stiffmat21,
    LINALG::TMatrix<double,6,6>& stiffmat22)
{
  CheckInitSetup();

  LINALG::TMatrix<double,6,1> disp_totlag_centerline;

  for (unsigned int i=0; i<3; ++i)
  {
    disp_totlag_centerline(i) = GetBindSpotPos1()(i);
    disp_totlag_centerline(3+i) = GetBindSpotPos2()(i);
  }

  std::vector<LINALG::TMatrix<double,4,1> > Qnode;
  Qnode.push_back(GetBindSpotQuaternion1());
  Qnode.push_back(GetBindSpotQuaternion2());


  Epetra_SerialDenseVector force(12);
  Epetra_SerialDenseMatrix stiffmat(12,12);

  linkele_->CalcInternalAndInertiaForcesAndStiff<2,2,1,double>(
      disp_totlag_centerline,
      Qnode,
      &stiffmat,
      NULL,
      &force,
      NULL);

  for (unsigned int i=0; i<6; ++i)
  {
    forcevec1(i) = force(i);
    forcevec2(i) = force(6+i);

    for (unsigned int j=0; j<6; ++j)
    {
      stiffmat11(i,j) = stiffmat(i,j);
      stiffmat12(i,j) = stiffmat(i,6+j);
      stiffmat21(i,j) = stiffmat(6+i,j);
      stiffmat22(i,j) = stiffmat(6+i,6+j);
    }
  }

  // *************************** DEBUG *******************************************
//  std::cout <<"\nBeam3rLin2Linkage::EvaluateForceStiff: f_ele=\n";
//  force.Print(std::cout);
//  std::cout <<"\nstiff_ele=\n";
//  stiffmat.Print(std::cout);

  return true;
}
