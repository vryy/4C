/*----------------------------------------------------------------------*/
/*!
\file beam_link_beam3r_lin2_pinjointed.cpp

\brief Wrapper for a linear Reissner beam element used as mechanical pin joint
       between two other beam elements

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------*/

#include "beam_link_beam3r_lin2_pinjointed.H"

#include "../drt_beam3/beam3r.H"

#include "../drt_fem_general/largerotations.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_factory.H"

#include <Teuchos_RCP.hpp>
#include "beam_link.H"


BEAMINTERACTION::BeamLinkBeam3rLin2PinJointedType BEAMINTERACTION::BeamLinkBeam3rLin2PinJointedType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::ParObject* BEAMINTERACTION::BeamLinkBeam3rLin2PinJointedType::Create( const std::vector<char> & data )
{
  BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed * my_beam3rlin2 = new BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed();
  my_beam3rlin2->Unpack(data);
  return my_beam3rlin2;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::BeamLinkBeam3rLin2PinJointed() :
    BeamLinkPinJointed(),
    triad_( true ),
    linkele_( Teuchos::null ),
    bspotforces_( 2, LINALG::SerialDenseVector(true) )
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::BeamLinkBeam3rLin2PinJointed(
    const BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed & old) :
    BEAMINTERACTION::BeamLinkPinJointed(old),
    triad_( old.triad_ ),
    bspotforces_( 2, LINALG::SerialDenseVector(true) )
{
  if ( linkele_ != Teuchos::null )
    linkele_ =  Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Beam3r>( Teuchos::rcp( old.linkele_->Clone(), true ) );
  else
    linkele_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLink> BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::Clone() const
{
  Teuchos::RCP<BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed> newlinker =
      Teuchos::rcp( new BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed( *this ) );
  return newlinker;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::Init(
    int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix<3,1> >& initpos,
    const std::vector<LINALG::Matrix<3,3> >& inittriad,
    double timelinkwasset)
{
  issetup_ = false;

  BeamLinkPinJointed::Init( id, eleids, initpos, inittriad, timelinkwasset);

  // *** initialization of the two triads of the connecting element ***
  /* they are determined such that:
   * - the first base vector points in the direction of the distance vector
   *   of the two connection sites; (will be axis of connecting element)
   * - second and third base vector are arbitrarily constructed from cross-product
   *   of first base vector with either first or second base vector of global
   *   coordinate system; this avoids any singularities */
  LINALG::Matrix<3,3> linkeletriad(true);
  LINALG::Matrix<3,1> distvec(true);

  distvec.Update( 1.0, GetBindSpotPos2(), -1.0, GetBindSpotPos1() );

  // feasibility check regarding coinciding connection sites
  if ( distvec.Norm2() < 1e-12 )
  {
    std::cout << "\nBeamLinkPinJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    initpos[0].Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    initpos[1].Print(std::cout);

    dserror("Initialization of BeamLinkPinJointed failed because the two given binding "
        "spot positions are almost identical, i.e. extremely short linker!");
  }

  // FIRST base vector
  distvec.Scale( 1.0 / distvec.Norm2() );

  std::copy( distvec.A(), distvec.A() + 3, &linkeletriad( 0, 0 ) );

  // SECOND base vector
  // check included angle of desired crosslinker axis (normalized distvec = first
  // base vector) and (1,0,0), i.e. scalar product which in this case simplifies to
  // first component of distvec
  LINALG::Matrix<3,1> unit_vector_global_x(true), unit_vector_global_y(true);
  unit_vector_global_x(0) = 1.0;
  unit_vector_global_y(1) = 1.0;

  const double scalarproduct = distvec(0);

  LINALG::Matrix<3,1> second_base_vecor_linkerele(true);

  // is included angle smaller than 45° ? then avoid singularity at angle=0° ...
  if ( std::abs(scalarproduct) > 0.5 * std::sqrt(2) )
  {
    second_base_vecor_linkerele.CrossProduct( distvec, unit_vector_global_y );
  }
  else
  {
    second_base_vecor_linkerele.CrossProduct( distvec, unit_vector_global_x );
  }

  // feasibility check
  if ( second_base_vecor_linkerele.Norm2() < 1e-12 )
  {
    std::cout << "\nBeamLinkPinJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    initpos[0].Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    initpos[1].Print(std::cout);

    std::cout << "\ndistvec = ";
    distvec.Print(std::cout);
    std::cout << "\nsecond_base_vecor_linkerele = ";
    second_base_vecor_linkerele.Print(std::cout);

    dserror("Initialization of BeamLinkPinJointed failed because the second base vector of the"
        "linker element's triad has almost length zero!");
  }
  else
  {
    second_base_vecor_linkerele.Scale( 1.0 / second_base_vecor_linkerele.Norm2() );
  }


  // THIRD base vector to complete orthonormal triad
  LINALG::Matrix<3,1> third_base_vecor_linkerele(true);
  third_base_vecor_linkerele.CrossProduct( distvec, second_base_vecor_linkerele );

  // feasibility check
  if ( std::abs( third_base_vecor_linkerele.Norm2() - 1.0 ) > 1e-12 )
  {
    std::cout << "\nBeamLinkPinJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    initpos[0].Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    initpos[1].Print(std::cout);

    std::cout << "\ndistvec = ";
    distvec.Print(std::cout);
    std::cout << "\nsecond_base_vecor_linkerele = ";
    second_base_vecor_linkerele.Print(std::cout);
    std::cout << "\nthird_base_vecor_linkerele = ";
    third_base_vecor_linkerele.Print(std::cout);

    dserror("Initialization of BeamLinkRigidJointed failed because the third base vector of the"
        "linker element's triad is no unit vector!");
  }


  /* store the initial triads as quaternions in class variables for the subsequent
   * use in setup of reference configuration of the connecting element */
  std::copy( second_base_vecor_linkerele.A(), second_base_vecor_linkerele.A() + 3, &linkeletriad( 0, 1 ) );
  std::copy( third_base_vecor_linkerele.A(), third_base_vecor_linkerele.A() + 3, &linkeletriad( 0, 2 ) );

  LARGEROTATIONS::triadtoquaternion( linkeletriad, triad_ );

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::Setup( const int matnum )
{
  CheckInit();

  // call setup of base class first
  BeamLinkPinJointed::Setup( matnum );

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

  // set material
  linkele_->SetMaterial(matnum);

  // Todo @grill: safety check for proper material type (done on element anyway, but do it here as well)?!

  linkele_->SetCenterlineHermite(false);

  // set dummy node Ids, in order to make NumNodes() method of element return the correct number of nodes
  int nodeids[2];
  for ( unsigned int i = 0; i < 2; ++i )
   nodeids[i] = -1;
  linkele_->SetNodeIds( 2, &nodeids[0] );


  // the triads at the two connection sites are chosen identical initially, so we only use the first one
  LINALG::Matrix< 3, 1 > linkelerotvec(true);
  LARGEROTATIONS::quaterniontoangle( triad_, linkelerotvec );

  std::vector<double> refpos(6,0.0);
  std::vector<double> refrotvec(6,0.0);

  for ( unsigned int i = 0; i < 3 ; ++i )
  {
    refpos[i] = GetBindSpotPos1()(i);
    refpos[3 + i] = GetBindSpotPos2()(i);

    refrotvec[i] = linkelerotvec(i);
    refrotvec[3 + i] = linkelerotvec(i);
  }

  linkele_->SetUpReferenceGeometry< 2, 2, 1 >( refpos, refrotvec );

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::Pack(DRT::PackBuffer& data) const
{
  CheckInitSetup();

  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack( data, type );
  // add base class
  BeamLinkPinJointed::Pack(data);

  // pack linker element
  if ( linkele_ != Teuchos::null )
    linkele_->Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  BeamLinkPinJointed::Unpack(basedata);

  // Unpack data of sub material (these lines are copied from drt_element.cpp)
  std::vector<char> dataele;
  ExtractfromPack( position, data, dataele );
  if ( dataele.size() > 0 )
  {
    DRT::ParObject* object = DRT::UTILS::Factory(dataele);  // Unpack is done here
    DRT::ELEMENTS::Beam3r* linkele = dynamic_cast<DRT::ELEMENTS::Beam3r*>(object);
    if (linkele==NULL)
      dserror("failed to unpack Beam3r object within BeamLinkBeam3rLin2PinJointed");
    linkele_ = Teuchos::rcp(linkele);
  }
  else linkele_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::EvaluateForce(
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2)
{
  CheckInitSetup();

  LINALG::TMatrix<double,6,1> disp_totlag_centerline;
  std::vector<LINALG::TMatrix<double,4,1> > Qnode;

  FillStateVariablesForElementEvaluation(
      disp_totlag_centerline,
      Qnode);

  LINALG::SerialDenseVector force(12,true);

  linkele_->CalcInternalAndInertiaForcesAndStiff< 2, 2, 1>(
      disp_totlag_centerline,
      Qnode,
      NULL,
      NULL,
      &force,
      NULL);

  // Todo maybe we can avoid this copy by setting up 'force' as a view on the
  //      two separate force vectors ?
  std::copy( &force(0), &force(0) + 3, &forcevec1(0) );
  std::copy( &force(0) + 6, &force(0) + 9, &forcevec2(0) );

  bspotforces_[0] = forcevec1;
  bspotforces_[1] = forcevec2;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::EvaluateStiff(
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  LINALG::TMatrix<double,6,1> disp_totlag_centerline;
  std::vector<LINALG::TMatrix<double,4,1> > Qnode;

  FillStateVariablesForElementEvaluation(
      disp_totlag_centerline,
      Qnode);

  LINALG::SerialDenseMatrix stiffmat(12,12,true);

  linkele_->CalcInternalAndInertiaForcesAndStiff<2,2,1>(
      disp_totlag_centerline,
      Qnode,
      &stiffmat,
      NULL,
      NULL,
      NULL);

  // Todo the linearization is incomplete yet. fix this or delete related code and
  // resort to truss linker element
  dserror("we miss stiffness contributions from rotation of the nodal triads here! "
      "implement the transformation matrices that describe the dependency of triad "
      "rotation on current position of binding spots, i.e. nodal positions and "
      "tangents!");

  // Todo can we use std::copy here or even set up 'stiffmat' as a view on the
  //      four individual sub-matrices ?
  for ( unsigned int i = 0; i < 3; ++i )
  {
    for ( unsigned int j = 0; j < 3; ++j )
    {
      stiffmat11( i, j) = stiffmat( i, j);
      stiffmat12( i, j) = stiffmat( i, 6 + j );
      stiffmat21( i, j) = stiffmat( 6 + i, j );
      stiffmat22( i ,j) = stiffmat( 6 + i, 6 + j );
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::EvaluateForceStiff(
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{
  CheckInitSetup();

  LINALG::TMatrix< double, 6, 1 > disp_totlag_centerline;
  std::vector<LINALG::TMatrix< double, 4, 1 > > Qnode;

  FillStateVariablesForElementEvaluation(
      disp_totlag_centerline,
      Qnode);

  LINALG::SerialDenseVector force( 12, true );
  LINALG::SerialDenseMatrix stiffmat( 12, 12, true );

  linkele_->CalcInternalAndInertiaForcesAndStiff< 2, 2, 1 >(
      disp_totlag_centerline,
      Qnode,
      &stiffmat,
      NULL,
      &force,
      NULL);

  std::copy( &force(0), &force(0) + 3, &forcevec1(0) );
  std::copy( &force(0) + 6, &force(0) + 9, &forcevec2(0) );

  // Todo the linearization is incomplete yet. fix this or delete related code and
  // resort to truss linker element
  dserror("we miss stiffness contributions from rotation of the nodal triads here! "
      "implement the transformation matrices that describe the dependency of triad "
      "rotation on current position of binding spots, i.e. nodal positions and "
      "tangents!");

  for ( unsigned int i = 0; i < 3; ++i )
  {
    for ( unsigned int j = 0; j < 3; ++j )
    {
      stiffmat11( i, j ) = stiffmat( i, j );
      stiffmat12( i, j ) = stiffmat( i, 6 + j );
      stiffmat21( i, j ) = stiffmat( 6 + i, j );
      stiffmat22( i, j ) = stiffmat( 6 + i, 6 + j );
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::ResetState(
    std::vector<LINALG::Matrix<3,1> >& bspotpos,
    std::vector<LINALG::Matrix<3,3> >& bspottriad)
{
  CheckInitSetup();

  BeamLinkPinJointed::ResetState( bspotpos, bspottriad );

  // *** re-initialization of the two triads of the connecting element ***

  /* the idea is that for this pin jointed linker element the underlying linear
   * Reissner beam element is used as a truss model only; that means, it only
   * reacts with an axial force - no moments and no transverse forces;
   * this shall be achieved by rotating the nodal triads according to the beam
   * axis in every new configuration given here from outside;
   * to keep the Reisner element shear-, bending- and torsion-free, we use the
   * same strategy to determine the nodal triads as for initialization of any
   * linker (see Init() ) */



  /* they are determined such that:
   * - the first base vector points in the direction of the distance vector
   *   of the two connection sites; (will be axis of connecting element)
   * - second and third base vector are arbitrarily constructed from cross-product
   *   of first base vector with either first or second base vector of global
   *   coordinate system; this avoids any singularities */
  LINALG::Matrix<3,3> linkeletriad(true);
  LINALG::Matrix<3,1> distvec(true);

  distvec.Update( 1.0, GetBindSpotPos2(), -1.0, GetBindSpotPos1() );

  // feasibility check regarding coinciding connection sites
  if (distvec.Norm2() < 1e-12)
  {
    std::cout << "\nBeamLinkRigidJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    GetBindSpotPos2().Print(std::cout);
    std::cout << "\ninitbspottriad1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspottriad2 =";
    GetBindSpotPos2().Print(std::cout);

    dserror("Initialization of BeamLinkRigidJointed failed because the two given binding "
        "spot positions are almost identical!");
  }

  // first base vector
  distvec.Scale( 1.0 / distvec.Norm2() );

  std::copy( distvec.A(), distvec.A() + 3, &linkeletriad( 0, 0 ) );


  // second base vector
  // check included angle of desired crosslinker axis (normalized distvec = first
  // base vector) and (1,0,0), i.e. scalar product which in this case simplifies to
  // first component of distvec
  LINALG::Matrix<3,1> unit_vector_global_x(true), unit_vector_global_y(true);
  unit_vector_global_x(0) = 1.0;
  unit_vector_global_y(1) = 1.0;

  const double scalarproduct = distvec(0);

  LINALG::Matrix<3,1> second_base_vecor_linkerele(true);

  // is included angle smaller than 45° ? then avoid singularity at angle=0° ...
  if ( std::abs(scalarproduct) > 0.5 * std::sqrt(2) )
  {
    second_base_vecor_linkerele.CrossProduct( distvec, unit_vector_global_y );
  }
  else
  {
    second_base_vecor_linkerele.CrossProduct( distvec, unit_vector_global_x );
  }

  // feasibility check
  if ( second_base_vecor_linkerele.Norm2() < 1e-12 )
  {
    std::cout << "\nBeamLinkRigidJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    GetBindSpotPos2().Print(std::cout);
    std::cout << "\ninitbspottriad1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspottriad2 =";
    GetBindSpotPos2().Print(std::cout);

    std::cout << "\ndistvec = ";
    distvec.Print(std::cout);
    std::cout << "\nsecond_base_vecor_linkerele = ";
    second_base_vecor_linkerele.Print(std::cout);

    dserror("Initialization of BeamLinkRigidJointed failed because the second base vector of the"
        "linker element's triad has almost length zero!");
  }
  else
  {
    second_base_vecor_linkerele.Scale( 1.0 / second_base_vecor_linkerele.Norm2() );
  }


  // third base vector to complete orthonormal triad
  LINALG::Matrix<3,1> third_base_vecor_linkerele(true);
  third_base_vecor_linkerele.CrossProduct( distvec, second_base_vecor_linkerele );

  // feasibility check
  if ( std::abs( third_base_vecor_linkerele.Norm2() - 1.0 ) > 1e-12 )
  {
    std::cout << "\nBeamLinkRigidJointed initialized with ...";
    std::cout << "\ninitbspotpos1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspotpos2 =";
    GetBindSpotPos2().Print(std::cout);
    std::cout << "\ninitbspottriad1 =";
    GetBindSpotPos1().Print(std::cout);
    std::cout << "\ninitbspottriad2 =";
    GetBindSpotPos2().Print(std::cout);

    std::cout << "\ndistvec = ";
    distvec.Print(std::cout);
    std::cout << "\nsecond_base_vecor_linkerele = ";
    second_base_vecor_linkerele.Print(std::cout);
    std::cout << "\nthird_base_vecor_linkerele = ";
    third_base_vecor_linkerele.Print(std::cout);

    dserror("Initialization of BeamLinkRigidJointed failed because the third base vector of the"
        "linker element's triad is no unit vector!");
  }

  /* store the initial triads as quaternions in class variables for the subsequent
   * use in setup of reference configuration of the connecting element */
  std::copy( second_base_vecor_linkerele.A(), second_base_vecor_linkerele.A()+3, &linkeletriad(0,1) );
  std::copy( third_base_vecor_linkerele.A(), third_base_vecor_linkerele.A()+3, &linkeletriad(0,2) );

  LARGEROTATIONS::triadtoquaternion( linkeletriad, triad_ );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::FillStateVariablesForElementEvaluation(
    LINALG::TMatrix<double,6,1>&               disp_totlag_centerline,
    std::vector<LINALG::TMatrix<double,4,1> >& Qnode
    ) const
{
  for ( unsigned int i = 0; i < 3; ++i )
  {
    disp_totlag_centerline(i) = GetBindSpotPos1()(i);
    disp_totlag_centerline(3+i) = GetBindSpotPos2()(i);
  }

  Qnode.push_back( triad_ );
  Qnode.push_back( triad_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed::GetBindingSpotForce(
    int bspotid,
    LINALG::SerialDenseVector & bspotforce ) const
{
  bspotforce = bspotforces_[bspotid];
}
