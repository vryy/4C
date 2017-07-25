/*----------------------------------------------------------------------*/
/*!
\file beam_link_beam3r_lin2_rigidjointed.H

\brief Wrapper for a linear Reissner beam element used as mechanical link between two other beam elements

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BeamLinkBeam3rLin2RigidJointed_H
#define BeamLinkBeam3rLin2RigidJointed_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "beam_link_rigidjointed.H"

// forward declarations
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}

namespace DRT {
namespace ELEMENTS {
  class Beam3r;
}
}

namespace BEAMINTERACTION
{

class BeamLinkBeam3rLin2RigidJointedType: public DRT::ParObjectType
{
public:

  std::string Name() const {return "BeamLinkBeam3rLin2RigidJointedType";};

  static BeamLinkBeam3rLin2RigidJointedType & Instance(){return instance_;};

  virtual DRT::ParObject* Create(const std::vector<char> & data);

private:

  static BeamLinkBeam3rLin2RigidJointedType instance_;
};


/*!
 \brief element for link between two 3D beam elements via a linear (2 noded) Simo-Reissner beam element
 */
class BeamLinkBeam3rLin2RigidJointed : public BeamLinkRigidJointed
{
public:
  //! @name Friends
  //no friend classes defined
  //@}

  //! @name Constructors and destructors and related methods
  /*!
  \brief Standard Constructor
  */
  BeamLinkBeam3rLin2RigidJointed();

  /*!
  \brief Destructor
  */
  virtual ~BeamLinkBeam3rLin2RigidJointed() {};

  //! Setup [derived]
  virtual void Setup( const int matnum );

  /*!
  \brief Return unique ParObject id [derived]

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return BeamLinkBeam3rLin2RigidJointedType::Instance().UniqueParObjectId(); };

  /*!
  \brief Pack this class so it can be communicated [derived]

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class [derived]

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  //@}


  //! @name Access methods

  //@}

  //! @name Public evaluation methods

  /*!
  \brief Evaluate forces and stiffness contribution [derived]
  */
  virtual bool EvaluateForce(
      LINALG::SerialDenseVector& forcevec1,
      LINALG::SerialDenseVector& forcevec2);

  /*!
  \brief Evaluate stiffness contribution [derived]
  */
  virtual bool EvaluateStiff(
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22);

  /*!
  \brief Evaluate forces and stiffness contribution [derived]
  */
  virtual bool EvaluateForceStiff(
      LINALG::SerialDenseVector& forcevec1,
      LINALG::SerialDenseVector& forcevec2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22);

  //@}

private:
  //! @name Private evaluation methods

  /*!
  \brief Fill absolute nodal positions and nodal quaternions with current values
  */
  void FillStateVariablesForElementEvaluation(
      LINALG::TMatrix<double,6,1>&               disp_totlag_centerline,
      std::vector<LINALG::TMatrix<double,4,1> >& Qnode
      ) const;

  //@}

private:
  //! @name member variables

  //! new connecting element
  Teuchos::RCP<DRT::ELEMENTS::Beam3r> linkele_;

  //@}

};

} // namespace BEAMINTERACTION

#endif
