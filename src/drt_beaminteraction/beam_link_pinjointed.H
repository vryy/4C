/*----------------------------------------------------------------------*/
/*!
\file beam_link_pinjointed.H

\brief connecting beam linked by pin joint

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------*/
#ifndef BEAMLINKPINJOINTED_H
#define BEAMLINKPINJOINTED_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "beam_link.H"

//#include <Teuchos_RCP.hpp>

// forward declarations
//class Epetra_Vector;
//class Epetra_SerialDenseMatrix;
namespace DRT {
  class Element;
  class PackBuffer;
namespace ELEMENTS {
  class Beam3Base;
  class Beam3r;
}
}
namespace LINALG {
  class SerialDenseVector;
  class SerialDenseMatrix;
}

namespace BEAMINTERACTION
{

class BeamLinkPinJointedType: public DRT::ParObjectType
{
public:

  std::string Name() const{return "BeamLinkPinJointedType";};

  static BeamLinkPinJointedType & Instance(){return instance_;};

private:

  static BeamLinkPinJointedType instance_;
};


/*!
 \brief element for interaction of two 3D beam elements via a mechanical linkage
 */
class BeamLinkPinJointed : public BeamLink
{
public:
  //! @name Constructors and destructors and related methods

  //! Constructor
  BeamLinkPinJointed();

  //! Destructor
  virtual ~BeamLinkPinJointed() {}

  //! Initialization
  virtual void Init(
      int id,
      const std::vector<std::pair<int, int> >& eleids,
      const std::vector<LINALG::Matrix<3,1> >& initpos,
      const std::vector<LINALG::Matrix<3,3> >& inittriad);

  //! Setup
  virtual void Setup( const int matnum);

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const = 0;

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Return global id
  */

  //@}

  //! @name Public evaluation methods

  /*!
  \brief Evaluate forces
  */
  virtual bool EvaluateForce(
      LINALG::SerialDenseVector& forcevec1,
      LINALG::SerialDenseVector& forcevec2) = 0;

  /*!
  \brief Evaluate stiffness contribution
  */
  virtual bool EvaluateStiff(
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) = 0;

  /*!
  \brief Evaluate forces and stiffness contribution
  */
  virtual bool EvaluateForceStiff(
      LINALG::SerialDenseVector& forcevec1,
      LINALG::SerialDenseVector& forcevec2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) = 0;

  /*
  \brief Update position and triad of both connection sites (a.k.a. binding spots)
  */
  void ResetState(
      std::vector<LINALG::Matrix<3,1> >& bspotpos,
      std::vector<LINALG::Matrix<3,3> >& bspottriad);

  //! return appropriate instance of the desired class (acts as a simple factory)
  static Teuchos::RCP<BeamLinkPinJointed> Create(
      INPAR::BEAMINTERACTION::JointType type);

  void Print(std::ostream& out) const;
  //@}

private:


  //@}
};

}

#endif
