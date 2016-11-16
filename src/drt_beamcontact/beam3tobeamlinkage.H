/*----------------------------------------------------------------------*/
/*!
\file beam3tobeamlinkage.H

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BEAMTOBEAMLINKAGE_H
#define BEAMTOBEAMLINKAGE_H

#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
class Epetra_Vector;
class Epetra_SerialDenseMatrix;
namespace DRT {
  class Element;
namespace ELEMENTS {
  class Beam3Base;
  class Beam3r;
}
}

namespace LINALG {
  class SparseMatrix;
}

namespace BEAMINTERACTION
{

/*!
 \brief element for interaction of two 3D beam elements via a mechanical linkage
 */
class BeamToBeamLinkage
{
public:
  //! @name Constructors and destructors and related methods

  //! Constructor
  BeamToBeamLinkage();

  //! Destructor
  virtual ~BeamToBeamLinkage() {}

  //! Initialization
  virtual void Init(
//      DRT::ELEMENTS::Beam3Base*  element1,
//      DRT::ELEMENTS::Beam3Base*  element2,
//      const double&  xi1,
//      const double&  xi2,
      const LINALG::Matrix<3,1>& initpos1,
      const LINALG::Matrix<3,1>& initpos2,
      const LINALG::Matrix<3,3>& inittriad1,
      const LINALG::Matrix<3,3>& inittriad2);

  //! Setup
  virtual void Setup();

  //@}

  //! @name Access methods

  //! return position of first connection site
  inline const LINALG::TMatrix<double,3,1>& GetBindSpotPos1() const
    { return bspotpos1_; }

  //! return position of second connection site
  inline const LINALG::TMatrix<double,3,1>& GetBindSpotPos2() const
    { return bspotpos2_; }

  //! return orientation of first connection site as quaternion
  inline const LINALG::TMatrix<double,4,1>& GetBindSpotQuaternion1() const
    { return bspottriad1_; }

  //! return orientation of first connection site as quaternion
  inline const LINALG::TMatrix<double,4,1>& GetBindSpotQuaternion2() const
    { return bspottriad2_; }

  //@}

  //! @name Public evaluation methods

  /*!
  \brief Evaluate forces
  */
  virtual bool EvaluateForce(
      LINALG::TMatrix<double,6,1>& forcevec1,
      LINALG::TMatrix<double,6,1>& forcevec2) = 0;

  /*!
  \brief Evaluate stiffness contribution
  */
  virtual bool EvaluateStiff(
      LINALG::TMatrix<double,6,6>& stiffmat11,
      LINALG::TMatrix<double,6,6>& stiffmat12,
      LINALG::TMatrix<double,6,6>& stiffmat21,
      LINALG::TMatrix<double,6,6>& stiffmat22) = 0;

  /*!
  \brief Evaluate forces and stiffness contribution
  */
  virtual bool EvaluateForceStiff(
      LINALG::TMatrix<double,6,1>& forcevec1,
      LINALG::TMatrix<double,6,1>& forcevec2,
      LINALG::TMatrix<double,6,6>& stiffmat11,
      LINALG::TMatrix<double,6,6>& stiffmat12,
      LINALG::TMatrix<double,6,6>& stiffmat21,
      LINALG::TMatrix<double,6,6>& stiffmat22) = 0;

  /*
  \brief Update position and triad of both connection sites (a.k.a. binding spots)
  */
  void ResetState(
      LINALG::Matrix<3,1>& bspotpos1,
      LINALG::Matrix<3,1>& bspotpos2,
      LINALG::Matrix<3,3>& bspottriad1,
      LINALG::Matrix<3,3>& bspottriad2);

  //! return appropriate instance of the desired class (acts as a simple factory)
  static Teuchos::RCP<BeamToBeamLinkage> Create();

  //@}

protected:
  //! returns init state
  inline const bool& IsInit() const
  { return isinit_; };

  //! returns setup state
  inline const bool& IsSetup() const
  { return issetup_; };

  //! Check the init state
  inline void CheckInit() const
  {
    if (not IsInit())
      dserror("Call Init() first!");
  }

  //! Check the init and setup state
  inline void CheckInitSetup() const
  {
    if (not IsInit() or not IsSetup())
      dserror("Call Init() and Setup() first!");
  }

protected:
  //! @name member variables

  //! indicates if the Init() function has been called
  bool isinit_;

  //! indicates if the Setup() function has been called
  bool issetup_;

private:
//  //! first element of pair
//  Teuchos::RCP<const DRT::ELEMENTS::Beam3Base> element1_;
//
//  //! second element of pair
//  Teuchos::RCP<const DRT::ELEMENTS::Beam3Base> element2_;

  //! position of link (in parameter space [-1,1]) on first/second element
//  double xi1_;
//  double xi2_;

  //! current position of the two connection sites (a.k.a. binding spots)
  LINALG::TMatrix<double,3,1> bspotpos1_;
  LINALG::TMatrix<double,3,1> bspotpos2_;

  //! current triad of the two connection sites as quaternions
  LINALG::TMatrix<double,4,1> bspottriad1_;
  LINALG::TMatrix<double,4,1> bspottriad2_;

  //! triads representing the (constant) relative rotation between the two nodal triads of the linker element
  // and cross-section orientation of its connection sites (a.k.a. binding spots)
  LINALG::TMatrix<double,3,3> Lambdarel1_;
  LINALG::TMatrix<double,3,3> Lambdarel2_;

  //@}
};


/*!
 \brief element for link between two 3D beam elements via a linear (2 noded) Simo-Reissner beam element
 */
class Beam3rLin2Linkage: public BEAMINTERACTION::BeamToBeamLinkage
{
public:

  //! @name Friends
  //no friend classes defined
  //@}

  //! @name Constructors and destructors and related methods
  /*!
  \brief Standard Constructor
  */
  Beam3rLin2Linkage();

  /*!
  \brief Destructor
  */
  virtual ~Beam3rLin2Linkage() {};

  //! Setup [derived]
  virtual void Setup();
  //@}


  //! @name Access methods

  //@}

  //! @name Public evaluation methods

  /*!
  \brief Evaluate forces and stiffness contribution [derived]
  */
  virtual bool EvaluateForce(
      LINALG::TMatrix<double,6,1>& forcevec1,
      LINALG::TMatrix<double,6,1>& forcevec2);

  /*!
  \brief Evaluate stiffness contribution [derived]
  */
  virtual bool EvaluateStiff(
      LINALG::TMatrix<double,6,6>& stiffmat11,
      LINALG::TMatrix<double,6,6>& stiffmat12,
      LINALG::TMatrix<double,6,6>& stiffmat21,
      LINALG::TMatrix<double,6,6>& stiffmat22);

  /*!
  \brief Evaluate forces and stiffness contribution [derived]
  */
  virtual bool EvaluateForceStiff(
      LINALG::TMatrix<double,6,1>& forcevec1,
      LINALG::TMatrix<double,6,1>& forcevec2,
      LINALG::TMatrix<double,6,6>& stiffmat11,
      LINALG::TMatrix<double,6,6>& stiffmat12,
      LINALG::TMatrix<double,6,6>& stiffmat21,
      LINALG::TMatrix<double,6,6>& stiffmat22);


  //@}
private:

  //! @name member variables

  //! new connecting element
  Teuchos::RCP<DRT::ELEMENTS::Beam3r> linkele_;

  //@}

  //! @name Private evaluation methods

  //@}

};

} // namespace BEAMINTERACTION

#endif
