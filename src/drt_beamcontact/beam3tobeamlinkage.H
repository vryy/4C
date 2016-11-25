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

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"

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

class BeamToBeamLinkageType: public DRT::ParObjectType
{
public:

  std::string Name() const{return "BeamToBeamLinkageType";};

  static BeamToBeamLinkageType & Instance(){return instance_;};

private:

  static BeamToBeamLinkageType instance_;
};


/*!
 \brief element for interaction of two 3D beam elements via a mechanical linkage
 */
class BeamToBeamLinkage : public DRT::ParObject
{
public:
  //! @name Constructors and destructors and related methods

  //! Constructor
  BeamToBeamLinkage();

  //! Destructor
  virtual ~BeamToBeamLinkage() {};

  //! Initialization
  virtual void Init(
      const int id,
      const std::vector<std::pair<int, int> >& eleids,
      const std::vector<LINALG::Matrix<3,1> >& initpos,
      const std::vector<LINALG::Matrix<3,3> >& inittriad);

  //! Setup
  virtual void Setup();

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
  virtual inline int Id() const { return id_; }

  /*!
  \brief Return gid of element one
  */
  virtual inline const int& GetEleGid(const int& elenum) const {
    return eleids_[elenum].first; }

  /*!
  \brief Return gid of element one
  */
  virtual inline const int& GetLocBSpotNum(const int& elenum) const {
    return eleids_[elenum].second; }

  /*!
  \brief Return gid of element one
  */
  virtual inline const std::set<int>& GetInvolvedProcs() const {
    return involvedprocs_; }

  /*!
  \brief Return gid of element one
  */
  virtual inline void SetInvolvedProcs(std::set<int>& involvedprocs) {
    involvedprocs_.clear();
    involvedprocs_ = involvedprocs;
    return; }

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

  void Print(std::ostream& out) const;

protected:
  //! @name member variables

  //! indicates if the Init() function has been called
  bool isinit_;

  //! indicates if the Setup() function has been called
  bool issetup_;

private:

  //! a unique global id
  int id_;

  //! elegid and local binding spot number of first [0] and second [1] element
  std::vector<std::pair<int, int> > eleids_;

  //! involved procs
  std::set<int> involvedprocs_;

  //! current position of the two connection sites (a.k.a. binding spots)
  LINALG::Matrix<3,1> bspotpos1_;
  LINALG::Matrix<3,1> bspotpos2_;

  //! current triad of the two connection sites as quaternions
  LINALG::Matrix<4,1> bspottriad1_;
  LINALG::Matrix<4,1> bspottriad2_;

  //! triads representing the (constant) relative rotation between the two nodal triads of the linker element
  // and cross-section orientation of its connection sites (a.k.a. binding spots)
  LINALG::Matrix<3,3> Lambdarel1_;
  LINALG::Matrix<3,3> Lambdarel2_;

  //@}
};

} // namespace BEAMINTERACTION

#endif
