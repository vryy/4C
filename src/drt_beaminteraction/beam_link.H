/*----------------------------------------------------------------------*/
/*!
\file beam_link.H

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BeamLink_H
#define BeamLink_H

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

class BeamLinkType: public DRT::ParObjectType
{
public:

  std::string Name() const{return "BeamLinkType";};

  static BeamLinkType & Instance(){return instance_;};

private:

  static BeamLinkType instance_;
};


/*!
 \brief element for interaction of two 3D beam elements via a mechanical linkage
 */
class BeamLink : public DRT::ParObject
{
public:
  //! @name Constructors and destructors and related methods

  //! Constructor
  BeamLink();

  //! Destructor
  virtual ~BeamLink() {}

  //! Initialization
  virtual void Init(
      const int id,
      const std::vector<std::pair<int, int> >& eleids,
      const std::vector<LINALG::Matrix<3,1> >& initpos);

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
  virtual inline int Id() const { return id_; }

  /*!
  \brief Return gid of first/second element (specified via given local element number 0/1)
  */
  virtual inline const int& GetEleGid(const int& elenum) const {
    return bspotIds_[elenum].first; }

  /*!
  \brief Return element-local binding spot number of first/second element
         (specified via given local element number 0/1)
  */
  virtual inline const int& GetLocBSpotNum(const int& elenum) const {
    return bspotIds_[elenum].second; }

  //! return position of first connection site
  inline const LINALG::TMatrix<double,3,1>& GetBindSpotPos1() const
    { return bspotpos1_; }

  //! return position of second connection site
  inline const LINALG::TMatrix<double,3,1>& GetBindSpotPos2() const
    { return bspotpos2_; }

  //@}

  //! @name Public evaluation methods

  /*
  \brief Update position and triad of both connection sites (a.k.a. binding spots)
  */
  void ResetState(
      std::vector<LINALG::Matrix<3,1> >& bspotpos);

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

  //! a unique global id
  int id_;

  //! unique identifiers for first [0] and second [1] binding spot:
  // each is a pair of element GID and local binding spot number
  std::vector<std::pair<int, int> > bspotIds_;

  //! current position of the two connection sites (a.k.a. binding spots)
  LINALG::Matrix<3,1> bspotpos1_;
  LINALG::Matrix<3,1> bspotpos2_;

  //@}
};

}

#endif
