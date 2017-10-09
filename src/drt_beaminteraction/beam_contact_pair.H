/*----------------------------------------------------------------------*/
/*!
\file beam_contact_pair.H

\brief one generic (beam-to-?) contact element pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BEAMCONTACTPAIR_H
#define BEAMCONTACTPAIR_H

#include <Teuchos_RCP.hpp>

//forward declaration ...
namespace LINALG {
  class SerialDenseVector;
  class SerialDenseMatrix;
}
namespace DRT {
class Element;
}

namespace BEAMINTERACTION
{

//forward declaration ...
class BeamContactParams;


/*!
 \brief
 */
class BeamContactPair
{
public:

  //! @name Friends
  //no friend classes defined
  //@}

  //! @name Constructors and destructors and related methods

  BeamContactPair();

  /*!
  \brief Destructor
  */
  virtual ~BeamContactPair() {}

  //! Initialization
  void Init(
      const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
      std::vector< DRT::Element const*> elements);

  //! Setup
  virtual void Setup();

  //@}


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair, return value indicates whether pair is active,
         i.e. non-zero values for force and stiffmat are returned
  */
  virtual bool Evaluate(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22)=0;

  //! return appropriate internal implementation class (acts as a simple factory)
  static Teuchos::RCP<BeamContactPair> Create(
      std::vector< DRT::Element const *> const& ele_ptrs);

  /*
  \brief Update state of translational nodal DoFs (absolute positions and tangents) of both elements
  */
  virtual void ResetState(
      const std::vector<double>& centerline_dofvec_ele1,
      const std::vector<double>& centerline_dofvec_ele2)=0;

  //@}

  //! @name Access methods

  inline Teuchos::RCP<BEAMINTERACTION::BeamContactParams> Params() const
  {
    return params_;
  }

  /*!
  \brief Get first element
  */
  inline const DRT::Element* Element1() const { return element1_;};

  /*!
  \brief Get second element
  */
  inline const DRT::Element* Element2() const { return element2_;};

//  /*!
//  \brief Get number of contact points on this element pair
//  */
//  virtual int GetNumCps()=0;
//  virtual int GetNumGps()=0;
//  virtual int GetNumEps()=0;
//
//  /*!
//  \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end point contact) of all contact pairs
//  */
//  virtual std::vector<int> GetContactType()=0;
//
//  /*!
//  \brief Get vector of all gaps of this contact pair
//  */
//  virtual std::vector<double> GetGap()=0;
//
//  /*!
//  \brief Get vector of all contact forces of this contact pair
//  */
//  virtual std::vector<double> GetContactForce()=0;
//
//  /*!
//  \brief Get vector of all contact angles of this contact pair
//  */
//  virtual std::vector<double> GetContactAngle()=0;
//
//  /*!
//  \brief Get vector of all closest points of this contact pair
//  */
//  virtual std::vector<std::pair<double,double> > GetClosestPoint()=0;
//
//  /*!
//  \brief Return number of individual contact segments on element pair
//  */
//  virtual std::pair<int,int> GetNumSegments()=0;
//
//  /*!
//  \brief Return ids of active segments
//  */
//  virtual std::vector<std::pair<int,int> > GetSegmentIds()=0;

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  virtual bool GetContactFlag() const =0;

//  /*!
//  \brief Get coordinates of contact point of element1 and element2
//  */
//  virtual const std::vector<LINALG::Matrix<3,1> > GetX1()=0;
//
//  virtual const std::vector<LINALG::Matrix<3,1> > GetX2()=0;
//
//  virtual const Epetra_SerialDenseVector GetNormal()=0;
//
//  virtual LINALG::TMatrix<TYPE, 3, 1>* GetNormalOld()=0;
//
//  /*!
//    \Check, if there is a difference between the result of the new and old gap definition, i.e. if the beams centerlines have already crossed or not.
//  */
//  virtual bool GetNewGapStatus()=0;
//
//  /*!
//    \Get energy of penalty contact.
//  */
//  virtual double GetEnergy()=0;
//
//  /*!
//    \Get energy of perp penalty contact without transition factor contribution.
//  */
//  virtual double GetUnscaledPerpEnergy()=0;
//
//  /*!
//    \Get energy of parallel penalty contact without transition factor contribution.
//  */
//  virtual double GetUnscaledParallelEnergy()=0;
//
//  virtual bool FirstTimeStep()=0;
//
//  /*!
//  \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
//  */
//  virtual bool GetShiftStatus()=0;
//  //@}

  /** \brief print this beam contact element pair to screen
   *
   *  \author grill
   *  \date 05/16 */
  virtual void Print(std::ostream& out) const=0;

  /** \brief print this beam contact element pair to screen
   *
   *  \author grill
   *  \date 05/16 */
  virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const=0;

protected:
  //! returns init state
  inline const bool& IsInit() const
  { return isinit_; };

  //! returns setup state
  inline const bool& IsSetup() const
  { return issetup_; };

  //! Check the init state
  void CheckInit() const;

  //! Check the init and setup state
  void CheckInitSetup() const;

protected:
  //! @name member variables

  //! indicates if the Init() function has been called
  bool isinit_;

  //! indicates if the Setup() function has been called
  bool issetup_;

private:

  //! beam contact parameter data container
  Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_;

  //! first element of interacting pair
  const DRT::Element* element1_;

  //! second element of interacting pair
  const DRT::Element* element2_;
  //@}

};
}

#endif
