/*----------------------------------------------------------------------*/
/*!

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BEAMLINKRIGIDJOINTED_H
#define BEAMLINKRIGIDJOINTED_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "beam_link.H"

//#include <Teuchos_RCP.hpp>

// forward declarations
namespace DRT
{
  class Element;
  class PackBuffer;
  namespace ELEMENTS
  {
    class Beam3Base;
    class Beam3r;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

namespace BEAMINTERACTION
{
  class BeamLinkRigidJointedType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "BeamLinkRigidJointedType"; };

    static BeamLinkRigidJointedType& Instance() { return instance_; };

   private:
    static BeamLinkRigidJointedType instance_;
  };


  /*!
   \brief element for interaction of two 3D beam elements via a mechanical linkage
   */
  class BeamLinkRigidJointed : public BeamLink
  {
   public:
    //! @name Constructors and destructors and related methods

    //! Constructor
    BeamLinkRigidJointed();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkRigidJointed(const BeamLinkRigidJointed& old);

    //! Destructor
    virtual ~BeamLinkRigidJointed() {}

    //! Initialization
    virtual void Init(const int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<LINALG::Matrix<3, 1>>& initpos,
        const std::vector<LINALG::Matrix<3, 3>>& inittriad,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset);

    //! Setup
    virtual void Setup(const int matnum);

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

    //! return orientation of first connection site as quaternion
    inline const LINALG::TMatrix<double, 4, 1>& GetBindSpotQuaternion1() const
    {
      return bspottriad1_;
    }

    //! return orientation of first connection site as quaternion
    inline const LINALG::TMatrix<double, 4, 1>& GetBindSpotQuaternion2() const
    {
      return bspottriad2_;
    }

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(int bspotid, LINALG::SerialDenseVector& bspotforce) const
    {
      dserror(" needs to be implemented in derived classes.");
    }


    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces
    */
    virtual bool EvaluateForce(
        LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseVector& forcevec2) = 0;

    /*!
    \brief Evaluate stiffness contribution
    */
    virtual bool EvaluateStiff(LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22) = 0;

    /*!
    \brief Evaluate forces and stiffness contribution
    */
    virtual bool EvaluateForceStiff(LINALG::SerialDenseVector& forcevec1,
        LINALG::SerialDenseVector& forcevec2, LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22) = 0;

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(
        std::vector<LINALG::Matrix<3, 1>>& bspotpos, std::vector<LINALG::Matrix<3, 3>>& bspottriad);

    //! return appropriate instance of the desired class (acts as a simple factory)
    static Teuchos::RCP<BeamLinkRigidJointed> Create();
    void Print(std::ostream& out) const;
    //@}

   private:
    //! current triad of the two connection sites as quaternions
    LINALG::Matrix<4, 1> bspottriad1_;
    LINALG::Matrix<4, 1> bspottriad2_;

    //! triads representing the (constant) relative rotation between the two nodal triads of the
    //! linker element
    // and cross-section orientation of its connection sites (a.k.a. binding spots)
    LINALG::Matrix<3, 3> Lambdarel1_;
    LINALG::Matrix<3, 3> Lambdarel2_;

    //@}
  };

}  // namespace BEAMINTERACTION

#endif
