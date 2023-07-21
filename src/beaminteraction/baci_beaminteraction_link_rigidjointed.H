/*----------------------------------------------------------------------*/
/*! \file

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_LINK_RIGIDJOINTED_H
#define BACI_BEAMINTERACTION_LINK_RIGIDJOINTED_H

#include "baci_linalg_fixedsizematrix.H"
#include "baci_lib_parobject.H"
#include "baci_lib_parobjectfactory.H"
#include "baci_beaminteraction_link.H"


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

namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

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
        const std::vector<CORE::LINALG::Matrix<3, 1>>& initpos,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& inittriad,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset);

    //! Setup
    virtual void Setup(const int matnum);

    /*!
    \brief Return unique ParObject id

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
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
    inline const CORE::LINALG::Matrix<4, 1, double>& GetBindSpotQuaternion1() const
    {
      return bspottriad1_;
    }

    //! return orientation of first connection site as quaternion
    inline const CORE::LINALG::Matrix<4, 1, double>& GetBindSpotQuaternion2() const
    {
      return bspottriad2_;
    }

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const
    {
      dserror(" needs to be implemented in derived classes.");
    }


    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces
    */
    virtual bool EvaluateForce(
        CORE::LINALG::SerialDenseVector& forcevec1, CORE::LINALG::SerialDenseVector& forcevec2) = 0;

    /*!
    \brief Evaluate stiffness contribution
    */
    virtual bool EvaluateStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) = 0;

    /*!
    \brief Evaluate forces and stiffness contribution
    */
    virtual bool EvaluateForceStiff(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) = 0;

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(std::vector<CORE::LINALG::Matrix<3, 1>>& bspotpos,
        std::vector<CORE::LINALG::Matrix<3, 3>>& bspottriad);

    //! return appropriate instance of the desired class (acts as a simple factory)
    static Teuchos::RCP<BeamLinkRigidJointed> Create();
    void Print(std::ostream& out) const;
    //@}

   private:
    //! current triad of the two connection sites as quaternions
    CORE::LINALG::Matrix<4, 1> bspottriad1_;
    CORE::LINALG::Matrix<4, 1> bspottriad2_;

    //! triads representing the (constant) relative rotation between the two nodal triads of the
    //! linker element
    // and cross-section orientation of its connection sites (a.k.a. binding spots)
    CORE::LINALG::Matrix<3, 3> Lambdarel1_;
    CORE::LINALG::Matrix<3, 3> Lambdarel2_;

    //@}
  };

}  // namespace BEAMINTERACTION

#endif
