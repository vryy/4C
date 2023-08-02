/*----------------------------------------------------------------------*/
/*! \file

\brief connecting beam linked by pin joint

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_LINK_PINJOINTED_H
#define BACI_BEAMINTERACTION_LINK_PINJOINTED_H

#include "baci_beaminteraction_link.H"
#include "baci_inpar_beaminteraction.H"
#include "baci_lib_parobject.H"
#include "baci_lib_parobjectfactory.H"
#include "baci_linalg_fixedsizematrix.H"


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
  class BeamLinkPinJointedType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "BeamLinkPinJointedType"; };

    static BeamLinkPinJointedType& Instance() { return instance_; };

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

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkPinJointed(const BeamLinkPinJointed& old);

    //! Destructor
    virtual ~BeamLinkPinJointed() {}

    //! Initialization
    virtual void Init(int id, const std::vector<std::pair<int, int>>& eleids,
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

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const
    {
      dserror(" needs to be implemented in derived classes.");
    }

    // get current length of linker
    virtual double GetCurrentLinkerLength() const = 0;

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
    static Teuchos::RCP<BeamLinkPinJointed> Create(INPAR::BEAMINTERACTION::JointType type);

    void Print(std::ostream& out) const;
    //@}

   private:
    //@}
  };

}  // namespace BEAMINTERACTION

#endif
