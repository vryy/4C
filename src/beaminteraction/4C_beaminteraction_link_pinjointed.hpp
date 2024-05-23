/*----------------------------------------------------------------------*/
/*! \file

\brief connecting beam linked by pin joint

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_LINK_PINJOINTED_HPP
#define FOUR_C_BEAMINTERACTION_LINK_PINJOINTED_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_link.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace CORE::COMM
{
  class PackBuffer;
}
namespace DRT
{
  class Element;
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
  class BeamLinkPinJointedType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "BeamLinkPinJointedType"; };

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

    //! Initialization
    void Init(int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& initpos,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& inittriad,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset) override;

    //! Setup
    void Setup(const int matnum) override;

    /*!
    \brief Return unique ParObject id

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override = 0;

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    //! get force in first or second binding spot
    void GetBindingSpotForce(
        int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const override
    {
      FOUR_C_THROW(" needs to be implemented in derived classes.");
    }

    // get current length of linker
    virtual double get_current_linker_length() const = 0;

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces
    */
    bool EvaluateForce(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2) override = 0;

    /*!
    \brief Evaluate stiffness contribution
    */
    bool EvaluateStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override = 0;

    /*!
    \brief Evaluate forces and stiffness contribution
    */
    bool EvaluateForceStiff(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override = 0;

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(std::vector<CORE::LINALG::Matrix<3, 1>>& bspotpos,
        std::vector<CORE::LINALG::Matrix<3, 3>>& bspottriad) override;

    //! return appropriate instance of the desired class (acts as a simple factory)
    static Teuchos::RCP<BeamLinkPinJointed> Create(INPAR::BEAMINTERACTION::JointType type);

    void Print(std::ostream& out) const;
    //@}

   private:
    //@}
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
