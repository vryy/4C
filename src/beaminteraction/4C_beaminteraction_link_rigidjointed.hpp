/*----------------------------------------------------------------------*/
/*! \file

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_LINK_RIGIDJOINTED_HPP
#define FOUR_C_BEAMINTERACTION_LINK_RIGIDJOINTED_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_link.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Core::Communication
{
  class PackBuffer;
}
namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3Base;
    class Beam3r;
  }  // namespace ELEMENTS
}  // namespace Discret

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace BEAMINTERACTION
{
  class BeamLinkRigidJointedType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "BeamLinkRigidJointedType"; };

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

    //! Initialization
    void init(const int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<Core::LinAlg::Matrix<3, 1>>& initpos,
        const std::vector<Core::LinAlg::Matrix<3, 3>>& inittriad,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset) override;

    //! Setup
    void setup(const int matnum) override;

    /*!
    \brief Return unique ParObject id

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override = 0;

    /*!
    \brief Pack this class so it can be communicated

    \ref pack and \ref unpack are used to communicate this element

    */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref pack and \ref unpack are used to communicate this element

    */
    void unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    //! return orientation of first connection site as quaternion
    inline const Core::LinAlg::Matrix<4, 1, double>& get_bind_spot_quaternion1() const
    {
      return bspottriad1_;
    }

    //! return orientation of first connection site as quaternion
    inline const Core::LinAlg::Matrix<4, 1, double>& get_bind_spot_quaternion2() const
    {
      return bspottriad2_;
    }

    //! get force in first or second binding spot
    void GetBindingSpotForce(
        int bspotid, Core::LinAlg::SerialDenseVector& bspotforce) const override
    {
      FOUR_C_THROW(" needs to be implemented in derived classes.");
    }


    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces
    */
    bool evaluate_force(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2) override = 0;

    /*!
    \brief Evaluate stiffness contribution
    */
    bool evaluate_stiff(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override = 0;

    /*!
    \brief Evaluate forces and stiffness contribution
    */
    bool evaluate_force_stiff(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2, Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override = 0;

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(std::vector<Core::LinAlg::Matrix<3, 1>>& bspotpos,
        std::vector<Core::LinAlg::Matrix<3, 3>>& bspottriad) override;

    //! return appropriate instance of the desired class (acts as a simple factory)
    static Teuchos::RCP<BeamLinkRigidJointed> Create();
    void Print(std::ostream& out) const;
    //@}

   private:
    //! current triad of the two connection sites as quaternions
    Core::LinAlg::Matrix<4, 1> bspottriad1_;
    Core::LinAlg::Matrix<4, 1> bspottriad2_;

    //! triads representing the (constant) relative rotation between the two nodal triads of the
    //! linker element
    // and cross-section orientation of its connection sites (a.k.a. binding spots)
    Core::LinAlg::Matrix<3, 3> lambdarel1_;
    Core::LinAlg::Matrix<3, 3> lambdarel2_;

    //@}
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
