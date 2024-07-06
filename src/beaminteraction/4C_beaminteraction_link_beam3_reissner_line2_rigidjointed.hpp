/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a linear Reissner beam element used as mechanical link between two other beam
elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_HPP
#define FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_link_rigidjointed.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SerialDenseMatrix;
}

namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3r;
  }
}  // namespace Discret

namespace BEAMINTERACTION
{
  class BeamLinkBeam3rLine2RigidJointedType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "BeamLinkBeam3rLine2RigidJointedType"; };

    static BeamLinkBeam3rLine2RigidJointedType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    static BeamLinkBeam3rLine2RigidJointedType instance_;
  };


  /*!
   \brief element for link between two 3D beam elements via a linear (2 noded) Simo-Reissner beam
   element
   */
  class BeamLinkBeam3rLine2RigidJointed : public BeamLinkRigidJointed
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamLinkBeam3rLine2RigidJointed();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkBeam3rLine2RigidJointed(const BeamLinkBeam3rLine2RigidJointed& old);



    //! Setup [derived]
    void setup(int matnum) override;

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int unique_par_object_id() const override
    {
      return BeamLinkBeam3rLine2RigidJointedType::instance().unique_par_object_id();
    };

    /*!
    \brief Pack this class so it can be communicated [derived]

    \ref pack and \ref unpack are used to communicate this element

    */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class [derived]

    \ref pack and \ref unpack are used to communicate this element

    */
    void unpack(const std::vector<char>& data) override;

    /// return copy of this linking object
    Teuchos::RCP<BeamLink> clone() const override;

    //@}


    //! @name Access methods

    //! get internal linker energy
    double get_internal_energy() const override;

    //! get kinetic linker energy
    double get_kinetic_energy() const override;

    //! scale linker element reference length
    void scale_linker_reference_length(double scalefac) override
    {
      FOUR_C_THROW(" not yet implemented for beam3r element.");
    }

    //! get force in first or second binding spot
    void get_binding_spot_force(
        int bspotid, Core::LinAlg::SerialDenseVector& bspotforce) const override
    {
      bspotforce = bspotforces_[bspotid];
    }

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2) override;

    /*!
    \brief Evaluate stiffness contribution [derived]
    */
    bool evaluate_stiff(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override;

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force_stiff(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2, Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override;

    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Fill absolute nodal positions and nodal quaternions with current values
    */
    void fill_state_variables_for_element_evaluation(
        Core::LinAlg::Matrix<6, 1, double>& disp_totlag_centerline,
        std::vector<Core::LinAlg::Matrix<4, 1, double>>& Qnode) const;

    //@}

   private:
    //! @name member variables

    //! new connecting element
    Teuchos::RCP<Discret::ELEMENTS::Beam3r> linkele_;


    //! the following variables are for output purposes only (no need to pack or unpack)
    std::vector<Core::LinAlg::SerialDenseVector> bspotforces_;

    //@}
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
