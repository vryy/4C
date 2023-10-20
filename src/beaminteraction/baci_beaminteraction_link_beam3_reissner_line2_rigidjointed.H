/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a linear Reissner beam element used as mechanical link between two other beam
elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_H
#define BACI_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_H

#include "baci_beaminteraction_link_rigidjointed.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensevector.H"

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseMatrix;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class Beam3r;
  }
}  // namespace DRT

namespace BEAMINTERACTION
{
  class BeamLinkBeam3rLine2RigidJointedType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "BeamLinkBeam3rLine2RigidJointedType"; };

    static BeamLinkBeam3rLine2RigidJointedType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

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

    /*!
    \brief Destructor
    */
    ~BeamLinkBeam3rLine2RigidJointed() override{};

    //! Setup [derived]
    void Setup(int matnum) override;

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override
    {
      return BeamLinkBeam3rLine2RigidJointedType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Pack this class so it can be communicated [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    /// return copy of this linking object
    Teuchos::RCP<BeamLink> Clone() const override;

    //@}


    //! @name Access methods

    //! get internal linker energy
    double GetInternalEnergy() const override;

    //! get kinetic linker energy
    double GetKineticEnergy() const override;

    //! scale linker element reference length
    void ScaleLinkerReferenceLength(double scalefac) override
    {
      dserror(" not yet implemented for beam3r element.");
    }

    //! get force in first or second binding spot
    void GetBindingSpotForce(
        int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const override
    {
      bspotforce = bspotforces_[bspotid];
    }

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool EvaluateForce(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2) override;

    /*!
    \brief Evaluate stiffness contribution [derived]
    */
    bool EvaluateStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override;

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool EvaluateForceStiff(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override;

    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Fill absolute nodal positions and nodal quaternions with current values
    */
    void FillStateVariablesForElementEvaluation(
        CORE::LINALG::Matrix<6, 1, double>& disp_totlag_centerline,
        std::vector<CORE::LINALG::Matrix<4, 1, double>>& Qnode) const;

    //@}

   private:
    //! @name member variables

    //! new connecting element
    Teuchos::RCP<DRT::ELEMENTS::Beam3r> linkele_;


    //! the following variables are for output purposes only (no need to pack or unpack)
    std::vector<CORE::LINALG::SerialDenseVector> bspotforces_;

    //@}
  };

}  // namespace BEAMINTERACTION

#endif
