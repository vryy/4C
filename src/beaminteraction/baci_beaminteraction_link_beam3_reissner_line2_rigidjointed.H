/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a linear Reissner beam element used as mechanical link between two other beam
elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_H
#define BACI_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_H

#include "baci_linalg_fixedsizematrix.H"
#include "baci_beaminteraction_link_rigidjointed.H"
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
    std::string Name() const { return "BeamLinkBeam3rLine2RigidJointedType"; };

    static BeamLinkBeam3rLine2RigidJointedType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
    virtual ~BeamLinkBeam3rLine2RigidJointed(){};

    //! Setup [derived]
    virtual void Setup(int matnum);

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    virtual int UniqueParObjectId() const
    {
      return BeamLinkBeam3rLine2RigidJointedType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Pack this class so it can be communicated [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
    \brief Unpack data from a char vector into this class [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Unpack(const std::vector<char>& data);

    /// return copy of this linking object
    virtual Teuchos::RCP<BeamLink> Clone() const;

    //@}


    //! @name Access methods

    //! get internal linker energy
    virtual double GetInternalEnergy() const override;

    //! get kinetic linker energy
    virtual double GetKineticEnergy() const override;

    //! scale linker element reference length
    virtual void ScaleLinkerReferenceLength(double scalefac)
    {
      dserror(" not yet implemented for beam3r element.");
    }

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const
    {
      bspotforce = bspotforces_[bspotid];
    }

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    virtual bool EvaluateForce(
        CORE::LINALG::SerialDenseVector& forcevec1, CORE::LINALG::SerialDenseVector& forcevec2);

    /*!
    \brief Evaluate stiffness contribution [derived]
    */
    virtual bool EvaluateStiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22);

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    virtual bool EvaluateForceStiff(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22);

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
