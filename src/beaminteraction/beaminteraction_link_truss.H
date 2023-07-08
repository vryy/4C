/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a truss element used as mechanical link
       between two beam elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BEAMINTERACTION_LINK_TRUSS_H
#define BEAMINTERACTION_LINK_TRUSS_H

#include "linalg_fixedsizematrix.H"
#include "beaminteraction_link_pinjointed.H"
#include "linalg_serialdensevector.H"

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  namespace ELEMENTS
  {
    class Truss3;
  }
}  // namespace DRT

namespace BEAMINTERACTION
{
  class BeamLinkTrussType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "BeamLinkTrussType"; };

    static BeamLinkTrussType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static BeamLinkTrussType instance_;
  };


  /*!
   \brief element for link between two 3D beam elements via a truss element
   */
  class BeamLinkTruss : public BeamLinkPinJointed
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamLinkTruss();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkTruss(const BeamLinkTruss& old);

    /*!
    \brief Destructor
    */
    virtual ~BeamLinkTruss(){};

    //! Initialization [derived]
    virtual void Init(int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& initpos,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& inittriad,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset);

    //! Setup [derived]
    virtual void Setup(const int matnum);

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    virtual int UniqueParObjectId() const
    {
      return BeamLinkTrussType::Instance().UniqueParObjectId();
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
    virtual void ScaleLinkerReferenceLength(double scalefac);

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(
        int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const;

    virtual double GetCurrentLinkerLength() const;

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

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(std::vector<CORE::LINALG::Matrix<3, 1>>& bspotpos,
        std::vector<CORE::LINALG::Matrix<3, 3>>& bspottriad);


    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Fill absolute nodal positions and nodal quaternions with current values
    */
    void FillStateVariablesForElementEvaluation(
        CORE::LINALG::Matrix<6, 1, double>& absolute_nodal_positions) const;

    void GetDispForElementEvaluation(std::map<std::string, std::vector<double>>& ele_state) const;

    //@}

   private:
    //! @name member variables

    //! new connecting element
    Teuchos::RCP<DRT::ELEMENTS::Truss3> linkele_;

    //! the following variables are for output purposes only (no need to pack or unpack)
    std::vector<CORE::LINALG::SerialDenseVector> bspotforces_;

    //@}
  };

}  // namespace BEAMINTERACTION

#endif
