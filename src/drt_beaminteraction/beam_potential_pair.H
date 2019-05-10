/*-----------------------------------------------------------------------     */
/*!
\file beam_potential_pair.H

\brief one generic (beam-to-?) element pair interacting via potentials

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------     */
#ifndef BEAMPOTENTIALPAIR_H
#define BEAMPOTENTIALPAIR_H

#include "../drt_fem_general/drt_utils_integration.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>

// forward declaration ...
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

namespace DRT
{
  class Condition;
  class Element;
}  // namespace DRT

namespace BEAMINTERACTION
{
  // forward declaration ...
  class BeamPotentialParams;


  /*!
   \brief
   */
  class BeamPotentialPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Constructor
    */
    BeamPotentialPair();

    /*!
    \brief Destructor
    */
    virtual ~BeamPotentialPair() {}

    //! Initialization
    void Init(const Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> params_ptr,
        DRT::Element const* element1, DRT::Element const* element2);

    //! Setup
    virtual void Setup();

    //! return appropriate derived (templated) class (acts as a simple factory)
    static Teuchos::RCP<BeamPotentialPair> Create(std::vector<DRT::Element const*> const& ele_ptrs,
        BEAMINTERACTION::BeamPotentialParams const& beam_potential_params);

    //@}


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    virtual bool Evaluate(LINALG::SerialDenseVector* forcevec1,
        LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
        LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
        LINALG::SerialDenseMatrix* stiffmat22, const std::vector<DRT::Condition*> linechargeconds,
        const double k, const double m) = 0;

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    virtual void ResetState(double time, std::vector<double> const& centerline_dofvec_ele1,
        std::vector<double> const& centerline_dofvec_ele2) = 0;

    //@}

    //! @name Access methods

    inline Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> Params() const
    {
      return beam_potential_params_;
    }

    /*!
    \brief Get first element
    */
    inline DRT::Element const* Element1() const { return element1_; };

    /*!
    \brief Get second element
    */
    inline DRT::Element const* Element2() const { return element2_; };

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    virtual void GetAllInteractingPointCoordsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const = 0;

    virtual void GetAllInteractingPointCoordsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const = 0;

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    virtual void GetForcesAtAllInteractingPointsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& forces) const = 0;

    virtual void GetForcesAtAllInteractingPointsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& forces) const = 0;

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    virtual void GetMomentsAtAllInteractingPointsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& moments) const = 0;

    virtual void GetMomentsAtAllInteractingPointsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& moments) const = 0;

    /*!
    \brief Get interaction free energy / potential
    */
    virtual double GetEnergy() const = 0;

    /** \brief print this beam potential-based element pair to screen
     *
     *  \author grill */
    virtual void Print(std::ostream& out) const = 0;

    /** \brief print this beam potential element pair to screen
     *
     *  \author grill */
    virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const = 0;

   protected:
    //! returns init state
    inline bool const& IsInit() const { return isinit_; };

    //! returns setup state
    inline bool const& IsSetup() const { return issetup_; };

    //! Check the init state
    void CheckInit() const;

    //! Check the init and setup state
    void CheckInitSetup() const;

    //! get Gauss rule to be used
    DRT::UTILS::GaussRule1D GetGaussRule() const;

    //@}

   protected:
    //! @name member variables

    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

   private:
    //! beam potential parameter data container
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> beam_potential_params_;

    //! first element of interacting pair
    DRT::Element const* element1_;

    //! second element of interacting pair
    DRT::Element const* element2_;
    //@}
  };
}  // namespace BEAMINTERACTION

#endif
