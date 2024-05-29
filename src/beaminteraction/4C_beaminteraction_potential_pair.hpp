/*-----------------------------------------------------------------------     */
/*! \file

\brief one generic (beam-to-?) element pair interacting via potentials

\level 3

*/
/*-----------------------------------------------------------------------     */
#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace CORE::Elements
{
  class Element;
}

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
    virtual ~BeamPotentialPair() = default;
    //! Initialization
    void Init(const Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> params_ptr,
        CORE::Elements::Element const* element1, CORE::Elements::Element const* element2);

    //! Setup
    virtual void Setup();

    //! return appropriate derived (templated) class (acts as a simple factory)
    static Teuchos::RCP<BeamPotentialPair> Create(
        std::vector<CORE::Elements::Element const*> const& ele_ptrs,
        BEAMINTERACTION::BeamPotentialParams const& beam_potential_params);

    //@}


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    virtual bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22,
        const std::vector<CORE::Conditions::Condition*> linechargeconds, const double k,
        const double m) = 0;

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
    inline CORE::Elements::Element const* Element1() const { return element1_; };

    /*!
    \brief Get second element
    */
    inline CORE::Elements::Element const* Element2() const { return element2_; };

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    virtual void get_all_interacting_point_coords_element1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const = 0;

    virtual void get_all_interacting_point_coords_element2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const = 0;

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    virtual void get_forces_at_all_interacting_points_element1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const = 0;

    virtual void get_forces_at_all_interacting_points_element2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const = 0;

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    virtual void get_moments_at_all_interacting_points_element1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const = 0;

    virtual void get_moments_at_all_interacting_points_element2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const = 0;

    /*!
    \brief Get interaction free energy / potential
    */
    virtual double get_energy() const = 0;

    /** \brief print this beam potential-based element pair to screen
     *
     *  \author grill */
    virtual void Print(std::ostream& out) const = 0;

    /** \brief print this beam potential element pair to screen
     *
     *  \author grill */
    virtual void print_summary_one_line_per_active_segment_pair(std::ostream& out) const = 0;

   protected:
    //! returns init state
    inline bool const& is_init() const { return isinit_; };

    //! returns setup state
    inline bool const& is_setup() const { return issetup_; };

    //! Check the init state
    void check_init() const;

    //! Check the init and setup state
    void check_init_setup() const;

    //! get Gauss rule to be used
    CORE::FE::GaussRule1D get_gauss_rule() const;

    /*!
    \brief Set first element
    */
    inline void set_element1(CORE::Elements::Element const* element1) { element1_ = element1; };

    /*!
    \brief Set second element
    */
    inline void set_element2(CORE::Elements::Element const* element2) { element2_ = element2; };

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
    CORE::Elements::Element const* element1_;

    //! second element of interacting pair
    CORE::Elements::Element const* element2_;
    //@}
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
