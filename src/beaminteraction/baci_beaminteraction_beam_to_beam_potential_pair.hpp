/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief One beam-to-beam potential-based interacting pair (two beam elements)

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_POTENTIAL_PAIR_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_POTENTIAL_PAIR_HPP

#include "baci_config.hpp"

#include "baci_beaminteraction_potential_pair.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace DRT

namespace BEAMINTERACTION
{
  /*!
   \brief class for potential-based interaction between two 3D beam elements
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  class BeamToBeamPotentialPair : public BeamPotentialPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToBeamPotentialPair();



    //! Setup
    void Setup() override;
    //@}

    //! @name Derived methods from base class
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22,
        const std::vector<DRT::Condition*> linechargeconds, const double k,
        const double m) override;

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    void ResetState(double time, const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) override;

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    void GetAllInteractingPointCoordsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override
    {
      coords = centerline_coords_gp_1_;
    }

    void GetAllInteractingPointCoordsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override
    {
      coords = centerline_coords_gp_2_;
    }

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    void GetForcesAtAllInteractingPointsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const override
    {
      forces = forces_pot_gp_1_;
    }

    void GetForcesAtAllInteractingPointsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const override
    {
      forces = forces_pot_gp_2_;
    }

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    void GetMomentsAtAllInteractingPointsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const override
    {
      moments = moments_pot_gp_1_;
    }

    void GetMomentsAtAllInteractingPointsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const override
    {
      moments = moments_pot_gp_2_;
    }

    /*!
    \brief Get interaction free energy / potential
    */
    double GetEnergy() const override { return interaction_potential_; }

    /** \brief print this beam potential-based element pair to screen
     *
     *  \author grill */
    void Print(std::ostream& out) const override;

    /** \brief print this beam potential element pair to screen
     *
     *  \author grill */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;
    //@}

    //! @name Access methods
    /*!
    \brief Get ptr to first beam element
    */
    inline const DRT::ELEMENTS::Beam3Base* BeamElement1() const { return beam_element1_; };

    /*!
    \brief Get ptr to second beam element
    */
    inline const DRT::ELEMENTS::Beam3Base* BeamElement2() const { return beam_element2_; };
    //@}

   private:
    //! @name Private evaluation methods

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using double length specific potential with approximation for large separations
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateFpotandStiffpot_LargeSepApprox(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
        CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22);

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
        CORE::LINALG::Matrix<3, 1, double> const& dist, double const& norm_dist,
        double const& norm_dist_exp1, double q1q2_JacFac_GaussWeights,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& norm_dist_exp1,
        double q1q2_JacFac_GaussWeights,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using double length specific potential with approximation for small separations
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
        CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22);


    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
        CORE::LINALG::Matrix<3, 1, double> const& dist, double const& norm_dist, double const& gap,
        double const& gap_regularized, double const& gap_exp1, double q1q2_JacFac_GaussWeights,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& gap,
        Sacado::Fad::DFad<double> const& gap_regularized, Sacado::Fad::DFad<double> const& gap_exp1,
        double q1q2_JacFac_GaussWeights,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using single length specific potential with approximation for small separations
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateFpotandStiffpot_SingleLengthSpecific_SmallSepApprox(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
        CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22);

    /** \brief Evaluate all quantities for the full disk-cylinder potential law
     *
     *  \author grill
     *  \date 03/19 */
    bool EvaluateFullDiskCylinderPotential(T& interaction_potential_GP,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_slave_GP,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_master_GP,
        CORE::LINALG::Matrix<3, 1, T> const& r_slave,
        CORE::LINALG::Matrix<3, 1, T> const& r_xi_slave,
        CORE::LINALG::Matrix<3, 1, T> const& t1_slave,
        CORE::LINALG::Matrix<3, 1, T> const& r_master,
        CORE::LINALG::Matrix<3, 1, T> const& r_xi_master,
        CORE::LINALG::Matrix<3, 1, T> const& r_xixi_master,
        CORE::LINALG::Matrix<3, 1, T> const& t1_master, T alpha, T cos_alpha,
        CORE::LINALG::Matrix<3, 1, T> const& dist_ul,
        CORE::LINALG::Matrix<1, 3, T> const& xi_master_partial_r_slave,
        CORE::LINALG::Matrix<1, 3, T> const& xi_master_partial_r_master,
        CORE::LINALG::Matrix<1, 3, T> const& xi_master_partial_r_xi_master,
        double prefactor_visualization_data,
        CORE::LINALG::Matrix<3, 1, double>& vtk_force_pot_slave_GP,
        CORE::LINALG::Matrix<3, 1, double>& vtk_force_pot_master_GP,
        CORE::LINALG::Matrix<3, 1, double>& vtk_moment_pot_slave_GP,
        CORE::LINALG::Matrix<3, 1, double>& vtk_moment_pot_master_GP,
        double rho1rho2_JacFac_GaussWeight,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T> const& N_i_master,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T> const& N_i_xi_master);

    /** \brief scale given stiffness matrices by given scalar factor
     *
     *  \author grill
     *  \date 10/17 */
    void ScaleStiffpotAnalyticContributionsIfRequired(double const& scalefactor,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief scale given stiffness matrices by given scalar factor
     *
     *  \author grill
     *  \date 10/17 */
    void ScaleStiffpotAnalyticContributionsIfRequired(Sacado::Fad::DFad<double> const& scalefactor,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void CalcStiffmatAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type double, no automatic differentiation is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void CalcStiffmatAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot2,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief add contributions from linearization of parameter coordinate on master beam xi_master
     *         if determined via point-to-curve projection
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void AddStiffmatContributionsXiMasterAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot2,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_slaveDofs,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_masterDofs,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type double, no automatic differentiation is required
    }

    /** \brief add contributions from linearization of parameter coordinate on master beam xi_master
     *         if determined via point-to-curve projection
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void AddStiffmatContributionsXiMasterAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot2,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     *  \author grill
     *  \date 02/19 */
    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     *  \author grill
     *  \date 02/19 */
    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     *  \author grill
     *  \date 02/19 */
    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot1,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot2,
        double const& interaction_potential,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_slaveDofs,
        CORE::LINALG::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_masterDofs) const
    {
      // this is a dummy since for type double, no automatic differentiation is available
    }

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     *  \param stiffmat11 d (Residuum_vec_1) / d (dof_vec_1)
     *
     *  \author grill
     *  \date 04/19 */
    void EvaluateStiffpotAnalyticContributions_SingleLengthSpecific_SmallSepApprox_Simple(
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_master,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_master,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xixi_master,
        double const& xi_master, CORE::LINALG::Matrix<3, 1, double> const& r_xi_slave,
        CORE::LINALG::Matrix<3, 1, double> const& r_xi_master,
        CORE::LINALG::Matrix<3, 1, double> const& r_xixi_master, double const& norm_dist_ul,
        CORE::LINALG::Matrix<3, 1, double> const& normal_ul, double const& pot_ia_deriv_gap_ul,
        double const& pot_ia_deriv_cos_alpha, double const& pot_ia_2ndderiv_gap_ul,
        double const& pot_ia_deriv_gap_ul_deriv_cos_alpha, double const& pot_ia_2ndderiv_cos_alpha,
        CORE::LINALG::Matrix<3, 1, double> const& gap_ul_deriv_r_slave,
        CORE::LINALG::Matrix<3, 1, double> const& gap_ul_deriv_r_master,
        CORE::LINALG::Matrix<3, 1, double> const& cos_alpha_deriv_r_slave,
        CORE::LINALG::Matrix<3, 1, double> const& cos_alpha_deriv_r_master,
        CORE::LINALG::Matrix<3, 1, double> const& cos_alpha_deriv_r_xi_slave,
        CORE::LINALG::Matrix<3, 1, double> const& cos_alpha_deriv_r_xi_master,
        CORE::LINALG::Matrix<1, 3, double> const& xi_master_partial_r_slave,
        CORE::LINALG::Matrix<1, 3, double> const& xi_master_partial_r_master,
        CORE::LINALG::Matrix<1, 3, double> const& xi_master_partial_r_xi_master,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     *  \author grill
     *  \date 04/19 */
    void EvaluateStiffpotAnalyticContributions_SingleLengthSpecific_SmallSepApprox_Simple(
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_master,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_xi_master,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_xixi_master,
        Sacado::Fad::DFad<double> const& xi_master,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xi_slave,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xi_master,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xixi_master,
        Sacado::Fad::DFad<double> const& norm_dist_ul,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& normal_ul,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_gap_ul,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_cos_alpha,
        Sacado::Fad::DFad<double> const& pot_ia_2ndderiv_gap_ul,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_gap_ul_deriv_cos_alpha,
        Sacado::Fad::DFad<double> const& pot_ia_2ndderiv_cos_alpha,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& gap_ul_deriv_r_slave,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& gap_ul_deriv_r_master,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_slave,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_master,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_xi_slave,
        CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_xi_master,
        CORE::LINALG::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_slave,
        CORE::LINALG::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_master,
        CORE::LINALG::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_xi_master,
        CORE::LINALG::SerialDenseMatrix& stiffmat11, CORE::LINALG::SerialDenseMatrix& stiffmat12,
        CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /*!
    \brief Calculate shape function values for given parameter values
    Todo call more general utils method
    */
    void GetShapeFunctions(
        std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double>>& N1_i,
        std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double>>& N2_i,
        std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double>>& N1_i_xi,
        std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues, double>>& N2_i_xi,
        CORE::FE::IntegrationPoints1D& gausspoints) const;

    /*!
    \brief Compute coordinates of centerline point
    */
    template <typename T2>
    void ComputeCenterlinePosition(CORE::LINALG::Matrix<3, 1, T>& r,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T2>& N_i,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T> eledofvec) const;

    /*!
    \brief Compute tangent vector at centerline point
    */
    template <typename T2>
    void ComputeCenterlineTangent(CORE::LINALG::Matrix<3, 1, T>& r_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T2>& N_i_xi,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T> eledofvec) const;

    /** \brief set primary variables for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele1centerlinedofvec,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele2centerlinedofvec)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele1centerlinedofvec,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele2centerlinedofvec);

    /** \brief set primary variables including xi_master for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele1centerlinedofvec,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele2centerlinedofvec,
        double& xi_master)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables including xi_master for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele1centerlinedofvec,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele2centerlinedofvec,
        Sacado::Fad::DFad<double>& xi_master);

    /** \brief estimate whether the elements' separation is much more than the cutoff distance
     *
     *  The idea is to get a conservative lower bound estimate for the minimal centerline separation
     *  of both elements, which doesn't require any interpolation, i.e., information on the
     *  deformation of the elements and is thus cheap to compute before we begin with the actual,
     *  expensive evaluation of the pair on Gauss point level.
     *
     *  \author grill
     *  \date 08/19 */
    bool AreElementsMuchMoreSeparatedThanCutoffDistance();

    //@}

   private:
    //! @name member variables

    //! first element of pair
    DRT::ELEMENTS::Beam3Base const* beam_element1_;

    //! second element of pair
    DRT::ELEMENTS::Beam3Base const* beam_element2_;

    //! line charge conditions
    std::vector<DRT::Condition*> linechargeconds_;

    //! current time
    double time_;

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T> ele1pos_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, T> ele2pos_;

    //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
    double k_;
    double m_;

    //! initial element lengths
    double ele1length_;
    double ele2length_;

    //! Cross-section radius of beam 1
    double radius1_;

    //! Cross-section radius of beam 2
    double radius2_;
    //@}

    //! @name data storage for visualization output

    //! centerline coordinate vector of interacting points of element 1 and 2
    std::vector<CORE::LINALG::Matrix<3, 1, double>> centerline_coords_gp_1_;
    std::vector<CORE::LINALG::Matrix<3, 1, double>> centerline_coords_gp_2_;

    //! resulting forces at interacting points of element 1 and 2
    std::vector<CORE::LINALG::Matrix<3, 1, double>> forces_pot_gp_1_;
    std::vector<CORE::LINALG::Matrix<3, 1, double>> forces_pot_gp_2_;

    //! resulting moments at interacting points of element 1 and 2
    std::vector<CORE::LINALG::Matrix<3, 1, double>> moments_pot_gp_1_;
    std::vector<CORE::LINALG::Matrix<3, 1, double>> moments_pot_gp_2_;

    //! total interaction potential of this pair
    double interaction_potential_;
    //@}
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
