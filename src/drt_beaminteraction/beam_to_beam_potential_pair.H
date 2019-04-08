/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_beam_potential_pair.H

\brief One beam-to-beam potential-based interacting pair (two beam elements)

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BEAMTOBEAMPOTENTIALPAIR_H
#define BEAMTOBEAMPOTENTIALPAIR_H

#include "../drt_beaminteraction/beam_potential_pair.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Sacado.hpp>


// forward declarations
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

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

    /*!
    \brief Destructor
    */
    virtual ~BeamToBeamPotentialPair(){};

    //! Setup
    virtual void Setup();
    //@}

    //! @name Derived methods from base class
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    virtual bool Evaluate(LINALG::SerialDenseVector* forcevec1,
        LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
        LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
        LINALG::SerialDenseMatrix* stiffmat22, const std::vector<DRT::Condition*> linechargeconds,
        const double k, const double m);

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    virtual void ResetState(double time, const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2);

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    virtual void GetAllInteractingPointCoordsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const
    {
      coords = centerline_coords_GP1_;
    }

    virtual void GetAllInteractingPointCoordsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const
    {
      coords = centerline_coords_GP2_;
    }

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    virtual void GetForcesAtAllInteractingPointsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& forces) const
    {
      forces = forces_pot_GP1_;
    }

    virtual void GetForcesAtAllInteractingPointsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& forces) const
    {
      forces = forces_pot_GP2_;
    }

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    virtual void GetMomentsAtAllInteractingPointsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& moments) const
    {
      moments = moments_pot_GP1_;
    }

    virtual void GetMomentsAtAllInteractingPointsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& moments) const
    {
      moments = moments_pot_GP2_;
    }

    /*!
    \brief Get interaction free energy / potential
    */
    double GetEnergy() const override { return interaction_potential_; }

    /** \brief print this beam potential-based element pair to screen
     *
     *  \author grill */
    virtual void Print(std::ostream& out) const;

    /** \brief print this beam potential element pair to screen
     *
     *  \author grill */
    virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const;
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
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22);

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
        LINALG::TMatrix<double, 3, 1> const& dist, double const& norm_dist,
        double const& norm_dist_exp1, double q1q2_JacFac_GaussWeights,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 1> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& norm_dist_exp1,
        double q1q2_JacFac_GaussWeights,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using double length specific potential with approximation for small separations
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22);


    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<double, 3, 1> const& dist, double const& norm_dist, double const& gap,
        double const& gap_regularized, double const& gap_exp1, double q1q2_JacFac_GaussWeights,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 1> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& gap,
        Sacado::Fad::DFad<double> const& gap_regularized, Sacado::Fad::DFad<double> const& gap_exp1,
        double q1q2_JacFac_GaussWeights,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using single length specific potential with approximation for small separations
     *
     *  \author grill
     *  \date 10/17 */
    void EvaluateFpotandStiffpot_SingleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22);

    /** \brief Evaluate all quantities for the full disk-cylinder potential law
     *
     *  \author grill
     *  \date 03/19 */
    bool EvaluateFullDiskCylinderPotential(T& interaction_potential_GP,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot_slave_GP,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot_master_GP,
        LINALG::TMatrix<T, 3, 1> const& r_slave, LINALG::TMatrix<T, 3, 1> const& r_xi_slave,
        LINALG::TMatrix<T, 3, 1> const& g1_slave, LINALG::TMatrix<T, 3, 1> const& r_master,
        LINALG::TMatrix<T, 3, 1> const& r_xi_master, LINALG::TMatrix<T, 3, 1> const& r_xixi_master,
        LINALG::TMatrix<T, 3, 1> const& g1_master, T alpha, T cos_alpha,
        LINALG::TMatrix<T, 3, 1> const& dist_ul,
        LINALG::TMatrix<T, 1, 3> const& xi_master_partial_r_slave,
        LINALG::TMatrix<T, 1, 3> const& xi_master_partial_r_master,
        LINALG::TMatrix<T, 1, 3> const& xi_master_partial_r_xi_master, double prefactor_vtk,
        LINALG::TMatrix<double, 3, 1>& vtk_force_pot_slave_GP,
        LINALG::TMatrix<double, 3, 1>& vtk_force_pot_master_GP,
        LINALG::TMatrix<double, 3, 1>& vtk_moment_pot_slave_GP,
        LINALG::TMatrix<double, 3, 1>& vtk_moment_pot_master_GP, double rho1rho2_JacFac_GaussWeight,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_slave,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_xi_slave,
        LINALG::TMatrix<T, 1, numnodes * numnodalvalues> const& N_i_master,
        LINALG::TMatrix<T, 1, numnodes * numnodalvalues> const& N_i_xi_master);

    /** \brief scale given stiffness matrices by given scalar factor
     *
     *  \author grill
     *  \date 10/17 */
    void ScaleStiffpotAnalyticContributionsIfRequired(double const& scalefactor,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief scale given stiffness matrices by given scalar factor
     *
     *  \author grill
     *  \date 10/17 */
    void ScaleStiffpotAnalyticContributionsIfRequired(Sacado::Fad::DFad<double> const& scalefactor,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void CalcStiffmatAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> const& force_pot1,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> const& force_pot2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type double, no automatic differentiation is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void CalcStiffmatAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot1,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief add contributions from linearization of parameter coordinate on master beam xi_master
     *         if determined via point-to-curve projection
     *         using automatic differentiation based on Sacado::Fad package
     *
     *  \author grill
     *  \date 10/17 */
    void AddStiffmatContributionsXiMasterAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> const& force_pot1,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> const& force_pot2,
        LINALG::TMatrix<double, 1, 3 * numnodes * numnodalvalues> const& lin_xi_master_slaveDofs,
        LINALG::TMatrix<double, 1, 3 * numnodes * numnodalvalues> const& lin_xi_master_masterDofs,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
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
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot1,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot2,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_slaveDofs,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_masterDofs,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_slaveDofs,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_masterDofs) const;

    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_slaveDofs,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_masterDofs) const;

    void CalcFpotGausspointAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        double const& interaction_potential,
        LINALG::TMatrix<double, 1, 3 * numnodes * numnodalvalues> const& lin_xi_master_slaveDofs,
        LINALG::TMatrix<double, 1, 3 * numnodes * numnodalvalues> const& lin_xi_master_masterDofs)
        const
    {
      // this is a dummy since for type double, no automatic differentiation is available
    }

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     *  \author grill
     *  \date 04/19 */
    void EvaluateStiffpotAnalyticContributions_SingleLengthSpecific_SmallSepApprox_Simple(
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_slave,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_xi_slave,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_master,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_xi_master,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const;

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     *  \author grill
     *  \date 04/19 */
    void EvaluateStiffpotAnalyticContributions_SingleLengthSpecific_SmallSepApprox_Simple(
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_slave,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N_i_xi_slave,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, numnodes * numnodalvalues> const& N_i_master,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, numnodes * numnodalvalues> const&
            N_i_xi_master,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /*!
    \brief Calculate shape function values for given parameter values
    Todo call more general utils method
    */
    void GetShapeFunctions(std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N1_i,
        std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N2_i,
        std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N1_i_xi,
        std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N2_i_xi,
        DRT::UTILS::IntegrationPoints1D& gausspoints) const;

    /*!
    \brief Compute coordinates of centerline point
    */
    template <typename T2>
    void ComputeCenterlinePosition(LINALG::TMatrix<T, 3, 1>& r,
        const LINALG::TMatrix<T2, 1, numnodes * numnodalvalues>& N_i,
        const LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> eledofvec) const;

    /*!
    \brief Compute tangent vector at centerline point
    */
    template <typename T2>
    void ComputeCenterlineTangent(LINALG::TMatrix<T, 3, 1>& r_xi,
        const LINALG::TMatrix<T2, 1, numnodes * numnodalvalues>& N_i_xi,
        const LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> eledofvec) const;

    /** \brief set primary variables for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& ele1centerlinedofvec,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& ele2centerlinedofvec)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
            ele1centerlinedofvec,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
            ele2centerlinedofvec);

    /** \brief set primary variables including xi_master for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& ele1centerlinedofvec,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& ele2centerlinedofvec,
        double& xi_master)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables including xi_master for FAD if required
     *
     *  \author grill
     *  \date 10/17 */
    void SetAutomaticDifferentiationVariablesIfRequired(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
            ele1centerlinedofvec,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
            ele2centerlinedofvec,
        Sacado::Fad::DFad<double>& xi_master);


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
    LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> ele1pos_;
    LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> ele2pos_;

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
    std::vector<LINALG::TMatrix<double, 3, 1>> centerline_coords_GP1_;
    std::vector<LINALG::TMatrix<double, 3, 1>> centerline_coords_GP2_;

    //! resulting forces at interacting points of element 1 and 2
    std::vector<LINALG::TMatrix<double, 3, 1>> forces_pot_GP1_;
    std::vector<LINALG::TMatrix<double, 3, 1>> forces_pot_GP2_;

    //! resulting moments at interacting points of element 1 and 2
    std::vector<LINALG::TMatrix<double, 3, 1>> moments_pot_GP1_;
    std::vector<LINALG::TMatrix<double, 3, 1>> moments_pot_GP2_;

    //! total interaction potential of this pair
    double interaction_potential_;
    //@}
  };
}  // namespace BEAMINTERACTION

#endif
