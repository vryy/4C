/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an active stress material

\level 2
*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_ANISOACTIVESTRESS_EVOLUTION_HPP
#define FOUR_C_MATELAST_ANISOACTIVESTRESS_EVOLUTION_HPP

#include "4C_config.hpp"

#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_mat_par_aniso.hpp"
#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * <h3>Input line</h3>
       * MAT 1 ELAST_AnisoActiveStress_Evolution SIGMA 100.0 TAUC0 0.0 MAX_ACTIVATION 30.0
       * MIN_ACTIVATION -20.0 SOURCE_ACTIVATION 1 ACTIVATION_THRES 0 [STRAIN_DEPENDENCY No]
       * [LAMBDA_LOWER 0.707] [LAMBDA_UPPER 1.414]
       */
      class AnisoActiveStressEvolution : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        explicit AnisoActiveStressEvolution(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double sigma_;
        /// initial condition
        double tauc0_;
        /// Maximal value for rescaling the activation curve
        double maxactiv_;
        /// Minimal value for rescaling the activation curve
        double minactiv_;
        /// Threshold for stress activation function
        double activationthreshold_;
        /// Where the activation comes from: 0=scatra , >0 Id for FUNCT
        int sourceactiv_;
        /// is there a strain dependency for the active tension?
        bool strain_dep_;
        /// lower stretch threshold
        double lambda_lower_;
        /// upper stretch threshold
        double lambda_upper_;
        /// fiber angle
        double gamma_;
        /// fiber initalization status
        int init_;
        /// adapt angle during remodeling
        bool adapt_angle_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<CORE::MAT::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class AnisoActiveStress_Evolution

    }  // namespace PAR

    /*!
     * This is a simplification of the muscle contraction law proposed in [1],[2],
     * resulting in the following first order ODE for the active stress tau, compare [3]:
     *
     * \f[
     *   \frac{d}{dt} \tau = n_0 \sigma_0 |u|_+ - \tau |u|, \quad \tau(0) = tau_0
     * \f]
     *
     * where \f$\sigma_0\f$ is the contractility (asymptotic value of \tau) and u is a control
     * variable either provided by a electrophysiology simulation or by a user-specified function
     * therein, n0 is a strain-dependent factor that may take into account the Frank-Starling
     * effect! n0 \in [0; 1] scales the contractility depending on a lower and upper fiber stretch
     * lambda using a flipped parabola, n0 = -(lambda - lambda_lower)*(lambda - lambda_upper) *
     * 4/(lambda_lower-lambda_upper)^2, which is an approximation to the function from Sainte-Marie
     * et al. 2006, Fig. 2(ii) other laws might be thought of here, since hardly any literature
     * provides a meaningful dependency...
     *
     * Due to the active stress approach, see [4], the active stress will be added along a given
     * fiber direction f_0 to the 2nd Piola-Kirchhoff stress:
     * \f[
     *   S_{active} = \tau(t) f_0 \otimes f_0
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] 2012 Chapelle, Le Tallec, Moireau, Sorine - An energy-preserving muscle tissue
     * model: formulation and compatible discretizations, Journal for Multiscale Computational
     * Engineering 10(2):189-211 (2012) <li> [2] 2001 Bestel, Clement, Sorine - A Biomechanical
     * Model of Muscle Contraction (2001), Medical Image Computing and Computer-Assisted
     * Intervention (MICCAI'01), vol. 2208, Springer-Verlag Berlin, 1159-1161 <li> [3] 2002
     * Sermesant, Coudier, Delingette, Ayache - Progress towards an electromechanical model of the
     * heart for cardiac image analysis. (2002) IEEE International Symposium on Biomedical Imaging,
     * 10-13 <li> [4] 1998 Hunter, McCulloch, ter Keurs - Modelling the mechanical properties of
     * cardiac muscle (1998), Progress in Biophysics and Molecular Biology
     * </ul>
     */
    class AnisoActiveStressEvolution : public Summand
    {
     public:
      /// constructor with given material parameters
      explicit AnisoActiveStressEvolution(MAT::ELASTIC::PAR::AnisoActiveStressEvolution* params);

      ///@name Packing and Unpacking
      //@{

      void PackSummand(CORE::COMM::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      //@}

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_anisoactivestress_evolution;
      }

      //@}

      /*!
       * \brief Register the internally used AnisotropyExtension
       *
       * \param anisotropy Reference to the anisotropy manager
       */
      void register_anisotropy_extensions(MAT::Anisotropy& anisotropy) override;

      /// Setup of summand
      void Setup(int numgp, INPUT::LineDefinition* linedef) override;

      /*!
       * \brief post_setup routine of the element
       *
       * Here potential nodal fibers were passed to the Anisotropy framework
       *
       * @param params Container that potentially contains nodal fibers
       */
      void post_setup(Teuchos::ParameterList& params) override;

      /// Add anisotropic principal stresses
      void add_stress_aniso_principal(
          const CORE::LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          CORE::LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Set fiber directions
      void SetFiberVecs(double newgamma,             ///< new angle
          const CORE::LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const CORE::LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

      /// Setup of patient-specific materials
      void SetupAAA(Teuchos::ParameterList& params, const int eleGID) override {}

      // update internal stress variables
      void Update() override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        anisoprinc = true;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::AnisoActiveStressEvolution* params_;

      /// Active stress at current time step
      double tauc_np_;
      /// Active stress at last time step
      double tauc_n_;

      /// Special anisotropic behavior
      DefaultAnisotropyExtension<1> anisotropy_extension_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
