/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric part of the anisotropic Blemker active skeletal
muscle material

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef ELAST_ISOMUSCLEBLEMKER_H_
#define ELAST_ISOMUSCLEBLEMKER_H_

#include "matelast_summand.H"
#include "mat_par_parameter.H"
#include "mat_anisotropy.H"
#include "mat_anisotropy_extension_default.H"
#include "mat_anisotropy_extension_provider.H"
#include <Teuchos_RCP.hpp>
#include "lib_parobjectfactory.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      class IsoMuscleBlemker : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoMuscleBlemker(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        //! @name muscle shear moduli
        //! @{
        double G1_;  ///< along fiber shear modulus
        double G2_;  ///< cross fiber shear modulus
        //! @}

        //! @name along fiber parameters
        //! @{
        double P1_;           ///< linear material parameter for passive along-fiber response
        double P2_;           ///< exponential material parameter for passive along-fiber response
        double sigma_max_;    ///< maximal active isometric stress
        double lambda_ofl_;   ///< optimal fiber stretch
        double lambda_star_;  ///< stretch where normalized passive fiber force becomes linear
        //! @}

        //! @name time-dependent activation level
        //! @{
        double alpha_;        ///< tetanised activation level, 0<alpha<1
        double beta_;         ///< constant scaling tanh-type activation function, >=0
        double t_act_start_;  ///< starting time of muscle activation
        //! @}

        /// Create material instance of matching type with my parameters
        Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };
      };  // class IsoMuscleBlemker
    }     // namespace PAR

    /*!
     * \brief Isochoric part of the Blemker muscle material
     *
     * This material represents an active hyperelastic muscle material using an active stress
     * approach as described by Blemker et. al.
     *
     * Reference for the material model: S. S. Blemker, P. M. Pinsky und S. L. Delp, 'A 3D model of
     * muscle reveals the causes of nonuniform strains in the biceps brachii', Journal of
     * biomechanics, vol. 38, no. 4, pp. 657-665, 2005. doi: 10.1016/j.jbiomech.2004.04.009
     *
     * To account for a time dependent activation, a tanh-type function scaling the activation level
     * alpha is added in addition to the original version in the paper.
     */
    class IsoMuscleBlemker : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoMuscleBlemker(MAT::ELASTIC::PAR::IsoMuscleBlemker* params);

      /// Pack anisotropy
      void PackSummand(DRT::PackBuffer& data) const override;

      /// Unpack anisotropy
      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      /// Provide the material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_isomuscleblemker;
      }

      void RegisterAnisotropyExtensions(Anisotropy& anisotropy) override;

      /*!
       * \brief Add isochoric anisotropic stress and elasticity tensor
       *
       * Computation of the 2nd PK-stress and elasticity tensor with respect to the modified strains
       */
      void AddStressAnisoModified(
          const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
                                            ///< in strain-like-voigt-Notation
          const LINALG::Matrix<6, 1>& icg,  ///< inverse of right Cauchy Green Tensor
                                            ///< in stress-like-voigt-notation
          LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          double I3,                        ///< third principal invariant
          int gp,                           ///< Gauss point
          int eleGID,                       ///< element GID
          Teuchos::ParameterList& params    ///< Container for additional information
          ) override;

      /// Specify the formulation as anisomod
      void SpecifyFormulation(bool& isoprinc, bool& isomod, bool& anisoprinc, bool& anisomod,
          bool& viscogeneral) override
      {
        anisomod = true;
        return;
      };

     protected:
      /// Blemker material parameters
      MAT::ELASTIC::PAR::IsoMuscleBlemker* params_;

     private:
      /// Anisotropy extension holder
      MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;

      /*!
       * \brief Evaluate total fiber cauchy stress and derivative w.r.t the fibre stretch
       *
       * The total fiber cauchy stress is computed as the normalized total fiber force (sum of
       * active and passive component) scaled by the current maximal isometric stress (i.e. the
       * maximal isometric stress scaled by the time-dependency ft) and the ratio between current
       * fiber length and optimal fiber length.
       *
       * In contrast to the original publication, we use lambda_ofl=1 to evaluate the passive
       * force-stretch dependency.
       *
       * \param[in] lambdaM fiber stretch
       * \param[in] sigma_max_ft maximal fiber cauchy stress at current time t (i.e. scaled by ft)
       * \param[in, out] sigma_fiber_total total fiber cauchy stress
       * \param[in, out] deriv_sigma_fiber_total derivative of total fiber cauchy stress w.r.t.
       * fibre stretch
       */
      void EvaluateTotalFiberCauchyStressAndDerivative(double lambdaM, double sigma_max_ft,
          double& sigma_fiber_total, double& deriv_sigma_fiber_total);
    };
  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ELAST_ISOMUSCLEBLEMKER_H_
