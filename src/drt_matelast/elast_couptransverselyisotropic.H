/*----------------------------------------------------------------------*/
/*! \file
\brief Summand for the transversely hyperelastic isotropic material model
\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef SRC_DRT_MATELAST_ELAST_COUPTRANSVERSELYISOTROPIC_H_
#define SRC_DRT_MATELAST_ELAST_COUPTRANSVERSELYISOTROPIC_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*--------------------------------------------------------------------------*/
      /** \brief parameter container for the simple orthotropic, transversely
       *  isotropic hyperelastic constitutive equation for large strain computations
       *
       *  These parameters belong to a simple orthotropic, transversely isotropic
       *  material model which is suitable for large strains. The main purpose is to
       *  combine it with a simple hyperelastic material model (e.g. St. Venant Kirchhoff or
       *  Neo-Hookean), such that a easy parameter identification becomes possible
       *  by consideration of the parameters for the small strain regime via
       *
       *  \f{eqnarray}{
       *  \lambda &=& \Lambda = \frac{\nu_{12} \nu_{21} + \nu_{22}}{(1-\nu_{22}-2\nu_{12}
       * \nu_{21})(1+\nu_{22})} E_{2}, \\
       *  \mu &=& G_{23} = \frac{E_{2}}{2(1+\nu_{23})}, \\
       *  \alpha &=& G_{23} - G_{12}, \\
       *  \beta &=& \nu_{12} (\lambda + G_{23}) - \frac{\lambda}{2}, \\
       *  \gamma &=& \frac{1}{8} \left(\frac{1 -\nu_{23}}{1-\nu_{23} - \nu_{12} \nu_{21}} E_{1} + 4
       * \alpha - 4 \beta - \lambda - 2 \mu \right). \f} where \f$ E_{2} = E_{3} \f$ is the Young
       * modulus in the plane normal to the fiber direction and \f$ E_{1} = E_{A} \f$ denotes the
       * Young's modulus in fiber direction. \f$ \nu_{12} \f$ is the Poisson ratio for tension in
       * fiber direction. Furthermore, \f$ G_{12}=G_{13}=G_{A} \f$ represents shear modulus in fiber
       * direction, while the shear modulus in the normal plane \f$ G_{23} \f$ dependents on the
       * already known isotropic parameters via \f$ G_{23} = (E_{2}) / (2(1+\nu_{23})) \f$.
       *  Furthermore, the relation \f$ \frac{\nu_{21}}{E_{2}} = \frac{\nu_{12}}{E_{1}} \f$
       *  holds.
       *
       *  The implementation follows partly the reference
       *
       *  J. Bonet, A. J. Burton, "A simple orthotropic, transversely isotropic
       *  hyperelastic constitutive equation for large strain computations",
       *  Comput. Methods Appl. Mech. Engrg., 162, pp. 151-164, (1998)
       *
       *  The parameter identification is based on the linear compliance matrix
       *  given in
       *
       *  H. Schuermann, "Konstruieren mit Faser-Kunststoff-Verbunden", Springer, 2005.
       *
       *  \author hiermeier \date 09/17 */
      class CoupTransverselyIsotropic : public MAT::PAR::ParameterAniso
      {
       public:
        /// no default constructor
        CoupTransverselyIsotropic() = delete;

        /// standard constructor
        explicit CoupTransverselyIsotropic(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// default destructor
        virtual ~CoupTransverselyIsotropic() = default;

        /// alpha parameter of the constitutive law
        double alpha_ = 0.0;

        /// beta parameter of the constitutive law
        double beta_ = 0.0;

        /// gamma parameter of the constitutive law
        double gamma_ = 0.0;

        /// angle of the fiber direction
        double angle_ = 0.0;

        /// global fiber id number (1,2,3,...) used later as FIBER1,FIBER2,FIBER3,...
        int fiber_gid_ = 0;

        /// fiber initalization status
        int init_ = 0;

        /// print parameters
        void Print() const;

       private:
        /// derived function
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; }
      };  // class CoupAnisoSimple
    }     // namespace PAR

    /*--------------------------------------------------------------------------*/
    /** \brief Implementation of the simple orthotropic, transversely
     *  isotropic hyperelastic constitutive equation for large strain computations
     *
     *  The implementation follows partly the reference
     *
     *  J. Bonet, A. J. Burton, "A simple orthotropic, transversely isotropic
     *  hyperelastic constitutive equation for large strain computations",
     *  Comput. Methods Appl. Mech. Engrg., 162, pp. 151-164, (1998).
     *
     *  \note The given reference contains many typos, such that almost all
     *  equations (beside the strain energy function) were derived independently.
     *
     *  \author hiermeier \date 09/17 */
    class CoupTransverselyIsotropic : public Summand
    {
      typedef MAT::ELASTIC::PAR::CoupTransverselyIsotropic my_params;

     public:
      /// no default constructor
      CoupTransverselyIsotropic() = delete;

      /// standard constructor
      explicit CoupTransverselyIsotropic(my_params* params);

      /// default destructor
      virtual ~CoupTransverselyIsotropic() = default;

      /// @name Packing and Unpacking
      /// @{

      virtual void PackSummand(DRT::PackBuffer& data) const;

      virtual void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position);

      /// @}

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_couptransverselyisotropic;
      }

      /// Setup of summand
      virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

      /** \brief add strain energy [derived]
       *
       *  \f[
       *  \Psi_{\mathrm{trn}} = \left( \alpha + \frac{\beta}{2} \ln(I_{3}) +
       *  \gamma (I_{4}-1) \right) (I_{4}-1) - \frac{\alpha}{2} ( I_{5} - 1)
       *  \f]
       *
       *  \note \f$ \ln(J) = \ln(\sqrt(I_{3})) = \frac{1}{2} \ln(I_{3}) \f$
       *  holds.
       *
       *  \param[out] psi      Updated value of the strain energy function
       *  \param[in]  prinv    principal invariants of right Cauchy-Green tensor
       *  \param[in]  modinv   modified invariants of right Cauchy-Green tensor
       *  \param[in]  glstrain Green-Lagrange strain in strain like Voigt notation
       *  \param[in]  eleGID   element global ID
       *
       *  \author hiermeier \date 09/17 */
      virtual void AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
          const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain,
          const int eleGID);

      /** \brief Add anisotropic principal stresses
       *
       *  \param[in]  rcg    right Cauchy Green Tensor
       *  \param[out] cmat   material stiffness matrix
       *  \param[out] stress 2nd PK-stress
       *  \param[in]  params additional parameters for comp. of material properties
       *  \param[in]  eleGID element GID
       *
       *  \author hiermeier \date 09/17 */
      virtual void AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
          LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
          const int eleGID);

      /** \brief Set fiber directions
       *
       *  \param[in] newangle  new angle for the fiber direction
       *  \param[in] locsys    local coordinate system
       *  \param[in] defgrd    deformation gradient
       * */
      virtual void SetFiberVecs(const double newangle, const LINALG::Matrix<3, 3>& locsys,
          const LINALG::Matrix<3, 3>& defgrd);

      /// Get fiber directions
      virtual void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
      );

      /** \brief Indicator for formulation
       *
       *  \param[out] isoprinc     global indicator for isotropic principal formulation
       *  \param[out] isomod       global indicator for isotropic splitted formulation
       *  \param[out] anisoprinc   global indicator for anisotropic principal formulation
       *  \param[out] anisomod     global indicator for anisotropic splitted formulation
       *  \param[out] viscogeneral global indicator, if one viscoelastic formulation is used
       */
      virtual void SpecifyFormulation(
          bool& isoprinc, bool& isomod, bool& anisoprinc, bool& anisomod, bool& viscogeneral)
      {
        anisoprinc = true;
        return;
      };

     private:
      /** \brief Reset the pseudo invariants and the jacobian determinant
       *
       *  \param[in] rcg   right cauchy green tensor in perturbed Voigt strain notation
       *  \param[in] param parameter list pointer (optional)
       *
       *  \author hiermeier \date 09/17 */
      int ResetInvariants(
          const LINALG::Matrix<6, 1>& rcg, const Teuchos::ParameterList* params = NULL);

      /** \brief Add the material contributions to the second Piola Kirchhoff stress
       *  tensor
       *
       *  \f{eqnarray}{
       *  \mathbf{S}_{\mathrm{trn}} =& 2 \left\{ \frac{\partial \Psi_{\mathrm{trn}} }{\partial
       * I_{3}} \frac{\partial I_{3}}{\partial \mathbf{C}}
       *  + \frac{\partial \Psi_{\mathrm{trn}} }{\partial I_{4}} \frac{\partial I_{4}}{\partial
       * \mathbf{C}}
       *  + \frac{\partial \Psi_{\mathrm{trn}} }{\partial I_{5}} \frac{\partial I_{5}}{\partial
       * \mathbf{C}} \right\} \nonumber \\
       *  =& 2 \left\{
       *  \frac{\beta}{2 I_{3}} (I_{4}-1) I_{3} \mathbf{C}^{-1} + \left[\alpha + \frac{\beta}{2}
       * \ln(I_{3}) + 2 \gamma (I_{4}-1)\right] \mathbf{A} \otimes \mathbf{A}
       *  - \frac{\alpha}{2} ( \mathbf{A} \otimes \mathbf{C}\, \mathbf{A} + \mathbf{A}\, \mathbf{C}
       * \otimes \mathbf{A}
       * )\right\} \nonumber \\
       *  =& \beta (I_{4} - 1) \mathbf{C}^{-1} + 2\left[\alpha + \frac{\beta}{2} \ln(I_{3}) + 2
       * \gamma (I_{4}-1)\right] \mathbf{A} \otimes \mathbf{A}
       *  - \alpha ( \mathbf{A} \otimes \mathbf{C}\, \mathbf{A} + \mathbf{A}\, \mathbf{C} \otimes
       * \mathbf{A} ). \f}
       *
       *  \param[out] stress     Updated second Piola Kirchhoff stress tensor
       *  \param[in]  rcg_s      right Cauchy Green strain tensor in perturbed
       *                         Voigt stress notation
       *  \param[out] rcg_inv_s  Inverse of the right Cauchy green strain tensor
       *                         in perturbed Voigt stress notation
       *
       *  \author hiermeier \date 09/17 */
      void UpdateSecondPiolaKirchhoffStress(LINALG::Matrix<6, 1>& stress,
          const LINALG::Matrix<6, 1>& rcg_s, LINALG::Matrix<6, 1>& rcg_inv_s) const;

      /** \brief Update the elasticity tensor
       *
       *  \f{eqnarray}{
       *  \mathbb{C}_{\mathrm{trn}} = [\mathbb{C}_{\mathrm{trn}}]_{ABCD} =& 2 \frac{\partial
       * [\mathbf{S}_{\mathrm{trn}}]_{AB}}{\partial [\mathbf{C}]_{CD}} =
       *  2 \{ \beta [\mathbf{C}^{-1}]_{AB} [\underline{A} \otimes \underline{A}]_{CD}  \nonumber \\
       *  & - \beta (I_{4} -1 ) \frac{1}{2} \left([\mathbf{C}^{-1}]_{AC}[\mathbf{C}^{-1}]_{BD} +
       * [\mathbf{C}^{-1}]_{AD} [\mathbf{C}^{-1}]_{BC} \right) \nonumber \\
       *  &+ \beta \frac{1}{I_{3}} I_{3} [\mathbf{C}^{-1}]_{CD} [\underline{A} \otimes
       * \underline{A}]_{AB} + 4
       * \gamma  [\underline{A} \otimes \underline{A} \otimes \underline{A} \otimes
       * \underline{A}]_{ABCD} \nonumber \\
       *  &-\frac{\alpha}{2} \left( [\underline{A}]_{A} [\underline{A}]_{D} [\mathbf{\delta}]_{BC} +
       *  [\underline{A}]_{A} [\underline{A}]_{C} [\mathbf{\delta}]_{BD} +
       *  [\underline{A}]_{C} [\underline{A}]_{B} [\mathbf{\delta}]_{AD} \right. \nonumber \\
       *  &\left. +
       *  [\underline{A}]_{D} [\underline{A}]_{B} [\mathbf{\delta}]_{AC} \right) \nonumber \\
       *  =& 2 \beta \left( \mathbf{C}^{-1} \otimes \underline{A} \otimes \underline{A} +
       * \underline{A} \otimes \underline{A}\otimes \mathbf{C}^{-1} \right) -2 \beta (I_{4}-1)
       * \mathbf{C}^{-1} \odot \mathbf{C}^{-1}
       * \nonumber \\
       *  &-2 \alpha \left( (\underline{A} \otimes \underline{A} ) \odot \mathbf{\delta}
       *  + \mathbf{\delta} \odot (\underline{A} \otimes \underline{A} ) \right)
       *  \f}
       *
       *  The \f$\odot\f$ symbol follows the definition given in Holzapfel's book
       *  "Nonlinear Solid Mechanics - A continuum approach for engineering" on
       *  page 254, equation (6.164) and (6.165).
       *
       *  \param[out] cmat       updated elasticity tensor
       *  \param[in]  rcg_inv_s  right cauchy green tensor in perturbed Voigt
       *                         stress notation
       *
       *  \author hiermeier \date 09/17 */
      void UpdateElasticityTensor(
          LINALG::Matrix<6, 6>& cmat, const LINALG::Matrix<6, 1>& rcg_inv_s) const;

      /// error handling in case of a negative deformation gradient determinant
      void ErrorHandling(const Teuchos::ParameterList* params, std::stringstream& msg) const;

     private:
      /// pointer to the fiber parameters
      my_params* params_ = NULL;

      /// fiber direction
      LINALG::Matrix<3, 1> A_;

      /** \brief outer product of the fiber directions
       *  \f$ \underline{A} \otimes \underline{A}\f$
       *
       *  \note We are following a perturbed Voigt notation:
       *  {11, 22, 33, 12, 23, 13}. */
      LINALG::Matrix<6, 1> AA_;

      /// pseudo invariant \f$ I_4 \f$ ( strain measure in fiber direction )
      double I4_ = 0.0;

      /// pseudo invariant \f$ I_5 \f$ ( quadratic strain measure in fiber direction )
      double I5_ = 0.0;

      /// determinant of the deformation gradient
      double j_ = 0.0;
    };  // class CoupAnisoSimple

  }  // namespace ELASTIC
}  // namespace MAT


#endif /* SRC_DRT_MATELAST_ELAST_COUPTRANSVERSELYISOTROPIC_H_ */
