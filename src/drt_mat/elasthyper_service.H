/*----------------------------------------------------------------------*/
/*! \file

\brief Contains the declaration of servce functions for hylerelastic materials

\level 1

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ELASTHYPER_EVALUATE_H
#define BACI_ELASTHYPER_EVALUATE_H

#include <NOX.H>
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_matelast/elast_summand.H"
namespace MAT
{
  // Forward declaration
  class SummandProperties;

  /*!
   * \brief Evaluate the stress response and the elasticity tensor of an hyperelastic material
   *
   * This static method evaluates the 2nd Piola-Kirchhoff stress tensor and its linearization with
   * respect to the Green-Lagrange strain.
   *
   * @param defgrd      (in)      : deformation gradient
   * @param glstrain    (in)      : Green lagrange strain in strain-like Voigt notation
   * @param params      (in/out)  : Container for additional information
   * @param stress      (out)     : 2nd Piola Kirchhoff stress
   * @param cmat        (out)     : Elasticity tensor
   * @param gp          (in)      : Gauss point
   * @param eleGID      (in)      : Element id
   * @param potsum      (in)      : Summands of the Free-energy function
   * @param properties    (in)      : Data class with flags of the type of the summands
   * @param checkpolyconvexity (in) : Flag, whether to check the polyconvexity
   */
  void ElastHyperEvaluate(const LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<6, 1>& glstrain,
      Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& stress, LINALG::Matrix<6, 6>& cmat,
      int gp, int eleGID, const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      const SummandProperties& properties, bool checkpolyconvexity = false);

  /*!
   * Evaluates the Right Cauchy-Green strain tensor in strain like Voigt notation
   *
   * C = 2.0 * E + I
   *
   * @param E_strain (in) : Green-Lagrange in strain-like Voigt notation
   * @param C_strain (out) : Cauchy-Green in strain-like Voigt notation
   */
  void EvaluateRightCauchyGreenStrainLikeVoigt(
      const LINALG::Matrix<6, 1>& E_strain, LINALG::Matrix<6, 1>& C_strain);

  /*!
   * \brief Evaluates the first and second derivatives w.r.t. principal invariants
   *
   * @param prinv Principle invariants
   * @param dPI First derivative of potsum
   * @param ddPII Second derivative of potsum
   * @param potsum List of summands
   * @param properties Properties of the summands
   * @param gp Gauss point
   * @param eleGID Global element id
   */
  void ElastHyperEvaluateInvariantDerivatives(const LINALG::Matrix<3, 1>& prinv,
      LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      const SummandProperties& properties, int gp, int eleGID);

  /*!
   * \brief Converts the derivatives with respect to the modified principal invariants in
   * derivatives with respect to the principal invariants.
   *
   * The used conversions are:
   * \f[
   *   \overline{I}_{\boldsymbol{C}} = J^{-2/3} I_{\boldsymbol{C}},
   * \f]
   * \f[
   *   \overline{II}_{\boldsymbol{C}} = J^{-4/3} II_{\boldsymbol{C}},
   * \f]
   * \f[
   *   J = \sqrt{III_{\boldsymbol{C}}}
   * \f]
   *
   * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
   * @param dPmodI (in) : First derivatives of the Free-energy function with respect to the
   * modified principal invariants
   * @param ddPmodII (in) : Second derivatives of the Free-energy function with respect to the
   * modified principal invariants
   * @param dPI (out) : First derivatives of the Free-energy function with respect to the
   * principal invariants
   * @param ddPII (out) : Second derivatives of the Free-energy function with respect to the
   * principal invariants
   */
  void ConvertModToPrinc(const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPmodI,
      const LINALG::Matrix<6, 1>& ddPmodII, LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII);

  /*!
   * \brief Evaluates the Isotropic stress response from the potsum and adds it to the global stress
   * vector
   *
   * @param S_stress 2nd Piola-Kirchhoff stress tensor in stress-like Voigt-Notation
   * @param cmat Elasticity tensor in stress-like Voigt notation (rows and columns)
   * @param C_strain Right Cauchy-Green deformation tensor in strain like Voigt notation
   * @param iC_strain Inverse Right Cauchy-Green deformation tensor in strain like Voigt notation
   * @param prinv Principal invariants of the Right Cauchy-Green strain tensor
   * @param dPI First derivatives of the Free-energy function with respect to the
   * principal invariants
   * @param ddPII Second derivatives of the Free-energy function with respect to the
   * principal invariants
   */
  void ElastHyperAddIsotropicStressCmat(LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
      const LINALG::Matrix<6, 1>& C_strain, const LINALG::Matrix<6, 1>& iC_strain,
      const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
      const LINALG::Matrix<6, 1>& ddPII);

  /**
   * \brief Determine PK2 stress response and material elasticity tensor
   *  due to energy densities (MAT::ELASTIC::Summand) described
   *  in (modified) principal stretches.
   *
   *  The stress response is achieved by collecting the coefficients
   *  \f$\gamma_\alpha\f$ and \f$\delta_{\alpha\beta}\f$ due to
   *  MAT::ELASTIC::Summand::AddCoefficientsStretchesPrincipal()
   *  (and/or \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$
   *  due to MAT::ELASTIC::Summand::AddCoefficientsStretchesModified()).
   *  The collected coefficients build the principal 2nd Piola--Kirchhoff
   *  stresses which are transformed into ordinary Cartesian co-ordinates
   *  applying the principal directions. Similarly, the elasticity
   *  4-tensor is obtained.
   *
   *  Please note, unlike suggested by Holzapfel, p 263-264, the modified
   *  coefficients are transformed to unmodified coefficients and than added
   *  onto the respective quantities. The transformation goes along the following
   *  lines. The first derivatives or \f$\gamma_\alpha\f$ coefficients are related
   *  to the modified coefficients \f$\bar{\gamma}_\alpha\f$ via the chain rule
   * \f[
   *    \gamma_\alpha
   *    = \frac{\partial \Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha}
   *    = \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
   *    \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
   *    = \bar{\gamma}_\eta \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
   * \f]
   *  utilising Holzapfel Eq (6.142):
   * \f[
   *    \frac{\partial\bar{\lambda}_\alpha}{\partial\lambda_\beta}
   *    = J^{-1/3} \Big( 1_{\alpha\beta} - \frac{1}{3} \lambda_\alpha \lambda_\beta^{-1} \Big)
   * \f]
   *  in which Kronecker's delta is denoted \f$1_{\alpha\beta}\f$.
   *  Likewise (and once without Holzapfel), the second derivatives \f$\delta_{\alpha\beta}\f$ can
   *  be recovered directly by knowledge of the modified coefficients
   *  \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$, i.e.
   * \f[
   *    \delta_{\alpha\beta}
   *    = \frac{\partial^2\Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha
   *    \partial\lambda_\beta} =
   *    \frac{\partial^2\Psi(\bar{\boldsymbol{\lambda}})}{\partial\bar{\lambda}_\eta\partial\bar{\lambda}_\epsilon}
   *      \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
   *      \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
   *    + \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
   *      \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
   *    = \bar{\delta}_{\eta\epsilon}
   *      \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
   *      \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
   *    + \bar{\gamma}_\eta
   *      \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
   * \f]
   *  with
   * \f[
   *    \frac{\partial^2 \bar{\lambda}_\alpha}{\partial\lambda_\beta\partial\lambda_\eta}
   *    = -\frac{1}{3} J^{-1/3} \Big(
   *      1_{\alpha\beta} \lambda_\eta^{-1}
   *      + 1_{\alpha\eta}\lambda_\beta^{-1}
   *      - \big( 1_{\beta\eta} + \frac{1}{3} \big) \lambda_\alpha \lambda_\beta^{-1}
   *      \lambda_\eta^{-1}
   *    \Big)
   * \f]
   *
   *  <h3>References</h3>
   *  See Holzapfel, p 245-246, p 257-259, p 263-264
   *
   * \author bborn
   * \date 08/09
   * @param cmat (out) : Material elasticity tensor in Voigt notation
   * @param S_stress (out) : 2n Piola-Kirchhoff stress tensor in stress-like Voigt notation
   * @param C_strain (in) : Right Cauchy-Green strain tensor in strain-like Voigt notation
   * @param potsum (in) : Summands of the Free-energy function
   * @param properties (in) : Properties of the summands
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   */
  void ElastHyperAddResponseStretches(LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& S_stress,
      const LINALG::Matrix<6, 1>& C_strain,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      const SummandProperties& properties, int gp, int eleGID);

  /*!
   * \brief Evaluates the anisotropic stress response from the potsum elements formulated in the
   * principle invariants and add it to the global stress tensor.
   *
   * @param S_stress 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation
   * @param cmat Elasticity tensor in stress-like Voigt notation (rows and columns)
   * @param C_strain Right Cauchy-Green deformation tensor in strain like Voigt notation
   * @param params Container for additional information
   * @param gp Gauss point
   * @param eleGID Global element id
   * @param potsum Summands of the Free-energy function
   */
  void ElastHyperAddAnisotropicPrinc(LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
      const LINALG::Matrix<6, 1>& C_strain, Teuchos::ParameterList& params, int gp, int eleGID,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

  /*!
   * \brief Evaluates the anisotropic stress response from the potsum elements formulated in the
   * modified principle invariants and add it to the global stress tensor.
   *
   * @param S_stress 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation
   * @param cmat Elasticity tensor in stress-like Voigt notation (rows and columns)
   * @param C_strain Right Cauchy-Green deformation tensor in strain like Voigt notation
   * @param iC_strain Inverse Right Cauchy-Green deformation tensor in strain like Voigt notation
   * @param prinv Principal invariants of the Right Cauchy-Green strain tensor
   * @param eleGID Global element id
   * @param potsum Summands of the Free-energy function
   */
  void ElastHyperAddAnisotropicMod(LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
      const LINALG::Matrix<6, 1>& C_strain, const LINALG::Matrix<6, 1>& iC_strain,
      const LINALG::Matrix<3, 1>& prinv, int gp, int eleGID,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

  /*!
   * \brief Calculate coefficients gamma and delta from partial derivatives w.r.t. invariants.
   *
   * The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
   *  on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$
   *  in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
   *  are the arguments, cf. Holzapfel [1],
   * \f[
   *   I_{\boldsymbol{C}} = C_{AA},
   *   \quad
   *   II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) -
   *   \mathrm{trace}(\boldsymbol{C}^2) \big), \quad III_{\boldsymbol{C}} = \det(\boldsymbol{C})
   * \f]
   *
   *  cf. Holzapfel [1], p. 216 (6.32) and p. 248
   * \f[
   *   \bold{S} = \gamma_1 \ \bold{Id} + \gamma_2 \ \bold{C} + \gamma_3 \ \bold{C}^{-1}
   * \f]
   * \f[
   *  \gamma_1 = 2\left( \frac{\partial \Psi}{\partial I_{\boldsymbol{C}}}
   *           + I_{\boldsymbol{C}}\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}} \right);
   * \f]
   * \f[
   *  \gamma_2 = -2\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}};
   * \f]
   * \f[
   *  \gamma_3 = 2III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}};
   * \f]
   *
   *  material constitutive tensor coefficients
   *  cf. Holzapfel [1], p. 261
   * \f[
   *   \mathbb{C} = \delta_1 \left( \bold{Id} \otimes \bold{Id} \right) + \delta_2 \left(
   *   \bold{Id} \otimes \bold{C} + \bold{C} \otimes \bold{Id} \right)
   *   + \delta_3 \left( \bold{Id} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{Id} \right)
   *   + \delta_4 \left( \bold{C} \otimes \bold{C} \right)
   *   + \delta_5 \left( \bold{C} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{C} \right) +
   *   \delta_6 \left( \bold{C}^{-1} \otimes \bold{C}^{-1} \right)
   *   + \delta_7 \mathbb{P}
   * \f]
   * \f[
   *  \delta_1 = 4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}^2}
   *           + 2 Ic \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial
   *           II_{\boldsymbol{C}}}
   *           + \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
   *           + I_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
   *           \right)
   * \f]
   * \f[
   *  \delta_2 = -4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial
   *  II_{\boldsymbol{C}}}
   *           + I_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
   * \f]
   * \f[
   *  \delta_3 = 4\left( III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}
   *  \partial III_{\boldsymbol{C}}}
   *           + I_{\boldsymbol{C}} III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial IIc
   *           \partial III_{\boldsymbol{C}}} \right)
   * \f]
   * \f[
   *  \delta_4 = 4\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
   * \f]
   * \f[
   *  \delta_5 = -4 III_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}
   *  \partial III_{\boldsymbol{C}}}
   * \f]
   * \f[
   *  \delta_6 = 4\left( III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
   *           + III_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial III_{\boldsymbol{C}}^2}
   *           \right)
   * \f]
   * \f[
   *  \delta_7 = -4 III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
   * \f]
   * \f[
   *  \delta_8 = -4 \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
   * \f]
   *
   * @param gamma (out) : Coefficients gamma, cf. Holzapfel [1]
   * @param delta (out) : Coefficients delta, cf. Holzapfel [1]
   * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
   * @param dPI (in) : First derivatives of the Free energy function w.r.t the principal
   * invariants of the Right Cauchy-Green strain tensor
   * @param ddPII (in) : Second derivatives of the Free energy function w.r.t the principal
   * invariants of the Right Cauchy-Green strain tensor
   */
  void CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma, LINALG::Matrix<8, 1>& delta,
      const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
      const LINALG::Matrix<6, 1>& ddPII);

  /**
   * \brief Evaluates the properties if the summands
   *
   * @param potsum List of the summands
   * @param properties Class holding the properties of the formulation of the summands
   */
  void ElastHyperProperties(const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      SummandProperties& properties);

  /**
   * \brief Check if material state is polyconvex
   *
   *  Polyconvexity of isotropic hyperelastic material
   *  dependent on principal or modified invariants)
   *  is tested with eq. (5.31) of Vera Ebbing - PHD-thesis (p. 79).
   *  Partial derivatives of SEF are used.
   *  \author birzle
   *  \date 04/2016
   *
   * @param defgrd (in) : Deformation gradient
   * @param prinv (in) : Principal invariants of the Right Cauchy-Green tensor
   * @param dPI (in) : First derivative of the Free-energy function with respect to the principal
   * invariants
   * @param ddPII (in) : Second derivative of the Free-energy function with respect to the
   * principal invariants
   * @param params (in/out) : Container for additional information
   * @param eleGID : Global element id
   * @param properties (in) : Class holding the properties of the formulation of the summands
   */
  void ElastHyperCheckPolyconvexity(const LINALG::Matrix<3, 3>& defgrd,
      const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
      const LINALG::Matrix<6, 1>& ddPII, Teuchos::ParameterList& params, int eleGID,
      const SummandProperties& properties);

  /**
   * \brief Class for holding the summand formulation properties
   */
  class SummandProperties
  {
   public:
    /// @name Flags to specify the elastic formulations (initialize with false)
    //@{
    ///< Indicator for isotropic principal formulation
    bool isoprinc = false;

    ///< Indicator for isotropic splitted formulation
    bool isomod = false;

    ///< Indicator for anisotropic principal formulation
    bool anisoprinc = false;

    ///< Indicator for anisotropic splitted formulation
    bool anisomod = false;

    ///< Indicator for coefficient stretches principal formulation
    bool coeffStretchesPrinc = false;

    ///< Indicator for coefficient stretches splitted formulation
    bool coeffStretchesMod = false;

    ///< Indicator for general viscoelastic behavior
    bool viscoGeneral = false;
    //@}

    /**
     * Pack all data to distribute it on other processors
     * @param data data where to store the values
     */
    void Pack(DRT::PackBuffer& data) const
    {
      DRT::ParObject::AddtoPack(data, isoprinc);
      DRT::ParObject::AddtoPack(data, isomod);
      DRT::ParObject::AddtoPack(data, anisoprinc);
      DRT::ParObject::AddtoPack(data, anisomod);
      DRT::ParObject::AddtoPack(data, coeffStretchesPrinc);
      DRT::ParObject::AddtoPack(data, coeffStretchesMod);
      DRT::ParObject::AddtoPack(data, viscoGeneral);
    }

    /**
     * Unpack all data received from another processor
     * @param position position where to start to unpack
     * @param data data holding the values
     */
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      isoprinc = (bool)DRT::ParObject::ExtractInt(position, data);
      isomod = (bool)DRT::ParObject::ExtractInt(position, data);
      anisoprinc = (bool)DRT::ParObject::ExtractInt(position, data);
      anisomod = (bool)DRT::ParObject::ExtractInt(position, data);
      coeffStretchesPrinc = (bool)DRT::ParObject::ExtractInt(position, data);
      coeffStretchesMod = (bool)DRT::ParObject::ExtractInt(position, data);
      viscoGeneral = (bool)DRT::ParObject::ExtractInt(position, data);
    }

    /**
     * Set all flags to false
     */
    void Clear()
    {
      isoprinc = false;
      isomod = false;
      anisoprinc = false;
      anisomod = false;
      coeffStretchesPrinc = false;
      coeffStretchesMod = false;
      viscoGeneral = false;
    }

    /**
     * This is an or operation: $this = $this || $other (for all items)
     * @param other (in) : Other Summand properties
     */
    void Merge(const SummandProperties& other)
    {
      isoprinc |= other.isoprinc;
      isomod |= other.isomod;
      anisoprinc |= other.anisoprinc;
      anisomod |= other.anisomod;
      coeffStretchesPrinc |= other.coeffStretchesPrinc;
      coeffStretchesMod |= other.coeffStretchesMod;
      viscoGeneral |= other.viscoGeneral;
    }

    /**
     * This is an or operation: $this = $other (for all items)
     * @param other (in) : Other Summand properties
     */
    void Update(const SummandProperties& other)
    {
      isoprinc = other.isoprinc;
      isomod = other.isomod;
      anisoprinc = other.anisoprinc;
      anisomod = other.anisomod;
      coeffStretchesPrinc = other.coeffStretchesPrinc;
      coeffStretchesMod = other.coeffStretchesMod;
      viscoGeneral = other.viscoGeneral;
    }
  };
}  // namespace MAT

#endif  // BACI_ELASTHYPER_EVALUATE_H
