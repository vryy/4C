/*----------------------------------------------------------------------*/
/*!
\file elast_summand.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELAST_SUMMAND_H
#define ELAST_SUMMAND_H



#include "../drt_inpar/inpar_material.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"

// forward declarations

namespace DRT
{
  class PackBuffer;

  namespace INPUT
  {
    class LineDefinition;
  }
}


namespace MAT
{

  /// This namespace contains hyperelastic summands
  /// and their derivatives with respect to invariants
  namespace ELASTIC
  {

    /// This namespace contains the input data -- parameters -- for
    /// hyperelastic summands
    namespace PAR { }

    /*----------------------------------------------------------------------*/
    /// Interface for hyperelastic potentials
    ///
    /// The interface defines the way how MAT::ElastHyper can access
    /// coefficients to build the actual stress response and
    /// elasticity tensor.
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
    /// <li> [2] C Sansour, "On the physical assumptions underlying the volumatric-isochoric split and the case of anisotropy", European Journal of Mechanics 2007
    /// </ul>
    ///
    /// \author rausch,tk,bborn
    /// \date 05/09
    class Summand : public DRT::ParObject
    {
    public:

      /// standard constructor
      Summand() { ; }

      /// destructor
      virtual ~Summand() { ; }

      ///@name Packing and Unpacking (dummy routines)
      //@{

      virtual int UniqueParObjectId() const;

      virtual void Pack(DRT::PackBuffer& data) const;

      virtual void Unpack(const std::vector<char>& data);

      virtual void PackSummand(DRT::PackBuffer& data) const { return; };

      virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position) { return; };

      //@}

      /// provide material type
      virtual INPAR::MAT::MaterialType MaterialType() const = 0;

      /// Create summand object by input parameter ID
      static Teuchos::RCP<Summand> Factory(
        int matnum  ///< material ID
        );

      /// Dummy routine for setup of summand (needed, e.g., for HU-dependent material
      /// parameters)
      virtual void Setup(DRT::INPUT::LineDefinition* linedef) { return; };

      /// Dummy routine for setup of patient-specific materials
      virtual void SetupAAA(Teuchos::ParameterList& params) { return; };

      /// Dummy routine for setup update of summand
      virtual void Update() { return; };

      // add strain energy
      virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy functions
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv ///< modified invariants of right Cauchy-Green tensor
      )
      { dserror("Summand does not support calculation of strain energy"); };

      /// add shear modulus equivalent
      virtual void AddShearMod(
        bool& haveshearmod,  ///< non-zero shear modulus was added
        double& shearmod  ///< variable to add upon
        ) const;

      /// retrieve coefficients of first and second derivative
      /// of summand with respect to principal invariants
      ///
      /// The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
      /// on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$
      /// in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
      /// are the arguments, cf. Holzapfel [1],
      ///\f[
      ///  I_{\boldsymbol{C}} = C_{AA},
      ///  \quad
      ///  II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) - \mathrm{trace}(\boldsymbol{C}^2) \big),
      ///  \quad
      ///  III_{\boldsymbol{C}} = \det(\boldsymbol{C})
      ///\f]
      ///
      /// cf. Holzapfel [1], p. 216 (6.32) and p. 248
      ///\f[
      ///  \bold{S} = \gamma_1 \ \bold{Id} + \gamma_2 \ \bold{C} + \gamma_3 \ \bold{C}^{-1}
      ///\f]
      ///\f[
      /// \gamma_1 = 2\left( \frac{\partial \Psi}{\partial I_{\boldsymbol{C}}}
      ///          + I_{\boldsymbol{C}}\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}} \right);
      ///\f]
      ///\f[
      /// \gamma_2 = -2\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}};
      ///\f]
      ///\f[
      /// \gamma_3 = 2III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}};
      ///\f]
      ///
      /// material constitutive tensor coefficients
      /// cf. Holzapfel [1], p. 261
      ///\f[
      ///  \mathbb{C} = \delta_1 \left( \bold{Id} \otimes \bold{Id} \right) + \delta_2 \left( \bold{Id} \otimes \bold{C} + \bold{C} \otimes \bold{Id} \right)
      ///  + \delta_3 \left( \bold{Id} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{Id} \right) + \delta_4 \left( \bold{C} \otimes \bold{C} \right)
      ///  + \delta_5 \left( \bold{C} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{C} \right) + \delta_6 \left( \bold{C}^{-1} \otimes \bold{C}^{-1} \right)
      ///  + \delta_7 \mathbb{P}
      ///\f]
      ///\f[
      /// \delta_1 = 4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}^2}
      ///          + 2 Ic \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial II_{\boldsymbol{C}}}
      ///          + \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
      ///          + I_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
      ///\f]
      ///\f[
      /// \delta_2 = -4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial II_{\boldsymbol{C}}}
      ///          + I_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
      ///\f]
      ///\f[
      /// \delta_3 = 4\left( III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial III_{\boldsymbol{C}}}
      ///          + I_{\boldsymbol{C}} III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial IIc \partial III_{\boldsymbol{C}}} \right)
      ///\f]
      ///\f[
      /// \delta_4 = 4\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
      ///\f]
      ///\f[
      /// \delta_5 = -4 III_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}} \partial III_{\boldsymbol{C}}}
      ///\f]
      ///\f[
      /// \delta_6 = 4\left( III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
      ///          + III_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial III_{\boldsymbol{C}}^2} \right)
      ///\f]
      ///\f[
      /// \delta_7 = -4 III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
      ///\f]
      ///\f[
      /// \delta_8 = -4 \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
      ///\f]
      virtual void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<8,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
        )
      {
        return;  // do nothing
      };

      /// retrieve coefficients of first and second derivative
      /// of summand (decoupled form) with respect to modified invariants
      ///
      /// The coefficients \f$\bar\gamma_i\f$ and \f$\bar\delta_j\f$ are based
      /// on the summand \f$\Psi(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}},J)=\Psi_{iso}(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}})+\Psi_{vol}(J)\f$
      /// in which the principal invariants of the modified right Cauchy-Green tensor \f$\overline{\boldsymbol{C}}=(\det\boldsymbol{C})^{-1/3} \boldsymbol{C}\f$
      /// are the arguments, cf. Holzapfel [1],
      ///\f[
      ///  \overline{I}_{\boldsymbol{C}} = J^{-2/3} \bold{C},
      ///  \quad
      ///  \overline{II}_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) -
      ///         \mathrm{trace}(\boldsymbol{C}^2) \big),
      ///  \quad
      ///  J = \sqrt{\det(\boldsymbol{C})}
      ///\f]
      ///
      /// cf. Holzapfel [1], p. 234 and p. 262
      ///
      ///\f[
      /// \bold{S} = \bold{S}_{iso} + \bold{S}_{vol} = J^{-\frac{2}{3}} \left( \mathbb{P} : \bar{\bold{S}} \right) + J p \bold{C}^{-1}
      ///\f]
      ///\f[
      ///  \bar{\bold{S}} = \bar{\gamma}_1 \ \bold{Id} + \bar{\gamma}_2 \ \bold{C}
      ///\f]
      ///\f[
      /// \bar\gamma_1 = 2\left( \frac{\partial \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}}}
      ///          + \overline{I}_{\boldsymbol{C}} \frac{\partial \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}} \right);
      ///\f]
      ///\f[
      /// \bar\gamma_2 = -2 \frac{\partial \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}};
      ///\f]
      ///\f[
      /// \bar\gamma_3 = p = \frac{\partial \Psi_{vol}}{\partial J};
      ///\f]
      ///
      /// material constitutive tensor coefficients
      /// cf. Holzapfel [1], p. 261
      ///\f[
      /// \mathbb{C} = \mathbb{C}_{iso} + \mathbb{C}_{vol} = \left[ \mathbb{P}: \bar{\mathbb{C}} : \mathbb{P}^T + \frac{2}{3} Tr(J^{-\frac{2}{3} } \bar{\bold{S}} ) \tilde{\mathbb{P}} - \frac{2}{3} \left( \bold{C}^{-1} \otimes J^{-\frac{2}{3}} \left( \mathbb{P} : \bar{\bold{S}} \right) + J^{-\frac{2}{3}} \left( \mathbb{P} : \bar{\bold{S}} \right) \otimes \bold{C}^{-1} \right) \right] + \left[ J \tilde{p} \bold{C}^{-1} \otimes \bold{C}^{-1} - 2 p J \bold{C}^{-1} \odot \bold{C}^{-1} \right]
      ///\f]
      ///\f[
      ///  \bar{\mathbb{C}} = \bar\delta_1 \left( \bold{Id} \otimes \bold{Id} \right) + \bar\delta_2 \left( \bold{Id} \otimes \bold{C} + \bold{C} \otimes \bold{Id} \right)
      ///  + \bar\delta_3 \left( \bold{C} \otimes \bold{C} \right) + \bar\delta_4 \mathbb{P}
      ///\f]
      ///\f[
      /// \bar\delta_1 = 4 \left( \frac{\partial^2 \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}}^2}
      ///          + 2 \overline{I}_{\boldsymbol{C}} \frac{\partial^2 \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}} \partial \overline{II}_{\boldsymbol{C}}}
      ///          + \frac{\partial \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}}
      ///          + \overline{I}_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi_{sio}}{\partial \overline{II}_{\boldsymbol{C}}^2} \right)
      ///\f]
      ///\f[
      /// \bar\delta_2 = -4 \left( \frac{\partial^2 \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}} \partial \overline{II}_{\boldsymbol{C}}}
      ///          + \overline{I}_{\boldsymbol{C}} \frac{\partial^2 \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}^2} \right)
      ///\f]
      ///\f[
      /// \bar\delta_3 = 4 \frac{\partial^2 \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}^2}
      ///\f]
      ///\f[
      /// \bar\delta_4 = -4 \frac{\partial \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}}
      ///\f]
      ///\f[
      /// \bar\delta_5 = \tilde{p} = J \frac{\partial^2 \Psi_{vol}}{\partial J^2} + p
      ///\f]
      virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<5,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
        )
      {
        return;  // do nothing
      };

      /// retrieve coefficients of first and second derivative
      /// of summand for fiber directions with respect to principal invariants
      ///
      /// The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
      /// on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}},IV_{\boldsymbol{C}},VI_{\boldsymbol{C}},VIII_{\boldsymbol{C}})\f$
      /// in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
      /// are the arguments,
      ///\f[
      ///  I_{\boldsymbol{C}} = \mathbf{C},
      ///  \quad
      ///  II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) - \mathrm{trace}(\boldsymbol{C}^2) \big),
      ///  \quad
      ///  III_{\boldsymbol{C}} = \det(\boldsymbol{C})
      ///\f]
      ///\f[
      ///  IV_{\boldsymbol{C}} = \boldsymbol{C} : \boldsymbol{A}_1
      ///  \quad
      ///  VI_{\boldsymbol{C}} = \boldsymbol{C} : \boldsymbol{A}_2
      ///  \quad
      ///  VIII_{\boldsymbol{C}} = \boldsymbol{a}_1*\boldsymbol{C}\boldsymbol{a}_2
      ///\f]
      ///
      ///where \f$ \boldsymbol{a}_i \f$ is one fiber direction and \f$ \boldsymbol A_i= \boldsymbol a_i \otimes \boldsymbol a_i\f$,
      ///\f$ \boldsymbol{A}_1\boldsymbol{A}_2 = \left( \boldsymbol a_1 \otimes \boldsymbol a_2 + \boldsymbol a_2 \otimes \boldsymbol a_1 \right) \f$, cf. Holzapfel [1], p. 274
      ///
      ///\f[
      /// \bold{S} = \gamma^{aniso}_1 \ \boldsymbol{A}_1 + \gamma^{aniso}_2 \ \boldsymbol{A}_2 + \gamma^{aniso}_3 \ \boldsymbol{A}_1\boldsymbol{A}_2
      ///\f]
      ///\f[
      /// \gamma^{aniso}_1 = 2\frac{\partial \Psi}{\partial IV_{\boldsymbol{C}}}
      ///\f]
      ///\f[
      /// \gamma^{aniso}_2 = 2\frac{\partial \Psi}{\partial VI_{\boldsymbol{C}}};
      ///\f]
      ///\f[
      /// \gamma^{aniso}_3 = \frac{\partial \Psi}{\partial VIII_{\boldsymbol{C}}};
      ///\f]
      ///
      /// material constitutive tensor coefficients
      /// cf. Holzapfel [1], p. 261
      ///\f[
      ///  \mathbb{C} = \delta^{aniso}_1 \left( \bold{A}_1 \otimes \bold{A}_1 \right) + \delta^{aniso}_2 \left( \bold{A}_2 \otimes \bold{A}_2 \right)
      ///  + \delta^{aniso}_3 \left( \bold{Id} \otimes \bold{A}_1 + \bold{A}_1 \otimes \bold{Id} \right) + \delta^{aniso}_4 \left( \bold{Id} \otimes \bold{A}_2 + \bold{A}_2 \otimes \bold{Id} \right)
      ///  + \delta^{aniso}_5 \left( \bold{C} \otimes \bold{A}_1 + \bold{A}_1 \otimes \bold{C} \right) + \delta^{aniso}_6 \left( \bold{C} \otimes \bold{A}_2 + \bold{A}_2 \otimes \bold{C} \right)
      ///  + \delta^{aniso}_7 \left( \bold{C}^{-1} \otimes \bold{A}_1 + \bold{A}_1 \otimes \bold{C}^{-1} \right)
      ///\f]
      ///\f[
      ///  + \delta^{aniso}_8 \left( \bold{C}^{-1} \otimes \bold{A}_2 + \bold{A}_2 \otimes \bold{C}^{-1} \right)
      ///  + \delta^{aniso}_9 \left( \bold{A}_1 \otimes \bold{A}_2 + \bold{A}_2 \otimes \bold{A}_1 \right)
      ///  + \delta^{aniso}_{10} \left( \bold{Id} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1 \bold{A}_2 \otimes \bold{Id} \right)
      ///  + \delta^{aniso}_{11} \left( \bold{C} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1 \bold{A}_2 \otimes \bold{C} \right)
      ///  + \delta^{aniso}_{12} \left( \bold{C}^{-1} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1 \bold{A}_2 \otimes \bold{C}^{-1} \right)
      ///  + \delta^{aniso}_{13} \left( \bold{A}_1 \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1 \bold{A}_2 \otimes \bold{A}_1 \right)
      ///  + \delta^{aniso}_{14} \left( \bold{A}_2 \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1 \bold{A}_2 \otimes \bold{A}_2 \right)
      ///  + \delta^{aniso}_{15} \left( \bold{A}_1 \bold{A}_2 \otimes \bold{A}_1 \bold{A}_2 \right)
      ///\f]
      ///\f[
      /// \delta^{aniso}_1 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_2 = 4\frac {\partial^2\Psi}{\partial VI_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_3 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial IV_{\boldsymbol C}}
      ///                  + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial IV_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_4 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial VI_{\boldsymbol C}}
      ///                  + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VI_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_5 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial IV_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_6 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_7 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VI_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_8 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VI_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_9 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}\partial VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{10} = 2\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///                     + 2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VIII_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{11} = -2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{12} = 2\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VIII_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{13} = 2\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{14} = 2\frac {\partial^2\Psi}{\partial VI_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{15} = \frac {\partial^2\Psi}{\partial VIII_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///\f]
      virtual void AddStressAnisoPrincipal(
        const LINALG::Matrix<6,1> rcg,  ///<right Cauchy Green Tensor
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        Teuchos::ParameterList& params
        )
      {
        return;  // do nothing
      };

      /// retrieve stress and cmat
      /// of summand for fiber directions with respect to modified strains
      virtual void AddStressAnisoModified(
        const LINALG::Matrix<6,1> rcg,  ///<right Cauchy Green Tensor
        const LINALG::Matrix<6,1> icg,  ///< inverse of right Cauchy Green Tensor
        LINALG::Matrix<6,6>& cmat,  ///< material stiffness matrix
        LINALG::Matrix<6,1>& stress,  ///< 2nd PK-stress
        double I3   ///< third principal invariant
        )
      {
        return;  // do nothing
      };

      /// Answer if coefficients with respect to principal stretches are provided
      virtual bool HaveCoefficientsStretchesPrincipal()
      {
        return false;
      }

      /// Add coefficients with respect to principal stretches (or zeros)
      ///
      /// The coefficients \f$\gamma_\alpha\f$ are based on  \f$\alpha=1,2,3\f$
      ///\f[
      /// \gamma_\alpha = \frac{\partial \Psi}{\partial \lambda_\alpha}
      ///\f]
      ///
      /// The coefficients \f$\delta_{\alpha\beta}\f$ are based on \f$\alpha,\beta=1,2,3\f$
      ///\f[
      /// \delta_{\alpha\beta} = \frac{\partial^2\Psi}{\partial\lambda_\alpha \partial\lambda_\beta}
      ///\f]
      /// <i>Note:</i> These parameters have \e nothing in common with Kronecker's delta.
      ///
      /// \date 08/09
      /// \author bborn
      virtual void AddCoefficientsStretchesPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
        LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
        const LINALG::Matrix<3,1>& prstr  ///< principal stretches, [lambda_1, lambda_2, lambda_3]
      )
      {
        return;  // do nothing
      }

      /// Answer if coefficients with respect to modified principal stretches are provided
      virtual bool HaveCoefficientsStretchesModified()
      {
        return false;
      }

      /// Add coefficients with respect to modified principal stretches (or zeros)
      ///
      /// The coefficients \f$\bar{\gamma}_\alpha\f$ are based on \f$\alpha=1,2,3\f$
      ///\f[
      /// \bar{\gamma}_\alpha
      /// = \frac{\partial \Psi}{\partial \bar{\lambda}_\alpha}
      ///\f]
      /// and \f$\bar{\lambda}_\alpha = J^{-1/3} \lambda_\alpha\f$.
      ///
      /// The coefficients \f$\bar{\delta}_{\alpha\beta}\f$ are based on \f$\alpha,\beta=1,2,3\f$
      ///\f[
      /// \bar{\delta}_{\alpha\beta}
      /// = \frac{\partial^2\Psi}{\partial\bar{\lambda}_\alpha \partial\bar{\lambda}_\beta}
      ///\f]
      ///
      /// \date 08/09
      /// \author bborn
      virtual void AddCoefficientsStretchesModified(
        LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
        LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
        const LINALG::Matrix<3,1>& modstr  ///< modified principal stretches, [lambda_1, lambda_2, lambda_3]
      )
      {
        return;  // do nothing
      }

      /// Set fiber directions
      virtual void SetFiberVecs(
          const double newgamma,   ///< new angle
          const LINALG::Matrix<3,3> locsys,   ///< local coordinate system
          const LINALG::Matrix<3,3> defgrd   ///< deformation gradient
      )
      {
        return;  // do nothing
      };

      /// Get fiber directions
      virtual void GetFiberVecs(
          std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
      )
      {
        return;  // do nothing
      };

      /// Indicator for the chosen formulations
      virtual void SpecifyFormulation(
          bool& isoprinc_, ///< global indicator for anisotropy
          bool& isomod_, ///< global indicator for anisotropy
          bool& anisoprinc_, ///< global indicator for anisotropy
          bool& anisomod_ ///< global indicator for anisotropy
        ) = 0;

    };  // class Summand

  }  // namespace ELASTIC

}

#endif

