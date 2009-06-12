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

#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"

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
    /// </ul>
    ///
    /// \author rausch,tk,bborn
    /// \date 05/09
    class Summand
    {
    public:

      /// standard constructor
      Summand() { ; }

      /// destructor
      virtual ~Summand() { ; }

      /// provide material type
      virtual INPAR::MAT::MaterialType MaterialType() const = 0;

      /// Create summand object by input parameter ID
      static Teuchos::RCP<Summand> Factory(
        int matnum  ///< material ID
        );

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
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<8,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
        ) = 0;

      /// retrieve coefficients of first and second derivative
      /// of summand (decoupled form) with respect to modified invariants
      ///
      /// The coefficients \f$\bar\gamma_i\f$ and \f$\bar\delta_j\f$ are based
      /// on the summand \f$\Psi(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}},J)=\Psi_{iso}(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}})+\Psi_{vol}(J)\f$
      /// in which the principal invariants of the modified right Cauchy-Green tensor \f$\overline{\boldsymbol{C}}=(\det\boldsymbol{C})^{-1/3} \boldsymbol{C}\f$
      /// are the arguments, cf. Holzapfel [1], 
      ///\f[
      ///  \overline{I}_{\boldsymbol{C}} = J^{-2/3} C_{AA},
      ///  \quad
      ///  \overline{II}_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) - 
      ///         \mathrm{trace}(\boldsymbol{C}^2) \big),
      ///  \quad
      ///  J = \sqrt{\det(\boldsymbol{C})}
      ///\f]
      /// 
      /// cf. Holzapfel [1], p. 234 and p. 262
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
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<5,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
        ) = 0;

      /// retrieve coefficients of first and second derivative
      /// of summand for fiber directions with respect to principal invariants
      ///
      /// The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
      /// on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}},IV_{\boldsymbol{C}},VI_{\boldsymbol{C}},VIII_{\boldsymbol{C}})\f$
      /// in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
      /// are the arguments,  
      ///\f[
      ///  I_{\boldsymbol{C}} = C_{AA},
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
      ///where \f$ \boldsymbol{a}_i \f$ is one fiber direction and \f$ \boldsymbol A_i= \boldsymbol a_i \otimes \boldsymbol a_i\f$.
      ///  
      ///
      /// cf. Holzapfel [1], p. 274
      ///\f[
      /// \gamma^{aniso}_1 = 2\frac{\partial \Psi}{\partial IV_{\boldsymbol{C}}} 
      ///\f]
      ///\f[
      /// \gamma^{aniso}_2 = 2\frac{\partial \Psi}{\partial VI_{\boldsymbol{C}}};
      ///\f]
      ///\f[
      /// \gamma^{aniso}_3 = 2\frac{\partial \Psi}{\partial VIII_{\boldsymbol{C}}};
      ///\f]
      ///
      /// material constitutive tensor coefficients
      /// cf. Holzapfel [1], p. 261
      ///\f[
      /// \delta^{aniso}_1 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_2 = 4\frac {\partial^2\Psi}{\partial VI_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_3 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial IV_{\boldsymbol C}} + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial IV_{\boldsymbol C}}I_\boldsymbol C
      ///\f]
      ///\f[
      /// \delta^{aniso}_4 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial VI_{\boldsymbol C}} + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VI_{\boldsymbol C}}I_\boldsymbol C
      ///\f]
      ///\f[
      /// \delta^{aniso}_5 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial IV_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_6 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_7 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VI_{\boldsymbol C}} III_\boldsymbol C
      ///\f]
      ///\f[
      /// \delta^{aniso}_8 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VI_{\boldsymbol C}} III_\boldsymbol C
      ///\f]      
      ///\f[
      /// \delta^{aniso}_9 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}\partial VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{10} = 2\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial VIII_{\boldsymbol C}} + 2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VIII_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{11} = -2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{12} = 2\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial VIII_{\boldsymbol C}} III_\boldsymbol C
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
      virtual void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<15,1>& delta,  ///< see above
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor         
      ) = 0;
      
      /// To be done and commented by Moritz :-)      
      virtual void AddCoefficientsModifiedAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<15,1>& delta,  ///< see above
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor         
      ) = 0;
      
    private:
      

    };  // class Summand

  }  // namespace ELASTIC

}

#endif

#endif
