/*----------------------------------------------------------------------*/
/*!
\file elast_summand.H

\brief
Interface class for materials of (visco)elasthyper toolbox.

\level 1

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELAST_SUMMAND_H
#define ELAST_SUMMAND_H



#include "../drt_inpar/inpar_material.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"
#include "../headers/FADmatrix_utils.H"

// forward declarations

namespace DRT
{
  class PackBuffer;

  namespace INPUT
  {
    class LineDefinition;
  }
}  // namespace DRT


namespace MAT
{
  /// This namespace contains hyperelastic summands
  /// and their derivatives with respect to invariants
  namespace ELASTIC
  {
    /// This namespace contains the input data -- parameters -- for
    /// hyperelastic summands
    namespace PAR
    {
      // !brief enum for mapping between material parameter and entry in the matparams_ vector
      // for all elasthyper summands
      enum matparamelastnames_
      {
        coupneohooke_c,
        coupneohooke_beta,
        coup1pow_c,
        coup1pow_d,
        first = coupneohooke_c,
        last = coup1pow_d
      };
    }  // namespace PAR

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
    /// <li> [2] C Sansour, "On the physical assumptions underlying the volumatric-isochoric split
    /// and the case of anisotropy", European Journal of Mechanics 2007
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

      virtual void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position)
      {
        return;
      };

      //@}

      /// provide material type
      virtual INPAR::MAT::MaterialType MaterialType() const = 0;

      /// Create summand object by input parameter ID
      static Teuchos::RCP<Summand> Factory(int matnum  ///< material ID
      );

      /// Dummy routine for setup of summand (needed, e.g., for HU-dependent material
      /// parameters)
      virtual void Setup(DRT::INPUT::LineDefinition* linedef) { return; };

      /// Dummy routine for setup of patient-specific materials
      virtual void SetupAAA(Teuchos::ParameterList& params, const int eleGID) { return; };

      /// Dummy routine for setup update of summand
      virtual void Update() { return; };

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy functions
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>&
              glstrain,     ///< Green-Lagrange strain in strain like Voigt notation
          const int eleGID  ///< element GID
      )
      {
        dserror("Summand does not support calculation of strain energy");
      };

      /// add shear modulus equivalent
      virtual void AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
          double& shearmod                          ///< variable to add upon
          ) const;

      /// add young's modulus equivalent
      virtual void AddYoungsMod(double& young, double& shear, double& bulk)
      {
        dserror("Summand does not support calculation of youngs modulus");
      };

      /// retrieve coefficients of first and second derivative
      /// of summand with respect to principal invariants
      ///
      /// The derivatives of the summand
      /// \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$ in which the
      /// principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$ are the
      /// arguments are defined as following:
      ///
      /// First derivatives:
      ///
      ///\f[
      /// dPI_0 = \frac{\partial \Psi}{\partial I_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dPI_1 = \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dPI_2 = \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}} ;
      ///\f]
      ///
      /// Second derivatives:
      ///\f[
      /// ddPII_0 = \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// ddPII_1 = \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// ddPII_2 = \frac{\partial^2 \Psi}{\partial III_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// ddPII_3 = \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}} \partial
      /// III_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// ddPII_4 = \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial
      /// III_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// ddPII_5 = \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial II_{\boldsymbol{C}}}
      /// ;
      ///\f]

      virtual void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const int eleGID                    ///< element GID
      )
      {
        return;  // do nothing
      };

      /// retrieve coefficients of third derivative
      /// of summand with respect to principal isotropic invariants
      ///
      /// The derivatives of the summand
      /// \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$ in which the
      /// principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$ are the
      /// arguments are defined as following:
      ///
      /// Third derivatives:
      ///
      ///\f[
      /// dddPIII_0 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}^3} ;
      ///\f]
      ///\f[
      /// dddPIII_1 = \frac{\partial^3 \Psi}{\partial II_{\boldsymbol{C}}^3} ;
      ///\f]
      ///\f[
      /// dddPIII_2 = \frac{\partial^3 \Psi}{\partial III_{\boldsymbol{C}}^3} ;
      ///\f]
      ///\f[
      /// dddPIII_3 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}
      /// \partial II_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// dddPIII_4 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}
      /// \partial III_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// dddPIII_5 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}^2
      /// \partial II_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dddPIII_6 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}^2
      /// \partial III_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dddPIII_7 = \frac{\partial^3 \Psi}{\partial II_{\boldsymbol{C}}^2
      /// \partial III_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dddPIII_8 = \frac{\partial^3 \Psi}{\partial III_{\boldsymbol{C}}^2
      /// \partial II_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dddPIII_9 = \frac{\partial^3 \Psi}{\partial I_{\boldsymbol{C}}
      /// \partial II_{\boldsymbol{C}} \partial III_{\boldsymbol{C}}} ;
      ///\f]
      virtual void AddThirdDerivativesPrincipalIso(
          LINALG::Matrix<10, 1>& dddPIII_iso,     ///< third derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv_iso,  ///< principal isotropic invariants
          const int eleGID)                       ///< element GID
      {
        dserror("not implemented");
        return;  // do nothing
      }

      /// retrieve coefficients of first and second derivative
      /// of summand (decoupled form) with respect to modified invariants
      ///
      /// The derivatives of the summand
      /// \f$\Psi(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}},J)=\Psi_{iso}(\overline{I}_{\boldsymbol{C}},\overline{II}_{\boldsymbol_{C}})+\Psi_{vol}(J)\f$
      /// in which the modified invariants of the modified right Cauchy-Green tensor
      /// \f$\overline{\boldsymbol{C}}=(\det\boldsymbol{C})^{-1/3} \boldsymbol{C}\f$ are the
      /// arguments, defined as following:
      ///
      /// First derivatives:
      ///
      ///\f[
      /// dmodPI_0 = \frac{\partial \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dmodPI_1 = \frac{\partial \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}} ;
      ///\f]
      ///\f[
      /// dmodPI_2 = \frac{\partial \Psi_{vol}}{\partial J} ;
      ///\f]
      ///
      /// Second derivatives:
      ///\f[
      /// ddmodPII_0 = \frac{\partial^2 \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// ddmodPII_1 = \frac{\partial^2 \Psi_{iso}}{\partial \overline{II}_{\boldsymbol{C}}^2} ;
      ///\f]
      ///\f[
      /// ddmodPII_2 = \frac{\partial^2 \Psi_{vol}}{\partial J^2} ;
      ///\f]
      ///\f[
      /// ddmodPII_3 = \frac{\partial^2 \Psi}{\partial \overline{II}_{\boldsymbol{C}} \partial J} =
      /// 0 ;
      ///\f]
      ///\f[
      /// ddmodPII_4 = \frac{\partial^2 \Psi}{\partial \overline{I}_{\boldsymbol{C}} \partial J} = 0
      /// ;
      ///\f]
      ///\f[
      /// ddmodPII_5 = \frac{\partial^2 \Psi_{iso}}{\partial \overline{I}_{\boldsymbol{C}} \partial
      /// \overline{II}_{\boldsymbol{C}}} ;
      ///\f]

      virtual void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddPmodII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const int eleGID                     ///< global ID of element
      )
      {
        return;  // do nothing
      };

      /// retrieve coefficients for the third derivative
      /// of volumetric summand with respect to modified invariants
      /// This is needed for TSI problems where
      ///\f[
      /// \hat{\mathbb{M}}(J,\Delta T)=-3\alpha_T\Delta T \frac{\partial
      /// \Psi_{\text{vol}}(J)}{\partial J}
      ///\f]
      /// and the PK2 stress
      ///\f[
      /// \mathbf{S}=-\frac{\partial \mathcal{M}}{\partial J}J\mathbf{C}^{-1} = -3\alpha_T\Delta T
      /// (\bar{\delta}_5-\bar{\gamma}_3)\mathbf{C}^{-1}
      ///\f]
      /// and the tangent
      ///\f[
      /// \mathbb{C} = 6\alpha_T\Delta T (\bar{\delta}_5-\bar{\gamma}_3)
      /// \mathbf{C}^{-1}\odot\mathbf{C}^{-1}
      ///             -\frac{3}{4}\alpha_T\Delta T  J^2  \frac{\partial^3 \Psi_{vol}}{\partial J^3}
      ///             \mathbf{C}^{-1}\otimes\mathbf{C}^{-1}
      ///             -\frac{3}{4}\alpha_T\Delta T \underbrace{J  \frac{\partial^2
      ///             \Psi_{vol}}{\partial J^2}}_{=\bar{\delta}_5-\bar{\gamma}_3}
      ///             \mathbf{C}^{-1}\otimes\mathbf{C}^{-1}
      ///\f]
      virtual void Add3rdVolDeriv(const LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
      {
        return;  // do nothing
      };

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation This is not to be called in "usual" calculation of stress/stiffness
      /// since the effects are already included in the derivatives w.r.t. to the invariants. For
      /// special applications, those derivatives may be interesting.
      virtual void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
      {
        return;  // do nothing
      }

      /// retrieve coefficients of first and second derivative
      /// of summand for fiber directions with respect to principal invariants
      ///
      /// The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
      /// on the summand
      /// \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}},IV_{\boldsymbol{C}},VI_{\boldsymbol{C}},VIII_{\boldsymbol{C}})\f$
      /// in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
      /// are the arguments,
      ///\f[
      ///  I_{\boldsymbol{C}} = \mathbf{C},
      ///  \quad
      ///  II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) -
      ///  \mathrm{trace}(\boldsymbol{C}^2) \big), \quad III_{\boldsymbol{C}} = \det(\boldsymbol{C})
      ///\f]
      ///\f[
      ///  IV_{\boldsymbol{C}} = \boldsymbol{C} : \boldsymbol{A}_1
      ///  \quad
      ///  VI_{\boldsymbol{C}} = \boldsymbol{C} : \boldsymbol{A}_2
      ///  \quad
      ///  VIII_{\boldsymbol{C}} = \boldsymbol{a}_1*\boldsymbol{C}\boldsymbol{a}_2
      ///\f]
      ///
      /// where \f$ \boldsymbol{a}_i \f$ is one fiber direction and \f$ \boldsymbol A_i= \boldsymbol
      /// a_i \otimes \boldsymbol a_i\f$, \f$ \boldsymbol{A}_1\boldsymbol{A}_2 = \left( \boldsymbol
      ///a_1 \otimes \boldsymbol a_2 + \boldsymbol a_2 \otimes \boldsymbol a_1 \right) \f$, cf.
      ///Holzapfel [1], p. 274
      ///
      ///\f[
      /// \bold{S} = \gamma^{aniso}_1 \ \boldsymbol{A}_1 + \gamma^{aniso}_2 \ \boldsymbol{A}_2 +
      /// \gamma^{aniso}_3 \ \boldsymbol{A}_1\boldsymbol{A}_2
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
      ///  \mathbb{C} = \delta^{aniso}_1 \left( \bold{A}_1 \otimes \bold{A}_1 \right) +
      ///  \delta^{aniso}_2 \left( \bold{A}_2 \otimes \bold{A}_2 \right)
      ///  + \delta^{aniso}_3 \left( \bold{Id} \otimes \bold{A}_1 + \bold{A}_1 \otimes \bold{Id}
      ///  \right) + \delta^{aniso}_4 \left( \bold{Id} \otimes \bold{A}_2 + \bold{A}_2 \otimes
      ///  \bold{Id} \right)
      ///  + \delta^{aniso}_5 \left( \bold{C} \otimes \bold{A}_1 + \bold{A}_1 \otimes \bold{C}
      ///  \right) + \delta^{aniso}_6 \left( \bold{C} \otimes \bold{A}_2 + \bold{A}_2 \otimes
      ///  \bold{C} \right)
      ///  + \delta^{aniso}_7 \left( \bold{C}^{-1} \otimes \bold{A}_1 + \bold{A}_1 \otimes
      ///  \bold{C}^{-1} \right)
      ///\f]
      ///\f[
      ///  + \delta^{aniso}_8 \left( \bold{C}^{-1} \otimes \bold{A}_2 + \bold{A}_2 \otimes
      ///  \bold{C}^{-1} \right)
      ///  + \delta^{aniso}_9 \left( \bold{A}_1 \otimes \bold{A}_2 + \bold{A}_2 \otimes \bold{A}_1
      ///  \right)
      ///  + \delta^{aniso}_{10} \left( \bold{Id} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1
      ///  \bold{A}_2 \otimes \bold{Id} \right)
      ///  + \delta^{aniso}_{11} \left( \bold{C} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1
      ///  \bold{A}_2 \otimes \bold{C} \right)
      ///  + \delta^{aniso}_{12} \left( \bold{C}^{-1} \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1
      ///  \bold{A}_2 \otimes \bold{C}^{-1} \right)
      ///  + \delta^{aniso}_{13} \left( \bold{A}_1 \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1
      ///  \bold{A}_2 \otimes \bold{A}_1 \right)
      ///  + \delta^{aniso}_{14} \left( \bold{A}_2 \otimes \bold{A}_1 \bold{A}_2 + \bold{A}_1
      ///  \bold{A}_2 \otimes \bold{A}_2 \right)
      ///  + \delta^{aniso}_{15} \left( \bold{A}_1 \bold{A}_2 \otimes \bold{A}_1 \bold{A}_2 \right)
      ///\f]
      ///\f[
      /// \delta^{aniso}_1 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_2 = 4\frac {\partial^2\Psi}{\partial VI_{\boldsymbol C}^2}
      ///\f]
      ///\f[
      /// \delta^{aniso}_3 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial
      /// IV_{\boldsymbol C}}
      ///                  + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      ///                  IV_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_4 = 4\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial
      /// VI_{\boldsymbol C}}
      ///                  + 4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      ///                  VI_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_5 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      /// IV_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_6 = -4\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      /// VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_7 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial
      /// VI_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_8 = 4\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial
      /// VI_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_9 = 4\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}\partial
      /// VI_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{10} = 2\frac {\partial^2\Psi}{\partial I_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}}
      ///                     + 2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      ///                     VIII_{\boldsymbol C}}I_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{11} = -2\frac {\partial^2\Psi}{\partial II_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{12} = 2\frac {\partial^2\Psi}{\partial III_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}} III_{\boldsymbol C}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{13} = 2\frac {\partial^2\Psi}{\partial IV_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{14} = 2\frac {\partial^2\Psi}{\partial VI_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}}
      ///\f]
      ///\f[
      /// \delta^{aniso}_{15} = \frac {\partial^2\Psi}{\partial VIII_{\boldsymbol C}\partial
      /// VIII_{\boldsymbol C}}
      ///\f]
      virtual void AddStressAnisoPrincipal(
          const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList& params, const int eleGID)
      {
        return;  // do nothing
      };

      virtual void AddCoefficientsViscoPrincipal(
          const LINALG::Matrix<3, 1>& inv,  ///< invariants of right Cauchy-Green tensor
          LINALG::Matrix<8, 1>& mu,         ///< see above
          LINALG::Matrix<33, 1>& xi,        ///< see above
          LINALG::Matrix<7, 1>& rateinv, Teuchos::ParameterList& params, const int eleGID)
      {
        return;  // do nothing
      };

      virtual void AddCoefficientsViscoModified(
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          LINALG::Matrix<8, 1>& modmu,         ///< see above
          LINALG::Matrix<33, 1>& modxi,        ///< see above
          LINALG::Matrix<7, 1>& modrateinv, Teuchos::ParameterList& params, const int eleGID)
      {
        return;  // do nothing
      };

      /// Read material parameters of viscogenmax or viscofract
      virtual void ReadMaterialParametersVisco(double& tau,  ///< relaxation parameter tau
          double& beta,   ///< emphasis of viscous to elastic part
          double& alpha,  ///< fractional order derivative (just for visoc_fract)
          std::string&
              solve  //!< solution variant for time evolution of viscous stress (just for genmax)
      )
      {
        return;  // do nothing
      };

      //! GeneralizedGenMax
      virtual void ReadMaterialParameters(int& numbranch,  //!< number of visco branches
          const std::vector<int>*& matids,                 //!< material ids of visco branches
          std::string& solve  //!< solution variant for time evolution of viscous stress
      )
      {
        return;  // not implemented in base class. May be overridden in subclass.
      };

      //! GeneralizedGenMax
      virtual void ReadMaterialParameters(double& nummat,  //!< number of visco branches
          const std::vector<int>*& matids                  //!< material ids of visco branches
      )
      {
        return;  // not implemented in base class. May be overridden in subclass.
      };

      /// retrieve stress and cmat
      /// of summand for fiber directions with respect to modified strains
      virtual void AddStressAnisoModified(
          const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          const LINALG::Matrix<6, 1>& icg,  ///< inverse of right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          double I3,                        ///< third principal invariant
          const int eleGID                  ///< element GID
      )
      {
        return;  // do nothing
      };

      /// Answer if coefficients with respect to principal stretches are provided
      virtual bool HaveCoefficientsStretchesPrincipal() { return false; }

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
          LINALG::Matrix<3, 1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
          LINALG::Matrix<6, 1>&
              delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
          const LINALG::Matrix<3, 1>&
              prstr  ///< principal stretches, [lambda_1, lambda_2, lambda_3]
      )
      {
        return;  // do nothing
      }

      /// Answer if coefficients with respect to modified principal stretches are provided
      virtual bool HaveCoefficientsStretchesModified() { return false; }

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
          LINALG::Matrix<3, 1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
          LINALG::Matrix<6, 1>&
              delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
          const LINALG::Matrix<3, 1>&
              modstr  ///< modified principal stretches, [lambda_1, lambda_2, lambda_3]
      )
      {
        return;  // do nothing
      }

      /// Set fiber directions
      virtual void SetFiberVecs(const double newgamma,  ///< new angle
          const LINALG::Matrix<3, 3>& locsys,           ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd            ///< deformation gradient
      )
      {
        return;  // do nothing
      };

      /// Set fiber directions
      virtual void SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec  ///< new fiber vector
      )
      {
        dserror("Not implemented yet for this type of anisotropic material;");
      };

      /// Get fiber directions
      virtual void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
      )
      {
        return;  // do nothing
      };

      /// Read FIBERn
      void ReadFiber(DRT::INPUT::LineDefinition* linedef, std::string kenner,
          LINALG::Matrix<3, 1>& fiber_vector);

      /// Read RAD-AXI-CIR
      void ReadRadAxiCir(DRT::INPUT::LineDefinition* linedef, LINALG::Matrix<3, 3>& locsys);

      /// Indicator for the chosen formulations
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) = 0;

      /// Indicator for the chosen viscoelastic formulations
      virtual void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according to the
                              ///< SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
      ){/* do nothing for non viscoelastic material models */};

      /// Copy matparams_ to summands for stat inverse analysis
      virtual void CopyStatInvAnaMatParams(std::vector<Teuchos::RCP<Epetra_Vector>> input){};

      /// Add parameters of elasthyper-summands for stat inverse analysis to matparams
      virtual void SetStatInvAnaSummandMatParams(){};

      /// Add paramters to optimize from summands for stat inverse analysis
      virtual void AddElastOptParams(std::map<std::string, int>* pnames)
      {
        dserror(
            "This summand does not support the stat inverse analysis. \n"
            "Have a look in CoupNeoHooke, how to add an additional summand.");
      };
      //@}

      //! @name Visualization methods

      /// Return names of visualization data
      virtual void VisNames(std::map<std::string, int>& names){
          /* do nothing for simple material models */};

      /// Return visualization data
      virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleId)
      { /* do nothing for simple material models */
        return false;
      };

      //@}

    };  // class Summand

  }  // namespace ELASTIC

}  // namespace MAT

#endif
