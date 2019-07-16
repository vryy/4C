/*----------------------------------------------------------------------*/
/*!
\brief evaluation of a generic material whose deformation gradient is modeled to be split
multiplicatively into elastic and inelastic parts

\level 3

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#ifndef MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_H
#define MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_H

#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*----------------------------------------------------------------------*/
/*! \class MultiplicativeSplitDefgrad_ElastHyper

    In this class the deformation gradient is modeled to be split multiplicatively
    in elastic and inelastic deformation gradients (F = F_{el} * F_{in}).
    The elastic contribution can be any elastic material law provided in elast_summand.
    Only elastic strains cause stresses!
    The inelastic deformation gradient itself can be a product of different inelastic
    deformation gradients, i.e. F_{in} = F_{in,1} * F_{in,2} * ... * F_{in,n}.
    The inverse inelastic deformation gradients and their derivative w.r.t. the primary variables
    that are needed to set up the system to be solved are evaluated in the derived classes
    of the interface class 'InelasticDefgradFactors'.
*/

// forward declaration
namespace MAT
{
  class InelasticDefgradFactors;

  namespace ELASTIC
  {
    class Summand;
  }

  namespace PAR
  {
    class MultiplicativeSplitDefgrad_ElastHyper : public Parameter
    {
     public:
      /// standard constructor
      MultiplicativeSplitDefgrad_ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~MultiplicativeSplitDefgrad_ElastHyper() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// length of elastic material list
      const int nummat_elast_;

      /// the list of elastic material IDs
      const std::vector<int>* matids_elast_;

      /// number of factors of inelastic deformation gradient F_{in} = F_{in,1} . F_{in,2}. ... .
      /// F_{in,n} (n factors)
      const int numfac_inel_;

      /// IDs of inelastic deformation gradient factors (i-th ID specifies calculation of F_{in,i})
      const std::vector<int>* inel_defgradfacids_;

      /// material mass density
      const double density_;

    };  // class MultiplicativeSplitDefgrad_ElastHyper

  }  // namespace PAR

  class MultiplicativeSplitDefgrad_ElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MultiplicativeSplitDefgrad_ElastHyperType"; }

    static MultiplicativeSplitDefgrad_ElastHyperType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MultiplicativeSplitDefgrad_ElastHyperType instance_;
  };  // class MultiplicativeSplitDefgrad_ElastHyperType


  /*----------------------------------------------------------------------*/
  class Material;

  class MultiplicativeSplitDefgrad_ElastHyper : public So3Material
  {
   public:
    /// construct empty material object
    MultiplicativeSplitDefgrad_ElastHyper();

    /// construct the material object given material parameters
    explicit MultiplicativeSplitDefgrad_ElastHyper(
        MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper* params);

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return MultiplicativeSplitDefgrad_ElastHyperType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    virtual void Unpack(const std::vector<char>& data);

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_multiplicative_split_defgrad_elasthyper;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new MultiplicativeSplitDefgrad_ElastHyper(*this));
    }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// evaluate stress response and elasticity tensor
    void Evaluate(const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,           ///< Green-Lagrange strain
        Teuchos::ParameterList& params,                 ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,                   ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                     ///< Constitutive matrix
        const int eleGID) final;                        ///< Element ID

    /// setup
    void Setup(const int numgp, DRT::INPUT::LineDefinition* linedef) final;

    /// update
    void Update() final;

    /*! \brief Calculate coefficients gamma and delta from partial derivatives w.r.t. invariants.
     *
     *  The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based on the the
     *  strain energy function \f$\Psi(I_{\bold{C}},II_{\bold{C}},III_{\bold{C}})\f$
     *  in which the principal invariants of the right Cauchy-Green tensor \f$\bold{C}\f$
     *  are the arguments, cf. Holzapfel [1].
     *  Here, the formulation is slightly adapted since only the elastic strains \f$\bold{C}_e\f$
     *  cause stresses, thus the invariants are also only composed of the elastic strains
     *  \f$\bold{C}_e\f$ The strain energy function is thus named \f$\Psi_e\f$
     *
     * \f[
     *   I_{\bold{C}_e} = {C_e}_{AA},
     *   \quad
     *   II_{\bold{C}_e} = 1/2 \left( \mathrm{trace}^2(\bold{C}_e) - \mathrm{trace}({\bold{C}_e}^2)
     *   \right), \quad III_{\bold{C_e}} = \det(\bold{C}_e)
     * \f]
     *
     * \f[
     *  \gamma_1 = 2 \left( \frac{\partial \Psi_e}{\partial I_{\bold{C}_e}}
     *           + I_{\bold{C}_e}\frac{\partial \Psi_e}{\partial II_{\bold{C}_e}} \right);
     * \f]
     * \f[
     *  \gamma_2 = -2 \frac{\partial \Psi_e}{\partial II_{\bold{C}_e}};
     * \f]
     * \f[
     *  \gamma_3 = 2 III_{\bold{C}_e} \frac{\partial \Psi_e}{\partial III_{\bold{C}_e}};
     * \f]
     *
     *  material constitutive tensor coefficients
     *  cf. Holzapfel [1], p. 261
     *  But formulation is slightly adapted as above since only the elastic strains \f$\bold{C}_e\f$
     *  cause stresses, thus the invariants are also only composed of the elastic strains
     *  \f$\bold{C}_e\f$
     * \f[
     *  \delta_1 = 4 \left( \frac{\partial^2 \Psi_e}{\partial I_{\bold{C}_e}^2}
     *           + 2 I_{\bold{C}_e} \frac{\partial^2 \Psi_e}{\partial I_{\bold{C}_e} \partial
     *           II_{\bold{C}_e}}
     *           + \frac{\partial \Psi_e}{\partial II_{\bold{C}_e}}
     *           + I_{\bold{C}_e}^2 \frac{\partial^2 \Psi_e}{\partial II_{\bold{C}_e}^2} \right)
     * \f]
     * \f[
     *  \delta_2 = -4 \left( \frac{\partial^2 \Psi_e}{\partial I_{\bold{C}_e} \partial
     *  II_{\bold{C}_e}}
     *           + I_{\bold{C}_e}\frac{\partial^2 \Psi_e}{\partial II_{\bold{C}_e}^2} \right)
     * \f]
     * \f[
     *  \delta_3 = 4 \left( III_{\bold{C}_e} \frac{\partial^2 \Psi_e}{\partial I_{\bold{C}_e}
     *  \partial III_{\bold{C}_e}}
     *           + I_{\bold{C}_e} III_{\bold{C}_e} \frac{\partial^2 \Psi_e}{\partial II_{\bold{C}_e}
     *           \partial III_{\bold{C}_e}} \right)
     * \f]
     * \f[
     *  \delta_4 = 4 \frac{\partial^2 \Psi_e}{\partial II_{\bold{C}_e}^2}
     * \f]
     * \f[
     *  \delta_5 = -4 III_{\bold{C}_e}\frac{\partial^2 \Psi_e}{\partial II_{\bold{C}_e} \partial
     *  III_{\bold{C}_e}}
     * \f]
     * \f[
     *  \delta_6 = 4 \left( III_{\bold{C}_e} \frac{\partial \Psi_e}{\partial III_{\bold{C}_e}}
     *           + III_{\bold{C}_e}^2 \frac{\partial^2 \Psi_e}{\partial III_{\bold{C}_e}^2} \right)
     * \f]
     * \f[
     *  \delta_7 = -4 III_{\bold{C}_e} \frac{\partial \Psi_e}{\partial III_{\bold{C}_e}}
     * \f]
     * \f[
     *  \delta_8 = -4 \frac{\partial \Psi_e}{\partial II_{\bold{C}_e}}
     * \f]
     */
    void CalculateGammaDelta(
        const LINALG::Matrix<3, 1>&
            prinv,  ///< principal invariants of the elastic right Cauchy-Green tensor
        const LINALG::Matrix<3, 1>& dPI,     ///< first derivative with respect to invariants
        const LINALG::Matrix<6, 1>& ddPII,   ///< second derivative with respect to invariants
        LINALG::Matrix<3, 1>& gamma,         ///< factors for stress calculation
        LINALG::Matrix<8, 1>& delta) const;  ///< factors for elasticity tensor calculation

    /// Evaluate off-diagonal stiffness matrix
    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const std::vector<LINALG::Matrix<3, 3>>&
            iFinjM,  ///< Vector that holds all inverse inelastic deformation gradient factors as
                     ///< 3x3 matrices
        const LINALG::Matrix<6, 9>& dSdiFin,  ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                                              ///< the inverse inelastic deformation gradient
        LINALG::Matrix<6, 1>& dstressdx);     ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                                              ///< primary variable of different field

    /*! \brief Evaluate additional terms of the elasticity tensor
     *
     * \f[
     *   cmatadd = \frac{\partial\bold{S}}{\partial \bold{F}_{in}^{-1}}
     *             :\frac{\partial\bold{F}_{in}^{-1}}{\partial \bold{C}})
     * \f]
     */
    void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const std::vector<LINALG::Matrix<3, 3>>&
            iFinjM,  ///< Vector that holds all inverse inelastic deformation gradient factors as
                     ///< 3x3 matrices
        const LINALG::Matrix<6, 1>& iCV,       ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>& dSdiFin,   ///< Derivative of 2nd Piola Kirchhoff stress w.r.t.
                                               ///< the inverse inelastic deformation gradient
        LINALG::Matrix<6, 6>& cmatadd) const;  ///< Additional elasticity tensor

    /// Evaluate derivative of 2nd Piola Kirchhoff stresses w.r.t. the inelastic deformation
    /// gradient
    void EvaluatedSdiFin(const LINALG::Matrix<3, 1>& gamma,  ///< Factors for stress calculation
        const LINALG::Matrix<8, 1>& delta,        ///< Factors for elasticity tensor calculation
        const LINALG::Matrix<3, 3>& iFinM,        ///< Inverse inelastic deformation gradient
        const LINALG::Matrix<3, 3>& iCinCM,       ///< C_{in}^{-1} * C
        const LINALG::Matrix<6, 1>& iCinV,        ///< Inverse inelastic right Cauchy-Green tensor
        const LINALG::Matrix<9, 1>& CiFin9x1,     ///< C * F_{in}^{-1}
        const LINALG::Matrix<9, 1>& CiFinCe9x1,   ///< C * F_{in}^{-1} * C_e
        const LINALG::Matrix<6, 1>& iCinCiCinV,   ///< C_{in}^{-1} * C * C_{in}^{-1}
        const LINALG::Matrix<9, 1>& CiFiniCe9x1,  ///< C * F_{in}^{-1} * C_e^{-1}
        const LINALG::Matrix<6, 1>& iCV,          ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<3, 3>& iFinCeM,      ///< F_{in}^{-1} * C_e
        LINALG::Matrix<6, 9>& dSdiFin)
        const;  ///< derivative of 2nd Piola Kirchhoff stresses w.r.t. inverse inelastic right
                ///< Cauchy-Green tensor

    /// evaluate inelastic deformation gradient
    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        std::vector<LINALG::Matrix<3, 3>>&
            iFinjM,  ///< Vector that holds all inverse inelastic deformation gradient factors as
                     ///< 3x3 matrices
        LINALG::Matrix<3, 3>& iFinM) const;  ///< Inverse inelastic deformation gradient

    /*! \brief Evaluates stress and cmat
     *
     * \f[
     *   cmat = 2 \frac{\partial \bold{S}}{\partial \bold{C}})
     * \f]
     *
     *  cf. Holzapfel [1], p. 216 (6.32) and p. 248
     *  The formulation in Holzapfel is slightly changed as only the elastic strains cause stresses
     *  with \f$\bold{C}_{in} = \bold{F}_{in}^T \bold{F}_{in}\f$
     *  and \f$\bold{F}_{in}\f$ detoting the inelastic deformation gradient
     * \f[
     *   \bold{S} = 2 \bold{F}_{in}^{-1} \frac{\partial \Psi_e}{\partial \bold{C}_e}
     *   \bold{F}_{in}^{-T}
     * \f]
     * \f[
     *   \bold{S} = \gamma_1 \ \bold{C}_{in}^{-1} + \gamma_2 \ \bold{C}_{in}^{-1} . \bold{C} .
     *   \bold{C}_{in}^{-1} + \gamma_3 \ \bold{C}^{-1}
     * \f]
     * with \f$\gamma_i\f$ as defined above for CalculateGammaDelta()
     *
     * \f[
     *   \mathbb{C} = \delta_1 \left( \bold{C}_{in}^{-1} \otimes \bold{C}_{in}^{-1} \right)
     *              + \delta_2 \left( \bold{C}_{in}^{-1} \otimes \bold{C}_{in}^{-1} \bold{C}
     *              \bold{C}_{in}^{-1}  + \bold{C}_{in}^{-1} \bold{C} \bold{C}_{in}^{-1} \otimes
     *              \bold{C}_{in}^{-1} \right)
     *              + \delta_3 \left( \bold{C}_{in}^{-1} \otimes \bold{C}^{-1} + \bold{C}^{-1}
     *              \otimes \bold{C}_{in}^{-1} \right)
     *              + \delta_4 \left( \bold{C}_{in}^{-1} \bold{C} \bold{C}_{in}^{-1} \otimes
     *              \bold{C}_{in}^{-1} \bold{C} \bold{C}_{in}^{-1} \right)
     *              + \delta_5 \left( \bold{C}_{in}^{-1} \bold{C} \bold{C}_{in}^{-1} \otimes
     *              \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{C}_{in}^{-1} \bold{C}
     *              \bold{C}_{in}^{-1} \right)
     *              + \delta_6 \left( \bold{C}^{-1} \otimes \bold{C}^{-1} \right)
     *              + \delta_7 \left( \bold{C}^{-1} \odot \bold{C}^{-1} \right)
     *              + \delta_8 \left( \bold{C}_{in}^{-1} \odot \bold{C}_{in}^{-1} \right)
     * \f]
     *
     * \f$\odot\f$ is defined as follows:
     * \f[
     *   - (A^{-1} \odot A^{-1})_{abcd} = 1/2 (A^{-1}_{ac} A^{-1}_{bd} +
     *   A^{-1}_{ad} A^{-1}_{bc}); \text{ for } A_{ab} = A_{ba}
     * \f] meaning a symmetric 2-tensor \f$A\f$
     */
    void EvaluateStressCmatIso(
        const LINALG::Matrix<6, 1>& iCV,         ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 1>& iCinV,       ///< Inverse inelastic right Cauchy-Green tensor
        const LINALG::Matrix<6, 1>& iCinCiCinV,  ///< C_{in}^{-1} * C * C_{in}^{-1}
        const LINALG::Matrix<3, 1>& gamma,       ///< Factors for stress calculation
        const LINALG::Matrix<8, 1>& delta,       ///< Factors for elasticity tensor calculation
        LINALG::Matrix<6, 1>& stress,            ///< Isotropic stress tensor
        LINALG::Matrix<6, 6>& cmatiso) const;    ///< Isotropic stiffness matrix

    /// Evaluates some kinematic quantities that are used in stress and elasticity tensor
    /// calculation
    void EvaluateKinQuantElast(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinM,  ///< Inverse inelastic deformation gradient
        LINALG::Matrix<6, 1>& iCinV,        ///< Inverse inelastic right Cauchy-Green tensor
        LINALG::Matrix<6, 1>& iCinCiCinV,   ///< C_{in}^{-1} * C * C_{in}^{-1}
        LINALG::Matrix<6, 1>& iCV,          ///< Inverse right Cauchy-Green tensor
        LINALG::Matrix<3, 3>& iCinCM,       ///< C_{in}^{-1} * C
        LINALG::Matrix<3, 3>& iFinCeM,      ///< F_{in}^{-1} * C_e
        LINALG::Matrix<9, 1>& CiFin9x1,     ///< C * F_{in}^{-1}
        LINALG::Matrix<9, 1>& CiFinCe9x1,   ///< C * F_{in}^{-1} * C_e
        LINALG::Matrix<9, 1>& CiFiniCe9x1,  ///< C * F_{in}^{-1} * C_e^{-1}
        LINALG::Matrix<3, 1>& prinv)
        const;  ///< Principal invariants of the elastic right Cauchy-Green tensor

    /// calculates the derivatives of the hyperelastic laws with respect to the invariants
    void EvaluateInvariantDerivatives(
        const LINALG::Matrix<3, 1>&
            prinv,         ///< Principal invariants of the elastic right Cauchy-Green tensor
        const int eleGID,  ///< Element ID
        LINALG::Matrix<3, 1>& dPI,           ///< First derivative with respect to invariants
        LINALG::Matrix<6, 1>& ddPII) const;  ///< Second derivative with respect to invariants

    /// pre-evaluation, intended to be used for stuff that have to be done only once
    void PreEvaluate(
        Teuchos::ParameterList& params) const;  ///< parameter list as handed in from the element

   private:
    /// My material parameters
    MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper* params_;

    /// map to elastic materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsumel_;

    /// map to inelastic deformation gradient factors
    std::vector<Teuchos::RCP<MAT::InelasticDefgradFactors>> facdefgradin_;
  };

}  // namespace MAT
#endif
