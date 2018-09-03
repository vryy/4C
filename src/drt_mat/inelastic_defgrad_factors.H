/*----------------------------------------------------------------------*/
/*!
\file inelastic_defgrad_factors.H

\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

<pre>
\maintainer Christoph Schmidt
</pre>
*/

#ifndef INELASTIC_DEFGRAD_FACTORS_H
#define INELASTIC_DEFGRAD_FACTORS_H

#include "matpar_parameter.H"
#include "../linalg/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*/
/*! \class InelasticDefgradFactors

    Provides the interface called by the class "MultiplicativeSplitDefgrad_ElastHyper"
    and is needed to evaluate the inelastic deformation gradient and
    their derivatives w.r.t. the primary variables.

    In the material "MultiplicativeSplitDefgrad_ElastHyper" the deformation gradient is split
    multiplicatively in elastic and inelastic deformation gradients (F = F_{el} * F_{in}).
    The inelastic deformation gradient itself can be a product of different inelastic
    deformation gradients, i.e. F_{in} = F_{in,1} * F_{in,2} * ... * F_{in,n}.
    The derived classes below are needed to evaluate the inverse of the j-th inelastic
    deformation gradient F_{in,j}^{-1} and its derivatives w.r.t. the primary variables.
*/

/* forward declarations */
namespace MAT
{
  class InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradFactors();

    /// construct material with specific material params
    explicit InelasticDefgradFactors(MAT::PAR::Parameter* params);

    /// destructor
    virtual ~InelasticDefgradFactors() { ; }

    /// create object by input parameter ID
    static Teuchos::RCP<InelasticDefgradFactors> Factory(int matnum);  ///< material ID

    /// provide material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    /// evaluate the deformation gradient and its inverse
    virtual void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
        LINALG::Matrix<3, 3>& iFinM) = 0;     ///< Inverse inelastic deformation gradient

    /// evaluate additional terms for the elasticity tensor
    virtual void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 6>& cmatadd) = 0;  ///< Additional elasticity tensor

    /// evaluate off-diagonal stiffness matrix for monolithic systems to get the
    /// cross-linearizations
    virtual void EvaluateODStiffMat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 1>& dstressdx) = 0;  ///< Derivative of 2nd Piola Kirchhoff stresses
                                               ///< w.r.t. primary variable of different field

    /// pre-evaluation, intended to be used for stuff that has to be done only once per Evaluate()
    virtual void PreEvaluate(
        Teuchos::ParameterList& params) = 0;  ///< parameter list as handed in from the element

    /// return material parameters
    MAT::PAR::Parameter* Parameter() { return params_; }

   private:
    /// material parameters
    MAT::PAR::Parameter* params_;

  };  // InelasticDefgradFactors



  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarIso

      This inelastic deformation gradient provides an isotropic growth law.
      Volumetric change due to this law is linearly dependent on a
      concentration mapped to the reference frame as follows:
      F_{in} = [1.0 + Scalar1growthfac_ * det(F) * (c - Scalar1refconc_)]^(1/3) * I
     */

  namespace PAR
  {
    class InelasticDefgradLinScalarIso : public MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      InelasticDefgradLinScalarIso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradLinScalarIso() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      /// return scalar that causes growth
      int Scalar1() { return Scalar1_; };
      /// return growth factor for scalar 1
      double Scalar1growthfac() { return Scalar1growthfac_; };
      /// return value of scalar 1 for that there is no growth
      double Scalar1refconc() { return Scalar1refconc_; };

     private:
      /// scalar that causes growth
      const int Scalar1_;
      /// growth factor for scalar 1
      const double Scalar1growthfac_;
      /// value of scalar 1 for that there is no growth
      const double Scalar1refconc_;

    };  // class InelasticDefgradLinScalarIso

  }  // namespace PAR


  class InelasticDefgradLinScalarIso : public InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarIso();

    /// construct material with specific material params
    explicit InelasticDefgradLinScalarIso(MAT::PAR::Parameter* params);

    /// destructor
    virtual ~InelasticDefgradLinScalarIso() { ; }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mfi_lin_scalar_iso; }

    /// evaluate the deformation gradient and its inverse
    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
        LINALG::Matrix<3, 3>& iFinM) final;   ///< Inverse inelastic deformation gradient

    /// evaluate additional terms of the elasticity tensor
    void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 6>& cmatadd) final;  ///< Additional elasticity tensor

    /// evaluate off-diagonal stiffness matrix for monolithic systems to get the
    /// cross-linearizations
    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 1>& dstressdc)
        final;  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. concentration

    /// pre-evaluation, intended to be used for stuff that has to be done only once per Evaluate()
    void PreEvaluate(
        Teuchos::ParameterList& params) final;  ///< parameter list as handed in from the element

    /// return the material parameters
    MAT::PAR::InelasticDefgradLinScalarIso* Parameter()
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradLinScalarIso*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;

  };  // InelasticDefgradLinScalarIso



  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarAniso

      This inelastic deformation gradient provides an anisotropic growth law.
      Volumetric change due to this law is linearly dependent on a
      concentration mapped to the reference frame as follows:
      F_{in} = I + [Scalar1growthfac_ * det(F) * (c - Scalar1refconc_)] * G,

      where G (growthdirmat_) is a matrix providing the information of the growth
      direction, that is constructed as follows:

      G = g \otimes g, where g is the growth direction vector given in the input file.

      g is normalized to length 1 before calculation of G.
     */

  namespace PAR
  {
    class InelasticDefgradLinScalarAniso : public Parameter
    {
     public:
      /// standard constructor
      InelasticDefgradLinScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradLinScalarAniso() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      /// return scalar that causes growth
      int Scalar1() { return Scalar1_; };
      /// return growth factor for scalar 1
      double Scalar1growthfac() { return Scalar1growthfac_; };
      /// return value of scalar 1 for that there is no growth
      double Scalar1refconc() { return Scalar1refconc_; };
      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() { return growthdirmat_; };

     private:
      /// scalar that causes growth
      const int Scalar1_;
      /// growth factor for scalar 1
      const double Scalar1growthfac_;
      /// value of scalar 1 for that there is no growth
      const double Scalar1refconc_;
      /// matrix that determines growth direction
      LINALG::Matrix<3, 3> growthdirmat_;

    };  // class InelasticDefgradLinScalarAniso

  }  // namespace PAR


  class InelasticDefgradLinScalarAniso : public InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarAniso();

    /// construct material with specific material params
    explicit InelasticDefgradLinScalarAniso(MAT::PAR::Parameter* params);

    /// destructor
    virtual ~InelasticDefgradLinScalarAniso() { ; }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mfi_lin_scalar_aniso; }

    /// evaluate the deformation gradient and its inverse
    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
        LINALG::Matrix<3, 3>& iFinM) final;   ///< Inverse inelastic deformation gradient

    /// evaluate additional terms for the elasticity tensor
    void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 6>& cmatadd) final;  ///< Additional elasticity tensor

    /// evaluate off-diagonal stiffness matrix for monolithic systems to get the
    /// cross-linearizations
    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 1>& dstressdc)
        final;  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. concentration

    /// pre-evaluation, intended to be used for stuff that has to be done only once per Evaluate()
    void PreEvaluate(
        Teuchos::ParameterList& params) final;  ///< parameter list as handed in from the element

    /// return the material parameters
    MAT::PAR::InelasticDefgradLinScalarAniso* Parameter()
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;

  };  // InelasticDefgradLinScalarAniso


}  // namespace MAT

#endif
