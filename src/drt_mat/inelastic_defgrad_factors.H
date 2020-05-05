/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

\maintainer Christoph Schmidt
*/

#ifndef INELASTIC_DEFGRAD_FACTORS_H
#define INELASTIC_DEFGRAD_FACTORS_H

#include "matpar_parameter.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT
{
  namespace PAR
  {
    enum class InelasticSource;

    /*----------------------------------------------------------------------*/
    /*! \class InelasticDeformationDirection
     Calculates and holds growth direction in matrix format for anisotropic growth
     */
    class InelasticDeformationDirection
    {
     public:
      /*!
       * @brief standard constructor
       * @param[in] growthdirection  direction of anisotropic growth
       */
      explicit InelasticDeformationDirection(std::vector<double> growthdirection);

      /// destructor
      virtual ~InelasticDeformationDirection() = default;

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() const { return growthdirmat_; }

     private:
      /// matrix that determines growth direction
      LINALG::Matrix<3, 3> growthdirmat_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradScalar
     This is a parameter class holding parameters for evaluation of inelastic deformation (incl.
     linearization) induced by a scalar. Special quantities are hold in derived classes
     */
    class InelasticDefgradScalar : public MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      explicit InelasticDefgradScalar(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradScalar() override = default;

      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; }

      /// return scalar that causes growth
      int Scalar1() const { return Scalar1_; }
      /// return value of scalar 1 for that there is no growth
      double Scalar1refconc() const { return Scalar1refconc_; }
      /// return saturation concentration of material
      double Cmax() const { return Cmax_; }

     private:
      /// scalar that causes growth
      const int Scalar1_;
      /// value of scalar 1 for that there is no growth
      const double Scalar1refconc_;
      /// material ID of corresponding scatra material
      const int Matid_;
      /// saturation concentration of material
      double Cmax_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradLinScalar : public InelasticDefgradScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalar(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradLinScalar() override = default;

      /// return growth factor for scalar 1
      double Scalar1growthfac() const { return Scalar1growthfac_; }

     private:
      /// growth factor for scalar 1
      const double Scalar1growthfac_;
    };
    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradLinScalarIso : public InelasticDefgradLinScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalarIso(Teuchos::RCP<MAT::PAR::Material> matdata)
          : InelasticDefgradLinScalar(matdata){};

      ~InelasticDefgradLinScalarIso() override = default;
    };  // end of class InelasticDefgradLinScalarIso

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradLinScalarAniso : public InelasticDefgradLinScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradLinScalarAniso() override = default;

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() { return growthdir_.Growthdirmat(); }

     private:
      /// object, that calculates and holds direction of inelastic deformation
      InelasticDeformationDirection growthdir_;
    };  // end of class InelasticDefgradLinScalarAniso

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradPolyScalar : public InelasticDefgradScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyScalar(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradPolyScalar() override = default;

      /// return vector containing polynomial coefficients
      std::vector<double> Polycoeffs() const { return Polycoeffs_; }
      /// return lower bound of validity of polynomial
      double Xmin() const { return Xmin_; }
      /// return upper bound of validity of polynomial
      double Xmax() const { return Xmax_; }
      /// return value of polynomial at reference intercalation fraction
      double MATPolynomReference() const { return MATPolynomReference_; }

     private:
      /// coefficients of polynomial
      const std::vector<double> Polycoeffs_;
      /// lower bound of validity of polynomial
      const double Xmin_;
      /// upper bound of validity of polynomial
      const double Xmax_;
      /// /// value of polynomial at reference intercalation fraction
      double MATPolynomReference_;
    };  // namespace PAR

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradPolyScalarIso : public InelasticDefgradPolyScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyScalarIso(Teuchos::RCP<MAT::PAR::Material> matdata)
          : InelasticDefgradPolyScalar(matdata){};

      ~InelasticDefgradPolyScalarIso() override = default;
    };  // end of class InelasticDefgradPolyScalarIso

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradPolyScalarAniso : public InelasticDefgradPolyScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradPolyScalarAniso() override = default;

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() const { return growthdir_.Growthdirmat(); }

     private:
      /// object, that calculates and holds direction of inelastic deformation
      InelasticDeformationDirection growthdir_;
    };  // end of class InelasticDefgradPolyScalarAniso

  }  // namespace PAR

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
  class InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradFactors();

    /// construct material with specific material params
    explicit InelasticDefgradFactors(MAT::PAR::Parameter* params);

    /// destructor
    virtual ~InelasticDefgradFactors() = default;

    /*!
     * @brief create object by input parameter ID
     *
     * @param[in] matnum  material ID
     * @return pointer to material that is defined by material ID
     */
    static Teuchos::RCP<InelasticDefgradFactors> Factory(int matnum);

    /// provide material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    /*!
     * @brief evaluate the deformation gradient and its inverse
     *
     * @param[in] defgrad  Deformation gradient
     * @param[out] iFinM   Inverse inelastic deformation gradient
     */
    virtual void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) = 0;

    /*!
     * @brief evaluate additional terms for the elasticity tensor
     *
     * @param[in] defgrad  Deformation gradient
     * @param[in] iFinjM   Inverse inelastic deformation gradient of current inelastic contribution
     *                     as 3x3 matrix
     * @param[in] iCV      Inverse right Cauchy-Green tensor
     * @param[in] dSdiFinj Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
     *                     deformation gradient of current inelastic contribution
     * @param[in/out] cmatadd  Additional elasticity tensor
     */
    virtual void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) = 0;


    /*!
     * @brief evaluate off-diagonal stiffness matrix for monolithic systems to get the
     *        cross-linearizations
     *
     * @param[in] defgrad Deformation gradient
     * @param[in] iFinjM  Inverse inelastic deformation gradient of current inelastic contribution
     *                    as 3x3 matrix
     * @param[in] dSdiFinj  Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
     *                      deformation gradient of current inelastic contribution
     * @param[in/out] dstressdx Derivative of 2nd Piola Kirchhoff stresses w.r.t. primary variable
     *                          of different field
     */
    virtual void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
        LINALG::Matrix<6, 1>& dstressdx) = 0;

    /*!
     * @brief pre-evaluation, intended to be used for stuff that has to be done only once per
     *        Evaluate()
     *
     * @param[in] params parameter list as handed in from the element
     */
    virtual void PreEvaluate(Teuchos::ParameterList& params) = 0;

    /// return material parameters
    virtual MAT::PAR::Parameter* Parameter() { return params_; }

    /// Get type of scalar, that leads to deformation
    virtual PAR::InelasticSource GetInelasticSource() = 0;

   private:
    /// material parameters
    MAT::PAR::Parameter* params_;

  };  // end of InelasticDefgradFactors

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinear

   This class evaluates delta growth which is required in various routines of subclasses for
   isotropic and anisotropic case. Delta growth is used to calculate the scaling factor which gives
   the magnitude of material growth.
   */
  class InelasticDefgradLinear : public InelasticDefgradFactors
  {
   public:
    /// constructor
    InelasticDefgradLinear();

    /// construct material with specific material params
    explicit InelasticDefgradLinear(MAT::PAR::Parameter* params);

    ~InelasticDefgradLinear() override = default;

    /*!
     * @brief evaluate the growth scaling factor
     *
     * @param[in] Growthfac      growth factor
     * @param[in] DetJacobian    determinant of the deformation gradient
     * @param[in] Referenceconc  reference concentration (has to equal the inital concentration)
     * @param[in] concentration  current concentration
     * @param[in] Cmax           maximum concentration
     * @return  growth scaling factor
     */
    static double EvaluateDeltaGrowth(double Growthfac, double DetJacobian, double Referenceconc,
        double concentration, double Cmax);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      dserror("not implemented here");
      return INPAR::MAT::m_none;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override
    {
      dserror("not implemented here");
    }

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* const defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override
    {
      dserror("not implemented here");
    }

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
        LINALG::Matrix<6, 1>& dstressdx) override
    {
      dserror("not implemented here");
    }

    void PreEvaluate(Teuchos::ParameterList& params) override { dserror("not implemented here"); }

    MAT::PAR::InelasticSource GetInelasticSource() override;
  };  // end of class InelasticDefgradLinear

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolynomial

   This class evaluates polynomial and its first derivative w.r.t. intercalation fraction which is
   required in various routines of subclasses for isotropic and anisotropic case. This polynomial
   describes the growth of material with respect to intercalation fraction and it is prescribed by
   user in input file by defining it coefficients.
   */
  class InelasticDefgradPolynomial : public InelasticDefgradFactors
  {
   public:
    /// constructor
    InelasticDefgradPolynomial() = default;

    /// construct material with specific material params
    explicit InelasticDefgradPolynomial(MAT::PAR::Parameter* params);

    ~InelasticDefgradPolynomial() override = default;

    /*!
     * @brief evaluate polynomial describing growth of material with regard to concentration
     *
     * @param[in] PolyCoeff     coefficients of the polynomial
     * @param[in] X_min         lower bound for lithiation
     * @param[in] X_max         upper bound for lithiation
     * @param[in] concentration current concentration
     * @param[in] Cmax          maximum concentration
     * @param[in] DetJacobian   determinant of the deformation gradient
     * @return value of polynomial describing the growth for current concentration
     */
    double EvaluateMatPolynomial(const std::vector<double>& PolyCoeff, double X_min, double X_max,
        double concentration, double Cmax, double DetJacobian);

    ///
    /*!
     * @brief evaluate the first derivative of the polynomial describing the growth
     *
     * @param[in] PolyCoeff     coefficients of the polynomial
     * @param[in] concentration current concentration
     * @param[in] Cmax          maximum concentration
     * @param[in] DetJacobian   determinant of the deformation gradient
     * @return value of the first derivative of the polynomial describing the growth for current
     *         concentration w.r.t. current lithiation
     */
    static double EvaluateMatPolynomialDerivative(const std::vector<double>& PolyCoeff,
        double concentration, double Cmax, double DetJacobian);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      dserror("not implemented here");
      return INPAR::MAT::m_none;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override
    {
      dserror("not implemented here");
    }

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* const defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override
    {
      dserror("not implemented here");
    }

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
        LINALG::Matrix<6, 1>& dstressdx) override
    {
      dserror("not implemented here");
    }

    void PreEvaluate(Teuchos::ParameterList& params) override { dserror("not implemented here"); }

    MAT::PAR::InelasticSource GetInelasticSource() override;

  };  // end of class InelasticDefgradPolynomial

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarIso

      This inelastic deformation gradient provides an isotropic growth law.
      Volumetric change due to this law is linearly dependent on a
      concentration mapped to the reference frame as follows:
      F_{in} = [1.0 + Scalar1growthfac_ * det(F) * (c - Scalar1refconc_) / C_MAX]^(1/3) * I
     */
  class InelasticDefgradLinScalarIso : public InelasticDefgradLinear
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarIso();

    /// construct material with specific material params
    explicit InelasticDefgradLinScalarIso(MAT::PAR::Parameter* params);

    ~InelasticDefgradLinScalarIso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_lin_scalar_iso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradLinScalarIso* Parameter() override
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

  };  // end of class InelasticDefgradLinScalarIso

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarAniso

      This inelastic deformation gradient provides an anisotropic growth law.
      Volumetric change due to this law is linearly dependent on a
      concentration mapped to the reference frame as follows:
      F_{in} = I + [Scalar1growthfac_ * det(F) * (c - Scalar1refconc_) / C_MAX] * G,
      where G (growthdirmat_) is a matrix providing the information of the growth
      direction, that is constructed as follows:

      G = g \otimes g, where g is the growth direction vector given in the input file.

      g is normalized to length 1 before calculation of G.
     */
  class InelasticDefgradLinScalarAniso : public InelasticDefgradLinear
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarAniso();

    /// construct material with specific material params
    explicit InelasticDefgradLinScalarAniso(MAT::PAR::Parameter* params);

    ~InelasticDefgradLinScalarAniso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_lin_scalar_aniso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradLinScalarAniso* Parameter() override
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

  };  // end of InelasticDefgradLinScalarAniso

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyScalarIso

   This inelastic deformation gradient provides an isotropic growth law.
   Volumetric change due to this law is non-linearly dependent on a
   concentration mapped to the reference frame as follows:
   F_{in} = [(1.0 + f(X)) / (1 + f(X_ref))]^(1/3) * I
   where f(X) is prescribed by user in input file.
   */
  class InelasticDefgradPolyScalarIso : public InelasticDefgradPolynomial
  {
   public:
    /// construct empty material
    InelasticDefgradPolyScalarIso();

    /// construct material with specific material params
    explicit InelasticDefgradPolyScalarIso(MAT::PAR::Parameter* params);

    ~InelasticDefgradPolyScalarIso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_scalar_iso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradPolyScalarIso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyScalarIso*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
  };  // end of class InelasticDefgradPolyScalarIso

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyScalarAniso

   This inelastic deformation gradient provides an anisotropic growth law.
   Volumetric change due to this law is nonlinearly dependent on a
   concentration mapped to the reference frame as follows:
   F_{in} = I + [(f(X) - f(X_ref)) / (f(X_ref) + 1)] * G,
   where G (growthdirmat_) is a matrix providing the information of the growth
   direction, that is constructed as follows:

   G = g \otimes g, where g is the growth direction vector given in the input file.
   g is normalized to length 1 before calculation of G. f(X) is defined by user in input file.
   */
  class InelasticDefgradPolyScalarAniso : public InelasticDefgradPolynomial
  {
   public:
    /// construct empty material
    InelasticDefgradPolyScalarAniso();

    /// construct material with specific material params
    explicit InelasticDefgradPolyScalarAniso(MAT::PAR::Parameter* params);

    ~InelasticDefgradPolyScalarAniso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_scalar_aniso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradPolyScalarAniso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyScalarAniso*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
  };  // end of class InelasticDefgradPolyScalarAniso

}  // namespace MAT

#endif