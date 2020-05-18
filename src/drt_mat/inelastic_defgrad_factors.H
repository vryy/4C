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
     *
     * Calculates and holds growth direction in matrix format for anisotropic growth
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

      /// return reference to matrix that determines growth direction
      const LINALG::Matrix<3, 3>& Growthdirmat() const { return growthdirmat_; }

     private:
      /// matrix that determines growth direction
      LINALG::Matrix<3, 3> growthdirmat_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradScalar
     *
     * This is a parameter class holding parameters for evaluation of inelastic deformation (incl.
     * linearization) induced by a scalar. Special quantities are hold in derived classes
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

     private:
      /// scalar that causes growth
      const int Scalar1_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradLinScalarAniso
     *
     * This is a specialized parameter class that can return the anisotropic growth direction
     * represented as a growth matrix
     */
    class InelasticDefgradLinScalarAniso : public InelasticDefgradScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradLinScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradLinScalarAniso() override = default;

      /// return reference to matrix that determines growth direction
      const LINALG::Matrix<3, 3>& Growthdirmat() { return growthdir_->Growthdirmat(); }

     private:
      /// pointer to object, that calculates and holds direction of inelastic deformation
      Teuchos::RCP<InelasticDeformationDirection> growthdir_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradIntercalFrac
     *
     * This parameter class provides all electrochemical quantities that are needed to calculate the
     * intercalation fraction from a given species concentration
     */
    class InelasticDefgradIntercalFrac : public InelasticDefgradScalar
    {
     public:
      /// standard constructor
      explicit InelasticDefgradIntercalFrac(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradIntercalFrac() override = default;

      /// return saturation concentration of material
      double Cmax() const { return Cmax_; }
      /// return intercalation fraction at saturation concentration of material
      double Chimax() const { return Chimax_; }

     private:
      /// saturation concentration of material
      double Cmax_;
      /// intercalation fraction at saturation concentration of material
      double Chimax_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradPolyIntercalFrac
     *
     * This parameter class provides the value of the polynomial that models the growth evaluated in
     * the reference configuration
     */
    class InelasticDefgradPolyIntercalFrac : public InelasticDefgradIntercalFrac
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyIntercalFrac(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradPolyIntercalFrac() override = default;

      /// return value of polynomial at reference intercalation fraction
      double GetPolynomReferenceValue() const { return PolynomReferenceValue_; }

      // set value of polynomial at reference intercalation fraction
      void SetPolynomReferenceValue(double PolynomReferenceValue)
      {
        PolynomReferenceValue_ = PolynomReferenceValue;
      }

     private:
      /// value of polynomial at reference intercalation fraction
      double PolynomReferenceValue_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradPolyIntercalFracAniso
     *
     * This is a specialized parameter class that can return the anisotropic growth direction
     * represented as a growth matrix
     */
    class InelasticDefgradPolyIntercalFracAniso : public InelasticDefgradPolyIntercalFrac
    {
     public:
      /// standard constructor
      explicit InelasticDefgradPolyIntercalFracAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      ~InelasticDefgradPolyIntercalFracAniso() override = default;

      /// return reference to matrix that determines growth direction
      const LINALG::Matrix<3, 3>& Growthdirmat() const { return growthdir_->Growthdirmat(); };

     private:
      /// pointer to object, that calculates and holds direction of inelastic deformation
      Teuchos::RCP<InelasticDeformationDirection> growthdir_;
    };
  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinearShape
   *
   * This class provides the functionality to be used if the growth law obeys a linear relation
   */
  class InelasticDefgradLinearShape
  {
   public:
    /*!
     * @brief constructor with required parameters
     *
     * @param[in] GrowthFac       linear growth factor (slope of linear function)
     * @param[in] ReferenceValue  reference value
     */
    explicit InelasticDefgradLinearShape(double GrowthFac, double ReferenceValue);

    /// destructor
    ~InelasticDefgradLinearShape() = default;

    /*!
     * @brief evaluation of the linear growth law
     *
     * @param[in] Value           value the linear relation shall be evaluated for
     * @return growth factor
     */
    double EvaluateLinearGrowth(double Value) const;

    // return the growth factor (needed for linearizations)
    double GrowthFac() const { return GrowthFac_; }

   private:
    /// growth factor
    const double GrowthFac_;
    /// reference value
    const double ReferenceValue_;
  };  // namespace MAT

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolynomialShape
   *
   * This class provides the functionality to be used if the growth law obeys a polynomial relation
   */
  class InelasticDefgradPolynomialShape
  {
   public:
    /*!
     * @brief  constructor with required parameters
     *
     * @param[in] PolyCoeffs  coefficients describing the polynomial to be evaluated
     * @param[in] X_min       lower bound of validity of the polynomial
     * @param[in] X_max       upper bound of validity of the polynomial
     */
    explicit InelasticDefgradPolynomialShape(
        std::vector<double> PolyCoeffs, double Xmin, double Xmax);

    /// destructor
    ~InelasticDefgradPolynomialShape() = default;

    /*!
     * @brief checks the bounds of validity of the polynomial and writes a warning to screen if
     * bounds are violated
     *
     * @param[in] X  value the polynomial is evaluated at
     */
    void CheckPolynomialBounds(double X) const;

    /*!
     * @brief Evaluate the polynomial defined by #PolyCoeffs_ at the current position X
     *
     * @param[in] X  value the polynomial is evaluated at
     * @return value of the polynomial evaluated at X
     */
    double ComputePolynomial(double X);

    /*!
     * @brief Evaluate the first derivative of the polynomial defined by #PolyCoeffs_ at the current
     * position X
     *
     * @param[in] X  value the polynomial is evaluated at
     * @return value the first derivative of the polynomial evaluated at X
     */
    double ComputePolynomialDerivative(double X);

   private:
    /// coefficients of the polynomial to be evaluated
    const std::vector<double> PolyCoeffs_;
    /// lower bound of validity of polynomial
    const double Xmin_;
    /// upper bound of validity of polynomial
    const double Xmax_;
  };

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
     * @brief evaluate the inelastic deformation gradient and its inverse
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
     * @param[in,out] cmatadd  Additional elasticity tensor
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
     * @param[in,out] dstressdx Derivative of 2nd Piola Kirchhoff stresses w.r.t. primary variable
     *                          of different field
     */
    virtual void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
        LINALG::Matrix<6, 1>& dstressdx) = 0;

    /*!
     * @brief pre-evaluation, intended to be used for stuff that has to be done only once per
     *        Evaluate()
     *
     * @param[in] params  parameter list as handed in from the element
     */
    virtual void PreEvaluate(Teuchos::ParameterList& params) = 0;

    /// return material parameters
    virtual MAT::PAR::Parameter* Parameter() { return params_; }

    /// Get type of scalar, that leads to deformation
    virtual PAR::InelasticSource GetInelasticSource() = 0;

   private:
    /// material parameters
    MAT::PAR::Parameter* params_;
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFrac

   This class evaluates polynomial and its first derivative w.r.t. intercalation fraction which is
   required in various routines of subclasses for isotropic and anisotropic case. This polynomial
   describes the growth of material with respect to intercalation fraction and it is prescribed by
   user in input file by defining it coefficients.
   */
  class InelasticDefgradPolyIntercalFrac : public InelasticDefgradFactors
  {
   public:
    /// constructor
    InelasticDefgradPolyIntercalFrac() = default;

    /*!
     * @brief construct material with reqired inputs
     *
     * @param[in] params           pointer to material specific parameters
     * @param[in] PolynomialGrowth pointer to object that evaluates the polynomial as prescribed in
     *            the input file
     */
    explicit InelasticDefgradPolyIntercalFrac(MAT::PAR::Parameter* params,
        const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth);

    ~InelasticDefgradPolyIntercalFrac() override = default;

    /*!
     * @brief evaluate polynomial describing growth of material with regard to intercalation
     * fraction based on the current concentration
     *
     * @param[in] concentration current concentration
     * @param[in] DetJacobian   determinant of the deformation gradient
     * @return value of polynomial describing the growth according to current intercalation fraction
     */
    double EvaluatePolynomial(double concentration, double DetJacobian);

    /*!
     * @brief evaluate the first derivative of the polynomial describing the growth
     *
     * @param[in] concentration current concentration
     * @param[in] DetJacobian   determinant of the deformation gradient
     * @return first derivative of the polynomial describing the growth
     */
    double EvaluatePolynomialDerivative(double concentration, double DetJacobian);

    INPAR::MAT::MaterialType MaterialType() const override = 0;

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override = 0;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override = 0;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdx) override = 0;

    void PreEvaluate(Teuchos::ParameterList& params) override = 0;

    MAT::PAR::InelasticSource GetInelasticSource() override;

    MAT::PAR::InelasticDefgradPolyIntercalFrac* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFrac*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// pointer to class that evaluates the polynomial growth law
    Teuchos::RCP<InelasticDefgradPolynomialShape> polynomial_growth_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarIso
        This inelastic deformation gradient provides an isotropic growth law. Volumetric change due
        to this law is dependent on the current concentration \f$ c \f$ as follows :
      \f[
      \boldsymbol{F} _\text{in} = \left[1 + \text { Scalar1MolarGrowthFac }
      \left(c \det \boldsymbol{F} - \text { Scalar1refconc } \right) \right] ^ { 1 / 3 }
      \boldsymbol{I}
      \f]
      */
  class InelasticDefgradLinScalarIso : public InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarIso();

    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] LinearGrowth    pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarIso(
        MAT::PAR::Parameter* params, const Teuchos::RCP<InelasticDefgradLinearShape>& LinearGrowth);

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

    MAT::PAR::InelasticSource GetInelasticSource() override;

    MAT::PAR::InelasticDefgradScalar* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradScalar*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
    /// pointer to class that evaluates the linear growth law
    Teuchos::RCP<InelasticDefgradLinearShape> linear_growth_;
  };  // namespace MAT

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarAniso

     This inelastic deformation gradient provides an anisotropic growth law.
     Volumetric change due to this law is dependent on the current concentration \f$ c \f$ as
     follows:
     \f[
     \bold{F}_\text{in} = \bold{I} + \left[ \text{Scalar1MolarGrowthFac}
     \left( c \det\bold{F}  - \text{Scalar1refconc} \right) \right] \bold{G},
     \f]
     where \f$ \bold{G} \f$ (#growthdirmat_) is a matrix providing the information of the
     growth direction, that is constructed as follows:
     \f$ \bold{G} = \bold{g} \otimes \bold{g} \f$,
     where \f$ \bold{g} \f$ is the growth direction vector given in the input file.
     \f$ \bold{g} \f$ is normalized to length 1 before calculation of \f$ \bold{G} \f$.
     \f$ f(\chi) \f$ is defined by the user in the input file.
     */
  class InelasticDefgradLinScalarAniso : public InelasticDefgradFactors
  {
   public:
    /// construct empty material
    InelasticDefgradLinScalarAniso();

    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] LinearGrowth    pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarAniso(
        MAT::PAR::Parameter* params, const Teuchos::RCP<InelasticDefgradLinearShape>& LinearGrowth);

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

    MAT::PAR::InelasticSource GetInelasticSource() override;

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
    /// store pointer to class that evaluates the linear growth law
    Teuchos::RCP<InelasticDefgradLinearShape> linear_growth_;
  };  // end of InelasticDefgradLinScalarAniso

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFracIso

   This inelastic deformation gradient provides an isotropic growth law.
   Volumetric change due to this law is non-linearly dependent on the intercalation fraction
   \f$ \chi \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} =
   \left[ \frac{f(\chi) + 1 }{f(\chi^0) + 1} \right]^{1/3} \boldsymbol{I},
   \f]
   where \f$ f(\chi) \f$ is defined by the user in the input file.
   */
  class InelasticDefgradPolyIntercalFracIso : public InelasticDefgradPolyIntercalFrac
  {
   public:
    /// construct empty material
    InelasticDefgradPolyIntercalFracIso();

    /*!
     * @brief construct material with reqired inputs
     *
     * @param[in] params           pointer to material specific parameters
     * @param[in] PolynomialGrowth pointer to object that evaluates the polynomial as prescribed in
     *                             the input file
     */
    explicit InelasticDefgradPolyIntercalFracIso(MAT::PAR::Parameter* params,
        const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth);


    ~InelasticDefgradPolyIntercalFracIso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_intercal_frac_iso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradPolyIntercalFrac* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFrac*>(
          MAT::InelasticDefgradPolyIntercalFrac::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFracAniso

   This inelastic deformation gradient provides an anisotropic growth law.
   Volumetric change due to this law is nonlinearly dependent on the intercalation fraction
   \f$ \chi \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} =
   \boldsymbol{I} + \left[ \frac{f(\chi) - f(\chi^0)}{f(\chi^0) + 1} \right] \boldsymbol{G},
   \f]
   where \f$ \boldsymbol{G} \f$ (#growthdirmat_) is a matrix providing the information of the growth
   direction, that is constructed as follows:
   \f$ \boldsymbol{G} = \boldsymbol{g} \otimes \boldsymbol{g} \f$, where \f$ \boldsymbol{g} \f$ is
   the growth direction vector given in the input file.
   \f$ \boldsymbol{g} \f$ is normalized to length 1 before calculation of \f$ \boldsymbol{G} \f$.
   \f$ f(\chi) \f$ is defined by the user in the input file.
   */
  class InelasticDefgradPolyIntercalFracAniso : public InelasticDefgradPolyIntercalFrac
  {
   public:
    /// construct empty material
    InelasticDefgradPolyIntercalFracAniso();

    /*!
     * @brief construct material with reqired inputs
     *
     * @param[in] params           pointer to material specific parameters
     * @param[in] PolynomialGrowth pointer to object that evaluates the polynomial as prescribed in
     *                             the input file
     */
    explicit InelasticDefgradPolyIntercalFracAniso(MAT::PAR::Parameter* params,
        const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth);

    ~InelasticDefgradPolyIntercalFracAniso() override = default;

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_intercal_frac_aniso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    void PreEvaluate(Teuchos::ParameterList& params) override;

    MAT::PAR::InelasticDefgradPolyIntercalFracAniso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFracAniso*>(
          MAT::InelasticDefgradPolyIntercalFrac::Parameter());
    }

   private:
    /// store current gauss point id
    int gp_;
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
  };
}  // namespace MAT

#endif