/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

*/

#ifndef INELASTIC_DEFGRAD_FACTORS_H
#define INELASTIC_DEFGRAD_FACTORS_H

#include "mat_par_parameter.H"
#include "linalg_fixedsizematrix.H"

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

      /// return reference to matrix that determines growth direction
      const LINALG::Matrix<3, 3>& Growthdirmat() const { return growthdirmat_; }

     private:
      /// matrix that determines growth direction
      LINALG::Matrix<3, 3> growthdirmat_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradNoGrowth
     *
     * This is a parameter class that is only needed to implement the pure virtual method
     * 'CreateMaterial()'.
     */
    class InelasticDefgradNoGrowth : public MAT::PAR::Parameter
    {
     public:
      //! standard constructor
      explicit InelasticDefgradNoGrowth(Teuchos::RCP<MAT::PAR::Material> matdata);

      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; }
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

      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; }

      /// return scalar that causes growth
      int Scalar1() const { return scalar1_; }

     private:
      /// scalar that causes growth
      const int scalar1_;
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

      /// return saturation concentration of material
      double Cmax() const { return c_max_; }
      /// return intercalation fraction at saturation concentration of material
      double Chimax() const { return chi_max_; }

     private:
      /// saturation concentration of material
      double c_max_;
      /// intercalation fraction at saturation concentration of material
      double chi_max_;
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

      /// return value of polynomial at reference intercalation fraction
      double GetPolynomReferenceValue() const { return polynomReferenceValue_; }

      // set value of polynomial at reference intercalation fraction
      void SetPolynomReferenceValue(double polynomReferenceValue)
      {
        polynomReferenceValue_ = polynomReferenceValue;
      }

     private:
      /// value of polynomial at reference intercalation fraction
      double polynomReferenceValue_;
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

      /// return reference to matrix that determines growth direction
      const LINALG::Matrix<3, 3>& Growthdirmat() const { return growthdir_->Growthdirmat(); };

     private:
      /// pointer to object, that calculates and holds direction of inelastic deformation
      Teuchos::RCP<InelasticDeformationDirection> growthdir_;
    };

    /*----------------------------------------------------------------------
    ----------------------------------------------------------------------*/
    /*! \class InelasticDefgradLinTempIso

    Parameter class of InelasticDefgradLinTempIso.
    */
    class InelasticDefgradLinTempIso : public MAT::PAR::Parameter
    {
     public:
      explicit InelasticDefgradLinTempIso(Teuchos::RCP<MAT::PAR::Material> matdata);

      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      /// return temperature related growth factor
      double GetTempGrowthFac() const { return tempgrowthfac_; };

      /// return value of temperature that causes no growth
      double RefTemp() const { return reftemp_; };

     private:
      /// value of temperature that causes no growth
      const double reftemp_;

      /// growth factor
      const double tempgrowthfac_;
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
     * @param[in] growthFac       linear growth factor (slope of linear function)
     * @param[in] referenceValue  reference value
     */
    explicit InelasticDefgradLinearShape(double growthFac, double referenceValue);

    /*!
     * @brief evaluation of the linear growth law
     *
     * @param[in] value           value the linear relation shall be evaluated for
     * @return growth factor
     */
    double EvaluateLinearGrowth(double value) const;

    // return the growth factor (needed for linearizations)
    double GrowthFac() const { return growthFac_; }

   private:
    /// growth factor
    const double growthFac_;
    /// reference value
    const double referenceValue_;
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
     * @param[in] poly_coeffs  coefficients describing the polynomial to be evaluated
     * @param[in] x_min        lower bound of validity of the polynomial
     * @param[in] x_max        upper bound of validity of the polynomial
     */
    explicit InelasticDefgradPolynomialShape(
        std::vector<double> poly_coeffs, double x_min, double x_max);

    /*!
     * @brief checks the bounds of validity of the polynomial and writes a warning to screen if
     * bounds are violated
     *
     * @param[in] x  value the polynomial is evaluated at
     */
    void CheckPolynomialBounds(double x) const;

    /*!
     * @brief Evaluate the polynomial defined by #PolyCoeffs_ at the current position X
     *
     * @param[in] x  value the polynomial is evaluated at
     * @return value of the polynomial evaluated at x
     */
    double ComputePolynomial(double x);

    /*!
     * @brief Evaluate the first derivative of the polynomial defined by #PolyCoeffs_ at the current
     * position x
     *
     * @param[in] x  value the polynomial is evaluated at
     * @return value the first derivative of the polynomial evaluated at x
     */
    double ComputePolynomialDerivative(double x);

   private:
    /// coefficients of the polynomial to be evaluated
    const std::vector<double> poly_coeffs_;
    /// lower bound of validity of polynomial
    const double x_min_;
    /// upper bound of validity of polynomial
    const double x_max_;
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
    /// construct material with specific material params
    explicit InelasticDefgradFactors(MAT::PAR::Parameter* params);

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
     * @brief calculate the derivative of the inelastic deformation gradient
     *
     * @param[in] detjacobian  determinant of the deformation gradient
     * @param[out] dFindx      derivative of inelastic deformation gradient w.r.t. primary variable
     *                         of different field
     */
    virtual void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) = 0;

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
     * @param[in] gp      Gauss point
     */
    virtual void PreEvaluate(Teuchos::ParameterList& params, int gp) = 0;

    /*!
     * @brief set gauss point concentration to parameter class
     *
     * @param[in] concentration  gauss point concentration to be set to internal member of parameter
     *                           class
     *
     * @note This method is used by methods called from the contact algorithm. Since the gauss point
     * ids do not match anyways (volume vs. surface element gauss point ids) and the id is not
     * relevant since the method is only called for one gauss point anyways, we set it to a dummy
     * gauss point id of 0 here
     */
    virtual void SetConcentrationGP(double concentration){};

    /// return material parameters
    virtual MAT::PAR::Parameter* Parameter() { return params_; }

    /// Get type of scalar, that leads to deformation
    virtual PAR::InelasticSource GetInelasticSource() = 0;

   protected:
    //! Get ID of current Gauss point
    int GetGP() const { return gp_; };

    //! Set ID of current Gauss point
    void SetGP(int gp) { gp_ = gp; };

   private:
    //! ID of current Gauss point
    int gp_;

    /// material parameters
    MAT::PAR::Parameter* params_;
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradNoGrowth

   This class models materials in combination with the multiplicative split material that feature
   no volume changes, i.e. the inelastic deformation gradient is always the identity tensor and
   contributions to the linearizations therefore vanish.
   */
  class InelasticDefgradNoGrowth : public InelasticDefgradFactors
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params           pointer to material specific parameters
     */
    explicit InelasticDefgradNoGrowth(MAT::PAR::Parameter* params);

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdx) override;

    PAR::InelasticSource GetInelasticSource() override;

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mfi_no_growth; }

    void PreEvaluate(Teuchos::ParameterList& params, int gp) override;

   private:
    // identity tensor
    LINALG::Matrix<3, 3> identity_;
  };

  class InelasticDefgradScalar : public InelasticDefgradFactors
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params           pointer to material specific parameters
     */
    explicit InelasticDefgradScalar(MAT::PAR::Parameter* params);

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override = 0;

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override = 0;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdx) override = 0;

    PAR::InelasticSource GetInelasticSource() override = 0;

    INPAR::MAT::MaterialType MaterialType() const override = 0;

    void PreEvaluate(Teuchos::ParameterList& params, int gp) override;

    void SetConcentrationGP(double concentration) override;

    MAT::PAR::InelasticDefgradScalar* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradScalar*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

   protected:
    //! Get vector of concentration at current Gauss point
    std::vector<double>& GetConcentrationGP() const { return concentrations_->at(GetGP()); };

   private:
    /// store vector of gauss point concentrations calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<std::vector<double>>> concentrations_;
  };

  /*--------------------------------------------------------------------*/
  /*! \class InelasticDefgradPolyIntercalFrac

   This class evaluates polynomial and its first derivative w.r.t. intercalation fraction which is
   required in various routines of subclasses for isotropic and anisotropic case. This polynomial
   describes the growth of material with respect to intercalation fraction and it is prescribed by
   user in input file by defining it coefficients.
   */
  class InelasticDefgradPolyIntercalFrac : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFrac(MAT::PAR::Parameter* params,
        Teuchos::RCP<InelasticDefgradPolynomialShape> polynomial_growth);

    /*!
     * @brief evaluate polynomial describing growth of material with regard to intercalation
     * fraction based on the current concentration
     *
     * @param[in] concentration current concentration
     * @param[in] detjacobian   determinant of the deformation gradient
     * @return value of polynomial describing the growth according to current intercalation fraction
     */
    double EvaluatePolynomial(double concentration, double detjacobian);

    /*!
     * @brief evaluate the first derivative of the polynomial describing the growth
     *
     * @param[in] concentration current concentration
     * @param[in] detjacobian   determinant of the deformation gradient
     * @return first derivative of the polynomial describing the growth
     */
    double EvaluatePolynomialDerivative(double concentration, double detjacobian);

    INPAR::MAT::MaterialType MaterialType() const override = 0;

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override = 0;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override = 0;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override = 0;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdx) override = 0;

    MAT::PAR::InelasticSource GetInelasticSource() override;

    MAT::PAR::InelasticDefgradPolyIntercalFrac* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFrac*>(
          MAT::InelasticDefgradScalar::Parameter());
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
  class InelasticDefgradLinScalarIso : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] linear_growth   pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarIso(
        MAT::PAR::Parameter* params, Teuchos::RCP<InelasticDefgradLinearShape> linear_growth);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_lin_scalar_iso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    MAT::PAR::InelasticSource GetInelasticSource() override;

   private:
    /// pointer to class that evaluates the linear growth law
    Teuchos::RCP<InelasticDefgradLinearShape> linear_growth_;
  };  // namespace MAT

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinScalarAniso

     This inelastic deformation gradient provides an anisotropic growth law.
     Volumetric change due to this law is dependent on the current concentration \f$ c \f$ as
     follows:
     \f[
     \mathbf{F}_\text{in} = \mathbf{I} + \left[ \text{Scalar1MolarGrowthFac}
     \left( c \det\mathbf{F}  - \text{Scalar1refconc} \right) \right] \mathbf{G},
     \f]
     where \f$ \mathbf{G} \f$ (#growthdirmat_) is a matrix providing the information of the
     growth direction, that is constructed as follows:
     \f$ \mathbf{G} = \mathbf{g} \otimes \mathbf{g} \f$,
     where \f$ \mathbf{g} \f$ is the growth direction vector given in the input file.
     \f$ \mathbf{g} \f$ is normalized to length 1 before calculation of \f$ \mathbf{G} \f$.
     \f$ f(\chi) \f$ is defined by the user in the input file.
     */
  class InelasticDefgradLinScalarAniso : public InelasticDefgradScalar
  {
   public:
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params          pointer to material specific parameters
     * @param[in] linear_growth   pointer to object that evaluates the linear relation as prescribed
     *                            in the input file
     */
    explicit InelasticDefgradLinScalarAniso(
        MAT::PAR::Parameter* params, Teuchos::RCP<InelasticDefgradLinearShape> linear_growth);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_lin_scalar_aniso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    MAT::PAR::InelasticSource GetInelasticSource() override;

    MAT::PAR::InelasticDefgradLinScalarAniso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(
          MAT::InelasticDefgradScalar::Parameter());
    }

   private:
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
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFracIso(MAT::PAR::Parameter* params,
        const Teuchos::RCP<InelasticDefgradPolynomialShape>& polynomial_growth);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_intercal_frac_iso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    MAT::PAR::InelasticDefgradPolyIntercalFrac* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFrac*>(
          MAT::InelasticDefgradPolyIntercalFrac::Parameter());
    }
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
    /*!
     * @brief construct material with required inputs
     *
     * @param[in] params             pointer to material specific parameters
     * @param[in] polynomial_growth  pointer to object that evaluates the polynomial as prescribed
     *                               in the input file
     */
    explicit InelasticDefgradPolyIntercalFracAniso(MAT::PAR::Parameter* params,
        const Teuchos::RCP<InelasticDefgradPolynomialShape>& polynomial_growth);

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::mfi_poly_intercal_frac_aniso;
    }

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc) override;

    MAT::PAR::InelasticDefgradPolyIntercalFracAniso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFracAniso*>(
          MAT::InelasticDefgradPolyIntercalFrac::Parameter());
    }
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticDefgradLinTempIso
   *Volumetric change due to this law is linearily dependent on the temperature
   \f$ T \f$ as follows:
   \f[
   \boldsymbol{F}_\text{in} = \boldsymbol{I} \left[ 1 + \beta \left( T - T_\text{ref} \right)
   \right]^\frac{1}{3}, \f]
   */
  class InelasticDefgradLinTempIso : public InelasticDefgradFactors
  {
   public:
    explicit InelasticDefgradLinTempIso(MAT::PAR::Parameter* params);

    void EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
        const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd) override;

    void EvaluateInelasticDefGradDerivative(
        double detjacobian, LINALG::Matrix<9, 1>& dFindx) override;

    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM) override;

    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<3, 3>& iFinjM,
        const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdT) override;

    MAT::PAR::InelasticSource GetInelasticSource() override;

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mfi_lin_temp_iso; };

    MAT::PAR::InelasticDefgradLinTempIso* Parameter() override
    {
      return dynamic_cast<MAT::PAR::InelasticDefgradLinTempIso*>(
          MAT::InelasticDefgradFactors::Parameter());
    }

    void PreEvaluate(Teuchos::ParameterList& params, int gp) override;

   private:
    //! Get temperature at current Gauss point
    double GetTemperatureGP() const { return temperatures_->at(GetGP()); };

    /// store vector of gauss point temperature calculated in the pre-evaluate of the so3_scatra
    /// element
    Teuchos::RCP<std::vector<double>> temperatures_;
  };
}  // namespace MAT

#endif
