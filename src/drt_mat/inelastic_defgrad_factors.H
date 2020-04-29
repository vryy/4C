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

/* forward declarations */
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class InelasticDeformationDirection
     Calculates and holds growth direction in matrix format for anisotropic growth
     */
    class InelasticDeformationDirection
    {
     public:
      /// standard constructor
      InelasticDeformationDirection(std::vector<double> growthdirection);

      /// destructor
      virtual ~InelasticDeformationDirection() { ; }

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() { return growthdirmat_; };

     private:
      /// matrix that determines growth direction
      LINALG::Matrix<3, 3> growthdirmat_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    /*! \class InelasticDefgradScalar
     This is a parameter class holding parameters for evaluation of inelastic deformation (incl.
     linearization) induced by a scalar.
     It's inheritance hierarchy is as follows
     InelasticDefgradScalar
       -InelasticDefGradLinScalar
         --InelasticDefGradLinScalarIso
         --InelasticDefGradLinScalarAniso
       -InelasticDefGradPolyScalar
         --InelasticDefGradPolyScalarIso
         --InelasticDefGradPolyScalarAniso
     */
    class InelasticDefgradScalar : public MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      InelasticDefgradScalar(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradScalar() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      /// return scalar that causes growth
      int Scalar1() { return Scalar1_; };
      /// return value of scalar 1 for that there is no growth
      double Scalar1refconc() { return Scalar1refconc_; };
      /// return saturation concentration of material
      double Cmax() { return Cmax_; };

     protected:
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
      InelasticDefgradLinScalar(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradLinScalar() { ; }

      /// return growth factor for scalar 1
      double Scalar1growthfac() { return Scalar1growthfac_; };

     protected:
      /// growth factor for scalar 1
      const double Scalar1growthfac_;
    };
    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradLinScalarIso : public InelasticDefgradLinScalar
    {
     public:
      /// standard constructor
      InelasticDefgradLinScalarIso(Teuchos::RCP<MAT::PAR::Material> matdata)
          : InelasticDefgradLinScalar(matdata){};

    };  // end of class InelasticDefgradLinScalarIso

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradLinScalarAniso : public InelasticDefgradLinScalar
    {
     public:
      /// standard constructor
      InelasticDefgradLinScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradLinScalarAniso() { ; }

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() { return growthdir_.Growthdirmat(); };

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
      InelasticDefgradPolyScalar(Teuchos::RCP<MAT::PAR::Material> matdata);
      /// destructor
      virtual ~InelasticDefgradPolyScalar() { ; }

      /// return number of polynomial coefficients
      unsigned Polyparanum() { return Polyparanum_; };
      /// return vector containing polynomial coefficients
      std::vector<double> Polyparams() { return Polyparams_; };
      /// return lower bound of validity of polynomial
      double Xmin() { return Xmin_; };
      /// return upper bound of validity of polynomial
      double Xmax() { return Xmax_; };
      /// return saturation concentration of material
      /// /// return value of polynomial at reference intercalation fraction
      double MATPolynomReference() { return MATPolynomReference_; };

     protected:
      /// number of parameters of the polynomial describing volume change
      const unsigned Polyparanum_;
      /// coefficients of polynomial
      const std::vector<double> Polyparams_;
      /// lower bound of validity of polynomial
      const double Xmin_;
      /// upper bound of validity of polynomial
      const double Xmax_;
      /// /// value of polynomial at reference intercalation fraction
      double MATPolynomReference_;
    };

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradPolyScalarIso : public InelasticDefgradPolyScalar
    {
     public:
      /// standard constructor
      InelasticDefgradPolyScalarIso(Teuchos::RCP<MAT::PAR::Material> matdata)
          : InelasticDefgradPolyScalar(matdata){};

      /// destructor
      virtual ~InelasticDefgradPolyScalarIso() { ; }
    };  // end of class InelasticDefgradPolyScalarIso

    /*----------------------------------------------------------------------
     *----------------------------------------------------------------------*/
    class InelasticDefgradPolyScalarAniso : public InelasticDefgradPolyScalar
    {
     public:
      /// standard constructor
      InelasticDefgradPolyScalarAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~InelasticDefgradPolyScalarAniso() { ; }

      /// return matrix that determines growth direction
      LINALG::Matrix<3, 3> Growthdirmat() { return growthdir_.Growthdirmat(); };

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

    /// destructor
    ~InelasticDefgradLinear() { ; }

    /// evaluate scaling factor
    double EvaluateDeltaGrowth(const double Growthfac, const double DetJacobian,
        const double Referenceconc, const double concentration, const double Cmax);

    /// provide material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      dserror("not implemented here");
      return INPAR::MAT::m_none;
    }

    /// evaluate the deformation gradient and its inverse
    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad,   ///< Deformation gradient
        LINALG::Matrix<3, 3>& iFinM) override  ///< Inverse inelastic deformation gradient
    {
      dserror("not implemented here");
    }

    /// evaluate additional terms for the elasticity tensor
    void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 6>& cmatadd) override  ///< Additional elasticity tensor
    {
      dserror("not implemented here");
    }

    /// evaluate off-diagonal stiffness matrix for monolithic systems to get the
    /// cross-linearizations
    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 1>& dstressdx) override  ///< Derivative of 2nd Piola Kirchhoff stresses
                                                   ///< w.r.t. primary variable of different field
    {
      dserror("not implemented here");
    }

    /// pre-evaluation, intended to be used for stuff that has to be done only once per Evaluate()
    void PreEvaluate(
        Teuchos::ParameterList& params) override  ///< parameter list as handed in from the element
    {
      dserror("not implemented here");
    }
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
    InelasticDefgradPolynomial();

    /// construct material with specific material params
    explicit InelasticDefgradPolynomial(MAT::PAR::Parameter* params);

    /// destructor
    ~InelasticDefgradPolynomial() { ; }

    /// evaluate polynomial describing growth of material with regard to concentration
    double EvaluateMatPolynomial(const double Scalar, const int NumberPolyparams,
        const std::vector<double> ParamsPoly, const double X_min, const double X_max,
        double concentration, const double Cmax, double DetJacobian);

    /// evaluate the first derivative of polynomial
    double EvaluateMatPolynomialDerivative(const double Scalar, const int NumberPolyparams,
        const std::vector<double> ParamsPoly, double concentration, const double Cmax,
        double DetJacobian);

    /// provide material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      dserror("not implemented here");
      return INPAR::MAT::m_none;
    }

    /// evaluate the deformation gradient and its inverse
    void EvaluateInverseInelasticDefGrad(
        const LINALG::Matrix<3, 3>* defgrad,   ///< Deformation gradient
        LINALG::Matrix<3, 3>& iFinM) override  ///< Inverse inelastic deformation gradient
    {
      dserror("not implemented here");
    }

    /// evaluate additional terms for the elasticity tensor
    void EvaluateAdditionalCmat(
        const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 6>& cmatadd) override  ///< Additional elasticity tensor
    {
      dserror("not implemented here");
    }

    /// evaluate off-diagonal stiffness matrix for monolithic systems to get the
    /// cross-linearizations
    void EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                             ///< inelastic contribution as 3x3 matrix
        const LINALG::Matrix<6, 9>&
            dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                       ///< deformation gradient of current inelastic contribution
        LINALG::Matrix<6, 1>& dstressdx) override  ///< Derivative of 2nd Piola Kirchhoff stresses
                                                   ///< w.r.t. primary variable of different field
    {
      dserror("not implemented here");
    }

    /// pre-evaluation, intended to be used for stuff that has to be done only once per Evaluate()
    void PreEvaluate(
        Teuchos::ParameterList& params) override  ///< parameter list as handed in from the element
    {
      dserror("not implemented here");
    }

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

    /// destructor
    virtual ~InelasticDefgradPolyScalarIso() { ; }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mfi_poly_scalar_iso; }

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
    MAT::PAR::InelasticDefgradPolyScalarIso* Parameter()
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

    /// destructor
    virtual ~InelasticDefgradPolyScalarAniso() { ; }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mfi_poly_scalar_aniso; }

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
    MAT::PAR::InelasticDefgradPolyScalarAniso* Parameter()
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
