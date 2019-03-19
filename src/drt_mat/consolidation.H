/*----------------------------------------------------------------------*/
/*!
 \file consolidation.H
 \brief material managing consolidation or other phase change

 \level 3
 <pre>
 \maintainer Sebastian Pröll
 </pre>
 */

/*----------------------------------------------------------------------*
 |  definitions                                                         |
 *----------------------------------------------------------------------*/
#ifndef CONSOLIDATION_H_
#define CONSOLIDATION_H_

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"

#include "../drt_lib/drt_parobjectfactory.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for consolidation material
    ///
    /// <h3>Input line</h3>
    /// MAT 2 MAT_Consolidation SOLIDUS 1615 LIQUIDUS 1674 DELTA 2 LATENTHEAT 5e6
    class Consolidation : public Parameter
    {
     public:
      /// standard constructor
      Consolidation(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Consolidation() { ; }

      /// @name material parameters
      //@{
      const double solidustemp_;
      const double liquidustemp_;
      /// smoothing width, if 0 no smoothing but linear interpolation
      const double delta_;
      const double latentheat_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// if true use heat integration instead of apparent capacity
      bool heatint_;
      /// enthalpy dependent parameters instead of temperature
      bool enthalpydep_;
    };

  }  // namespace PAR

  class ConsolidationType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ConsolidationType"; }

    static ConsolidationType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ConsolidationType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Consolidation material
  ///
  /// Manager for consolidation process which provides a common storage and evaluation of
  /// consolidation related stuff for SLM.
  /// This material needs to be referenced by material definitions that use consolidation.
  ///
  /// \author proell
  /// \date 05/18
  class Consolidation : public Material
  {
   public:
    /// interface of a general interpolant
    class Interpolant
    {
     public:
      /// evaluate interpolation at point x
      virtual double valueAt(double x) = 0;
      /// derivative of interpoaltion at point x
      virtual double derivativeAt(double x) = 0;
      /// inverse of interpolation, i.e., find
      virtual double inverseAt(double x) = 0;
    };

    /// linear interpolant
    class LinearInterpolant : public Interpolant
    {
     public:
      /// ctor
      LinearInterpolant();
      /// ctor
      LinearInterpolant(double x1, double x2, double y1, double y2)
          : x1_(x1), x2_(x2), y1_(y1), y2_(y2)
      {
      }

      double valueAt(double x)
      {
        if (x <= x1_)
          return y1_;
        else if (x >= x2_)
          return y2_;
        else
        {
          const double t = (x - x1_) / (x2_ - x1_);
          return t * y2_ + (1 - t) * y1_;
        }
      }

      double derivativeAt(double x)
      {
        if (x <= x1_)
          return 0;
        else if (x >= x2_)
          return 0;
        else
        {
          return (y2_ - y1_) / (x2_ - x1_);
        }
      }

      double inverseAt(double y)
      {
        dserror("Inverse not implemented.");
        return -1;
      }

     private:
      // interpolation points
      const double x1_, x2_, y1_, y2_;
    };

    /// Smoothed linear interpolation
    ///
    /// Construct and manage a quadratic-linear-quadratic approximation to a linear curve.
    /// Value and derivatives at x1,x2 are continuous. Outside of [x1,x2] linear with given slope.
    class QuadLinQuadInterpolant : public Interpolant
    {
     public:
      QuadLinQuadInterpolant();

      // create a new interpolation from two points with slope
      QuadLinQuadInterpolant(
          double x1, double x2, double y1, double y2, double k1, double k2, double delta)
          : x1_(x1), x2_(x2), y1_(y1), y2_(y2), k1_(k1), k2_(k2), delta_(delta)
      {
        if (delta <= 0) dserror("delta=%f must be positive.", delta);
        if (x1_ > x2_) dserror("x2 must be greater x1.");
        if (x2_ - x1_ < 2 * delta) dserror("x1 and x2 within 2*delta.");

        if (abs(y1_ - y2_) < 1e-8)
        {
          // connect with constant
          a1_ = 0;
          a2_ = 0;
          a3_ = 0;
          b3_ = y1;
          b1_ = 0;
          b2_ = 0;
          c1_ = y1;
          c2_ = y1;
        }
        else
        {
          // coefficients generated with Matlab
          a1_ = (2 * y1 - 2 * y2 + delta * (k1 + k2)) / (4 * delta * (delta + x1 - x2)) -
                k1 / (2 * delta);
          a2_ = k2 / (2 * delta) -
                (2 * y1 - 2 * y2 + delta * (k1 + k2)) / (4 * delta * (delta + x1 - x2));
          //        sx1_ = x1 - (2*delta*k1*(delta + x1 - x2))/(2*y1 - 2*y2 - delta*k1 + delta*k2 -
          //        2*k1*x1 + 2*k1*x2); sx2_ = x2 + (2*delta*k2*(delta + x1 - x2))/(2*y1 - 2*y2 +
          //        delta*k1 - delta*k2 - 2*k2*x1 + 2*k2*x2); sy1_ = y1 - (delta*k1*k1*(delta + x1 -
          //        x2))/(2*y1 - 2*y2 - delta*k1 + delta*k2 - 2*k1*x1 + 2*k1*x2); sy2_ = y2 +
          //        (delta*k2*k2*(delta + x1 - x2))/(2*y1 - 2*y2 + delta*k1 - delta*k2 - 2*k2*x1 +
          //        2*k2*x2);
          a3_ = (y1 - y2 + (delta * (k1 + k2)) / 2) / (delta + x1 - x2);
          b3_ = y1 + (delta * k1) / 2 -
                ((delta + 2 * x1) * (2 * y1 - 2 * y2 + delta * k1 + delta * k2)) /
                    (4 * (delta + x1 - x2));

          // new coeffs
          const double d = delta;
          c1_ = ((x1 * x1 * (2 * y1 - 2 * y2)) / 4 + (delta * x1 * x1 * (k1 + k2)) / 4) /
                    (delta * (delta + x1 - x2)) -
                (k1 * x1 * x1 + 2 * delta * k1 * x1 - 2 * delta * y1) / (2 * delta);
          c2_ = (k2 * x2 * x2 - 2 * d * k2 * x2 + 2 * d * y2) / (2 * d) -
                ((x2 * x2 * (2 * y1 - 2 * y2)) / 4 + (d * x2 * x2 * (k1 + k2)) / 4) /
                    (d * (d + x1 - x2));
          b1_ = (k1 * (d + x1)) / d -
                ((x1 * (2 * y1 - 2 * y2)) / 2 + (d * x1 * (k1 + k2)) / 2) / (d * (d + x1 - x2));
          b2_ = (k2 * (d - x2)) / d +
                ((x2 * (2 * y1 - 2 * y2)) / 2 + (d * x2 * (k1 + k2)) / 2) / (d * (d + x1 - x2));
        }
      }

      // evaluate interpolator at given coordinate
      double valueAt(double x)
      {
        double val = -1;
        if ((x >= x1_ + delta_) && (x <= x2_ - delta_))
          val = a3_ * x + b3_;
        else if ((x >= x1_) && (x < x1_ + delta_))
          //        val = a1_*(x-sx1_)*(x-sx1_)+sy1_;
          val = a1_ * x * x + b1_ * x + c1_;
        else if ((x > x2_ - delta_) && (x <= x2_))
          //        val = a2_*(x-sx2_)*(x-sx2_)+sy2_;
          val = a2_ * x * x + b2_ * x + c2_;
        else if (x < x1_)
          val = k1_ * (x - x1_) + y1_;
        else if (x > x2_)
          val = k2_ * (x - x2_) + y2_;
        else
          dserror("implementation error");
        //      if((val <= y2_+TOL && val>=y1_+TOL) or(val <= y1_+TOL && val >= y2_+TOL)){
        return val;
        //      }
        //      else
        //        dserror("x=%f: y=%f out of range [%f,%f]",x,val,y1_,y2_);
      }

      // evaluate derivative of interpolator at given coordinate
      double derivativeAt(double x)
      {
        if ((x >= x1_ + delta_) && (x <= x2_ - delta_))
          return a3_;
        else if ((x >= x1_) && (x < x1_ + delta_))
          //          return 2*a1_*(x-sx1_);
          return 2 * a1_ * x + b1_;
        else if ((x > x2_ - delta_) && (x <= x2_))
          //        return 2*a2_*(x-sx2_);
          return 2 * a2_ * x + b2_;
        else if (x < x1_)
          return k1_;
        else if (x > x2_)
          return k2_;
        else
          dsassert(false, "implementation error");
        return -1;
      }

      // inverse of interpolation, only works reliably for k1=k2=0
      double inverseAt(double y)
      {
        const double yc1 = a3_ * (x1_ + delta_) + b3_;
        const double yc2 = a3_ * (x2_ - delta_) + b3_;
        // safety to avoid getting imaginary parts
        const double tol = 1e-6;
        if (abs(y - y1_) < tol)
          return x1_;
        else if (abs(y - y2_) < tol)
          return x2_;
        else if (yc1 <= y && y <= yc2)
          return (y - b3_) / a3_;
        else if (y >= y1_ && y < yc1)
          //        return sqrt((y-sy1_)/a1_)+sx1_;
          return (-b1_ + sqrt(b1_ * b1_ - 4 * a1_ * (c1_ - y))) / (2 * a1_);
        else if (y > yc2 && y <= y2_)
          //        return -sqrt((y-sy2_)/a2_)+sx2_;
          return (-b2_ + sqrt(b2_ * b2_ - 4 * a2_ * (c2_ - y))) / (2 * a2_);
        else
          dserror("y out of range [%f,%f]", y1_, y2_);
        return -1;
      }

     private:
      // input points
      const double x1_, x2_;
      const double y1_, y2_;
      const double k1_, k2_;
      const double delta_;
      // coefficients
      double a1_, a2_, a3_;
      double c1_, c2_, b1_, b2_, b3_;

      const double TOL = 1e-6;
    };

    using InterpolantPtr = Teuchos::RCP<MAT::Consolidation::Interpolant>;
    using InterpolantCache = std::map<const std::vector<int>, InterpolantPtr>;

    /// empty constructor
    Consolidation();

    /// constructor with given material parameters
    Consolidation(MAT::PAR::Consolidation* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return ConsolidationType::Instance().UniqueParObjectId();
    }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    virtual void Pack(DRT::PackBuffer& data  ///< (in/out): char vector to store class information
        ) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    virtual void Unpack(
        const std::vector<char>& data  ///< vector storing all data to be unpacked into this
    );

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_consolidation; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Consolidation(*this)); }

    //! setup routine to create history variables based on number of Gauss points
    void Setup(const int numgp);

    //! update history variables
    void Update();

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// @name evaluation of functions and their derivatives based on phase and temperature
    //! @{

    /// \brief evaluation routine for function evaluation based on temperature and consolidation
    ///
    /// This is the main method to be used by classes holding an instance of Consolidation.
    /// Based on temperature and GP history the two functions are interpolated.
    double EvaluateTempDependentFunction(const double temperature,  ///< temperature at GP
        const int gp,  ///< current GP to select internal phase information
        const std::vector<int>
            functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
    );

    /// \brief evaluation routine for derivative evaluation based on temperature and consolidation
    ///
    /// This is the main method to be used by classes holding an instance of Consolidation.
    /// Based on temperature and GP history the correct derivative of the interpolation is
    /// calculated.
    double EvaluateTempDependentDerivative(const double temperature,  ///< temperature at GP
        const int gp,  ///< current GP to select internal phase information
        const std::vector<int>
            functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
    );

    /// @}

    /**
     * Set value of shapefunctions
     * @param shapefunct
     */
    template <unsigned int nen>
    void SetFunct(const LINALG::Matrix<nen, 1>& shapefunct)
    {
      enth_ = params_->latentheat_;
      for (unsigned i = 0; i < nen; i++) enth_ += shapefunct(i) * lhremnp_->at(i);
    }

    // new wrapper until refactored
    double EvaluateFunction(const double temperature,  ///< temperature at GP
        const int gp,  ///< current GP to select internal phase information
        const std::vector<int>
            functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
    )
    {
      // TODO this should be solved via templates (preferably) or polymorhpism
      if (params_->enthalpydep_)
      {
        return EvaluateTempDependentFunction(enth_, gp, functions);
      }
      else
      {
        return EvaluateTempDependentFunction(temperature, gp, functions);
      }
    }

    // new wrapper
    double EvaluateDerivative(const double temperature,  ///< temperature at GP
        const int gp,  ///< current GP to select internal phase information
        const std::vector<int>
            functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
    )
    {
      if (params_->enthalpydep_)
      {
        // for now just skip derivative evaluation, not clear how temp dependency should be done
        // with enthalpy
        return 0;
      }
      else
      {
        return EvaluateTempDependentDerivative(temperature, gp, functions);
      }
    }

    //! get the middle temperature between solidus and liquidus
    double MeltingTemperature() { return (params_->solidustemp_ + params_->liquidustemp_) / 2; }

    //! get the value of latent heat distributed on the capacity
    double ApparentCapacityPeak(double temperature)
    {
      if (params_->heatint_)
        // currently apparent capacity not allowed in this scenario
        return 0;
      if (temperature <= Ts_mod_ or temperature >= Tl_mod_)
        return 0;
      else if (temperature > Ts_mod_ and temperature <= (Ts_mod_ + Tl_mod_) / 2)
      {
        return leftInterpAppCapa_->valueAt(temperature);
      }
      else if (temperature > (Ts_mod_ + Tl_mod_) / 2 and temperature < Tl_mod_)
      {
        return rightInterpAppCapa_->valueAt(temperature);
      }
      else
      {
        dserror("Temperature %f out of range", temperature);
        return -1;
      }
    }
    //! linearization contribution of latent heat distributed on the capacity
    double ApparentCapacityDerivative(double temperature)
    {
      if (params_->heatint_)
        // currently apparent capacity not allowed in this scenario
        return 0;
      if (temperature <= Ts_mod_ or temperature >= Tl_mod_)
        return 0;
      else if (temperature > Ts_mod_ and temperature <= (Ts_mod_ + Tl_mod_) / 2)
      {
        return leftInterpAppCapa_->derivativeAt(temperature);
      }
      else if (temperature > (Ts_mod_ + Tl_mod_) / 2 and temperature < Tl_mod_)
      {
        return rightInterpAppCapa_->derivativeAt(temperature);
      }
      else
      {
        dserror("Temperature %f out of range", temperature);
        return -1;
      }
    }

    /// \brief limit the latent heat source term based on history at nodes
    ///
    /// \note this method is only implemented for melting
    template <unsigned int nen>
    void LimitLatentHeatSourceTerm(
        LINALG::Matrix<nen, 1>&
            nodalsourceinc,  //< raw increment without history, (return) limited increment
        LINALG::Matrix<nen, 1>& nodaltemp  // nodal temperature, (return) reset if source term
    )
    {
      for (unsigned i = 0; i < nen; i++)
      {
        // take the full increment or the available rest
        nodalsourceinc(i) = std::max(lhremnp_->at(i), nodalsourceinc(i));
        // limit by total
        nodalsourceinc(i) = std::min(lhremnp_->at(i) + params_->latentheat_, nodalsourceinc(i));
        lhremnp_->at(i) -= nodalsourceinc(i);
        // finally scale by lumped factor
        nodalsourceinc(i) *= nodallump_->at(i);
        // perform nodal reset if there is a source term
        if (abs(nodalsourceinc(i)) > 1e-8) nodaltemp(i) = MeltingTemperature();
      }
    }

    //! save the nodally lumped factors and reset the latent heat
    template <unsigned int nen>
    void SetTotalLatentHeatAvailable(LINALG::Matrix<nen, 1>& nodalfactor)
    {
      // nodalfactor.Scale(params_->latentheat_);
      for (unsigned i = 0; i < nen; i++)
      {
        nodallump_->at(i) = nodalfactor(i);
        lhrem_->at(i) = -params_->latentheat_;
        lhremnp_->at(i) = -params_->latentheat_;
      }
    }

    //! reset memory of
    void ResetLatentHeatLastStep() { *lhremnp_ = *lhrem_; }

    /// service function for postprocessing, evaluate "phase"
    /// at nodes based on history and temperature
    ///
    /// prerequisite: Update() called in this step
    template <unsigned int nen>
    void NodalPhase(LINALG::Matrix<nen, 1>& nodaltemp,  ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalphase              ///< (return)
    )
    {
      for (unsigned int i = 0; i < nen; i++)
      {
        // overwrite the nodal history and pass as output
        cfracnode_->at(i) = CFracNext(nodaltemp(i), cfracnode_->at(i));
        nodalphase(i) = cfracnode_->at(i);
      }
    }

    /// service function for postprocessing, arbitrary functions
    /// at nodes based on history and temperature
    /// prerequisite: Update() called in this step
    template <unsigned int nen>
    void EvaluateTempDependentFunctionAtNodes(
        LINALG::Matrix<nen, 1>& nodaltemp,   ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalvalue,  ///< (return) function interpolated at nodes
        std::vector<int>
            functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
    )
    {
      if (functions.size() != NUMFUNCT)
        dserror("Required %d functions but only got %d.", NUMFUNCT, functions.size());

      for (unsigned int i = 0; i < nen; i++)
      {
        // update the nodal history
        cfracnode_->at(i) = CFracNext(nodaltemp(i), cfracnode_->at(i));
        // fixme should be checked only once and called polymorphically?
        if (functions[0] == functions[1] && functions[0] == functions[2])
        {
          // this call currently ensures that volumetric heat capacity can have a peak in the middle
          nodalvalue(i) = FunctionValue(nodaltemp(i), functions[0]);
        }
        else if (nodaltemp(i) > Tl_mod_ || nodaltemp(i) < Ts_mod_)
        {
          EvaluateFractions(nodaltemp(i), cfracnode_->at(i));
          const double powderVal = FunctionValue(nodaltemp(i), functions[0]);
          const double meltVal = FunctionValue(nodaltemp(i), functions[1]);
          const double solidVal = FunctionValue(nodaltemp(i), functions[2]);
          nodalvalue(i) =
              powderVal * frac_->at(0) + meltVal * frac_->at(1) + solidVal * frac_->at(2);
        }
        else
        {
          if (HasSmoothing())
          {
            const double rc = cfracnode_->at(i);
            nodalvalue(i) = GetCurrentInterpolant(rc, functions)->valueAt(nodaltemp(i));
          }
          else
          {
            EvaluateFractions(nodaltemp(i), cfracnode_->at(i));
            const double powderVal = FunctionValue(nodaltemp(i), functions[0]);
            const double meltVal = FunctionValue(nodaltemp(i), functions[1]);
            const double solidVal = FunctionValue(nodaltemp(i), functions[2]);
            nodalvalue(i) =
                powderVal * frac_->at(0) + meltVal * frac_->at(1) + solidVal * frac_->at(2);
          }
        }
      }
    }

   private:
    /// number of functions used for interpolation
    ///
    /// indexing is 0: powder, 1: melt, 2: solid
    static const int NUMFUNCT = 3;

    /// internal evaluation for phase fractions based on current temperature and phase
    void EvaluateFractions(const double temperature, const double cfrac);

    /// internal evaluation of fraction derivatives
    void EvaluateFractionDerivatives(
        const double temperature, const double cfrac, const double cfracn);

    void EvaluateCFracnpAtGp(const double temperature, const int gp);

    /// evaluate derivative of consolidated fraction with respect to temperature, for linearization
    double CFrac_T(const double temperature,  ///< temperature at GP
        const int gp);                        ///< current GP

    /// helper evaluating the definition of cfrac, alsos ets fractional temperature
    double CFracNext(const double temperature, const double oldcfrac);

    /// evaluate between Ts and Tl using a smoothed interpolation
    double EvaluateSmoothTransition(
        const double temperature, const int gp, const std::vector<int> functions);

    /// evaluate derivative between Ts and Tl of smoothed interpolation
    double EvaluateSmoothDerivative(
        const double temperature, const int gp, const std::vector<int> functions);

    //! evaluate relative temperature between T1 and T1 and store
    void EvaluateTFrac(const double temperature);

    //! evaluate global functions by IDs for given temperature
    double FunctionValue(const double temperature, const int function);

    //! evaluate derivatives of global functions by IDs for given temperature
    double FunctionDerivative(const double temperature, const int function);

    //! \return vec1*vec2
    double ScalarProduct(std::vector<double>* vec1, std::vector<double>* vec2);

    /// internal temporary variables initialization
    void SetupFunctionEvaluateVariables();

    //! \return a if 0 <= a <= 1, 0 if a < 0, 1 if a > 1
    inline double limiter(const double a) const { return std::min(std::max(a, 0.0), 1.0); }

    /// select the interpolant based on the current consolidated fraction
    /// prerequisite: EvaluateCFracnpAtGp() called
    InterpolantPtr GetCurrentInterpolant(const double rc, const std::vector<int> functions);

    /// return the (cached) interpolant smoothly connecting all-powder with all-melt
    InterpolantPtr GetOrCreateBaseInterpolant(const std::vector<int> functions);

    /// create the interpolant smoothly connecting all-powder with all-melt
    InterpolantPtr CreateBaseInterpolant(const std::vector<int>& functions);

    /// Check if interpolation is linear or smoothed.
    ///
    /// \todo this switching should be done with polymorphism and different evaluation strategies
    bool HasSmoothing() { return (params_->delta_ > 1e-12); }

    /// my material parameters
    MAT::PAR::Consolidation* params_;

    /// consolidated phase fraction last step at GP
    Teuchos::RCP<std::vector<double>> cfracn_;
    /// save consolidated phase fraction at GP for update after convergence
    Teuchos::RCP<std::vector<double>> cfracnp_;

    /// nodally lumped N . N^T matrix
    std::vector<double>* nodallump_;
    /// remaining latent heat
    std::vector<double>* lhrem_;
    std::vector<double>* lhremnp_;

    /// consolidated fraction at nodes in last step, only for postproc
    Teuchos::RCP<std::vector<double>> cfracnode_;

    /// temperature fraction between Ts_mod and Tl_mod
    double Tfrac_;
    /// modified solidus and liquidus temperatures
    double Ts_mod_, Tl_mod_;

    /// (non-history) variables for fractions, functions and their derivatives
    /// avoids passing around vectors all the time
    std::vector<double>* frac_;
    std::vector<double>* frac_T_;

    /// initial interpolation curve from powder to melt, not dependent on GP
    /// since we might call this class with different functions, we cache all calculated
    /// interpolants
    InterpolantCache baseInterpolantCache;

    /// interpolates from 0 to 1 between [Ts-delta/2, Tl+delta/2]
    InterpolantPtr fractionInterpolant;

    /// apparent capacity interpolants
    InterpolantPtr leftInterpAppCapa_;
    InterpolantPtr rightInterpAppCapa_;

    /// enthalpy at current gp
    double enth_;


  };  // Consolidation

}  // namespace MAT
   // namespace MAT
   // namespace MAT
/*----------------------------------------------------------------------*/
#endif  // CONSOLIDATION_H_
