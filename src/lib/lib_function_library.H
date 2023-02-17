/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef LIB_FUNCTION_LIBRARY_H
#define LIB_FUNCTION_LIBRARY_H

#include "lib_function.H"
#include "lib_function_of_scalar.H"

namespace DRT::UTILS
{
  // forward declaration
  class CubicSplineInterpolation;

  /// add valid function lines
  void AddValidLibraryFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

  /// try to create functions from a given line definition
  Teuchos::RCP<FunctionOfScalar> TryCreateLibraryFunctionScalar(
      const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);


  /**
   * @brief special implementation of a 1D polynomial function
   *
   * This class defines a generic polynomial that that takes one argument. The polynomial can be
   * evaluated and the first derivative is also provided.
   *
   * This class does circumvent using automatic differentiation and may be used in
   * performance-sensitive applications (e.g. operations on Gauss-point level).
   * In particular, evaluating the derivative is faster for this function compared to
   * a SymbolicFunctionOfAnything representing the same polynomial.
   */
  class FastPolynomialFunction : public FunctionOfScalar
  {
   public:
    /**
     * ctor
     *
     * @param coefficients the coefficients of the monomials in ascending order
     */
    FastPolynomialFunction(std::vector<double> coefficients);

    /**
     * @brief Evaluate the polynomial.
     *
     * This is the only supported evaluation call for this class
     *
     * @param argument point to evaluate
     * @return value of the polynomial at `argument`
     */
    [[nodiscard]] double Evaluate(const double argument) const override;

    /**
     * Evaluate first derivative of polynomial.
     *
     * @param argument point to evaluate
     * @return value of first derivative of polynomial at `argument`
     */
    [[nodiscard]] double EvaluateDerivative(const double argument) const override;

   private:
    //! internal polynomial representation
    const DRT::UTILS::Polynomial mypoly_;
  };

  /**
   * @brief Cubic spline interpolation from csv file
   *
   * This class defines a cubic spline interpolation connecting \f$ n \f$ points
   * \f$ (x_i,y_i), i \in [0,1,...,n-1] \f$ expressed as splines \f$S_i(x)\f$ of the form:
   * \f[
   * S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3, j \in [0,1,...,n-1],
   * \f]
   * where \f$ S_i(x) \f$ is a third order polynomial on the interval \f$ [x_i, x_{i+1}] \f$.
   *
   * @note The spline is calculated such that the first and second derivatives match at the
   * intersection of two segments. Furthermore, so-called natural boundary conditions, i.e. the
   * second derivatives at the two boundaries are set to zero, are applied.
   */
  class CubicSplineFromCSV : public FunctionOfScalar
  {
   public:
    /*!
     * @brief Constructor of function defining a cubic spline interpolation created based on data
     * from csv-file
     *
     * @param[in] csv_file  absolute path to csv file
     */
    CubicSplineFromCSV(const std::string& csv_file);

    [[nodiscard]] double Evaluate(double scalar) const override;

    [[nodiscard]] double EvaluateDerivative(double scalar) const override;

   private:
    std::unique_ptr<CubicSplineInterpolation> cubic_spline_;
  };
}  // namespace DRT::UTILS

#endif