/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "../drt_lib/drt_function.H"
#include "../drt_mat/stvenantkirchhoff.H"


#ifndef STR_FUNCTION_H
#define STR_FUNCTION_H

namespace DRT
{
  class Discretization;
  namespace INPUT
  {
    class LineDefinition;
  }  // namespace INPUT
}  // namespace DRT


namespace STR
{
  /// add valid structure-specific function lines
  void StructureValidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);

  /// try to create structure-specific functions from a given line definition
  Teuchos::RCP<DRT::UTILS::Function> StructureTryCreateFunction(
      Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def);

  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIStructureFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIStructureFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    /*!
     *  \brief evaluate function at given position in space
     *
     * \param index (i) index defines the function-component which will be evaluated
     * \param x     (i) The point in space in which the function will be evaluated
     */
    double Evaluate(const int index, const double* x, double t) override;

    /*!
     *  \brief evaluate first or second time derivative
     *
     * \param index (i) index defines the function-component which will be evaluated
     * \param x     (i) The point in space in which the time derivative will be evaluated
     * \param t     (i) The point in time in which the time derivative will be evaluated
     * \param deg   (i) The degree of the time derivative
     */
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!
     *  \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components
     */
    std::size_t NumberComponents() override { return (2); };

   private:
    double poissonratio_;
  };

  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIStructureForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIStructureForceFunction(const MAT::PAR::StVenantKirchhoff& fparams);

    /*!
     *  \brief evaluate function at given position in space
     *
     * \param index (i) index defines the function-component which will be evaluated
     * \param x     (i) The point in space in which the function will be evaluated
     */
    double Evaluate(const int index, const double* x, double t) override;

    /*!
     *  \brief evaluate first or second time derivative
     *
     * \param index (i) index defines the function-component which will be evaluated
     * \param x     (i) The point in space in which the time derivative will be evaluated
     * \param t     (i) The point in time in which the time derivative will be evaluated
     * \param deg   (i) The degree of the time derivative
     */
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, const double t, const unsigned deg) override;

    /*!
     *  \brief Return the number of components of this spatial function (This is a vector-valued
     * function)
     *
     * \return number of components
     */
    std::size_t NumberComponents() override { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

}  // namespace STR

#endif
