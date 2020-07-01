/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

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
  /// special implementation for weakly compressible flow - Etienne FSI problem
  class WeaklyCompressibleEtienneFSIStructureFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIStructureFunction(int mat_id_struc);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    virtual int NumberComponents() { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };


  /// special implementation for weakly compressible flow - Etienne FSI problem (force)
  class WeaklyCompressibleEtienneFSIStructureForceFunction : public DRT::UTILS::Function
  {
   public:
    WeaklyCompressibleEtienneFSIStructureForceFunction(int mat_id_struc);
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(const int index, const double* x, double t);

    /*!

    \brief evaluate first or second time derivative

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the time derivative will
                     be evaluated
    \param t     (i) The point in time in which the time derivative will
                     be evaluated
    \param deg   (i) The degree of the time derivative

    */
    virtual std::vector<double> EvaluateTimeDerivative(const int index,  //< component of function
        const double* x,                                                 //< coordinates
        const double t,                                                  //< time
        const unsigned deg);  //< degree of maximal derivative

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components

    */
    virtual int NumberComponents() { return (2); };

   private:
    double youngmodulus_;
    double poissonratio_;
    double strucdensity_;
  };

}  // namespace STR

#endif
