/*----------------------------------------------------------------------*/
/*!
\file drt_function.H

\brief Declaration Managing and evaluating of spatial functions
<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_FUNCTION_H
#define DRT_FUNCTION_H


#include <iostream>
#include <vector>
#include <complex>
#include <Teuchos_RCP.hpp>

#include "drt_dserror.H"
#include "drt_parser.H"
//#include "drt_inputreader.H"

namespace MAT {
class Material;
}

namespace DRT {

  class Discretization;
  class Condition;
  class Problem;


namespace INPUT {

  class Lines;
  class DatFileReader;

}

namespace UTILS {


  /*!

  \brief Declaration abstract interface for spatial functions.

  Spatial functions can be defined either scalar-valued or vector-valued.

  There exist function definitions using explicitly defined spatial functions
  for test problems (beltrami, kim-moin etc).

  Spatial function based on parsed expressions are provided as well.

  */
  class Function
  {
  public:
    virtual ~Function() {}

    /*!

    \brief Declaration evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functionsb, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    virtual double Evaluate(int index, const double* x, double t, DRT::Discretization* dis) = 0;

    /*!

    \brief Declaration evaluate derivatives function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functionsb, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated
    \param t     (i) Absolut time in which the
                     function will be evaluated
    \param dis   (i) discretization

    */
    virtual std::vector<std::vector<double> > FctDer(int index, const double* x, const double t, DRT::Discretization* dis)
      {
        return std::vector<std::vector<double> >(0);
      };

    /*!

    \brief Declaration Return the number of components of this spatial function
    (1 for scalar functions, dim for vector-valued functions)

    \return number of components

    */
    virtual int NumberComponents()
      {
        return(-1);
      };

  };


  /// class that manages all defined spatial functions
  class FunctionManager
  {
  public:

    FunctionManager() {}

    Teuchos::RCP<DRT::INPUT::Lines> ValidFunctionLines();

    /// read the baci input file and set up all functions
    void ReadInput(DRT::INPUT::DatFileReader& reader);

    /// return a specified function
    Function& Funct(int num);

  private:

    /// no copying
    FunctionManager(const FunctionManager&);

    /// no copying
    FunctionManager& operator=(const FunctionManager&);

    /// all functions defined. The functions may be vector valued!
    std::vector<Teuchos::RCP<Function> > functions_;
  };


  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public Function
  {
  public:

    BeltramiUP(int mat_id);

    BeltramiUP( Teuchos::RCP<MAT::Material> & mat );


    /*!

      \brief Declaration evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Declaration Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents()
    {
      return(4);
    };

  private:
    double density_;
    double kinviscosity_;
  };



  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public Function
  {
  public:

    BeltramiGradU(int mat_id );

    BeltramiGradU( Teuchos::RCP<MAT::Material> & mat );


    /*!

      \brief Declaration evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Declaration Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity, pressure)
  class KimMoinUP : public Function
  {
  public:

    KimMoinUP(int mat_id, bool is_stationary);

    KimMoinUP( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief Declaration evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Declaration Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents()
    {
      return(4);
    };

  private:
    double density_;
    double kinviscosity_;
    bool   is_stationary_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of velocity)
  class KimMoinGradU : public Function
  {
  public:

    KimMoinGradU(int mat_id, bool is_stationary );

    KimMoinGradU( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief Declaration evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Declaration Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
    bool   is_stationary_;
  };


}

}


void PrintFunctionDatHeader();


#endif

