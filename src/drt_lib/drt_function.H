/*----------------------------------------------------------------------*/
/*!
\file drt_function.H

\brief Managing and evaluating of spatial functions

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
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

  \brief abstract interface for spatial functions.

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

    \brief evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functionsb, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    virtual double Evaluate(int index, const double* x, double t, DRT::Discretization* dis) = 0;

    /*!

    \brief Return the number of components of this spatial function
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

    BeltramiUP(int mat_id, bool is_stationary);

    BeltramiUP( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
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



  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public Function
  {
  public:

    BeltramiGradU(int mat_id, bool is_stationary );

    BeltramiGradU( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
    bool is_stationary_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity, pressure)
  class KimMoinUP : public Function
  {
  public:

    KimMoinUP(int mat_id, bool is_stationary);

    KimMoinUP( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
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

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
    bool is_stationary_;
  };


}

}


void PrintFunctionDatHeader();


#endif

