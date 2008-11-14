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

#ifdef CCADISCRET

#ifndef DRT_FUNCTION_H
#define DRT_FUNCTION_H

#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_dserror.H"
#ifdef PARALLEL
#include "mpi.h"
#endif
extern "C"
{
#include "../pss_full/pss_parser.h"
}


namespace DRT {
namespace UTILS {

  /// abstract interface for spatial functions
  class Function
  {
  public:
    virtual ~Function() {}

    /// evaluate function at given position in space
    virtual double Evaluate(int index, const double* x) = 0;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const Function& funct);
  };


  /// spatial function based on parsed expression
  class ExprFunction : public Function
  {
  public:
    /// construct function from expression with given origin
    /// (one function for all dimensions)
    ExprFunction(char* buf, double x, double y, double z);


    /// Default constructor creating empty object. Expressions are added with add function
    ExprFunction();

    /// clean up parse tree
    ~ExprFunction();

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x);

    /// add expression
    void AddExpr(char* buf, double x, double y, double z);

  private:

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;

    /// expression syntax tree
    std::vector<_ST_NODE*> expr_;
  };


  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public Function
  {
  public:
    /// evaluate function at given position in space
    double Evaluate(int index, const double* x);
  };


  /// special implementation for 2d Kim-Moin flow
  class KimMoinFunction : public Function
  {
  public:
    /// evaluate function at given position in space
    double Evaluate(int index, const double* x);
  };


  /// Singleton class that manages all defined spatial functions
  class FunctionManager
  {
  public:

    /// return the single instance we have of this class
    static FunctionManager& Instance() { return instance_; }

    /// read the baci input file and set up all functions
    void ReadInput();

    /// return a specified function
    Function& Funct(int num);



  private:

    /// private default constructor to disallow creation of instances
    FunctionManager() {}

    /// no copying
    FunctionManager(const FunctionManager&);

    /// no copying
    FunctionManager& operator=(const FunctionManager&);

    /// all functions defined. The functions may be vector valued!
    std::vector<Teuchos::RefCountPtr<Function> > functions_;

    /// the single instance
    static FunctionManager instance_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const FunctionManager& manager);
  };
}
}

#endif

#endif
