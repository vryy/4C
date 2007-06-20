/*----------------------------------------------------------------------*/
/*!
\file drt_function.H

\brief Managing and evaluating of spatial functions

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_dserror.H"

extern "C"
{
#include "../pss_full/pss_parser.h"
}


namespace DRT {

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

    /// construct spatial function from expression with given origin
    ExprFunction(char* buf, double x, double y, double z);

    /// clean up parse tree
    ~ExprFunction();

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x);

  private:

    double x_;
    double y_;
    double z_;

    /// expression syntax tree
    _ST_NODE* expr_;
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

    /// read the baci input file and set up all curves
    void ReadInput();

    /// return a specified function
    Function& Funct(int num) { return *functions_[num]; }

  private:

    /// private default constructor to disallow creation of instances
    FunctionManager() {}

    /// no copying
    FunctionManager(const FunctionManager&);

    /// no copying
    FunctionManager& operator=(const FunctionManager&);

    /// all curves defined
    std::vector<Teuchos::RefCountPtr<Function> > functions_;

    /// the single instance
    static FunctionManager instance_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const FunctionManager& manager);
  };
}

#endif
#endif
