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

#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_dserror.H"

extern "C"
{
#include "../pss_full/pss_parser.h"
}


namespace DRT {
namespace Utils {

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
    Function& Funct(int num) { return *(functions_[num]); }

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
