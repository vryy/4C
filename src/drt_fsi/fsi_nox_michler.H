#ifndef FSI_NOX_MICHLER_H
#define FSI_NOX_MICHLER_H

#ifdef TRILINOS_PACKAGE

#include <NOX_Direction_Generic.H> // base class
#include <NOX_Utils.H>
#include <Teuchos_RefCountPtr.hpp>
#include <vector>

namespace NOX {
  namespace FSI {

/// FSI coupling algorithm due to Michler
/*!
  The FSI coupling algorithm as suggested by Michler et al.
  This beast does not work. You need to specify a relaxation
  parameter. The algorithm is very sensitive with respect to that
  one. And furthermore the linear system of equations that needs to be
  solved after that is very badly conditioned.

  You are better of by far to use a simple fixed point scheme with
  relaxation. Aitken would be a good choice.
 */
class Michler : public NOX::Direction::Generic
{
public:

  //! Constructor
  Michler(const Teuchos::RefCountPtr<NOX::Utils>& utils,
              Teuchos::ParameterList& params);

  //! Destructor
  ~Michler();

  // derived
  bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
	       const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
		       const NOX::Solver::LineSearchBased& solver);

private:

  //! Print error message and throw error
  void throwError(const string& functionName, const string& errorMsg);

private:

  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > r_;
  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > x_;

  //! Printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils_;

  Teuchos::ParameterList& params_;

  double nu_;
};

  }
}

#endif
#endif
