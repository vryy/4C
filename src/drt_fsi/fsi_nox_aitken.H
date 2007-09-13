#ifndef FSI_NOX_AITKEN_H
#define FSI_NOX_AITKEN_H

#ifdef TRILINOS_PACKAGE

#include <NOX_LineSearch_Generic.H> // base class
#include <NOX_LineSearch_UserDefinedFactory.H>
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <Teuchos_RefCountPtr.hpp>

namespace NOX {
  namespace FSI {

//! Aikten line search - the simple relaxation.
/*!
  This line search can be called via NOX::LineSearch::Manager.

  The working horse in FSI.

 */
class AitkenRelaxation : public NOX::LineSearch::Generic
{
public:

  //! Constructor
  AitkenRelaxation(const Teuchos::RefCountPtr<NOX::Utils>& utils,
                   Teuchos::ParameterList& params);

  //! Destructor
  ~AitkenRelaxation();

  // derived
  bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);

private:

  //! difference of last two solutions
  Teuchos::RefCountPtr<NOX::Abstract::Vector> del_;

  //! difference of difference of last two pair of solutions
  Teuchos::RefCountPtr<NOX::Abstract::Vector> del2_;

  //! aitken factor
  double nu_;

  //! Printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils_;
};


/// simple factory that creates aitken relaxation class
class AitkenFactory : public NOX::LineSearch::UserDefinedFactory
{
public:
  Teuchos::RCP<NOX::LineSearch::Generic>
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
		  Teuchos::ParameterList& params) const
    {
      return Teuchos::rcp(new AitkenRelaxation(gd->getUtils(),params));
    }
};

  }
}

#endif
#endif
