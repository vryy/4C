#ifndef FSI_NOX_EXTRAPOLATE_H
#define FSI_NOX_EXTRAPOLATE_H

#ifdef TRILINOS_PACKAGE

#include <NOX_LineSearch_Generic.H> // base class
#include <NOX_Utils.H>
#include <vector>

//! Extrapolate line search - not even on a straight line.
/*!
  This line search can be called via NOX::LineSearch::Manager.

  Do the extrapolation that both Deparis and Michler apply (at
  different points in the FSI iteration, though.) This is not a real
  line search as we leave the line and extrapolate the direction.

  The code does more harm than good, it seems. However it might be a
  starting point to implement the Deparis version more closely.

 */
class Extrapolate : public NOX::LineSearch::Generic
{
public:

  //! Constructor
  Extrapolate(const Teuchos::RefCountPtr<NOX::Utils>& utils,
              Teuchos::ParameterList& params);

  //! Destructor
  ~Extrapolate();

  // derived
  bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);

private:

  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > r_;
  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > x_;

  //! Printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils_;
};

#endif
#endif
