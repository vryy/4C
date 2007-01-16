#ifndef FSI_AITKEN_NOX_H
#define FSI_AITKEN_NOX_H

#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#include <NOX_Direction_Generic.H>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class Aitken : public NOX::Direction::Generic
{
  Teuchos::RefCountPtr<NOX::Abstract::Vector> del_;
  Teuchos::RefCountPtr<NOX::Abstract::Vector> del2_;
  double nu_;

public:
  Aitken() : nu_(0) {}

  bool reset(const Teuchos::RefCountPtr< NOX::GlobalData > &gd, Teuchos::ParameterList &params);
  bool compute(NOX::Abstract::Vector &dir, NOX::Abstract::Group &grp, const NOX::Solver::Generic &solver);
};

#endif
#endif
#endif
