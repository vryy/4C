/*======================================================================*/
/*!
\file timeadaptivity.H

\class TimeAdaptivity

\brief Time step size adaptivity

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef TIMEADAPTIVITY_H
#define TIMEADAPTIVITY_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../io/io_drt.H"
#include "../drt_contact/drt_contact_manager.H"
#include "../drt_lib/linalg_mfoperator.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

/*----------------------------------------------------------------------*/
/* namespaces */
using namespace std;
using namespace Teuchos;

/*======================================================================*/
/*!
\class TimeAdaptivity

\brief Time step size adaptivity

\author bborn
\date 10/07
*/
class TimeAdaptivity {

public:

   /*!
     \brief Standard Constructor
   */
   explicit TimeAdaptivity(
      ParameterList& params,
      DRT::Discretization& dis,
      LINALG::Solver& solver,
      IO::DiscretizationWriter& output
   );

   /*!
     \brief Destructor
   */
   virtual ~TimeAdaptivity() {};

  /*!
  \brief Create default parameter list
  */
  static void SetDefaults(ParameterList& params);

  /*!
  \brief Do time integration
  */
//  void Integrate();

  /*!
  \brief do constant predictor step
  */
//  void ConstantPredictor();

  /*!
  \brief do matrix free constant predictor
  */
//  void MatrixFreeConstantPredictor();

  /*!
  \brief do consistent constant predictor step
  */
//  void ConsistentPredictor();

  /*!
  \brief do full Newton iteration
  */
//  void FullNewton();

  /*!
  \brief do modified Newton iteration (makes sense with direct solvers only)
  */
//  void ModifiedNewton();

  /*!
  \brief do matrix free Newton iteration (only with Aztec)
  */
//  void MatrixFreeNewton();

  /*!
  \brief do nonlinear cg iteration using nox
  */
//  void NonlinearCG();

  /*!
  \brief Compute the residual to a given state x

  This method is used by the nox user interface class and by the linalg matrix free operator

  \param x (in/out) : Current solution state (dirichlet dofs will be blanked)
  \param F (out): residual matching current solution state (dirichlet dofs will be blanked)
  */
//  virtual void computeF(const Epetra_Vector& x, Epetra_Vector& F);


  /*!
  \brief Compute the Jacobian to a given state x

  This method is used by the nox user interface

  \param x (in/out) : Current solution state (dirichlet dofs will be blanked)
  */
//  virtual void computeJacobian(const Epetra_Vector& x);



protected:

    // not wanted: = operator
    TimeAdaptivity operator = (const TimeAdaptivity& old);
    // not wanted: copy constructor
    TimeAdaptivity(const TimeAdaptivity& old);

  ParameterList&                params_;
  DRT::Discretization&          discret_;
  LINALG::Solver&               solver_;
  IO::DiscretizationWriter&     output_;
  int                           myrank_;
  int                           maxentriesperrow_;
  double                        norm_;

  bool                          havecontact_;
  RCP<CONTACT::Manager>         contactmanager_;

  RefCountPtr<Epetra_CrsMatrix> stiff_;
  RefCountPtr<Epetra_CrsMatrix> mass_;
  RefCountPtr<Epetra_CrsMatrix> damp_;

  RefCountPtr<Epetra_Vector>    dirichtoggle_;
  RefCountPtr<Epetra_Vector>    invtoggle_;
  RefCountPtr<Epetra_Vector>    zeros_;
  RefCountPtr<Epetra_Vector>    dis_;
  RefCountPtr<Epetra_Vector>    vel_;
  RefCountPtr<Epetra_Vector>    acc_;
  RefCountPtr<Epetra_Vector>    disn_;
  RefCountPtr<Epetra_Vector>    veln_;
  RefCountPtr<Epetra_Vector>    accn_;
  RefCountPtr<Epetra_Vector>    dism_;
  RefCountPtr<Epetra_Vector>    velm_;
  RefCountPtr<Epetra_Vector>    accm_;
  RefCountPtr<Epetra_Vector>    disi_;
  RefCountPtr<Epetra_Vector>    fint_;
  RefCountPtr<Epetra_Vector>    fext_;
  RefCountPtr<Epetra_Vector>    fextm_;
  RefCountPtr<Epetra_Vector>    fextn_;
  RefCountPtr<Epetra_Vector>    fresm_;

}; // class TimeAdaptivity



#endif  // #ifndef TIMEADAPTIVITY_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
