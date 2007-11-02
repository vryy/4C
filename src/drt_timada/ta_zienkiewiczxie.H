/*======================================================================*/
/*!
\file ta_zienkiewiczxie.H

\class ZienkiewiczXie

\brief ZienkiewiczXie time step indicator for time adaptivity

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
#ifndef ZIENKIEWICZXIE_H
#define ZIENKIEWICZXIE_H

#include "timeadaptivity.H"

#include "../drt_structure/strugenalpha.H"

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
\class ZienkiewiczXie

\brief Time step size adaptivity with Zienkiewicz-Xie error indicator

\author bborn
\date 10/07
*/
class ZienkiewiczXie : public TimeAdaptivity
{

public:

   /*\!
     \brief Ordinary Constructor with parameters
   */
   explicit ZienkiewiczXie
   (
      double timeinitial,
      double timefinal,
      int timestepinitial,
      int timestepfinal,
      double stepsizeinitial,
      //
      double stepsizemax,
      double stepsizemin,
      double sizeratiomax,
      double sizeratiomin,
      double sizeratioscale,
      TAErrNorm errnorm,
      double errtol,
      int errorder,
      int adaptstepmax,
      //
      DRT::Discretization& dis,
      LINALG::Solver& solver,
      IO::DiscretizationWriter& output
   );

   /*!
     \brief Destructor
   */
   virtual ~ZienkiewiczXie();

   /*!
     \brief Integrate in time
   */
   void Integrate(StruGenAlpha::StruGenAlpha& timint);

protected:

   // not wanted: = operator
   ZienkiewiczXie operator = (const ZienkiewiczXie& old);
   // not wanted: copy constructor
   ZienkiewiczXie(const ZienkiewiczXie& old);

}; // class ZienkiewiczXie



#endif  // #ifndef ZIENKIEWICZXIE_H
#endif  // #ifdef CCADISCRET
