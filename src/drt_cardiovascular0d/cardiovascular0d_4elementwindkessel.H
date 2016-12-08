/*!----------------------------------------------------------------------
\file cardiovascular0d_4elementwindkessel.H

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

A four-element windkessel model only (no valves!) (DESIGN SURF CARDIOVASCULAR 0D WINDKESSEL ONLY CONDITIONS):
orig. equation is C * dp/dt + (p - p_ref)/R_p - (1 + Z_c/R_p) * Q - (C Z_c  + L/R_p) * dQ/dt - L C * d2Q/dt2 = 0
where Q := -dV/dt

-> we reformulate the ODE to three first-order ODEs with variables p, q, s:

      [C * dp/dt + (p - p_ref)/R_p + (1 + Z_c/R_p) * q + (C Z_c  + L/R_p) * s + L C * ds/dt]   [ 0 ]
Res = [dV/dt - q                                                                           ] = [ 0 ]
      [dq/dt - s                                                                           ]   [ 0 ]

The classical 3- or 2-element windkessel models are reproduced by setting L, or L and Z_c to zero, respectively

*----------------------------------------------------------------------*/

#ifndef CARDIOVASCULAR0D4ELEMENTWINDKESSEL_H
#define CARDIOVASCULAR0D4ELEMENTWINDKESSEL_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>
#include "Epetra_FECrsMatrix.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_cardiovascular0d.H"
// class which is inherited from
#include "cardiovascular0d.H"


//forward declarations
namespace DRT
{
  class Condition;
  class Discretization;
}
namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}


namespace UTILS
{

  class Cardiovascular0D4ElementWindkessel : public Cardiovascular0D

  {
  public:

    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0D4ElementWindkessel
    (
      Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization where Cardiovascular0D lives on
      const std::string& CondName, ///< Name of condition to create Cardiovascular0D from
      std::vector<int>& curID ///< current ID
    );

    /*!
        \brief Destructor

     */
    virtual ~Cardiovascular0D4ElementWindkessel()
    {
      return;
    };

    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    virtual void Initialize
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1,  ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2  ///< distributed vector that may be filled by assembly of element contributions
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateCardiovascular0D routine is called
    virtual void Evaluate
    (
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseMatrix> sysmat1, ///< Cardiovascular0D stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat2, ///< Cardiovascular0D offdiagonal matrix dV/dd
      Teuchos::RCP<LINALG::SparseOperator> sysmat3, ///< Cardiovascular0D offdiagonal matrix dfext/dp
      Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> sysvec2,
      Teuchos::RCP<Epetra_Vector> sysvec3,
      Teuchos::RCP<Epetra_Vector> sysvec4,
      Teuchos::RCP<Epetra_Vector> sysvec5
    );

  private:

    // don't want = operator, cctor and destructor

    Cardiovascular0D4ElementWindkessel operator = (const Cardiovascular0D4ElementWindkessel& old);
    Cardiovascular0D4ElementWindkessel(const Cardiovascular0D4ElementWindkessel& old);



  }; //class
} //namespace

#endif /*CARDIOVASCULAR0D4ELEMENTWINDKESSEL_H*/
