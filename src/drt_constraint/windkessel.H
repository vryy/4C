/*!----------------------------------------------------------------------
\file windkessel.cpp

\brief Basic Windkessel class, dealing with Windkessel functions living on boundaries

**************************************************************************************************************************
Monolithic coupling of structure and a three-element Windkessel governed by
c dp/dt - c r2 dQ/dt + (p-p_0)/r1 - (1 + r2/r1) Q(d) = 0
[c: compliance, r1: first resistance, r2: second resistance, Q = -dV/dt: flux, p: pressure variable]
						                ____
											 ____|R_1 |___
   			    _____	   	|		 |____|		|
----Q >----|_R_2_|----|							|--->----| p_0
    									|______|C|____|
												     | |

There are three different versions:
a) the standard model,

b) a special trimodular one where every parameter is pressure-dependent in terms of
PARAM(p) = (PARAM_b - PARAM_c)*0.5*(1.0 - tanh[(p-p_bc)/k_p] ) + PARAM_c + (PARAM_a - PARAM_b)*(1.0 - tanh[(p-p_ab)/k_p] )

c) a heart-specific model governing the arterial pressure with a three-element Windkessel with an additional valve law
(resistive Windkessel) infront of it, opening when arterial pressure is reached in the ventricels and closing when flux is
beginning to be reversed (cf. Sainte-Marie et. al. "Modeling and estimation of the cardiac electromechanical activity",
Comp. & Struct. 84 (2006) 1743-1759)

The Windkessel is momolithically coupled with the standard structural dynamics governing equation

M a + C v + f_int(d) - f_ext(d,p) = 0,

with q being a function of the displacement vector d and f_ext additionally being a function of the Windkessel pressure p.
**************************************************************************************************************************

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#ifndef WINDKESSEL_H
#define WINDKESSEL_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>
#include "Epetra_FECrsMatrix.h"



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

  class Windkessel

  {
  public:

    //! Windkessel types
    enum WindkesselType
    {
      none,
      wk_std,
      wk_trimodular,
      wk_heartvalvearterial,
    };

    /*!
    \brief Constructor of a Windkessel based on conditions with a given name. It also
    takes care of the Windkessel IDs.
    */

    Windkessel
    (
      Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization where Windkessel lives on
      const std::string& CondName, ///< Name of condition to creat Windkessel from
      int& minID,  ///< minimum Windkessel ID so far
      int& maxID, ///< maximum Windkessel ID so far
      std::vector<int>& curID ///< current ID
    );

    /*!
    \brief Constructor of a Windkessel based on a conditions with a given name.
    */

    Windkessel
    (
      Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization where Windkessel funtion lives on
      const std::string& CondName ///< Name of condition to create Windkessel functions from
    );

    /*!
        \brief Destructor

     */
    ~Windkessel()
    {
        return;
    };

    /*!
     \brief Return if there are Windkessel functions
    */
    bool HaveWindkessel()
    {
        return windkesseltype_!=none;
    };

    /// Set state of the underlying discretization
    void SetState
    (
      const std::string& state,  ///< name of state to set
      Teuchos::RCP<Epetra_Vector> V  ///< values to set
    );


    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    void Initialize
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1,  ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2  ///< distributed vector that may be filled by assembly of element contributions
    );

    /// initialization routine called at restart to activate the right conditions
    void Initialize
    (
      const double& time  ///< current time
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateWindkessel routine is called
    void Evaluate
    (
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseMatrix> sysmat1, ///< Windkessel stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat2, ///< Windkessel offdiagonal matrix dV/dd
      Teuchos::RCP<LINALG::SparseOperator> sysmat3, ///< Windkessel offdiagonal matrix dfext/dp
      Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> sysvec2,
      Teuchos::RCP<Epetra_Vector> sysvec3,
      Teuchos::RCP<Epetra_Vector> sysvec4,
      Teuchos::RCP<Epetra_Vector> sysvec5,
      Teuchos::RCP<Epetra_Vector> sysvec6,
      Teuchos::RCP<Epetra_Vector> sysvec7,
      Teuchos::RCP<Epetra_Vector> sysvec8,
      Teuchos::RCP<Epetra_Vector> sysvec9,
      Teuchos::RCP<Epetra_Vector> sysvec10
    );

    ///Return type of Windkessel function
    WindkesselType Type()
    {
      return windkesseltype_;
    }

    ///Return vector with IDs of active conditions
    std::vector<int> GetActiveCondID();

  protected:

    Teuchos::RCP<DRT::Discretization> actdisc_; ///< standard discretization
    std::vector<DRT::Condition*> windkesselcond_; ///< conditions, that define the Windkessel (all of the same kind)
    WindkesselType windkesseltype_; ///< Windkessel type
    std::map<int,double>  inittimes_; ///< map with times at which Windkessel is supposed to become active
    std::map<int,bool> activecons_; ///< map with indicator if Windkessels are active

  private:

    // don't want = operator, cctor and destructor

    Windkessel operator = (const Windkessel& old);
    Windkessel(const Windkessel& old);

    //! Return the WindkesselType based on the condition name
    WindkesselType GetWindkesselType
    (
      const std::string& Name///<condition name
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateStdWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseMatrix> sysmat1, ///< Windkessel stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat2, ///< Windkessel offdiagonal matrix dV/dd
      Teuchos::RCP<LINALG::SparseOperator> sysmat3, ///< Windkessel offdiagonal matrix dfext/dp
      Teuchos::RCP<Epetra_Vector>    sysvec1,  ///< distributed vectors that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2,
      Teuchos::RCP<Epetra_Vector>    sysvec3,
      Teuchos::RCP<Epetra_Vector>    sysvec4,
      Teuchos::RCP<Epetra_Vector>    sysvec5,
      Teuchos::RCP<Epetra_Vector>    sysvec6
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateTrimodularWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseMatrix> sysmat1, ///< Windkessel stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat2, ///< Windkessel offdiagonal matrix dV/dd
      Teuchos::RCP<LINALG::SparseOperator> sysmat3, ///< Windkessel offdiagonal matrix dfext/dp
      Teuchos::RCP<Epetra_Vector>    sysvec1,  ///< distributed vectors that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2,
      Teuchos::RCP<Epetra_Vector>    sysvec3,
      Teuchos::RCP<Epetra_Vector>    sysvec4,
      Teuchos::RCP<Epetra_Vector>    sysvec5,
      Teuchos::RCP<Epetra_Vector>    sysvec6,
      Teuchos::RCP<Epetra_Vector>    sysvec7,
      Teuchos::RCP<Epetra_Vector>    sysvec8,
      Teuchos::RCP<Epetra_Vector>    sysvec9,
      Teuchos::RCP<Epetra_Vector>    sysvec10
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateHeartValveArterialWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseMatrix> sysmat1, ///< Windkessel stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat2, ///< Windkessel offdiagonal matrix dV/dd
      Teuchos::RCP<LINALG::SparseOperator> sysmat3, ///< Windkessel offdiagonal matrix dfext/dp
      Teuchos::RCP<Epetra_Vector>    sysvec1,  ///< distributed vectors that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2,
      Teuchos::RCP<Epetra_Vector>    sysvec3,
      Teuchos::RCP<Epetra_Vector>    sysvec4,
      Teuchos::RCP<Epetra_Vector>    sysvec5,
      Teuchos::RCP<Epetra_Vector>    sysvec6,
      Teuchos::RCP<Epetra_Vector>    sysvec7
    );

    //! Compute and assemble initial Windkessel values
    void InitializeStdAndTrimodularWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2 ///< distributed vector that may be filled by assembly of element contributions
    );

    //! Compute and assemble initial Windkessel values
    void InitializeHeartValveArterialWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2 ///< distributed vector that may be filled by assembly of element contributions
    );
  }; //class
} //namespace

#endif /*WINDKESSEL_H*/
