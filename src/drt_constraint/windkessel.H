/*!----------------------------------------------------------------------
\file windkessel.H

\brief Monolithic coupling of 3D structure 0D Windkessel models

\level 2

<pre>
\maintainer Marc Hirschvogel
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
#include "../drt_fem_general/drt_utils_integration.H"



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
      wk_heartvalvearterial_proxdist,
      wk_heartvalvecardiovascular_full
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
      Teuchos::RCP<Epetra_Vector>    sysvec2,  ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec3  ///< distributed vector that may be filled by assembly of element contributions
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
      Teuchos::RCP<Epetra_Vector> sysvec9
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
    std::vector<DRT::Condition*> windkesselcond_; ///< windkessel conditions
    std::vector<DRT::Condition*> windkesselstructcoupcond_; ///< windkessel structure coupling conditions
    WindkesselType windkesseltype_; ///< Windkessel type
    std::map<int,double>  inittimes_; ///< map with times at which Windkessel is supposed to become active
    std::map<int,bool> activecons_; ///< map with indicator if Windkessels are active
    //! gaussian integration to be used
    DRT::UTILS::GaussRule2D gaussrule_;

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
      Teuchos::RCP<Epetra_Vector>    sysvec6,
      Teuchos::RCP<Epetra_Vector>    sysvec7
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateHeartValveArterialProxDistWindkessel
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
      Teuchos::RCP<Epetra_Vector>    sysvec8
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateDStructDp
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseOperator> sysmat ///< Windkessel offdiagonal matrix dfext/dp
    );

    //! Evaluate Windkessel conditions and assemble the results
    void EvaluateHeartValveCardiovascularFullWindkessel
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
      Teuchos::RCP<Epetra_Vector>    sysvec9
    );

    //! Compute and assemble initial Windkessel values
    void InitializeStdWindkessel
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

    //! Compute and assemble initial Windkessel values
    void InitializeHeartValveArterialProxDistWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2 ///< distributed vector that may be filled by assembly of element contributions
    );

    //! Compute and assemble initial Windkessel values
    void InitializeHeartValveCardiovascularFullWindkessel
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    sysvec1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec2, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    sysvec3 ///< distributed vector that may be filled by assembly of element contributions
    );

  }; //class
} //namespace

#endif /*WINDKESSEL_H*/
