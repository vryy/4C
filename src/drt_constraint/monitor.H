/*!----------------------------------------------------------------------
\file monitor.H

\brief Basic monitor class, dealing with monitor boundary conditions

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef MONITOR_H
#define MONITOR_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_sparsematrix.H"


using namespace std;
using namespace Teuchos;

/*!
\brief UTILS: namespace constaining general utils for BACI, as constraints, monitors, ...

*/
namespace UTILS
{
  /*!
  \brief Basic monitor class, dealing with monitor boundary conditions. 
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class Monitor
  
  {
  public:
  
    //! Monitor types
    enum MoniType
    {
      none,
      volmonitor3d,
      areamonitor3d,
      areamonitor2d
    };
    
    /*!
    \brief Constructor of a monitor based on a conditions with a given name. It also
    takes care of the monitor IDs.
    */
  
    Monitor
    (
      RCP<DRT::Discretization> discr,  ///< Discretization monitor lives on
      const string& CondName, ///< Name of condition to creat monitor from
      int& minID,  ///< minimum monitor ID so far
      int& maxID  ///< maximum monitor ID so far   
    );  
  
    /*!
        \brief Destructor
  
     */
    ~Monitor()
    {
        return;
    };
  
    /*!
     \brief Return if there are monitors
    */
    bool HaveMonitor()
    {
        return montype_!=none;
    };
  
    /// Set state of the underlying discretization
    void SetState
    (
      const string& state,  ///< name of state to set
      RCP<Epetra_Vector> V  ///< values to set
    )
    {
      actdisc_->SetState(state,V);
    };
    
    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateConstraint routine is called 
    void Evaluate
    ( 
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector> systemvector1 ///< distributed vector that may be filled by assembly of element contributions
    ); 
    
    
    ///Return type of monitor
    MoniType Type()
    {
      return montype_;
    }
    

  protected:
    
    RCP<DRT::Discretization> actdisc_; /// standard discretization
    vector<DRT::Condition*> moncond_; /// conditions, that define the monitor (all of the same kind)
    MoniType montype_; /// monitor type
    map<int,double>  inittimes_; /// map with times at which monitor is supposed to become active
    map<int,bool> activemons_; /// map with indicator if monitors are active
      
  private:
    
    // don't want = operator, cctor and destructor
  
    Monitor operator = (const Monitor& old);
    Monitor(const Monitor& old);
  
    //! Return the MoniType based on the condition name
    MoniType GetMoniType
    (
      const string& Name///<condition name
    );
    
    
    //! Evaluate monitor values and assemble the results
    void EvaluateMonitor
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    systemvector ///< distributed vector that may be filled by aasembly of element contributions
    ); 
  }; //class
} //namespace

#endif /*CONSTRAINT_H*/
