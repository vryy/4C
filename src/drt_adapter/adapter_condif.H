/*----------------------------------------------------------------------*/
/*!
\file adapter_condif.H

\brief convection-diffusion field adapter

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_CONDIF_H
#define ADAPTER_CONDIF_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/condifimplicitintegration.H"
#include "../drt_fluid/condifresulttest.H"


/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general fluid field interface for multiphysics (FSI, ELCH, ...)
  /*!

  \author gjb
  \date 11/07
  */
  class ConDifTimeIntegration
  {
  public:

    /// virtual destructor to get polymorph destruction
    virtual ~ConDifTimeIntegration();
    
    //! @name Time step helpers

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;
    
    virtual void Integrate() = 0;
    virtual void Solve() = 0;
    virtual void SetVelocityField(int veltype, int velfuncno) = 0;
    virtual void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel) = 0;

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };
  
#if 1
  class ConDifImplicit : public ConDifTimeIntegration
    {
    public:
      ConDifImplicit(Teuchos::RCP<DRT::Discretization> dis,
                   Teuchos::RCP<LINALG::Solver> solver,
                   Teuchos::RCP<ParameterList> params,
                   Teuchos::RCP<IO::DiscretizationWriter> output
                   );

      virtual void PrepareTimeStep();
      virtual void Update();
      virtual void Output();
      virtual void ReadRestart(int step);
      virtual void Integrate();
      virtual void Solve();
      virtual void SetVelocityField(int veltype, int velfuncno);
      virtual void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel);

      virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();
    
    private:

      /// the actual convection diffusion algorithm
      CondifImplicitTimeInt condif_;

      //! @name local copies of input parameters
      Teuchos::RCP<DRT::Discretization>      dis_;
      Teuchos::RCP<LINALG::Solver>           solver_;
      Teuchos::RCP<Teuchos::ParameterList>   params_;
      Teuchos::RCP<IO::DiscretizationWriter> output_;
      //@}


 }; // class ConDifImplicit
 
#endif
  
} // namespace ADAPTER

#endif
#endif
