/*======================================================================*/
/*!
\file strtimada_joint.H

\brief Indicator based on general STR::TimInt object

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET
#ifndef STRTIMADA_JOINT_H
#define STRTIMADA_JOINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_contact/drt_contact_manager.H"

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

#include "strtimada.H"
#include "strtimint.H"
#include "strtimint_ab2.H"

/*----------------------------------------------------------------------*/
/* belongs to structure namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Time step size adaptivity with general STR::TimInt scheme
   *
   * \author bborn
   * \date 10/07
   */
  template <typename T> class TimAdaJoint : public TimAda
  {

  public:

    //! Map STR::TimInt::NameEnum to STR::TimAda::NameEnum
    enum TimAda::NameEnum MapNameTimIntToTimAda
    (
      const enum TimInt::NameEnum term  //!< input enum term
    )
    {
      switch (term)
      {
      case TimInt::name_ab2 :
        return TimAda::name_ab2; 
        break; 
      default : 
        dserror("Cannot handle requested time integrator");
        return TimAda::name_vague; 
        break;
      }
    }

    //! @name Life
    //@{

    //! Full-fledged constructor
    TimAdaJoint
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
      Teuchos::RCP<TimInt>& sti  //!< marching time integrator
    )
    : TimAda
      (
        sdynparams,
        adaparams,
        sti
      ),
      sta_(Teuchos::null)
    {
      // allocate auxiliar integrator
      sta_ = Teuchos::rcp(new T(ioparams, sdynparams, xparams,
                                sti->Discretization(), 
                                sti->GetSolver(),
                                sti->GetDiscretizationWriter()));
      // check explicitness
      if (sta_->MethodImplicit())
        dserror("Implicit might work, but please check carefully");

      // merge 
      sta_->Merge(*sti);
      
      // resize multi-step quantities
      sta_->ResizeMStep();
      
      // I am lost
      return;
    }

    //! Destructor
    ~TimAdaJoint() { ; }

    //@}

    //! @name Actions
    //@{

    //! Provide local discretisation error
    void IntegrateStepAuxiliar()
    {

      // get state vectors of marching integrator
      //const Teuchos::RCP<Epetra_Vector> dis = sti_->Disp();  // D_{n}^{A2}
      const Teuchos::RCP<Epetra_Vector> disn = sti_->Dispn();  // D_{n+1}^{A2}
      //const Teuchos::RCP<Epetra_Vector> vel = sti_->Vel();  // V_{n}^{A2}
      //const Teuchos::RCP<Epetra_Vector> acc = sti_->Acc();  // A_{n}^{A2}
      //const Teuchos::RCP<Epetra_Vector> accn = sti_->Accn();  // A_{n+1}^{A2}

      // integrate
      sta_->IntegrateStep();
      sta_->ResetStep();
      // copy onto target
      locerrdisn_->Update(1.0, *(sta_->disn_), 0.0);
      
      // provide local discretisation error vector
      // l_{n+1}^{A2} = D_{n+1}^{Joint} - D_{n+1}^{A1}
      //locerrdisn_->Update(-1.0, *disn, 1.0);

      // see you
      return;
    }

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    enum TimAda::NameEnum MethodName()
    {
      return MapNameTimIntToTimAda(sti_->MethodName());
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracy()
    {
      return 2;
    }

    //@}

  protected:

    //! not wanted: = operator
    TimAdaJoint operator = (const TimAdaJoint& old);
  
    //! not wanted: copy constructor
    TimAdaJoint(const TimAdaJoint& old);

    //! The auxiliar integrator 
    Teuchos::RCP<T> sta_;

  };

  /*====================================================================*/
  /*! \brief Explicitly declare #TimIntAB2 instantiated class 
   *         for joint auxiliar schemes
   */
  template class TimAdaJoint<TimIntAB2>;

}  // namespace STR



/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_JOINT_H
#endif  // #ifdef CCADISCRET
