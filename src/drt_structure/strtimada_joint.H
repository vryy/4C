/*======================================================================*/
/*!
\file strtimada_joint.H

\brief Indicator based on general STR::TimInt object

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef STRTIMADA_JOINT_H
#define STRTIMADA_JOINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
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


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class Solver;
}

namespace IO
{
  class DiscretizationWriter;
}

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
    enum INPAR::STR::TimAdaKind MapNameTimIntToTimAda
    (
      const enum INPAR::STR::DynamicType term  //!< input enum term
    ) const
    {
      switch (term)
      {
      case INPAR::STR::dyna_ab2 :
        return INPAR::STR::timada_kind_ab2;
        break;
      default :
        dserror("Cannot handle requested time integrator");
        return INPAR::STR::timada_kind_none;
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
      ada_(ada_vague),
      sta_(Teuchos::null)
    {
      // allocate auxiliary integrator
      sta_ = Teuchos::rcp(new T(ioparams, sdynparams, xparams,
                                sti->Discretization(),
                                sti->Solver(),
                                Teuchos::null, // no contact solver
                                sti->DiscWriter()));
      // check explicitness
      if (sta_->MethodImplicit())
      {
        dserror("Implicit might work, but please check carefully");
      }

      // check order
      if (sta_->MethodOrderOfAccuracyDis() > sti_->MethodOrderOfAccuracyDis())
      {
        ada_ = ada_upward;
      }
      else if (sta_->MethodOrderOfAccuracyDis() < sti_->MethodOrderOfAccuracyDis())
      {
        ada_ = ada_downward;
      }
      else if (sta_->MethodName() == sti_->MethodName())
      {
        ada_ = ada_ident;
      }
      else
      {
        ada_ = ada_orderequal;
      }

      // merge
      sta_->Merge(*sti);

      // resize multi-step quantities
      sta_->ResizeMStep();

      // I am lost
      return;
    }

    //! Destructor
    virtual ~TimAdaJoint() { ; }

    //@}

    //! @name Actions
    //@{

    //! Provide local discretisation error
    void IntegrateStepAuxiliar()
    {

      // get state vectors of marching integrator
      //const Teuchos::RCP<Epetra_Vector> dis = sti_->Dis();  // D_{n}^{A2}
      const Teuchos::RCP<Epetra_Vector> disn = sti_->DisNew();  // D_{n+1}^{A2}
      //const Teuchos::RCP<Epetra_Vector> vel = sti_->Vel();  // V_{n}^{A2}
      //const Teuchos::RCP<Epetra_Vector> acc = sti_->Acc();  // A_{n}^{A2}
      //const Teuchos::RCP<Epetra_Vector> accn = sti_->AccNew();  // A_{n+1}^{A2}

      // integrate
      sta_->IntegrateStep();
      // copy onto target
      locerrdisn_->Update(1.0, *(sta_->disn_), 0.0);
      // reset
      // remember: sta_ and sti_ are merged and work on the same vectors
      sta_->ResetStep();

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
    enum INPAR::STR::TimAdaKind MethodName() const
    {
      return MapNameTimIntToTimAda(sti_->MethodName());
    }

    //! Provide local order of accuracy of displacements
    int MethodOrderOfAccuracyDis() const
    {
      return sta_->MethodOrderOfAccuracyDis();
    }

    //! Provide local order of accuracy of velocities
    int MethodOrderOfAccuracyVel() const
    {
      return sta_->MethodOrderOfAccuracyVel();
    }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const
    {
      return sta_->MethodLinErrCoeffDis();
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const
    {
      return sta_->MethodLinErrCoeffVel();
    }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const
    {
      return ada_;
    }

    //@}

  protected:

    //! not wanted: = operator
    TimAdaJoint operator=(const TimAdaJoint& old);

    //! not wanted: copy constructor
    TimAdaJoint(const TimAdaJoint& old);

    //! type of adaptivity algorithm
    enum AdaEnum ada_;

    //! The auxiliary integrator
    Teuchos::RCP<T> sta_;

  };

}  // namespace STR



/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_JOINT_H
