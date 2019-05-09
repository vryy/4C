/*======================================================================*/
/*!
\brief Indicator based on general STR::TimInt object

\level 1

\maintainer Alexander Popp
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
  template <typename T>
  class TimAdaJoint : public TimAda
  {
   public:
    //! Map STR::TimInt::NameEnum to STR::TimAda::NameEnum
    enum INPAR::STR::TimAdaKind MapNameTimIntToTimAda(
        const enum INPAR::STR::DynamicType term  //!< input enum term
        ) const
    {
      switch (term)
      {
        case INPAR::STR::dyna_ab2:
          return INPAR::STR::timada_kind_ab2;
          break;
        case INPAR::STR::dyna_expleuler:
          return INPAR::STR::timada_kind_expleuler;
          break;
        case INPAR::STR::dyna_centrdiff:
          return INPAR::STR::timada_kind_centraldiff;
          break;
        default:
          dserror("Cannot handle requested time integrator");
          return INPAR::STR::timada_kind_none;
          break;
      }
    }

    //! @name Life
    //@{

    //! Full-fledged constructor
    TimAdaJoint(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& timeparams,        //!< TIS input parameters
        const Teuchos::ParameterList& sdyn,              //!< structural dynamic
        const Teuchos::ParameterList& xparams,           //!< extra flags
        const Teuchos::ParameterList& adaparams,         //!< adaptive input flags
        Teuchos::RCP<TimInt>& sti                        //!< marching time integrator
        )
        : TimAda(timeparams, adaparams, sti), ada_(ada_vague), sta_(Teuchos::null)
    {
      // allocate auxiliary integrator
      sta_ = Teuchos::rcp(
          new T(timeparams, ioparams, sdyn, xparams, sti->Discretization(), sti->Solver(),
              Teuchos::null,  // no contact solver
              sti->DiscWriter()));
      sta_->Init(timeparams, sdyn, xparams, sti->Discretization(), sti->Solver());

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

      // setup
      sta_->Setup();

      // Actually, we would like to call Merge() and ResizeMStep() now and in the
      // past we have done this right here. However, this requires that both Init()
      // and Setup() have already been called on both(!) the marching time integrator
      // and the auxiliary time integrator. Since Setup() has not yet been called
      // for the marching time integrator, we must postpone Merge() and ResizeMStep()
      // to a later point. They can now be found in the new Init() function below,
      // which is called at the beginning of Integrate(). [popp 01/2017]
      // sta_->Merge(*sti);
      // sta_->ResizeMStep();

      return;
    }

    //! Destructor
    virtual ~TimAdaJoint() { ; }

    //@}

    //! @name Actions
    //@{

    /*! Finalize the class initialization
     * Merge() and ResizeMStep() need to be called after(!) both Init()
     * and Setup() have been called on both the marching time integrator
     * and the auxiliary time integrator (popp 01/2017).
     */
    void Init(Teuchos::RCP<TimInt>& sti)
    {
      // merge
      sta_->Merge(*sti);

      // resize multi-step quantities
      sta_->ResizeMStep();

      return;
    }

    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    void IntegrateStepAuxiliar()
    {
      // integrate the auxiliary time integrator one step in time
      sta_->IntegrateStep();

      // copy onto target
      locerrdisn_->Update(1.0, *(sta_->disn_), 0.0);

      // reset
      // remember: sta_ and sti_ are merged and work on the same vectors
      sta_->ResetStep();

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
    int MethodOrderOfAccuracyDis() const { return sta_->MethodOrderOfAccuracyDis(); }

    //! Provide local order of accuracy of velocities
    int MethodOrderOfAccuracyVel() const { return sta_->MethodOrderOfAccuracyVel(); }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const { return sta_->MethodLinErrCoeffDis(); }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const { return sta_->MethodLinErrCoeffVel(); }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const { return ada_; }

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
