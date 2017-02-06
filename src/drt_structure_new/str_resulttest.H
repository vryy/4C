/*-----------------------------------------------------------*/
/*!
\file str_resulttest.H

\maintainer Michael Hiermeier

\date Dec 21, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_RESULTTEST_H_
#define SRC_DRT_STRUCTURE_NEW_STR_RESULTTEST_H_

#include <Epetra_Vector.h>
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_dserror.H"


// forward declarations
namespace DRT {
  class DiscretizationInterface;
} // namespace DRT

namespace LINALG {
  class Solver;
} // namespace Solver

namespace IO {
  class DiscretizationWriter;
} // namespace IO

namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
} // namespace TIMINT

/*! \brief Structure specific result test class */
class ResultTest : public DRT::ResultTest
{
public:
  //! Constructor for time integrators of general kind
  //! \author bborn \date 06/08 (originally)
  ResultTest();

  //! destructor
  virtual ~ResultTest() {};

  //! initialization of class variables
  virtual void Init(const STR::TIMINT::BaseDataGlobalState& gstate);

  //! setup of class variables
  virtual void Setup();

  //! \brief structure version of nodal value tests
  //!
  //! Possible position flags are "dispx", "dispy", "dispz",
  //!                             "velx", "vely", "velz",
  //!                             "accx", "accy", "accz"
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

protected:
  /// get the indicator state
  inline const bool& IsInit() const
  { return isinit_; };

  /// get the indicator state
  inline const bool& IsSetup() const
  { return issetup_; };

  /// Check if Init() and Setup() have been called
  inline void CheckInitSetup() const
  {
    if (!IsInit() or !IsSetup())
      dserror("Call Init() and Setup() first!");
  }

  /// Check if Init() has been called
  inline void CheckInit() const
  {
    if (not IsInit())
      dserror("Call Init() first!");
  }

protected:
  //! flag which indicates if the Init() routine has already been called
  bool isinit_;

  //! flag which indicates if the Setup() routine has already been called
  bool issetup_;
private:

  //! our discretisation
  Teuchos::RCP<const DRT::DiscretizationInterface> strudisc_;
  // our solution
  //! global displacement DOFs
  Teuchos::RCP<const Epetra_Vector> disn_;
  //! global material displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dismatn_;
  //! global velocity DOFs
  Teuchos::RCP<const Epetra_Vector> veln_;
  //! global acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> accn_;
  /* NOTE: these have to be present explicitly
   * as they are not part of the problem instance like in fluid3
   */
}; // class ResultTest
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_RESULTTEST_H_ */
