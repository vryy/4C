/*----------------------------------------------------------------------------*/
/**
\file xcontact_algorithm_partitioned.H

\brief partitioned algorithm for the inequality level-set approach for contact
       problems a.k.a. xcontact or extended contact

\maintainer Matthias Mayr

\date Jun 15, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_PARTITIONED_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_PARTITIONED_H_

#include "../drt_contact_xcontact/xcontact_algorithm_base.H"

namespace XCONTACT
{
  namespace ALGORITHM
  {
    class Partitioned : public Base
    {
     public:
      /// constructor
      Partitioned();

      /// destructor
      virtual ~Partitioned(){};

      /// Setup new class variables (derived)
      virtual void Setup();

     protected:
      //! Call the update routines in the different fields.
      virtual void OuterLoop();

      //! Prepares values and variables needed in the outer iteration
      void PrepareOuterIteration();

      //! Solve the level-set problem
      void DoScaTraField();

      //! Solve the structural contact problem
      void DoStructureField();

      //! convergence check of the partitioned outer loop
      bool ConvergenceCheck(const int& itnum);
    };
  }  // namespace ALGORITHM
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_PARTITIONED_H_ */
