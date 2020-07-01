/*----------------------------------------------------------------------------*/
/** \file
\brief xcontact level-set generic reinitialization algorithm


\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_GENERIC_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_GENERIC_H_

#include "../drt_lib/drt_dserror.H"

class Epetra_Vector;

namespace XCONTACT
{
  namespace LEVELSET
  {
    class Algorithm;
    namespace REINIT
    {
      class Generic
      {
       public:
        /// constructor
        Generic();

        /// destructor
        virtual ~Generic(){};

        void Init(XCONTACT::LEVELSET::Algorithm* const algorithm);

        virtual void Setup() = 0;

        void Compute(const Epetra_Vector& phinp);

       protected:
        /// access the calling xcontact levelset algorithm ( call-back function )
        XCONTACT::LEVELSET::Algorithm& Algorithm();
        const XCONTACT::LEVELSET::Algorithm& Algorithm() const;

        /// run at the beginning of the nonlinear solve call
        virtual void PreSolve(){/* do nothing in the default case */};

        virtual void Solve(const Epetra_Vector& phinp) = 0;

        /// run after the nonlinear solver stopped
        virtual void PostSolve(){/* do nothing in the default case */};

        /// check if the Init() routine has been called
        inline void CheckInit() const
        {
          if (not isinit_) dserror("Call Init() first!");
        }

        /// check if the Init() and Setup() routines have been called
        inline void CheckInitSetup() const
        {
          if (not isinit_ or not issetup_) dserror("Call Init() first!");
        }

       protected:
        bool isinit_;

        bool issetup_;

       private:
        XCONTACT::LEVELSET::Algorithm* algorithm_;
      };  // class Generic
    }     // namespace REINIT
  }       // namespace LEVELSET
}  // namespace XCONTACT



#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_GENERIC_H_ */
