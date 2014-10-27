/*----------------------------------------------------------------------------*/
/*!
\file fas_hierarchy.H

\brief Hierarchy of multigrid levels based on MueLu

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef FAS_HIERARCHY_H
#define FAS_HIERARCHY_H

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <vector>

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace NLNSOL
{
  class NlnProblem;

  namespace FAS
  {
    class NlnLevel;
  } // namespace FAS
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  namespace FAS
  {
    /*! \brief Hierarchy of multigrid levels based on MueLu
     *
     *  Here, we create and store the hierarchy of multigrid levels. For
     *  creation of the multigrid hierarchy, Setup() utilizes the MueLu package
     *  from Trilinos. We then extract level transfer operators and coarse level
     *  matrix representations and store them in our own NLNSOL::FAS::NlnLevel
     *  containers.
     *
     *  Level transfer (ProlongateToFineLevel(), RestrictToCurrentLevel()) is
     *  managed by this class.
     *
     *  \sa NLNSOL::NlnOperatorFas
     *  \sa NLNSOL::FAS::NlnLevel
     *
     *  \author mayr.mt \date 06/2014
     */
    class AMGHierarchy
    {

    public:

      //! @name Construction / Destruction
      //@{

      //! Constructor (empty)
      AMGHierarchy();

      //! Destructor
      virtual ~AMGHierarchy(){}

      //@}

      //! @name Setup
      //@{

      //! Initialization of member variables
      virtual void Init(const Epetra_Comm& comm, ///< communicator
          const Teuchos::ParameterList& params, ///< parameter list
          Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem ///< nonlinear problem to solve
          );

      //! Setup of multigrid hierarchy
      virtual void Setup();

      //@}

      //! @name Multigrid helpers
      //@{

      //! Restrict fine level vector to given coarse level
      Teuchos::RCP<Epetra_MultiVector> RestrictToCoarseLevel(
          const Epetra_MultiVector& vec, ///< vector
          const int targetlevel ///< target level (some coarse level)
      ) const;

      //! Prolongate vector from a given coarse level to fine level
      Teuchos::RCP<Epetra_MultiVector> ProlongateToFineLevel(
          const Epetra_MultiVector& veccoarse, ///< vector
          const int sourcelevel ///< source level (some coarse level)
      ) const;

      //@}

      //! @name Attribute access functions
      //@{

      //! Returns a character string describing the operator
      virtual const char* Label() const { return "NLNSOL::FAS::AMGHierarchy"; }

      //! Access to a certain nonlinear level object
      virtual Teuchos::RCP<NLNSOL::FAS::NlnLevel> NlnLevel(const int i) const;

      /*! \brief Get total number of levels
       *
       *  \return Total number of levels
       */
      virtual const int NumLevels() const { return nlnlevels_.size(); }

      //@}

      //! @name Utilities
      //@{

      //! Print all levels
      virtual void PrintNlnLevels(std::ostream& os) const;

      /*! \brief Check validity of mutligrid hierarchy
       *
       *  Check important stuff like if maps of matrices and transfer operator
       *  match.
       *
       *  \return Boolean flag indicating whether hierarchy has been set up
       *  correctly
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual const bool CheckValidity() const;

      /*! \brief Check if level ID is admissible
       *
       *  Level ID has to be in admissible range, i.e. between zero and the
       *  number of levels. Throw dserror() if level ID is not admissible.
       *
       *  \return Boolean flag indicating whether level ID is admissible
       *
       *  \sa NumLevels()
       *
       *  \author mayr.mt \date 11/2013
       */
      virtual const bool CheckLevelID(const int level) const;

      //@}

    protected:

    private:

      //! @name Setup
      //@{

      /*! \brief Add new level at the end of level container
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual void AddNlnLevel(Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel)
      {
        nlnlevels_.push_back(newlevel);
        return;
      }

      //! Init() has been called
      virtual void SetIsInit() { isinit_ = true; return; }

      //! Setup() has been called
      virtual void SetIsSetup() { issetup_ = true; return; }

      //@}

      //! Access to communicator
      virtual const Epetra_Comm& Comm() const;

      //! Access to parameter list
      virtual const Teuchos::ParameterList& Params() const;

      //! Access to the nonlinear problem
      virtual Teuchos::RCP<NLNSOL::NlnProblem> NlnProblem() const;

      //! @name Sanity and validity checks
      //@{

      //! Has Init() already been called?
      virtual const bool IsInit() const { return isinit_; }

      //! Has Setup() already been called?
      virtual const bool IsSetup() const { return issetup_; }

      //@}

      //! Init() has been called?
      bool isinit_;

      //! Setup() has been called?
      bool issetup_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! parameter list
      Teuchos::RCP<const Teuchos::ParameterList> params_;

      //! Nonlinear problem to be evaluated and solved
      Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem_;

      //! vector containing all nonlinear levels
      Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel> > nlnlevels_;

    }; // class AMGHierarchy
  } // namespace FAS
} // namespace NLNSOL


#endif /* FAS_HIERARCHY_H */
