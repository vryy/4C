/*----------------------------------------------------------------------------*/
/*!
\file fas_hierarchy.H

\level 2

\brief Declaration of hierarchy of multigrid levels based on MueLu

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef FAS_HIERARCHY_H
#define FAS_HIERARCHY_H

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <vector>

// MueLu
#ifdef HAVE_MueLu
//Xpetra
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>

// MueLu
//#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

/* MueLu typedefs: header files for default types, must be included after all
 * other MueLu/Xpetra headers */
#include <MueLu_UseDefaultTypes.hpp> // => Scalar = double, LocalOrdinal = GlobalOrdinal = int

#endif

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace NLNSOL
{
  class NlnProblemBase;

  namespace FAS
  {
    class NlnLevel;
  } // namespace FAS

  namespace UTILS
  {
    class NlnConfig;
  } // namespace UTILS
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  //! NLNSOL::FAS: A Full Approximation Scheme multigrid within the NLNSOL framework
  namespace FAS
  {
    /*! \class AMGHierarchy
     *  \brief Hierarchy of multigrid levels based on MueLu
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
    class AMGHierarchy : public Teuchos::VerboseObject<AMGHierarchy>
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
          Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config, ///< configuration
          const std::string listname, ///< name of sublist with configuration of this operator
          Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem, ///< nonlinear problem to solve
          const Teuchos::ParameterList& mlparams ///< ML parameter list from baci linear solver (just for null space)
          );

      //! Setup of multigrid hierarchy
      virtual void Setup();

      //! Refresh RAPs
      virtual void RefreshRAPs();

      //@}

      //! @name Multigrid helpers
      //@{

      /*! \brief Restrict fine level vector to given coarse level \c targetlevel
       *
       *  \return Coarse level vector on level \c targetlevel
       */
      Teuchos::RCP<Epetra_MultiVector> RestrictToCoarseLevel(
          const Epetra_MultiVector& vec, ///< vector
          const int targetlevel ///< target level (some coarse level)
      ) const;

      /*! \brief Prolongate vector from a given coarse level \c sourcelevel to fine level
       *
       *  \return Fine level vector
       */
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
      virtual int NumLevels() const { return nlnlevels_.size(); }

      //@}

      //! @name Utilities
      //@{

      //! Print all levels
      virtual void PrintNlnLevels(std::ostream& os) const;

      /*! \brief Check validity of multigrid hierarchy
       *
       *  Check important stuff like if maps of matrices and transfer operator
       *  match.
       *
       *  \return bool: true if setup was correct, false otherwise
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual bool CheckValidity() const;

      /*! \brief Check if level ID is admissible
       *
       *  Level ID has to be in admissible range, i.e. between zero and the
       *  number of levels. Throw dserror() if level ID is not admissible.
       *
       *  \return bool: true if level ID is admissible, false otherwise
       *
       *  \sa NumLevels()
       *
       *  \author mayr.mt \date 11/2013
       */
      virtual bool CheckLevelID(const int level) const;

      /*! \brief Check for stagnation of level smoothers
       *
       *  Loop over all levels and check, whether at least one level smoother
       *  got stuck in stagnation.
       *
       *  \return bool: true if stagnation has been detected, false otherwise
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual bool CheckAllLevelStagnation() const;

      //@}

    protected:

    private:

      //! @name Setup
      //@{

      /*! \brief Setup MueLu hierarchy
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual void SetupMueLuHierarchy();

      /*! \brief Setup NLNSOL hierarchy and fill it with MueLu content
       *
       *  \return bool: true if setup was successful, false otherwise
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual bool SetupNlnSolHierarchy();

      /*! \brief Add new level at the end of level container
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual void AddNlnLevel(
          Teuchos::RCP<NLNSOL::FAS::NlnLevel> newlevel ///< level to be added to level container #nlnlevels_
          );

      //! Init() has been called
      virtual void SetIsInit() { isinit_ = true; return; }

      //! Setup() has been called
      virtual void SetIsSetup() { issetup_ = true; return; }

#ifdef HAVE_MueLu
      /*! \brief Access to fine level matrix as Xpetra::Matrix
       *
       *  Ask the #nlnproblem_ for the fine level Epetra_Operator and wrap it
       *  into Xpetra::Matrix.
       *
       *  \return Fine level matrix as Xpetra::Matrix
       *
       *  \author mayr.mt \date 03/2015
       */
      Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> >
      GetXpetraFineLevelMatrix() const;

      /*! \brief Access to pre-computed null space as Xpetra::MultiVector
       *
       *  Baci pre-computes the null space and stores it in a
       *  Teuchos::ParameterList. Take this null space and transform it into a
       *  Xpetra::MultiVector.
       *
       *  \note Use NLNSOL::UTILS::GetXpetraNullSpaceFromBaci() for actual
       *  computation
       *
       *  \return Teuchos::RCP<Xpetra::MultiVector> with \f$N\f$ vectors, that
       *  contains the \f$N\f$ null space vectors
       *
       *  \sa NLNSOL::UTILS::GetXpetraNullSpaceFromBaci()
       *
       *  \author mayr.mt \date 06/2015
       */
      Teuchos::RCP<Xpetra::MultiVector<double,int,int,Node> >
      GetXpetraNullSpaceFromBaci() const;

      /*! \brief Debug-Tool: Check null space conservation properties
       *
       *  The null space is supposed to satisfy the following properties:
       *  \f[
       *  A_iVi=0
       *  A_0\left(P_{i0}V_i\right)=0
       *  \f]
       *  with the level index \f$i\f$, the matrix \f$A_i\f$ on level \f$i\f$,
       *  the prolongator \f$P_{i0}\f$ from level \f$i\f$ to level \f$0\f$ and
       *  the null space \f$V_i\f$ on level \f$i\f$.
       *
       *  Here, this can be computed for levels 0, 1, and 2.
       *
       *  \return bool: true if check was successful, false otherwise
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual bool CheckNullSpaceProperties() const;
#endif
      //@}

      //! Access to communicator
      virtual const Epetra_Comm& Comm() const;

      //! Access to global configuration
      virtual Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> Configuration() const;

      //! Get name of this operator's parameter list in the global configuration
      virtual const std::string& MyListName() const { return listname_; }

      //! Access to the nonlinear problem
      virtual Teuchos::RCP<NLNSOL::NlnProblemBase> NlnProblem() const;

      //! @name Sanity and validity checks
      //@{

      //! Has Init() already been called?
      virtual bool IsInit() const { return isinit_; }

      //! Has Setup() already been called?
      virtual bool IsSetup() const { return issetup_; }

      //@}

      //! Init() has been called?
      bool isinit_;

      //! Setup() has been called?
      bool issetup_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! @name Configuration
      //@{

      //! Configuration of this line search operator
      Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

      //! Name of sublist in #config_ for this operator
      std::string listname_;

      //@}

      //! Nonlinear problem to be evaluated and solved
      Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem_;

      //! ML parameter list of baci linear solver (just for nullspace)
      Teuchos::RCP<const Teuchos::ParameterList> mlparams_;

      //! vector containing all nonlinear levels
      Teuchos::Array<Teuchos::RCP<NLNSOL::FAS::NlnLevel> > nlnlevels_;

#ifdef HAVE_MueLu
      //! MueLu Hierarchy Factory to create #mueLuHierarchy_
      Teuchos::RCP<MueLu::HierarchyManager<double,int,int,Node> > mueLuFactory_;

      //! MueLu Hierarchy
      Teuchos::RCP<MueLu::Hierarchy<double,int,int,Node> > mueLuHierarchy_;
#endif

    }; // class AMGHierarchy
  } // namespace FAS
} // namespace NLNSOL


#endif /* FAS_HIERARCHY_H */
