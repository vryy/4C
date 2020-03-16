/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class for MueLu
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_SELECTIVESAPFACTORY_DECL_HPP_
#define MUELU_SELECTIVESAPFACTORY_DECL_HPP_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

#include <string>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_ParameterListAcceptor.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SelectiveSaPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu
{
  /*!
    @class SelectiveSaPFactory class.
    @brief Factory for building Smoothed Aggregation prolongators.
           This is an extension to the SaPFactory class which can selectively smooth some
           prolongator/restrictor basis functions and keep unsmoothed basis functions for
           user-given aggregates.
    @ingroup MueLuTransferClasses
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class SelectiveSaPFactory : public PFactory
  {
#undef MUELU_SELECTIVESAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

   public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
    */
    SelectiveSaPFactory() {}

    //! Destructor.
    virtual ~SelectiveSaPFactory() {}

    Teuchos::RCP<const ParameterList> GetValidParameterList(
        const ParameterList &paramList = ParameterList()) const;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      //FIXME what does the return code mean (unclear in MueMat)?
      */
    void Build(Level &fineLevel, Level &coarseLevel) const;

    void BuildP(Level &fineLevel, Level &coarseLevel) const;  // Build()

    //@}

    //! @name Set methods.
    //@{

    //! Deprecated: Set prolongator smoother damping factor.
    void SetDampingFactor(Scalar dampingFactor);

    //! Deprecated: Change view of diagonal.
    void SetDiagonalView(std::string const &diagView);
    //@}

    //! @name Get methods.
    //@{

    //! Deprecated: Returns prolongator smoother damping factor.
    Scalar GetDampingFactor();

    //! Deprecated: Returns current view of diagonal.
    std::string GetDiagonalView();

    //@}

   private:
    Teuchos::RCP<Xpetra::Matrix<double, LocalOrdinal, GlobalOrdinal>> MyTranspose(
        Teuchos::RCP<Xpetra::Matrix<double, LocalOrdinal, GlobalOrdinal>> const &Op,
        bool const &optimizeTranspose) const;
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> FixAPproduct(
        Level &fineLevel, Level &coarseLevel, Teuchos::RCP<Matrix> &A,
        Teuchos::RCP<Matrix> &DinvAP) const;
  };  // class SelectiveSaPFactory

}  // namespace MueLu

#define MUELU_SELECTIVESAPFACTORY_SHORT

#endif  // HAVE_MueLuContact

#endif /* MUELU_SELECTIVESAPFACTORY_DECL_HPP_ */
