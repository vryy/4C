/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu smoother interface
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_MYTRILINOSSMOOTHER_DECL_HPP_
#define MUELU_MYTRILINOSSMOOTHER_DECL_HPP_

#ifdef HAVE_MueLuContact

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_IfpackSmoother_fwd.hpp"
#include "MueLu_Ifpack2Smoother_fwd.hpp"

// Note: TrilinosSmoother is a SmootherPrototype that cannot be turned into a smoother using
// Setup().
//       When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu
{
  /*!
    @class TrilinosSmoother
    @brief Class that encapsulates external library smoothers. Autoselection of Ifpack or Ifpack2
    according to the underlying linear algebra library.
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MyTrilinosSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
#undef MUELU_MYTRILINOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

   public:
    //! @name Constructors / destructors
    //@{

    //! @brief Constructor
    MyTrilinosSmoother(std::string const &mapName, const Teuchos::RCP<const FactoryBase> &mapFact,
        std::string const &type = "",
        Teuchos::ParameterList const &paramList = Teuchos::ParameterList(), LO const &overlap = 0,
        Teuchos::RCP<const FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~MyTrilinosSmoother() {}

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! TrilinosSmoother cannot be turned into a smoother using Setup(). Setup() always returns a
    //! RuntimeError exception.
    void Setup(Level &currentLevel);

    //! TrilinosSmoother cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool InitialGuessIsZero = false) const;
    // void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
    Teuchos::RCP<SmootherPrototype> Copy() const;


    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    // using MueLu::Describable::describe; // overloading, not hiding
    // void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

   private:
    //! variable name for projected map
    std::string mapName_;

    //! Factory for projected variable
    Teuchos::RCP<const FactoryBase> mapFact_;

    //! map which has prescribed solution/rhs
    Teuchos::RCP<const Map> map_;

    //! ifpack1/2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Ifpack/Ifpack2 internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! A Factory
    Teuchos::RCP<const FactoryBase> AFact_;

    //
    // Underlying Smoother
    //

    //! Smoother
    Teuchos::RCP<SmootherPrototype> s_;  // TrilinosSmoother object

  };  // class MyTrilinosSmoother

}  // namespace MueLu

#endif  // HAVE_MueLuContact

#define MUELU_MYTRILINOSSMOOTHER_SHORT
#endif /* MUELU_MYTRILINOSSMOOTHER_DECL_HPP_ */
