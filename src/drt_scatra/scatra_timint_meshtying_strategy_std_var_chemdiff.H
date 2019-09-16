/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for variational chemical diffusion
     problems (without meshtying)

\level 3

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_VAR_CHEMDIFF_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_STD_VAR_CHEMDIFF_H

#include "scatra_timint_variational.H"

#include "scatra_timint_meshtying_strategy_std.H"

namespace SCATRA
{
  class MeshtyingStrategyStdVar : public MeshtyingStrategyStd
  {
   public:
    //! constructor
    explicit MeshtyingStrategyStdVar(SCATRA::TimIntVariational* vartimint);

    // destructor
    virtual ~MeshtyingStrategyStdVar() { return; };

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

    //! setup meshtying objects
    void SetupMeshtying();

    //! return map extractor associated with blocks of global system matrix
    const LINALG::MultiMapExtractor& BlockMaps() const { return *blockmaps_; };

   protected:
    //! build null spaces associated with blocks of global system matrix
    virtual void BuildBlockNullSpaces() const;

    //! map extractors associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_;

   private:
    //! copy constructor
    MeshtyingStrategyStdVar(const MeshtyingStrategyStdVar& old);

    //! return pointer to variational time integrator after cast
    SCATRA::TimIntVariational* VarTimInt() const
    {
      return dynamic_cast<SCATRA::TimIntVariational*>(scatratimint_);
    };

  };  // class MeshtyingStrategyStdVar
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_VAR_CHEMDIFF_H
