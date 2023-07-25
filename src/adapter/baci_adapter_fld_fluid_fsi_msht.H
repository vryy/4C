/*--------------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi with internal mesh tying or mesh sliding


\level 3
*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_ADAPTER_FLD_FLUID_FSI_MSHT_H
#define BACI_ADAPTER_FLD_FLUID_FSI_MSHT_H


#include "baci_adapter_fld_fluid_fsi.H"

namespace FLD
{
  namespace UTILS
  {
    class FsiMapExtractor;
  }
}  // namespace FLD

namespace ADAPTER
{
  /*! \brief Fluid field adapter for fsi with internal mesh tying or mesh sliding
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
  class FluidFSIMsht : public FluidFSI
  {
   public:
    /// Constructor
    FluidFSIMsht(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    virtual void Init();

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::FsiMapExtractor> const& FsiInterface() const
    {
      return fsiinterface_;
    }


   protected:
    /// create conditioned dof-map extractor for the fluid
    virtual void SetupFsiInterface();

    //! \brief interface map setup for fsi interface and other
    //!
    //! Note: full map contains velocity AND pressure DOFs
    Teuchos::RCP<FLD::UTILS::FsiMapExtractor> fsiinterface_;
  };
}  // namespace ADAPTER


#endif  // ADAPTER_FLD_FLUID_FSI_MSHT_H
