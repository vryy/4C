/*--------------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi with internal mesh tying or mesh sliding


\level 3
*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_FLD_FLUID_FSI_MSHT_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_FSI_MSHT_HPP


#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  namespace UTILS
  {
    class FsiMapExtractor;
  }
}  // namespace FLD

namespace Adapter
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
    FluidFSIMsht(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    void init() override;

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::FsiMapExtractor> const& fsi_interface() const
    {
      return fsiinterface_;
    }


   protected:
    /// create conditioned dof-map extractor for the fluid
    virtual void setup_fsi_interface();

    //! \brief interface map setup for fsi interface and other
    //!
    //! Note: full map contains velocity AND pressure DOFs
    Teuchos::RCP<FLD::UTILS::FsiMapExtractor> fsiinterface_;
  };
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
