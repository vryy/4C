/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fpsi. Can only be used in conjunction with #FluidImplicitTimeInt

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_ADAPTER_FLD_FLUID_FPSI_H
#define BACI_ADAPTER_FLD_FLUID_FPSI_H

#include "baci_adapter_fld_fluid_fsi.H"

namespace FPSI
{
  namespace UTILS
  {
    class MapExtractor;
  }
}  // namespace FPSI

namespace ADAPTER
{
  /*! \brief Fluid field adapter for fpsi
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
  class FluidFPSI : public FluidFSI
  {
   public:
    /// Constructor
    FluidFPSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    void Init() override;

    /// communication object at the interface without pressure dofs for FPSI problems
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const override
    {
      return fpsiinterface_;
    }

    void SetupInterface(const int nds_master = 0) override;

    void UseBlockMatrix(bool splitmatrix) override;
    virtual void UseBlockMatrix(
        bool splitmatrix, Teuchos::RCP<FPSI::UTILS::MapExtractor> const& shapederivSplitter);

   protected:
    /// the interface map setup for fpsi interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> fpsiinterface_;

  };  // class FluidFPSI
}  // namespace ADAPTER

#endif
