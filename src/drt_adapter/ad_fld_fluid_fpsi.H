/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fpsi. Can only be used in conjunction with #FluidImplicitTimeInt

\level 2

\maintainer  Johannes Kremheller
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FPSI_H
#define ADAPTER_FLUID_FPSI_H

#include "ad_fld_fluid_fsi.H"

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
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    virtual void Init();

    /// communication object at the interface without pressure dofs for FPSI problems
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const
    {
      return fpsiinterface_;
    }

    virtual void SetupInterface(const int nds_master = 0);

    virtual void UseBlockMatrix(bool splitmatrix);
    virtual void UseBlockMatrix(
        bool splitmatrix, Teuchos::RCP<FPSI::UTILS::MapExtractor> const& shapederivSplitter);

   protected:
    /// the interface map setup for fpsi interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> fpsiinterface_;

  };  // class FluidFPSI
}  // namespace ADAPTER

#endif
