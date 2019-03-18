/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fluid_xfsi.H

\brief Fluid field adapter for XFSI allowing multiple fluid discretizations

Can only be used in conjunction with XFluidFluid!

\level 2

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 </pre>
*/
/*----------------------------------------------------------------------*/

#ifndef AD_FLD_FLUID_FLUID_XFSI_H_
#define AD_FLD_FLUID_FLUID_XFSI_H_

#include "ad_fld_fluid_xfsi.H"

namespace FLD
{
  class XFluidFluid;
}

namespace ADAPTER
{
  class FluidFluidXFSI : public XFluidFSI
  {
   public:
    /// Constructor
    FluidFluidXFSI(Teuchos::RCP<Fluid> fluid,
        const std::string coupling_name_xfsi,  // name of the FSI coupling condition
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize algorithm
    virtual void Init();

   protected:
    /// A casted pointer to a fluid with multiple discretizations
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;
  };
}  // namespace ADAPTER


#endif /* AD_FLD_FLUID_FLUID_XFSI_H_ */
