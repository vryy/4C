/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_xfem.H

\brief Setting of specific XFEM based fluid parameter for element evaluation

<pre>
Maintainers: Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_XFEM_H
#define FLUID_ELE_PARAMETER_XFEM_H

#include "fluid_ele_parameter_std.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"


namespace DRT
{
namespace ELEMENTS
{

class FluidEleParameterXFEM : public FluidEleParameterStd
{

public:

  virtual ~FluidEleParameterXFEM() {}

  /// Singleton access method
  static FluidEleParameterXFEM* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /// set all the XFEM specific parameters
  void SetElementXFEMParameter(
      Teuchos::ParameterList& params   ///< parameter list
  );

  /*----------------------------------------------------*/
  //! @name access routines for integration on cut elements
  /*----------------------------------------------------*/

  //! get the volumecell integration method used for integrating cut elements?
  const INPAR::CUT::VCellGaussPts VolumeCellGaussPoints(){ return vcellgausspts_; };

  //! get the boundarycell integration method used for integrating the surface in cut elements?
  const INPAR::CUT::BCellGaussPts BoundaryCellGaussPoints(){ return bcellgausspts_; };


  /*----------------------------------------------------*/
  //! @name access routines for the viscous stabilization in Nitsche's method and MixedHybrid_LM methods
  /*----------------------------------------------------*/

  //! get the type of how to estimate the scaling of the trace inequality used for the viscous part of Nitsche's method?
  INPAR::XFEM::ViscStab_TraceEstimate ViscStabTracEstimate(){ return visc_stab_trace_estimate_; };

  //! get the element length definition used for viscous part of the penalty term in Nitsche's method
  INPAR::XFEM::ViscStab_hk ViscStabHK(){ return visc_stab_hk_; };

  //! get the dimensionless user defined scaling for the viscous part of the penalty term in Nitsche's method and
  //! scaling factor for the MHVS method (then gamma = 1/n, see publications)
  const double ViscStabGamma(){ return visc_stab_gamma_; };

private:

  /*----------------------------------------------------*/
  //! @name parameters for integration on cut elements
  /*----------------------------------------------------*/

  //! which volumecell integration is used for integrating cut elements?
  INPAR::CUT::VCellGaussPts vcellgausspts_;

  //! which boundarycell integration is used for integrating the surface in cut elements?
  INPAR::CUT::BCellGaussPts bcellgausspts_;


  /*----------------------------------------------------*/
  //! @name parameters for the viscous stabilization in Nitsche's method and MixedHybrid_LM methods
  /*----------------------------------------------------*/

  //! how to estimate the scaling of the trace inequality used for the viscous part of Nitsche's method?
  INPAR::XFEM::ViscStab_TraceEstimate visc_stab_trace_estimate_;

  //! element length definition used for viscous part of the penalty term in Nitsche's method
  INPAR::XFEM::ViscStab_hk visc_stab_hk_;

  //! dimensionless user defined scaling for the viscous part of the penalty term in Nitsche's method and
  //! scaling factor for the MHVS method (then gamma = 1/n, see publications),
  //! this parameter shall be independent of element type, shape and polynomial degree
  double visc_stab_gamma_;

protected:

  /// protected Constructor since we are a Singleton.
  FluidEleParameterXFEM();

};

}
}

#endif
