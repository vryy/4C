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

  void SetElementXFEMParameter(
      Teuchos::ParameterList& params   //> parameter list
  );


  INPAR::CUT::VCellGaussPts VolumeCellGaussPoints(){ return vcellgausspts_; };
  INPAR::CUT::BCellGaussPts BoundaryCellGaussPoints(){ return bcellgausspts_; };

private:

  INPAR::CUT::VCellGaussPts vcellgausspts_;
  INPAR::CUT::BCellGaussPts bcellgausspts_;

protected:

  /// protected Constructor since we are a Singleton.
  FluidEleParameterXFEM();

};

}
}

#endif
