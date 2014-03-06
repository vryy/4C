/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_poro.H

\brief Evaluation of general fluid parameter for porofluid

<pre>
Maintainers: Anh-Tu Vuong and Andreas Rauch
             {vuong,rauch}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15251/-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_PORO_H
#define FLUID_ELE_PARAMETER_PORO_H

#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_parameter.H"
#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
namespace ELEMENTS
{
class FluidEleParameterPoro : public FluidEleParameter
{

public:

  virtual ~FluidEleParameterPoro() {}

  /// Singleton access method
  static FluidEleParameterPoro* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  void SetElementPoroParameter( Teuchos::ParameterList& params, int myrank);

  //! print parameter to screen
  void PrintFluidParameterPoro();

  //! partial integration of porosity gradient in continuity equation
  bool PoroContiPartInt(){return poro_conti_partint_;};


private:

  //! Flag SetGeneralParameter was called
  bool set_fluid_parameter_poro_;

  /// parameter to switch the stabilization
  INPAR::FLUID::StabType stabtype_;

  //! partial integration of porosity gradient in continuity equation
  bool poro_conti_partint_;

  //! pointer to general fluid time parameter
  DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;

  /// private Constructor since we are a Singleton.
  FluidEleParameterPoro();

};

}
}

#endif
