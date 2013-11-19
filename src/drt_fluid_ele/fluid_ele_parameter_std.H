/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_loma.H

\brief Evaluation of general fluid parameter for loma

       As it is currently not feasible to split std and loma terms
       this class holds the loma specific paramters. Though this
       class provides the construction methods of the stdfluid it
       is nice to keep the loma parameters separated, since they are
       definitely not needed for poro / fpsi problems.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_STD_H
#define FLUID_ELE_PARAMETER_STD_H

#include "fluid_ele_parameter.H"


namespace DRT
{
namespace ELEMENTS
{

class FluidEleParameterStd : public FluidEleParameter
{

public:

  virtual ~FluidEleParameterStd() {}

  /// Singleton access method
  static FluidEleParameterStd* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /// set loma parameters
  void SetElementLomaParameter(Teuchos::ParameterList& params);

  //! some loma specific parameters
  //! flag for material update
  virtual bool UpdateMat(){return update_mat_;};
  //! flag to (de)activate continuity SUPG term
  virtual bool ContiSUPG(){return conti_supg_;};
  //! flag to (de)activate continuity cross-stress term -> residual-based VMM
  virtual INPAR::FLUID::CrossStress ContiCross(){return conti_cross_;};
  //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
  virtual INPAR::FLUID::ReynoldsStress ContiReynolds(){return conti_reynolds_;};
  //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
  virtual bool MultiFracLomaConti(){return multifrac_loma_conti_;};


private:

  //! flag for material update
  bool update_mat_;
  //! flag to (de)activate continuity SUPG term
  bool conti_supg_;
  //! flag to (de)activate continuity cross-stress term -> residual-based VMM
  INPAR::FLUID::CrossStress conti_cross_;
  //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
  INPAR::FLUID::ReynoldsStress conti_reynolds_;
  //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
  bool multifrac_loma_conti_;

protected:

  /// protected Constructor since we are a Singleton.
  FluidEleParameterStd();

};

}
}

#endif
