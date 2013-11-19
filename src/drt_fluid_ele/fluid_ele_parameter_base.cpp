/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl_parameter.cpp

\brief Evaluation of general fluid parameter

FluidEleParameter::SetParameter(Teuchos::ParameterList& params)
set all general fluid parameter once for all elements.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include "fluid_ele_parameter.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_io/io_pstream.H"

//----------------------------------------------------------------------*/
// private constructor of FluidEleParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterBase::FluidEleParameterBase()
  :
  set_general_fluid_parameter_(false),
  physicaltype_(INPAR::FLUID::incompressible),
  stabtype_(INPAR::FLUID::stabtype_nostab), // stabilization parameters
  is_conservative_(false),
  is_newton_(false),
  is_inconsistent_(false),
  reaction_(false),
  darcy_(false),
  reaction_topopt_(false)
{
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
}

//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterBase::SetElementGeneralFluidParameter( Teuchos::ParameterList& params, int myrank )
{
  if(set_general_fluid_parameter_ == false)
    set_general_fluid_parameter_ = true;
  // For turbulent inflow generation,
  // this function is indeed two times called.
  // In this sepcial case, calling this function twice
  // is ok!
  else
  {
    if (myrank == 0)
      std::cout << std::endl <<
      (" Warning: general fluid parameters should be set only once!!\n "
       " If you run a turbulent inflow generation, calling this function twice is ok!\n ")
      << std::endl << std::endl;
  }

  // set flag for type of linearization (fixed-point-like or Newton)
   //std::string newtonstr   = params.get<std::string>("Linearisation");
   if (DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(params, "Linearisation")==INPAR::FLUID::Newton)
     is_newton_       = true;

   // set flags for formuation of the convective velocity term (conservative or convective)
   std::string convformstr = params.get<std::string>("form of convective term");
   if (convformstr =="conservative")
   {
     is_conservative_ = true;
     if(myrank==0)
     {
       std::cout << std::endl << "Warning: \n"
         "a) Using PSPG stabilization yields a conservative formulation (Hughes & Wells 2005)\n"
         "b) Instablities may occur for complex flow situations" << std::endl;
     }
   }

   // set flag for physical type of fluid flow
   physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params, "Physical Type");
   if (
        (
             (physicaltype_ == INPAR::FLUID::loma)
          or (physicaltype_ == INPAR::FLUID::varying_density)
        )
       and (fldparatimint_->IsStationary() == true)
      )
     dserror("physical type is not supported in stationary FLUID implementation.");

} // set general parameters
