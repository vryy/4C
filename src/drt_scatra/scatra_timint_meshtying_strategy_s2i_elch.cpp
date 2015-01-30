/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.cpp

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_dserror.H"

#include "scatra_timint_meshtying_strategy_s2i_elch.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch* elchtimint
    ) :
MeshtyingStrategyS2I(elchtimint)
{
  return;
} // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*--------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions (electrochemistry)   fang 10/14 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying() const
{
  // safety check
  if(DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()),"BLOCKPRECOND"))
    dserror("Block preconditioning doesn't work for scatra-scatra interface coupling yet!");

  SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying();

  return;
} // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying
