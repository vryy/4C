/*!----------------------------------------------------------------------
\file inpar_pasi.H

\brief Input parameters for particle structure interaction problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/
#ifndef INPAR_PASI_H
#define INPAR_PASI_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace DRT
{

namespace INPUT
{

class ConditionDefinition;

} // namespace INPUT

} // namespace DRT

/*----------------------------------------------------------------------*
 | Input parameters for particle structure interaction   sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace INPAR
{

namespace PASI

{

/*----------------------------------------------------------------------*
 | Coupling Methods                                                     |
 *----------------------------------------------------------------------*/
enum _PASI_COUPLING
{
  partitioned_onewaycoup
};

//! set valid parameters for particle structure interaction
void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

//! set valid conditions for particle structure interaction
void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

} // namespace PASI

} // namespace INPAR

#endif // #ifndef INPAR_PASI_H
