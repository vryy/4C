/*---------------------------------------------------------------------*/
/*! \file

\brief Registration class for all singletons used within baci

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_singletondestruction.H"
#include "drt_globalproblem.H"

DRT::SingletonDestruction::SingletonDestruction() { DRT::Problem::Instance()->Register(this); }
