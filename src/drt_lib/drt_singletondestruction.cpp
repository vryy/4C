
#include "drt_singletondestruction.H"
#include "drt_globalproblem.H"

DRT::SingletonDestruction::SingletonDestruction()
{
  DRT::Problem::Instance()->Register(this);
}
