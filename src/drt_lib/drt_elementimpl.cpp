
#include "drt_elementimpl.H"
#include "drt_globalproblem.H"

DRT::ELEMENTS::ElementImpl::ElementImpl()
{
  DRT::Problem::Instance()->Register(this);
}
