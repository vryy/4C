/*!----------------------------------------------------------------------
\file drt_dofset_base.cpp
\brief A set of degrees of freedom

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset_base.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "../drt_io/io_pstream.H"

#include "../linalg/linalg_utils.H"


// list of all dof sets
std::list<DRT::DofSetInterface*> DRT::DofSetBase::static_dofsets_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::DofSetBase()
    : DofSetInterface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::~DofSetBase()
{
  // disconnect from proxies
  for (std::list<DofSetInterface*>::iterator i=registered_dofsets_.begin(); i!=registered_dofsets_.end(); ++i)
  {
    (*i)->Disconnect(this);
  }
  // remove dofset from static list if necessary
  static_dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::AddDofSettoList()
{
  if (std::find(static_dofsets_.begin(),static_dofsets_.end(),this)==static_dofsets_.end())
  {
    static_dofsets_.push_back(this);
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::ReplaceInStaticDofsets(Teuchos::RCP<DofSetInterface> olddofset)
{
  std::list<DRT::DofSetInterface*>::iterator iterold = std::find(static_dofsets_.begin(),static_dofsets_.end(),&(*olddofset));
  if (iterold == static_dofsets_.end())
  {
    static_dofsets_.push_back(this);
  }
  else
  {
    static_dofsets_.insert(iterold, this);
    static_dofsets_.remove(&(*olddofset));
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxGIDinList(const Epetra_Comm& comm) const
{
  int count = -1;
  for (std::list<DofSetInterface*>::const_iterator i=static_dofsets_.begin();
       i!=static_dofsets_.end();
       ++i)
  {
    if (*i==this)
      break;

    // ignore empty (no yet initialized) dof row maps
    if ((*i)->Initialized())
    {
      count = std::max((*i)->MaxAllGID(),count);
    }
  }
  int max;
  comm.MaxAll(&count,&max,1);
  return max;
}

void DRT::DofSetBase::PrintAllDofsets(const Epetra_Comm& comm) const
{
  if (comm.MyPID() == 0)
  {
    std::vector<int> min;
    std::vector<int> max;
    for (std::list<DofSetInterface*>::const_iterator i=static_dofsets_.begin(); i!=static_dofsets_.end(); ++i)
    {
      min.push_back((*i)->MinAllGID());
      max.push_back((*i)->MaxAllGID());
    }
    if (min.size() < 1)
      return;

    std::vector<int>::const_iterator largest = max_element( max.begin(), max.end() );
    int allmax = *largest;
    int availspace = 80;
    int arrowlen = availspace + 3;

    if (min[0] > 0)
      availspace -= 2;

    for (unsigned int i = 0; i < min.size(); ++i)
    {
      // left bar
      availspace--;
      // number
      availspace -= 12;
      //right bar
      if (i+1 < min.size())
      {
        if (max[i]+1 < min[i+1])
          availspace -= 2;
      }
      else
        availspace--;
    }

    if (availspace < 0)
    {
      arrowlen -= availspace;
      availspace = 0;
    }

    IO::cout << "All registered dofsets:" << IO::endl;

    if (min[0] > 0)
      IO::cout << "| ";
    for (unsigned int i = 0; i < min.size(); ++i)
    {
      // left bar
      if (i > 0 and max[i-1] > min[i])
        IO::cout << "X";
      else
        IO::cout << "|";

      // number
      IO::cout << std::left << std::setw(6) << min[i];
      for (int j = 0; j < (int)(max[i]-min[i])*1.0*availspace/allmax; ++j)
        IO::cout << " ";
      IO::cout << std::right << std::setw(6) << max[i];

      //right bar
      if (i+1 < min.size())
      {
        if (max[i]+1 < min[i+1])
          IO::cout << "| ";
      }
      else
      {
        IO::cout << "|";
      }
    }

    IO::cout << IO::endl << "+";
    for (int i = 0; i < arrowlen; ++i)
      IO::cout << "-";
    IO::cout << "> DofGID" << IO::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::Register(DofSetInterface* dofset)
{
  registered_dofsets_.push_back(dofset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::Unregister(DofSetInterface* dofset)
{
  registered_dofsets_.remove(dofset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::NotifyAssigned()
{
  for (std::list<DofSetInterface*>::iterator i=registered_dofsets_.begin(); i!=registered_dofsets_.end(); ++i)
    (*i)->NotifyAssigned();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::NotifyReset()
{
  for (std::list<DofSetInterface*>::iterator i=registered_dofsets_.begin(); i!=registered_dofsets_.end(); ++i)
    (*i)->NotifyReset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::Print(std::ostream& os) const
{
  dserror("Print() is not implemented in base class. Override Print() in subclass");
}
