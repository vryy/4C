/*!----------------------------------------------------------------------
\file drt_dofset_base.cpp
\brief A set of degrees of freedom

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
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
std::list<DRT::DofSetBase*> DRT::DofSetBase::static_dofsets_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::DofSetBase()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::~DofSetBase()
{
  static_dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::NumGlobalElements() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::NumGlobalElements(): dofrowmap_ not initialized, yet");
  return dofrowmap_->NumGlobalElements();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxAllGID() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::MaxAllGID(): dofrowmap_ not initialized, yet");
  return dofrowmap_->MaxAllGID();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MinAllGID() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::MinAllGID(): dofrowmap_ not initialized, yet");
  return dofrowmap_->MinAllGID();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSetBase::DofRowMap() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::DofRowMap(): dofrowmap_ not initialized, yet");
  return dofrowmap_.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSetBase::DofColMap() const
{
  if (dofcolmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::DofColMap(): dofcolmap_ not initialized, yet");
  return dofcolmap_.get();
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
void DRT::DofSetBase::ReplaceInStaticDofsets(Teuchos::RCP<DofSetBase> olddofset)
{
  std::list<DRT::DofSetBase*>::iterator iterold = std::find(static_dofsets_.begin(),static_dofsets_.end(),&(*olddofset));
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
bool DRT::DofSetBase::Initialized() const
{
  if (dofcolmap_ == Teuchos::null or dofrowmap_ == Teuchos::null)
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxGIDinList(const Epetra_Comm& comm) const
{
  int count = -1;
  for (std::list<DofSetBase*>::const_iterator i=static_dofsets_.begin();
       i!=static_dofsets_.end();
       ++i)
  {
    if (*i==this)
      break;

    // ignore empty (no yet initialized) dof row maps
    if ((*i)->NumGlobalElements()>0)
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
    for (std::list<DofSetBase*>::const_iterator i=static_dofsets_.begin(); i!=static_dofsets_.end(); ++i)
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
    	IO::cout << "|";
    }

    IO::cout << IO::endl << "+";
    for (int i = 0; i < arrowlen; ++i)
      IO::cout << "-";
    IO::cout << "> DofGID" << IO::endl;
  }
  return;
}

