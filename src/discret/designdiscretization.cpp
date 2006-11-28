/*!----------------------------------------------------------------------
\file designdiscretization.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "designdiscretization.H"
#include "dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignDiscretization::DesignDiscretization(RefCountPtr<Epetra_Comm> comm) :
Discretization(comm)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignDiscretization::DesignDiscretization(const DRT::DesignDiscretization& old) :
Discretization(old.comm_)
{
  dserror("Do not have copy-ctor");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignDiscretization::~DesignDiscretization()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::DesignDiscretization::FillComplete(
                                       DRT::DesignDiscretization* higher, 
                                       DRT::DesignDiscretization* lower)
{
  // build element to node pointers and maps (there might not be nodes)
  int err = Discretization::FillComplete();
  if (err) return err;
  if (!higher && !lower) return 0;

  // loop elements in this and set pointers to lower elements
  // and set pointers to higher elements
  for (int i=0; i<NumMyColElements(); ++i)
  {
    // we have to cast this to DesignElement
    DRT::DesignElement* element = 
                    dynamic_cast<DRT::DesignElement*>(lColElement(i));
    if (!element) return -1;
    if (lower)
    {
      bool success = element->BuildLowerElementPointers(*lower);
      if (!success) return -1;
    } // if (lower)
    if (higher)
    {
      if (!higher->Filled())
      {
        int err = higher->FillComplete();
        if (err) return err;
      }
      int eleid = element->Id();
      element->hentityid_.resize(0);
      element->horientation_.resize(0);
      element->hentity_.resize(0);
      for (int j=0; j<higher->NumMyColElements(); ++j)
      {
        DesignElement* highele = dynamic_cast<DesignElement*>(higher->lColElement(j));
        if (!highele) return -1;
        int        nents = highele->NumLowerEntityIds();
        const int* ents  = highele->LowerEntityIds();
        for (int k=0; k<nents; ++k)
          if (ents[k] == eleid)
          {
            const int size = element->hentityid_.size();
            element->hentityid_.resize(size+1);
            element->horientation_.resize(size+1);
            element->hentity_.resize(size+1);
            element->hentityid_[size] = highele->Id();
            element->horientation_[size] = 0;
            element->hentity_[size] = highele;
          }
      }
    } // if (higher)
  } // for (int i=0; i<NumMyElements(); ++i)
  filled_ = true;  
  return 0;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
