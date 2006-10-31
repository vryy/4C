/*!----------------------------------------------------------------------
\file element.cpp
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

#include "element.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::Element(int id, ElementType etype) :
ParObject(),
id_(id),
etype_(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::Element(const CCADISCRETIZATION::Element& old) :
ParObject(old),
id_(old.id_),
etype_(old.etype_)
{

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::~Element()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CCADISCRETIZATION::Element& element)
{
  element.Print(); 
  return os;
}


/*----------------------------------------------------------------------*
 |  read element input dummy (public)                        mwgee 11/06|
 |  this is a base class dummy for elements that do not need            |
 |  a reading method
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Element::ReadElement()
{
  cout << "CCADISCRETIZATION::Element::ReadElement:\n"
       << "Base class dummy routine Element::ReadElement() called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return false;
}















#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
