/*!
\file field_enriched.cpp

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "field_enriched.H"

/*----------------------------------------------------------------------*
 |  transform to a string                                       ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::FieldEnr::toString() const
{
  std::stringstream s;
  s << "(" << PHYSICS::physVarToString(field_) << ", " << enr_.toString() << ")";
  return s.str();
}



#endif  // #ifdef CCADISCRET
