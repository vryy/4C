/*!
\file field_enriched.cpp

\brief provides a class that represents an enriched physical scalar field

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 */


#include "enrichment.H"
#include "field_enriched.H"


//! constructor
XFEM::FieldEnr::FieldEnr(
    const PHYSICS::Field  physvar,  ///< physical variable
    const XFEM::Enrichment&     enr       ///< enrichment
) :
enr_(enr),
field_(physvar)
{
  return;
}


//! copy constructor
XFEM::FieldEnr::FieldEnr(
    const XFEM::FieldEnr& other           ///< source
) :
enr_(other.enr_),
field_(other.field_)
{
  assert(&other != this);
  return;
}


//! default constructor
XFEM::FieldEnr::FieldEnr(
) :
enr_(XFEM::Enrichment()),
field_(XFEM::PHYSICS::undefinedField)
{
  return;
}


/*----------------------------------------------------------------------*
 |  assignment operator                                      u.may 04/09|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr& XFEM::FieldEnr::operator = (const XFEM::FieldEnr& old)
{
  field_ = old.field_;
  enr_  = old.enr_;
  return *this;
}




/*----------------------------------------------------------------------*
 |  transform to a std::string                                       ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::FieldEnr::toString() const
{
  std::stringstream s;
  s << "(" << PHYSICS::physVarToString(field_) << ", " << enr_.toString() << ")";
  return s.str();
}


//! return physical field
XFEM::PHYSICS::Field XFEM::FieldEnr::getField() const
{
  return field_;
}


//! return enrichement for the field
const XFEM::Enrichment& XFEM::FieldEnr::getEnrichment() const
{
  return enr_;
}


//! order is given first by field, then by the enrichment
bool XFEM::FieldEnr::operator <(const FieldEnr& rhs) const
{
  if (field_ < rhs.field_)
    return true;
  else if (field_ > rhs.field_)
    return false;
  else
  {
    if (enr_ < rhs.enr_)
      return true;
    else
      return false;
  }
}


//! equality, if everything is equal
bool XFEM::FieldEnr::operator ==(const FieldEnr& rhs) const
{
  if (field_ == rhs.field_ and getEnrichment() == rhs.getEnrichment())
    return true;
  else
    return false;
}


//! inequality if any member differs
bool XFEM::FieldEnr::operator !=(const FieldEnr& rhs) const
{
  if (field_ != rhs.field_ or getEnrichment() != rhs.getEnrichment())
    return true;
  else
    return false;
}

