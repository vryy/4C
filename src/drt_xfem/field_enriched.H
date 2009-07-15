/*!
\file field_enriched.H

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef FIELD_ENRICHED_H
#define FIELD_ENRICHED_H

#include "physics.H"
#include "enrichment.H"


namespace XFEM
{

/*!
 \brief Corresponds to a specific enriched physical field variable

 Haskell-esque:
 type FieldEnr = (Field, Enrichment)

 */
class FieldEnr
{
  public:

    //! constructor
    explicit FieldEnr(
        const PHYSICS::Field  physvar,  ///< physical variable
        const Enrichment&     enr       ///< enrichment
    ) :
      enr_(enr),
      field_(physvar)
    {
      return;
    }

    //! copy constructor
    FieldEnr(
        const FieldEnr& other           ///< source
        ) :
          enr_(other.enr_),
          field_(other.field_)
    {
      assert(&other != this);
      return;
    }

    //! default constructor
    explicit FieldEnr(
    ) :
      enr_(XFEM::Enrichment()),
      field_(XFEM::PHYSICS::undefinedField)
    {
      return;
    }

    //! assignment operator
    FieldEnr& operator = (const FieldEnr& old);

    //! return string representation
    std::string toString() const;

    //! return physical field
    PHYSICS::Field getField() const {return field_;};

    //! return enrichement for the field
    const Enrichment& getEnrichment() const {return enr_;};

    //! order is given first by field, then by the enrichment
    bool operator <(const FieldEnr& rhs) const
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
    bool operator ==(const FieldEnr& rhs) const
    {
      if (field_ == rhs.field_ and getEnrichment() == rhs.getEnrichment())
        return true;
      else
        return false;
    }
    //! inequality if any member differs
    bool operator !=(const FieldEnr& rhs) const
    {
      if (field_ != rhs.field_ or getEnrichment() != rhs.getEnrichment())
        return true;
      else
        return false;
    }


  private:

    //! specific enrichment for this physical field
    Enrichment        enr_;

    //! physical field variable to be enriched
    PHYSICS::Field    field_;

}; // class FieldEnr

} // namespace XFEM

#endif  // #ifndef FIELD_ENRICHED_H
#endif  // #ifdef CCADISCRET
