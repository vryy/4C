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

#include "xfem.H"
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
      field_(physvar),
      enr_(enr)
    {
      return;
    }
        
    //! copy constructor
    FieldEnr(
            const FieldEnr& other           ///< source
            ) :
              field_(other.field_), 
              enr_(other.enr_)
    {
      assert(&other != this);
      return;
    }

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
//        return (!(*this == rhs));
        if (field_ != rhs.field_ or getEnrichment() != rhs.getEnrichment())
            return true;
        else
            return false;
    }
    

  private:
    
    //! hidden default constructor
    explicit FieldEnr();
  
    //! physical field variable to be enriched
    const PHYSICS::Field field_;
    //! specific enrichment for this physical field
    const Enrichment enr_;

}; // class FieldEnr

} // namespace XFEM

#endif  // #ifndef FIELD_ENRICHED_H
#endif  // #ifdef CCADISCRET
