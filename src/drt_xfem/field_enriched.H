/*!
\file field_enriched.H

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


#ifndef FIELD_ENRICHED_H
#define FIELD_ENRICHED_H


#include "enrichment.H"
#include "physics.H"


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
    );

    //! copy constructor
    FieldEnr(
        const FieldEnr& other           ///< source
    );

    //! default constructor
    explicit FieldEnr();

    //! assignment operator
    FieldEnr& operator = (const FieldEnr& old);

    //! return std::string representation
    std::string toString() const;

    //! return physical field
    PHYSICS::Field getField() const;

    //! return enrichement for the field
    const Enrichment& getEnrichment() const;

    //! order is given first by field, then by the enrichment
    bool operator <(const FieldEnr& rhs) const;

    //! equality, if everything is equal
    bool operator ==(const FieldEnr& rhs) const;

    //! inequality if any member differs
    bool operator !=(const FieldEnr& rhs) const;

  private:

    //! specific enrichment for this physical field
    Enrichment        enr_;

    //! physical field variable to be enriched
    PHYSICS::Field    field_;

}; // class FieldEnr

} // namespace XFEM

#endif  // #ifndef FIELD_ENRICHED_H

