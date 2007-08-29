/*!
\file xfem_enrichment.H

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_XFEM_ENRICHMENT_H
#define DRT_XFEM_ENRICHMENT_H


// namespace for all enrichment related stuff
namespace Enrichments
{

// Various kinds of enrichments, which differ in ther way,
// how the enrichment function is evaluated
enum EnrType
{
    enrtypeStandard,
    enrtypeJump,
    enrtypeVoid
};

// Special Id that groups degrees of freedom by cause, 
// e.g. a jump resulting from interface 2 -> EnrId 5
typedef int EnrId;

/*!
 \class Enrichment

 \brief Corresponds to a specific enrichment with unique id and type

 */
class Enrichment
{
public:

    Enrichment(
            const int enr_id,
            const EnrType enr_type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    void Print(
            std::ostream& os) const;

    inline EnrId getId() const
    {
        return enr_id_;
    };
    
    inline EnrType getType() const
    {
        return enr_type_;
    };

private:
    // a global id used to distuinguish several surfaces within one element
    EnrId enr_id_;
    // the specific type of enrichment
    EnrType enr_type_;
}; // class Enrichment

}; // namespace Enrichments

#endif  // #ifndef DRT_XFEM_ENRICHMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef XFEM
