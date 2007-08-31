/*!
\file enrichment.H

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;


namespace XFEM
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
            const EnrId enr_id,
            const EnrType enr_type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    void Print(
            ostream& os) const;

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
    // might become or contain reference to a certain condition
    EnrId enr_id_;
    // the specific type of enrichment
    EnrType enr_type_;
}; // class Enrichment

};


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM
