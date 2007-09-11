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
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;


namespace XFEM
{

/*!
 \class Enrichment

 \brief Corresponds to a specific enrichment with unique id and type

 */
class Enrichment
{
public:

    // Various kinds of enrichments, which differ in ther way,
    // how the enrichment function is evaluated
    enum Type
    {
        typeStandard,
        typeJump,
        typeVoid
    };
    
    Enrichment(
            const int id,
            const Type type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    void Print(
            ostream& os) const;

    inline int getId() const
    {
        return id_;
    };
    
    inline Type getType() const
    {
        return type_;
    };

private:
    // a global id used to distuinguish several surfaces within one element
    // might become or contain reference to a certain condition
    int id_;
    // the specific type of enrichment
    Type type_;
}; // class Enrichment

};


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
