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
#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include <blitz/array.h>
#include "Teuchos_RefCountPtr.hpp"

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

    // Various kinds of enrichments, which differ in the way,
    // the enrichment function is evaluated
    enum EnrType
    {
        typeStandard,
        typeJump,
        typeVoid
    };
    
    Enrichment(
            const int id,
            const EnrType type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    
    std::string toString() const;

    std::string enrTypeToString(const EnrType type) const;
    
    inline int Id() const
    {
        return id_;
    };
    
    inline EnrType Type() const
    {
        return type_;
    };

    inline bool operator ==(const Enrichment& rhs) const
    {
        if (Id() == rhs.Id() and Type() == rhs.Type())
            return true;
        else
            return false;
    };
    
    inline bool operator <(const Enrichment& rhs) const
    {
        if (Id() < rhs.Id())
            return true;
        else if (Id() > rhs.Id())
            return false;
        else
        {
            if (Type() < rhs.Type())
                return true;
            else
                return false;
        }
    }

    double enrValue(
            const blitz::Array<double,1>& actpos,
            const blitz::Array<double,1>& nodalpos,
            const blitz::Array<double,1>& cellcenterpos
            ) const;
    
private:
    // a global id used to distuinguish several surfaces within one element
    // might become or contain reference to a certain condition
    int id_;
    // the specific type of enrichment
    EnrType type_;
}; // class Enrichment

};


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef CCADISCRET
