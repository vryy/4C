/*!
\file enrichment.H

\brief describes the enrichment class

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
#include "../drt_lib/drt_discret.H"

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
            const int xfemconditionlabel,
            const EnrType type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    
    std::string toString() const;

    std::string enrTypeToString(const EnrType type) const;
    
    inline int XFEMConditionLabel() const
    {
        return xfemconditionlabel_;
    };
    
    inline EnrType Type() const
    {
        return type_;
    };

    inline bool operator ==(const Enrichment& rhs) const
    {
        if (XFEMConditionLabel() == rhs.XFEMConditionLabel() and Type() == rhs.Type())
            return true;
        else
            return false;
    };
    
    inline bool operator <(const Enrichment& rhs) const
    {
        if (XFEMConditionLabel() < rhs.XFEMConditionLabel())
            return true;
        else if (XFEMConditionLabel() > rhs.XFEMConditionLabel())
            return false;
        else
        {
            if (Type() < rhs.Type())
                return true;
            else
                return false;
        }
    }

    double ModifiedEnrValue(
            const blitz::Array<double,1>& actpos,
            const blitz::Array<double,1>& nodalpos,
            const RCP<DRT::Discretization>& cutterdis
            ) const;

    double EnrValue(
            const blitz::Array<double,1>& actpos,
            const RCP<DRT::Discretization>& cutterdis
            ) const;
    
private:
    // a global id used to distuinguish several surfaces within one element
    // might become or contain reference to a certain condition
    int xfemconditionlabel_;
    // the specific type of enrichment
    EnrType type_;
}; // class Enrichment

}


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef CCADISCRET
