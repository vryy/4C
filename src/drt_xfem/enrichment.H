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

    /*! Various kinds of enrichments, which differ in the way,
     *  the enrichment function is evaluated
     */
    enum EnrType
    {
        typeUndefined,
        typeStandard,
        typeJump,
        typeVoid
    };
    
    /*!
     * we want to be able to say, wether we approach froma particular side towards the discontinuity
     * this makes sense for jump or kink enrichments,
     * boundary layers or similar things might not need this
     */
    enum ApproachFrom
    {
        approachUnknown,
        approachFromPlus,
        approachFromMinus
    };
    
    explicit Enrichment();
    
    explicit Enrichment(
            const int xfemconditionlabel,
            const EnrType type);
    Enrichment(
            const Enrichment& other);
    ~Enrichment();
    
    std::string toString() const;

    std::string enrTypeToString(const EnrType type) const;
    
    int XFEMConditionLabel() const
    {
        return xfemconditionlabel_;
    };
    
    EnrType Type() const
    {
        return type_;
    };

    bool operator ==(const Enrichment& rhs) const
    {
        if (XFEMConditionLabel() == rhs.XFEMConditionLabel() and Type() == rhs.Type())
            return true;
        else
            return false;
    };
    
    bool operator <(const Enrichment& rhs) const
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
            const RCP<DRT::Discretization>& cutterdis,
            const XFEM::Enrichment::ApproachFrom approachdirection
            ) const;

    double EnrValue(
            const blitz::Array<double,1>& actpos,
            const RCP<DRT::Discretization>& cutterdis,
            const XFEM::Enrichment::ApproachFrom approachdirection
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
