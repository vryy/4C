/*!
\file dof_management.H

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef DOF_MANAGEMENT_H
#define DOF_MANAGEMENT_H

#include "physics.H"
#include "enrichment.H"
#include "integrationcell.H"
#include "../drt_lib/drt_discret.H"

using namespace Physics;

// namespace for all xfem related stuff
namespace XFEM
{


/*!
 \class EnrPhysVar

 \brief Corresponds to a specific enriched physical field variable

 type EnrPhysVar = (PhysVar, Enrichment)

 */
class EnrPhysVar
{
public:

    EnrPhysVar(
            const PhysVar physvar,
            const Enrichment enr);
    EnrPhysVar(
            const EnrPhysVar& other);
    ~EnrPhysVar();

    string toString() const;

    inline PhysVar getPhysVar() const
    {
        return physvar_;
    };

    inline Enrichment getEnrichment() const
    {
        return enr_;
    };
    
    inline bool operator <(const EnrPhysVar& rhs) const
    {
        if (getPhysVar() < rhs.getPhysVar())
            return true;
        else if (getPhysVar() > rhs.getPhysVar())
            return false;
        else
        {
            if (getEnrichment() < rhs.getEnrichment())
                return true;
            else
                return false;
        }
    } 
    

private:
    // physical field variable to be enriched
    PhysVar physvar_;
    // specific enrichment for this physical field
    Enrichment enr_;

}; // class EnrPhysVar

const map<int, const set<EnrPhysVar> > createNodalDofMap(
        RefCountPtr<DRT::Discretization> xfemdis,
        const map<int, vector<IntCell> >&  elementIntCellMap);

}; // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
