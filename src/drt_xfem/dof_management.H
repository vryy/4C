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
#ifdef TRILINOS_PACKAGE
#ifndef DOF_MANAGEMENT_H
#define DOF_MANAGEMENT_H

#include <physics.H>
#include <enrichment.H>

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
            const PhysVar var,
            const Enrichment enr);
    EnrPhysVar(
            const EnrPhysVar& other);
    ~EnrPhysVar();
    void Print(
            std::ostream& os) const;

    inline PhysVar getPhysVar() const
    {
        return var_;
    };
    
    inline Enrichment getEnrichment() const
    {
        return enr_;
    };

private:
    // physical field variable to be enriched
    PhysVar var_;
    // specific enrichment for this physical field
    Enrichment enr_;
    
}; // class EnrPhysVar

}; // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
