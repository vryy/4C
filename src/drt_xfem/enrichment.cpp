/*!
\file enrichment.cpp

\brief describes the enrichment class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "enrichment.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "intersection_service.H"
#include "interface.H"
#include <string>
#include <sstream>

using namespace XFEM;


/*----------------------------------------------------------------------*
 |  create string                                                       |
 *----------------------------------------------------------------------*/
std::string Enrichment::toString() const
{
    std::stringstream s;
    s << "Enrichment XFEMConditionLabel: " << this->xfemconditionlabel_ << ", type: " << enrTypeToString(this->type_);
    return s.str();
}

std::string Enrichment::enrTypeToString(const EnrType type) const
{
    std::string typetext;
    switch (type){
        case typeStandard:  typetext = "Standard"; break;
        case typeJump:      typetext = "Jump    "; break;
        case typeVoid:      typetext = "Void    "; break;
        default: dserror("no string defined for EnrType");
    };
    return typetext;
}


/*----------------------------------------------------------------------*
 |  get enrichment value                                        ag 11/07|
 *----------------------------------------------------------------------*/
double Enrichment::EnrValue(
        const BlitzVec3&                actpos,
        const XFEM::InterfaceHandle&    ih,
        const XFEM::Enrichment::ApproachFrom approachdirection
        ) const
{
    // return value
    double enrval = 1.0;
    map<int,bool> posInCondition;
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
    case XFEM::Enrichment::typeVoid:
    {
        switch (approachdirection)
        {
            case approachFromPlus:
            {
                enrval = 1.0;
                break;
            }
            case approachFromMinus:
            {
                enrval = 0.0;
                break;
            }
            case approachUnknown:
            {
                const int xfemcondition_label = this->XFEMConditionLabel();
                PositionWithinCondition(actpos,ih,posInCondition);
                double actpos_enr_val = 0.0;
                if (posInCondition.find(xfemcondition_label)->second) {
                    actpos_enr_val = 0.0;
                } else {
                    actpos_enr_val = 1.0;
                }
                enrval = actpos_enr_val;
                break;
            }
        }
        
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        switch (approachdirection)
        {
            case approachFromPlus:
            {
                enrval = 1.0;
                break;
            }
            case approachFromMinus:
            {
                enrval = -1.0;
                break;
            }
            case approachUnknown:
            {
                const int xfemcondition_label = this->XFEMConditionLabel();
                PositionWithinCondition(actpos,ih,posInCondition);
                double actpos_enr_val = 0.0;
                if (posInCondition.find(xfemcondition_label)->second) {
                    actpos_enr_val = -1.0;
                } else {
                    actpos_enr_val = 1.0;
                }
        
                enrval = actpos_enr_val;
                break;
            }
        }
        break;
    }
    default:
        dserror("unsupported enrichment!");
    }
    return enrval;
}


/*
 *  get modified enrichment value (satisfied interpolation property
 *                                                              ag 11/07
 */
double Enrichment::ModifiedEnrValue(
        const BlitzVec3&                actpos,
        const BlitzVec3&                nodalpos,
        const XFEM::InterfaceHandle&    ih,
        const XFEM::Enrichment::ApproachFrom approachdirection
        ) const
{
    dserror("needs update for the approach variable");
    // return value
    double enrval = 1.0;
    map<int,bool> posInCondition;
    
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
    case XFEM::Enrichment::typeVoid:
    {
        const int xfemcondition_label = this->XFEMConditionLabel();
        
        PositionWithinCondition(actpos,ih,posInCondition);
        double actpos_enr_val = 0.0;
        if (posInCondition.find(xfemcondition_label)->second) {
            actpos_enr_val = 0.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        PositionWithinCondition(nodalpos,ih,posInCondition);
        double nodepos_enr_val = 0.0;
        if (posInCondition.find(xfemcondition_label)->second) {
            nodepos_enr_val = 0.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        const int xfemcondition_label = this->XFEMConditionLabel();
        
        PositionWithinCondition(actpos,ih,posInCondition);
        double actpos_enr_val = 0.0;
        if (posInCondition.find(xfemcondition_label)->second) {
            actpos_enr_val = -1.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        PositionWithinCondition(nodalpos,ih,posInCondition);
        double nodepos_enr_val = 0.0;
        if (posInCondition.find(xfemcondition_label)->second) {
            nodepos_enr_val = -1.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    
    default:
        dserror("unsupported enrichment (modified)!");
    }
    return enrval;
}


#endif  // #ifdef CCADISCRET
