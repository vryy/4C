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

#include "enrichment.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_geometry/intersection_service.H"
#include "interfacexfsi.H"
#include <string>
#include <sstream>



/*----------------------------------------------------------------------*
 |  create string                                                       |
 *----------------------------------------------------------------------*/
std::string XFEM::Enrichment::toString() const
{
    std::stringstream s;
    int width = 1;
    if (xfemconditionlabel_ > 9)
      width = 2;
      
    s << "Enr(" << setw(width) << xfemconditionlabel_ << ", " << enrTypeToString(type_) << ")";
    return s.str();
}

std::string XFEM::Enrichment::enrTypeToString(const EnrType type) const
{
    std::string typetext;
    switch (type){
        case typeStandard:  typetext = "Stnd"; break;
        case typeJump:      typetext = "Jump"; break;
        case typeVoid:      typetext = "Void"; break;
        default: dserror("no string defined for EnrType");
    };
    return typetext;
}


/*----------------------------------------------------------------------*
 |  get enrichment value                                        ag 11/07|
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::EnrValue(
        const LINALG::Matrix<3,1>&            actpos,
        const XFEM::InterfaceHandle&          ih,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        ) const
{
    // return value
    double enrval = 1.0;
    switch (Type())
    {
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
//                double actpos_enr_val = 0.0;
                if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
                  enrval = 0.0;
                } else {
                  enrval = 1.0;
                }
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
                if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
                  enrval = -1.0;
                } else {
                  enrval = 1.0;
                }
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
double XFEM::Enrichment::ModifiedEnrValue(
        const LINALG::Matrix<3,1>&            actpos,
        const LINALG::Matrix<3,1>&            nodalpos,
        const XFEM::InterfaceHandle&          ih,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        ) const
{
    dserror("needs update for the approach variable");
    // return value
    double enrval = 1.0;
    
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
    case XFEM::Enrichment::typeVoid:
    {
        double actpos_enr_val = 0.0;
        if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
            actpos_enr_val = 0.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (ih.PositionWithinConditionNP(nodalpos) == this->XFEMConditionLabel()) {
            nodepos_enr_val = 0.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        double actpos_enr_val = 0.0;
        if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
            actpos_enr_val = -1.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (ih.PositionWithinConditionNP(nodalpos) == this->XFEMConditionLabel()) {
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
