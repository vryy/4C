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
#include "../drt_lib/drt_utils.H" //URSULA
#include "../drt_fem_general/drt_utils_gder2.H" //URSULA
#include <string>
#include <sstream>




/*----------------------------------------------------------------------*
 | ASSIGNMENT OPERATOR                                       u.may 04/09|
 *----------------------------------------------------------------------*/
XFEM::Enrichment& XFEM::Enrichment::operator = (const XFEM::Enrichment& old) 
{
  xfemconditionlabel_ = old.xfemconditionlabel_;
  type_ = old.type_;
  return *this;
}



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
        case typeKink:      typetext = "Kink"; break;
        default: dserror("no string defined for EnrType");
    };
    return typetext;
}


/*----------------------------------------------------------------------*
 | get enrichment value                                        ag 11/07 |
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
        // standard Heaviside function
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
        dserror("Use enrichment function for level set EnrValueLevelSet() instead!");
        // Heaviside function (jump height is 2!)
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
    case XFEM::Enrichment::typeKink:
    {
        dserror("Use enrichment functions based on level set function instead!");
    }
    default:
        dserror("unsupported type of enrichment!");
    }
    return enrval;
}


/*----------------------------------------------------------------------*
 | get enrichment value for integration cell (constant!)    henke 07/09 |
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::EnrValueIntCell(const GEO::DomainIntCell& cell) const
{
    //std::cout << "EnrValueIntCell()" << std::endl;
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
        if (cell.getDomainPlus()) {
            enrval = 1.0;
        } else {
            enrval = 0.0;
        }
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        if (cell.getDomainPlus()) {
            enrval = 1.0;
        } else {
            enrval = -1.0;
        }
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
        dserror("enrichment value not constant within integration cell for kink enrichment!");
    }
    default:
        dserror("unsupported type of enrichment!");
    }
    return enrval;
}


/*----------------------------------------------------------------------*
 | get modified enrichment value                               ag 11/07 |
 | remark: the enrichment function is 0 at nodes and hence the usual    |
 | interpolation property of the standard FEM is satisfied              |
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::ModifiedEnrValue(
        const LINALG::Matrix<3,1>&            actpos,
        const LINALG::Matrix<3,1>&            nodalpos,
        const XFEM::InterfaceHandle&          ih,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        ) const
{
    dserror("function not in use; needs update for the approach variable");
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
        dserror("Do not evaluate jump enrichment values using ModifiedEnrValue()!");
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
        dserror("Do not evaluate kink enrichment values using ModifiedEnrValue()!");
        break;
    }
    default:
        dserror("unsupported type of (modified) enrichment!");
    }
    return enrval;
}


#endif  // #ifdef CCADISCRET
