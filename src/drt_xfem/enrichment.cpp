/*!
\file enrichment.cpp

\brief describes the enrichment class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#include "enrichment.H"
#include "../drt_geometry/integrationcell.H"


/*----------------------------------------------------------------------*
 | assignment operator                                       u.may 04/09|
 *----------------------------------------------------------------------*/
XFEM::Enrichment& XFEM::Enrichment::operator = (const XFEM::Enrichment& old)
{
  xfemconditionlabel_ = old.xfemconditionlabel_;
  type_ = old.type_;
  return *this;
}



/*----------------------------------------------------------------------*
 |  create std::string                                                       |
 *----------------------------------------------------------------------*/
std::string XFEM::Enrichment::toString() const
{
    std::stringstream s;
    int width = 1;
    if (xfemconditionlabel_ > 9)
      width = 2;

    s << "Enr(" << enrTypeToString(type_) << ", " << std::setw(width) << xfemconditionlabel_ << ")";
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
        default: dserror("no std::string defined for EnrType");
    };
    return typetext;
}


/*----------------------------------------------------------------------*
 | get enrichment value                                        ag 11/07 |
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::EnrValue(
        const LINALG::Matrix<3,1>&             actpos,
        const COMBUST::InterfaceHandleCombust& ih,
        const XFEM::Enrichment::ApproachFrom   approachdirection
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
    case XFEM::Enrichment::typeJump:
    {
        dserror("Use enrichment functions based on level set function instead!");
        // Heaviside function (jump height is 2!)
//        switch (approachdirection)
//        {
//            case approachFromPlus:
//            {
//                enrval = 1.0;
//                break;
//            }
//            case approachFromMinus:
//            {
//                enrval = -1.0;
//                break;
//            }
//            case approachUnknown:
//            {
//                if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
//                  enrval = -1.0;
//                } else {
//                  enrval = 1.0;
//                }
//                break;
//            }
//        }
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
double XFEM::Enrichment::EnrValueIntCell(const GEO::IntCell& cell) const
{
    double enrval = -777.777;

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
 | get enrichment value at an embedded interface            henke 09/09 |
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::EnrValueAtInterface(
        const XFEM::Enrichment::ApproachFrom  approachdirection
        ) const
{
    // return value
    double enrval = -777.7;

    switch (Type())
    {
    case XFEM::Enrichment::typeStandard:
    {
        //dserror("there is no embedded interface for a standard enrichment!");
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
                dserror("specify side of embedded interface!");
                break;
            }
        }

        break;
    }
    case XFEM::Enrichment::typeJump:
    {
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
                dserror("specify side of embedded interface!");
                break;
            }
        }
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
        dserror("think before you use this for kink enrichments!");
        enrval = 1.0;
        break;
    }
    default:
        dserror("unsupported type of enrichment!");
    }
    return enrval;
}

