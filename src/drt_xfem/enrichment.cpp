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
#include <string>
#include <sstream>

using namespace XFEM;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const int xfemconditionlabel,
        const EnrType type) :
            xfemconditionlabel_(xfemconditionlabel),
            type_(type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const Enrichment& other) :
            xfemconditionlabel_(other.xfemconditionlabel_),
            type_(other.type_)
{
    assert(&other != this);
    return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::~Enrichment()
{
    return;
}

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
        const blitz::Array<double,1>& actpos,
        const RCP<DRT::Discretization>& cutterdis
        ) const
{
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
        const int xfemcondition_label = this->XFEMConditionLabel();
        
        double actpos_enr_val = 0.0;
        if (PositionWithinCondition(actpos, xfemcondition_label,cutterdis)) {
            actpos_enr_val = 0.0;
        } else {
            actpos_enr_val = 1.0;
        }

        enrval = actpos_enr_val;
        
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
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& nodalpos,
        const RCP<DRT::Discretization>& cutterdis
        ) const
{
    // return value
    double enrval = 1.0;
    
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
//    case XFEM::Enrichment::typeVoid:
//    {
//        const int xfemcondition_label = this->XFEMConditionLabel();
//        
//        double actpos_enr_val = 0.0;
//        if (PositionWithinCondition(actpos, xfemcondition_label,cutterdis)) {
//            actpos_enr_val = 0.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (PositionWithinCondition(nodalpos, xfemcondition_label,cutterdis)) {
//            nodepos_enr_val = 0.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        
//        break;
//    }
    case XFEM::Enrichment::typeJump:
    {
        const int xfemcondition_label = this->XFEMConditionLabel();
        
        double actpos_enr_val = 0.0;
        if (PositionWithinCondition(actpos, xfemcondition_label,cutterdis)) {
            actpos_enr_val = -1.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (PositionWithinCondition(nodalpos, xfemcondition_label,cutterdis)) {
            nodepos_enr_val = -1.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    
    default:
        dserror("unsupported enrichment!");
    }
    return enrval;
}


#endif  // #ifdef CCADISCRET
