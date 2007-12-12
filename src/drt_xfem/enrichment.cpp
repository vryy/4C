/*!
\file enrichment.cpp

\brief describes the enrichment types and classes

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
        const int id,
        const EnrType type) :
    id_(id), type_(type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const Enrichment& other) :
    id_(other.id_), type_(other.type_)
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
    s << "Enrichment id: " << this->id_ << ", type: " << enrTypeToString(this->type_);
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


bool inCircleCylinder(
        const blitz::Array<double,1>& pos,
        const blitz::Array<double,1>& center,
        const double cylinder_radius
        )
{
    blitz::Range _  = blitz::Range::all();
    const blitz::Array<double,1> origincircle(pos(_) - center(_));
    
    const double circle_radius = sqrt(origincircle(0)*origincircle(0) + origincircle(1)*origincircle(1));
    
    bool in_circle = false;
    if (circle_radius <= cylinder_radius){
        in_circle = true;
    } else {
        in_circle = false;
    }
    return in_circle;
}

/*----------------------------------------------------------------------*
 |  get enrichment value                                        ag 11/07|
 *----------------------------------------------------------------------*/
double Enrichment::enrValue(
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& nodalpos,
        const blitz::Array<double,1>& cellcenterpos,
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
//    case XFEM::Enrichment::typeJump:
//    {
//        // TODO: generalize
//        blitz::Array<double,1> center(3);
//        center(0) = 0.6; center(1) = 0.5; center(2) = 0.0;
//        const double cylinder_radius = 0.2;
//       
//        double actpos_enr_val = 0.0;
//        if (inCircleCylinder(cellcenterpos, center, cylinder_radius)) {
//            actpos_enr_val = -1.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (inCircleCylinder(nodalpos, center, cylinder_radius)) {
//            nodepos_enr_val = -1.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        //enrval = actpos_enr_val;
//        
//        break;
//    }
//    case XFEM::Enrichment::typeVoid:
//    {
//        // TODO: generalize
//        double actpos_enr_val = 0.0;
//        if (cellcenterpos(0) > 1.525) {
//            actpos_enr_val = 0.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (nodalpos(0) > 1.525) {
//            nodepos_enr_val = 0.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        //enrval = actpos_enr_val - nodepos_enr_val;
//        enrval = actpos_enr_val;
////        dserror("not yet");
//        
//        break;
//    }
//    case XFEM::Enrichment::typeJump:
//    {
//        // TODO: generalize
//        double actpos_enr_val = 0.0;
//        if (actpos(0) > 1.525) {
//            actpos_enr_val = -1.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (nodalpos(0) > 1.525) {
//            nodepos_enr_val = -1.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        //enrval = actpos_enr_val;
//        
//        break;
//    }
    case XFEM::Enrichment::typeJump:
    {
        const int xfemcondition_label = this->Id();
        
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
        //enrval = actpos_enr_val;
        
        break;
    }
    
    default:
        dserror("unsupported enrichment!");
    }
    return enrval;
}


#endif  // #ifdef CCADISCRET
