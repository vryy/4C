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
        // for the time being the standard jump enrichment function is not available   henke 05/09
        dserror("Use modified enrichment function instead!");
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
      // for the time being the standard kink enrichment function is not available   henke 05/09
      dserror("Use modified enrichment function instead!");
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
    // TODO @ Axel: What does that mean?
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
        /* literature (p. 1006, penultimate line):
         * Belytschko, T., Moës, N., Usui, S. and Parimi, C.
         * Arbitrary discontinuities in finite elements:
         * "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.
         */
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
    case XFEM::Enrichment::typeKink:
    {
        /* literature:
         * Moës, N., Cloirec, M.,Cartraud, P. and Remacle, J. F.
         * A computational approach to handle complex microstructure geometries:
         * "Computer Methods in Applied Mechanics and Engineering", 192:3163--3177, 2003.
         */
        dserror("kink enrichment function not implemented yet");
        // enrval = something involving the standard shape functions
    }
    default:
        dserror("unsupported enrichment (modified)!");
    }
    return enrval;
}



/*----------------------------------------------------------------------*
 | get modified enrichment value                               ag 11/07 |
 | remark: the enrichment function is 0 at nodes and hence the usual    |
 | interpolation property of the standard FEM is satisfied              |
 *----------------------------------------------------------------------*/
double XFEM::Enrichment::ModifiedKinkEnrValue(
        const LINALG::Matrix<3,1>&             actpos,
        //const LINALG::Matrix<3,1>&             nodalpos,
        const int                              nodelid,
        const XFEM::InterfaceHandle&           ih,
        const XFEM::Enrichment::ApproachFrom   approachdirection,
        const DRT::Element&                    ele
        //const LINALG::Matrix<1,8>&             phi
        ) const
{
/* für diese und die folgenden Funktionen gilt, dass ih.PostitionWithinConditionNP == this->XFEMCoditionlabel()
 * eventuell noch für Zweiphasenströmungen und das entsprechende InterfacehandeleCombust angepasst werden muss
 * es ist zu prüfen auf welcher Seite des Interfaces (+) oder (-) man sich befindet
 *  braucht man für ih.PostitionWithinConditionNP globale oder lokale Koordinaten
 * ich glaube global im Moment werden die lokalen Koordinaten des aktuellen Gausspunkts übergeben
 * ist für Kink-Enrichment ok
 * eventuell also noch umrechnen für den Rest

 * approachdirection ist hier nicht notwendig, verwende es doch, da Intergartionszelle weiß, wo sie liegt
 * 
 * henke (rasthofer) 06/09
 */
	//std::cout << "ModifiedEnrValue " << std::endl;
	const DRT::Element::DiscretizationType distype = ele.Shape();
	const int numnode = ele.NumNode();//8
	Epetra_SerialDenseVector  funct(numnode);
    DRT::UTILS::shape_function_3D(funct,actpos(0),actpos(1),actpos(2),distype);
    
    //phi an den Elementknoten wird noch benötigt
    //ih wird bisher nicht gebraucht, allerdings möglicherweise dann für die Knotenwerte von Phi
    //Dummy
    Epetra_SerialDenseVector  phi(numnode);
    for (int inode=0; inode<numnode; inode++)
    {
    	phi(inode) = 0.0;
    }
    
//    double phiactpos =0.0;
//    for (int inode=0; inode<numnode; inode++)
//    {
//    	phiactpos += phi(inode) * funct(inode);
//    }
    
    // TODO @ Axel: What does that mean?
    //dserror("needs update for the approach variable");
    // return value
    if (approachdirection==XFEM::Enrichment::approachUnknown)
    	dserror("approachUnknown not possible in ModifiedEnrValue");
    
    double enrval = 1.0;
    
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
    //was passiert mit der Bestimmung auf welcher Seite man sich befindet
    //if (phiactpos<0)
    //if (phi(nodelid)<0)
    //dann brauche ich approachdirection nicht und nodepos wird durch nodelid ersetzt
//    case XFEM::Enrichment::typeVoid:
//    {
//        double actpos_enr_val = 0.0;
//        if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
//            actpos_enr_val = 0.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (ih.PositionWithinConditionNP(nodalpos) == this->XFEMConditionLabel()) {
//            nodepos_enr_val = 0.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        
//        break;
//    }
//    case XFEM::Enrichment::typeJump:
//    {
//        /* literature (p. 1006, penultimate line):
//         * Belytschko, T., Moës, N., Usui, S. and Parimi, C.
//         * Arbitrary discontinuities in finite elements:
//         * "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.
//         */
//        double actpos_enr_val = 0.0;
//        if (ih.PositionWithinConditionNP(actpos) == this->XFEMConditionLabel()) {
//            actpos_enr_val = -1.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (ih.PositionWithinConditionNP(nodalpos) == this->XFEMConditionLabel()) {
//            nodepos_enr_val = -1.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        
//        break;
//    }
    
    //ist ok, falls das Interface das Element nur einmal schneidet
    //schneidet das Interface das Element mehrmals ist approachdirection in Bezug auf
    //das zugehörige Interfacestück zu wählen
    case XFEM::Enrichment::typeVoid:
    {
        double actpos_enr_val = 0.0;
        if (approachdirection == XFEM::Enrichment::approachFromMinus) {
            actpos_enr_val = 0.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (phi(nodelid)<0) {
            nodepos_enr_val = 0.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        /* literature (p. 1006, penultimate line):
         * Belytschko, T., Moës, N., Usui, S. and Parimi, C.
         * Arbitrary discontinuities in finite elements:
         * "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.
         */
        double actpos_enr_val = 0.0;
        if (approachdirection == XFEM::Enrichment::approachFromMinus) {
            actpos_enr_val = -1.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (phi(nodelid)<0) {
            nodepos_enr_val = -1.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
        /* literature:
         * Moës, N., Cloirec, M.,Cartraud, P. and Remacle, J. F.
         * A computational approach to handle complex microstructure geometries:
         * "Computer Methods in Applied Mechanics and Engineering", 192:3163--3177, 2003.
         */
        //dserror("kink enrichment function not implemented yet");
        // enrval = something involving the standard shape functions
    	
        /*
         * psi = sum(abs(phi(i))*N(i))-abs(sum(phi(i)*N(i)), i=0...numnode-1
         */
        //ICH
//    	const DRT::Element::DiscretizationType distype = ele.Shape();
//    	const int numnode = ele.NumNode();//8
//    	Epetra_SerialDenseVector  funct(numnode);
//        DRT::UTILS::shape_function_3D(funct,actpos(0),actpos(1),actpos(2),distype);
//        
//        //phi an den Elementknoten wird noch benötigt
//        //Dummy
//        Epetra_SerialDenseVector  phi(numnode);
//        for (int inode=0; inode<numnode; inode++)
//        {
//        	phi(inode,1) = 0.0;
//        }
        
        double sum_1 = 0.0;
        double sum_2 = 0.0;
        for (int inode=0; inode<numnode; inode++)
        {
        	sum_1 += fabs(phi(inode)) * funct(inode); //Betrag in C++: fabs() und include math.h
        	sum_2 += phi(inode) * funct(inode);
        }
        enrval = sum_1 - fabs(sum_2);
    }
    default:
        dserror("unsupported enrichment (modified)!");
    }
    return enrval;
}

std::vector<double> XFEM::Enrichment::ModifiedEnrValue_derxy(
        const LINALG::Matrix<3,1>&             actpos,
        const XFEM::InterfaceHandle&           ih,
        const DRT::Element&                    ele
        //const LINALG::Matrix<1,8>&             phi
        ) const
{
	std::vector<double> enrvalderxy(3);
	
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrvalderxy[0] = 0.0;
        enrvalderxy[1] = 0.0;
        enrvalderxy[2] = 0.0;
        break;
    }
    case XFEM::Enrichment::typeVoid:
    {
        enrvalderxy[0] = 0.0;
        enrvalderxy[1] = 0.0;
        enrvalderxy[2] = 0.0;
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
        enrvalderxy[0] = 0.0;
        enrvalderxy[1] = 0.0;
        enrvalderxy[2] = 0.0;
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
    	//Hier muss noch was hin
        enrvalderxy[0] = 0.0;
        enrvalderxy[1] = 0.0;
        enrvalderxy[2] = 0.0;
    	const DRT::Element::DiscretizationType distype = ele.Shape();
    	const int numnode = ele.NumNode();//8
    	Epetra_SerialDenseVector  funct(numnode);
    	Epetra_SerialDenseMatrix  deriv(3,numnode);
        DRT::UTILS::shape_function_3D(funct,actpos(0),actpos(1),actpos(2),distype);
        //local derivatives
        DRT::UTILS::shape_function_3D_deriv1(deriv,actpos(0),actpos(1),actpos(2),distype);
//      //global derivatives -> Jacobi-Matrix
//      //xyze global coordinates of nodes
        Epetra_SerialDenseMatrix  xyze(3,numnode);
          
//           // aus combust3_evaluate
//           Epetra_SerialDenseVector  funct(iel);
//           Epetra_SerialDenseMatrix  xjm(3,3);
//           Epetra_SerialDenseMatrix  deriv(3,iel);
//
//           // get node coordinates of element
//           Epetra_SerialDenseMatrix xyze(3,iel);
//           for(int inode=0;inode<iel;inode++)
//           {
//             xyze(0,inode)=Nodes()[inode]->X()[0];
//             xyze(1,inode)=Nodes()[inode]->X()[1];
//             xyze(2,inode)=Nodes()[inode]->X()[2];
//           } 
        
        for(int inode=0;inode<numnode;inode++)
        {
          xyze(0,inode) = ele.Nodes()[inode]->X()[0];
          xyze(1,inode) = ele.Nodes()[inode]->X()[1];
          xyze(2,inode) = ele.Nodes()[inode]->X()[2];
        }
        //Jacobi-Matix and inverse 
        Epetra_SerialDenseMatrix  xjm(3,3);
        xjm.Multiply('N','T',1,deriv,xyze,0); //-> xjm = deriv * xyze^T
        Epetra_SerialDenseMatrix  xji(3,3);
        xji.ApplyInverse(xjm,xji);
        // compute global derivates
        Epetra_SerialDenseMatrix  derxy(3,numnode);
        derxy.Multiply('N','N',1,xji,deriv,0);
        
        //phi an den Elementknoten wird noch benötigt
        //Dummy
        Epetra_SerialDenseVector  phi(numnode);
        for (int inode=0; inode<numnode; inode++)
        {
        	phi(inode) = 0.0;
        }
        
        for (int isd=0; isd<3; isd++)
        {
        	double sum_1_der = 0.0;
        	double sum_2_der = 0.0;
        	double sum_2 = 0.0;
        	for (int inode=0; inode<numnode; inode++)
        	{
        		sum_1_der += fabs(phi(inode)) * derxy(isd,inode);
        		sum_2_der += phi(inode) * derxy(isd,inode);
        		sum_2 += phi(inode) * funct(inode);
        	}
        	
        	if (sum_2<0)
        	{
        		enrvalderxy[isd] = sum_1_der + sum_2_der;
        	}
        	else
        	{
        		enrvalderxy[isd] = sum_1_der - sum_2_der;
        	}    	
        }
               
        break;
    }
    default:
        dserror("unsupported enrichment (modified)!");
    }
	
	return enrvalderxy;
}

//template <DRT::Element::DiscretizationType DISTYPE>
std::vector<double> XFEM::Enrichment::ModifiedEnrValue_derxy2(
        const LINALG::Matrix<3,1>&             actpos,
        const XFEM::InterfaceHandle&           ih,
        const DRT::Element&                    ele
        //const LINALG::Matrix<1,8>&             phi
        ) const
{
	std::vector<double> enrvalderxy2;
	
    switch (Type()){
    case XFEM::Enrichment::typeStandard:
    {
    	for (int isd=0; isd<6; isd++)
    	{
    		enrvalderxy2[isd] = 0.0;
    	}
        break;
    }
    case XFEM::Enrichment::typeVoid:
    {
    	for (int isd=0; isd<6; isd++)
    	{
    		enrvalderxy2[isd] = 0.0;
    	}
        break;
    }
    case XFEM::Enrichment::typeJump:
    {
    	for (int isd=0; isd<6; isd++)
    	{
    		enrvalderxy2[isd] = 0.0;
    	}
        break;
    }
    case XFEM::Enrichment::typeKink:
    {
        // number of nodes for element
    	const DRT::Element::DiscretizationType distype = ele.Shape();
    	//gibt Fehlermeldung aus, falls es sich um kein hex8 Element handelt
    	//besser wäre template
    	if (distype != DRT::Element::hex8)
    	{
    		dserror("Element not templated in kink enrichment second derivatives yet");
    	}
    	const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    	//const int numnode = ele.NumNode();//8
        static LINALG::Matrix<numnode,1> funct;
        static LINALG::Matrix<3,numnode> deriv;
        static LINALG::Matrix<6,numnode> deriv2;
//    	Epetra_SerialDenseMatrix funct(numnode,1);
//    	Epetra_SerialDenseMatrix deriv(3, numnode);
//    	Epetra_SerialDenseMatrix deriv2(6,numnode);
        DRT::UTILS::shape_function_3D(funct,actpos(0),actpos(1),actpos(2),distype);
        //local derivatives
        DRT::UTILS::shape_function_3D_deriv1(deriv,actpos(0),actpos(1),actpos(2),distype);
        DRT::UTILS::shape_function_3D_deriv2(deriv2,actpos(0),actpos(1),actpos(2),distype);
        static LINALG::Matrix<3,numnode> xyze;
//        Epetra_SerialDenseMatrix  xyze(3,numnode);
        for(int inode=0;inode<numnode;inode++)
        {
          xyze(0,inode) = ele.Nodes()[inode]->X()[0];
          xyze(1,inode) = ele.Nodes()[inode]->X()[1];
          xyze(2,inode) = ele.Nodes()[inode]->X()[2];
        }
//        //Jacobi-Matix and inverse 
//        Epetra_SerialDenseMatrix  xjm(3,3);
//        xjm.Multiply('N','T',1,deriv,xyze,0); //-> xjm = deriv * xyze^T
//        Epetra_SerialDenseMatrix  xji(3,3);
//        xji.ApplyInverse(xjm,xji);
//        // compute global derivates
//        Epetra_SerialDenseMatrix  derxy(3,numnode);
//        derxy.Multiply('N','N',1,xji,deriv,0);
//        Epetra_SerialDenseMatrix  derxy2(6,numnode);
        static LINALG::Matrix<3,3> xjm;
        xjm.MultiplyNT(deriv,xyze);
        static LINALG::Matrix<3,3> xji;
        xji.Invert(xjm);
        static LINALG::Matrix<3,numnode> derxy;
        derxy.Multiply(xji,deriv);
        
        static LINALG::Matrix<6,numnode> derxy2;
        
        //wie mit template umgehen?
        //Funktion benötigt Input als LINALG::Matrix
        DRT::UTILS::gder2<DRT::Element::hex8>(xjm, derxy, deriv2, xyze, derxy2);
        
        
        //phi an den Elementknoten wird noch benötigt
        //Dummy
        //Epetra_SerialDenseVector  phi(numnode);
        static LINALG::Matrix<numnode,1> phi;
        for (int inode=0; inode<numnode; inode++)
        {
        	phi(inode) = 0.0;
        }
        
        for (int isd=0; isd<6; isd++)
        {
        	double sum_1_der2 = 0.0;
        	double sum_2_der2 = 0.0;
        	double sum_2 = 0.0;
        	for (int inode=0; inode<numnode; inode++)
        	{
        		sum_1_der2 += fabs(phi(inode)) * derxy2(isd,inode);
        		sum_2_der2 += phi(inode) * derxy2(isd,inode);
        		sum_2 += phi(inode) * funct(inode);
        	}
        	
        	if (sum_2<0)
        	{
        		enrvalderxy2[isd] = sum_1_der2 + sum_2_der2;
        	}
        	else
        	{
        		enrvalderxy2[isd] = sum_1_der2 - sum_2_der2;
        	}    	
        }
    	
        break;
    }
    default:
        dserror("unsupported enrichment (modified)!");
    }
	
	return enrvalderxy2;
}


#endif  // #ifdef CCADISCRET
