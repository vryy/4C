/*----------------------------------------------------------------------*/
/*!
\file xfluid3_sysmat.H

\brief element formulations for 3d XFEM fluid element

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_SYSMAT_H
#define XFLUID3_SYSMAT_H

#include <math.h>
#include "xfluid3.H"
#include "xfluid3_utils.H"
#include <blitz/array.h>
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfem.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_xfem/coordinate_transformation.H"
#include "../drt_lib/drt_utils.H"
#include "fluid3_stabilization.H"
#include "xfluid3_local_assembler.H"
#include "xfluid3_interpolation.H"

#include <Epetra_SerialDenseSolver.h>

using namespace XFEM::PHYSICS;

namespace XFLUID
{
  
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void Sysmat(
        const DRT::ELEMENTS::XFluid3*     ele,           ///< the element those matrix is calculated
        const RCP<XFEM::InterfaceHandle>  ih,            ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const std::vector<double>&        locval,        ///< nodal unknowns at n+1, i
        const std::vector<double>&        locval_hist,   ///< nodal unknowns at n
        BlitzMat&                         estif,         ///< element matrix to calculate
        BlitzVec&                         eforce,        ///< element rhs to calculate
        const struct _MATERIAL*           material,      ///< fluid material
        const double                      time,          ///< current time (pseudotime for stationary formulation)
        const double                      timefac,       ///< One-step-Theta: theta*dt, BDF2: 2/3 * dt, stationary: 1.0
        const bool                        newton,        ///< full Newton or fixed-point-like
        const bool                        pstab,         ///< flag for stabilisation
        const bool                        supg,          ///< flag for stabilisation
        const bool                        cstab,         ///< flag for stabilisation
        const bool                        instationary   ///< switch between stationary and instationary formulation
        )
{
    // initialize arrays
    estif = 0.0;
    eforce = 0.0;
    
    // number of nodes for element
    const int numnode = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // get node coordinates of the current element
    static BlitzMat xyze(nsd, numnode);
    DRT::UTILS::fillPositionArray<DISTYPE>(ele, xyze);

    // dead load in element nodes
    //////////////////////////////////////////////////// , BlitzMat edeadng_(BodyForce(ele->Nodes(),time));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double visc = material->m.fluid->viscosity;

    // flag for higher order elements
    const bool higher_order_ele = ele->isHigherOrderElement<DISTYPE>();
    //const bool higher_order_ele = secondDerivativesAvailable<DISTYPE>();
    
    const DRT::UTILS::GaussRule3D gaussrule = getXFEMGaussrule<DISTYPE,ASSTYPE>(ih->ElementIntersected(ele->Id()));
    
    const LocalAssembler<DISTYPE, ASSTYPE> assembler(dofman, estif, eforce);
    
    const blitz::Range _  = blitz::Range::all();
    
    // We define the variables i,j,k to be indices to blitz arrays.
    // These are used for array expressions, that is matrix-vector
    // products in the following.
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index

    const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::stressdistype;
    const DRT::Element::DiscretizationType discpresdistype = XFLUID::getDiscPressureInterpolationType3D<DISTYPE>();
    
    // figure out whether we have stress unknowns at all
    const bool stress_unknowns_present = (getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
    const bool discpres_unknowns_present = (getNumParam<ASSTYPE>(dofman, DiscPres, 0) > 0);
//    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
//    const bool pressure_unknowns_present = (getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
//    cout << endl;
//    if (ASSTYPE == XFEM::standard_assembly)
//        cout << "standard assembly" << endl;
//    else
//        cout << "xfem assembly" << endl;
//            
//    cout << "stress unknowns present  : " << stress_unknowns_present << endl;
//    cout << "velocity unknowns present: " << velocity_unknowns_present << endl;
//    cout << "pressure unknowns present: " << pressure_unknowns_present << endl;
    
    
    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
    const int numparamvely = getNumParam<ASSTYPE>(dofman, Vely, numnode);
    const int numparamvelz = getNumParam<ASSTYPE>(dofman, Velz, numnode);
    const int numparampres = getNumParam<ASSTYPE>(dofman, Pres, numnode);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const int numparamsigmaxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
    const int numparamdiscpres = getNumParam<ASSTYPE>(dofman, DiscPres, 1);

    
    // split velocity and pressure (and stress)
    BlitzVec eprenp(numparampres);
    BlitzVec ediscpreele(numparamdiscpres);
    BlitzMat evelnp(3,numparamvelx,blitz::ColumnMajorArray<2>());
    BlitzMat evelnp_hist(3,numparamvelx,blitz::ColumnMajorArray<2>());
    BlitzMat esigma(6,numparamsigmaxx,blitz::ColumnMajorArray<2>());
    
    const vector<int>& velxdof(dofman.LocalDofPosPerField(Velx));
    const vector<int>& velydof(dofman.LocalDofPosPerField(Vely));
    const vector<int>& velzdof(dofman.LocalDofPosPerField(Velz));
    const vector<int>& presdof(dofman.LocalDofPosPerField(Pres));
    
    for (int iparam=0; iparam<numparamvelx; ++iparam)
    {
        evelnp(     0,iparam) = locval[     velxdof[iparam]];
        evelnp_hist(0,iparam) = locval_hist[velxdof[iparam]];
    }
    for (int iparam=0; iparam<numparamvely; ++iparam)
    {
        evelnp(     1,iparam) = locval[     velydof[iparam]];
        evelnp_hist(1,iparam) = locval_hist[velydof[iparam]];
    }
    for (int iparam=0; iparam<numparamvelz; ++iparam)
    {
        evelnp(     2,iparam) = locval[     velzdof[iparam]];
        evelnp_hist(2,iparam) = locval_hist[velzdof[iparam]];
    }
    for (int iparam=0; iparam<numparampres; ++iparam)
        eprenp(iparam) = locval[presdof[iparam]];
    if (stress_unknowns_present)
    {
        const int numparamsigmayy = getNumParam<ASSTYPE>(dofman, Sigmayy, 1);
        const int numparamsigmazz = getNumParam<ASSTYPE>(dofman, Sigmazz, 1);
        const int numparamsigmaxy = getNumParam<ASSTYPE>(dofman, Sigmaxy, 1);
        const int numparamsigmaxz = getNumParam<ASSTYPE>(dofman, Sigmaxz, 1);
        const int numparamsigmayz = getNumParam<ASSTYPE>(dofman, Sigmayz, 1);
        const vector<int>& sigmaxxdof(dofman.LocalDofPosPerField(Sigmaxx));
        const vector<int>& sigmayydof(dofman.LocalDofPosPerField(Sigmayy));
        const vector<int>& sigmazzdof(dofman.LocalDofPosPerField(Sigmazz));
        const vector<int>& sigmaxydof(dofman.LocalDofPosPerField(Sigmaxy));
        const vector<int>& sigmaxzdof(dofman.LocalDofPosPerField(Sigmaxz));
        const vector<int>& sigmayzdof(dofman.LocalDofPosPerField(Sigmayz));
        for (int iparam=0; iparam<numparamsigmaxx; ++iparam)   esigma(0,iparam) = locval[sigmaxxdof[iparam]];
        for (int iparam=0; iparam<numparamsigmayy; ++iparam)   esigma(1,iparam) = locval[sigmayydof[iparam]];
        for (int iparam=0; iparam<numparamsigmazz; ++iparam)   esigma(2,iparam) = locval[sigmazzdof[iparam]];
        for (int iparam=0; iparam<numparamsigmaxy; ++iparam)   esigma(3,iparam) = locval[sigmaxydof[iparam]];
        for (int iparam=0; iparam<numparamsigmaxz; ++iparam)   esigma(4,iparam) = locval[sigmaxzdof[iparam]];
        for (int iparam=0; iparam<numparamsigmayz; ++iparam)   esigma(5,iparam) = locval[sigmayzdof[iparam]];
    }
    if (discpres_unknowns_present)
    {
        const vector<int>& discpresdof(dofman.LocalDofPosPerField(DiscPres));
        for (int iparam=0; iparam<numparamdiscpres; ++iparam)   ediscpreele(iparam) = locval[discpresdof[iparam]];
    }
    
    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,visc,xyze);
    const double mk = XFLUID::MK<DISTYPE>();
    
    // temporary arrays
    BlitzVec enr_funct(numparamvelx);
    BlitzMat enr_derxy(3,numparamvelx,blitz::ColumnMajorArray<2>());
    BlitzMat enr_derxy2(6,numparamvelx,blitz::ColumnMajorArray<2>());

    BlitzVec enr_funct_stress(numparamsigmaxx);
    //BlitzMat enr_derxy_stress(3,numparamsigmaxx,blitz::ColumnMajorArray<2>());
    
    BlitzVec enr_funct_discpres(numparamdiscpres);
    //BlitzMat enr_derxy_discpres(3,numparamdiscpres,blitz::ColumnMajorArray<2>());
    
    // information about domain integration cells
    const XFEM::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele->Id(),DISTYPE));
    // loop over integration cells
    for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            static XFEM::PosEtaDomain pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static XFEM::PosXiDomain posXiDomain;
            XFEM::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = XFEM::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static BlitzVec funct(numnode);
            static BlitzMat deriv(3, numnode, blitz::ColumnMajorArray<2>());
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
            DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            // discontinouos stress shape functions
            static BlitzVec funct_stress(DRT::UTILS::getNumberOfElementNodes(stressdistype));
            static BlitzMat deriv_stress(3, DRT::UTILS::getNumberOfElementNodes(stressdistype),blitz::ColumnMajorArray<2>());
            if (stress_unknowns_present)
            {
                DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
                DRT::UTILS::shape_function_3D_deriv1(deriv_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
            }
            // discontinouos pressure shape functions
            static BlitzVec funct_discpres(DRT::UTILS::getNumberOfElementNodes(discpresdistype));
            static BlitzMat deriv_discpres(3, DRT::UTILS::getNumberOfElementNodes(discpresdistype),blitz::ColumnMajorArray<2>());
            if (discpres_unknowns_present)
            {
                DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);
                DRT::UTILS::shape_function_3D_deriv1(deriv_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);
            }

            // position of the gausspoint in physical coordinates
            static BlitzVec gauss_pos_xyz(3);
            gauss_pos_xyz = blitz::sum(funct(j)*xyze(i,j),j);
      
            // get transposed of the jacobian matrix d x / d \xi
            static BlitzMat xjm(3,3);
            xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);
            const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                               xjm(0,1)*xjm(1,2)*xjm(2,0)+
                               xjm(0,2)*xjm(1,0)*xjm(2,1)-
                               xjm(0,2)*xjm(1,1)*xjm(2,0)-
                               xjm(0,0)*xjm(1,2)*xjm(2,1)-
                               xjm(0,1)*xjm(1,0)*xjm(2,2);
            const double fac = intpoints.qwgt[iquad]*det*detcell;

            if (det < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
            }

            // inverse of jacobian
            static BlitzMat xji(3,3);
            xji = Inverse3x3(xjm, det);

            // compute global derivates
            static BlitzMat derxy(3, DRT::UTILS::getNumberOfElementNodes(DISTYPE));
            static BlitzMat derxy_stress(3, DRT::UTILS::getNumberOfElementNodes(stressdistype),blitz::ColumnMajorArray<2>());
            static BlitzMat derxy_discpres(3, DRT::UTILS::getNumberOfElementNodes(discpresdistype),blitz::ColumnMajorArray<2>());
            derxy          = blitz::sum(xji(i,k)*deriv(k,j),k);
            derxy_stress   = blitz::sum(xji(i,k)*deriv_stress(k,j),k);
            derxy_discpres = blitz::sum(xji(i,k)*deriv_discpres(k,j),k);

            // compute second global derivative
            static BlitzMat derxy2(6,numnode,blitz::ColumnMajorArray<2>());
            if (higher_order_ele)
            {
                static BlitzMat deriv2(6,numnode,blitz::ColumnMajorArray<2>());
                DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
                derxy2 = gder2<DISTYPE>(xjm, derxy, deriv2, xyze);
            }
            else
            {
                derxy2 = 0.;
            }

            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // shape function for nodal dofs
                XFEM::ComputeEnrichedNodalShapefunction(
                        *ele,
                        ih,
                        dofman,
                        Velx,
                        gauss_pos_xyz,
                        XFEM::Enrichment::approachUnknown,
                        funct,
                        derxy,
                        derxy2, 
                        enr_funct,
                        enr_derxy,
                        enr_derxy2);
          
                if (stress_unknowns_present)
                {
                    // shape functions for element dofs
                    XFEM::ComputeEnrichedElementShapefunction(
                            *ele,
                            ih,
                            dofman,
                            Sigmaxx,
                            gauss_pos_xyz,
                            XFEM::Enrichment::approachUnknown,
                            funct_stress,
                            enr_funct_stress);
                }
                
                if (discpres_unknowns_present)
                {
                    // shape functions for element dofs
                    XFEM::ComputeEnrichedElementShapefunction(
                            *ele,
                            ih,
                            dofman,
                            Sigmaxx,
                            gauss_pos_xyz,
                            XFEM::Enrichment::approachUnknown,
                            funct_discpres,
                            enr_funct_discpres);
                }
            }
            else
            {
                enr_funct.reference(funct);
                enr_derxy.reference(derxy);
                enr_derxy2.reference(derxy2);
          
                enr_funct_stress.reference(funct_stress);
                //enr_derxy_stress.reference(derxy_stress);
                
                enr_funct_discpres.reference(funct_discpres);
                //enr_derxy_discpres.reference(derxy_discpres);
            }
      
            // create views on the shape function arrays for easy handling in the assembly process
            const BlitzVec shp(enr_funct);
            const BlitzVec shp_dx(enr_derxy(0,_));
            const BlitzVec shp_dy(enr_derxy(1,_));
            const BlitzVec shp_dz(enr_derxy(2,_));
            const BlitzVec shp_dxdx(enr_derxy2(0,_)); const BlitzVec shp_dxdy(enr_derxy2(3,_)); const BlitzVec shp_dxdz(enr_derxy2(4,_));
            const BlitzVec shp_dydx(shp_dxdy);        const BlitzVec shp_dydy(enr_derxy2(1,_)); const BlitzVec shp_dydz(enr_derxy2(5,_));
            const BlitzVec shp_dzdx(shp_dxdz);        const BlitzVec shp_dzdy(shp_dydz);        const BlitzVec shp_dzdz(enr_derxy2(2,_));
            
            const BlitzVec shp_sigma(enr_funct_stress);
            //const BlitzVec shp_sigma_dx(enr_derxy_stress(0,_));
            //const BlitzVec shp_sigma_dy(enr_derxy_stress(1,_));
            //const BlitzVec shp_sigma_dz(enr_derxy_stress(2,_));
            
            const BlitzVec shp_discpres(enr_funct_discpres);
            //const BlitzVec shp_discpres_dx(enr_derxy_discpres(0,_));
            //const BlitzVec shp_discpres_dy(enr_derxy_discpres(1,_));
            //const BlitzVec shp_discpres_dz(enr_derxy_discpres(2,_));
      
            // get velocities (n+g,i) at integration point
            static BlitzVec velint(3);
            velint = blitz::sum(enr_funct(j)*evelnp(i,j),j);

            // get history data (n) at integration point
            static BlitzVec histvec(3);
            histvec = blitz::sum(enr_funct(j)*evelnp_hist(i,j),j);
            
            // get velocity (np,i) derivatives at integration point
            static BlitzMat vderxy(3,3);
            vderxy = blitz::sum(enr_derxy(j,k)*evelnp(i,k),k);
      
            // calculate 2nd velocity derivatives at integration point
            static BlitzMat vderxy2(3,6);
            if (higher_order_ele)
            {
                vderxy2 = blitz::sum(enr_derxy2(j,k)*evelnp(i,k),k);
            }
            else
            {
                vderxy2 = 0.;
            }

            // get pressure gradients
            static BlitzVec gradp(3);
            gradp = blitz::sum(enr_derxy(i,j)*eprenp(j),j);
    
            // get pressure
            const double press(blitz::sum(enr_funct*eprenp));
    
            // get viscous stress unknowns
            static BlitzMat sigma(3,3);
            if (stress_unknowns_present)
            {
                sigma(0,0) = blitz::sum(shp_sigma(_)*esigma(0,_)); sigma(0,1) = blitz::sum(shp_sigma(_)*esigma(3,_)); sigma(0,2) = blitz::sum(shp_sigma(_)*esigma(4,_));
                sigma(1,0) = sigma(0,1);                           sigma(1,1) = blitz::sum(shp_sigma(_)*esigma(1,_)); sigma(1,2) = blitz::sum(shp_sigma(_)*esigma(5,_));
                sigma(2,0) = sigma(0,2);                           sigma(2,1) = sigma(1,2);                           sigma(2,2) = blitz::sum(shp_sigma(_)*esigma(2,_));
            }
          
//            BlitzVec nabla_dot_tau(3);
//            if (stress_unknowns_present)
//            {
//                nabla_dot_tau(0) = blitz::sum(shp_tau_dx(_)*etau(0,_)) + blitz::sum(shp_tau_dy(_)*etau(3,_)) + blitz::sum(shp_tau_dz(_)*etau(4,_));
//                nabla_dot_tau(1) = blitz::sum(shp_tau_dx(_)*etau(3,_)) + blitz::sum(shp_tau_dy(_)*etau(1,_)) + blitz::sum(shp_tau_dz(_)*etau(5,_));
//                nabla_dot_tau(2) = blitz::sum(shp_tau_dx(_)*etau(4,_)) + blitz::sum(shp_tau_dy(_)*etau(5,_)) + blitz::sum(shp_tau_dz(_)*etau(2,_));
//            }
            
            
      
      
//            const BlitzMat eps(0.5*(vderxy_(i,j) + vderxy_(j,i));
//            const BlitzMat epstau(tau*reciproke_visc);
//      
//            if (ele->Id() == 0)
//            {
//                cout << endl;
//                cout << "eps^tau: " << epstau << endl;
//                cout << "eps^u:   " << eps << endl;
//                cout << "pressure = " << press << endl;
//            }

            
            // get bodyforce in gausspoint
            static BlitzVec bodyforce(3);
            bodyforce = 0.0;
            //////////////////////////////////////////BlitzVec bodyforce_(blitz::sum(enr_edeadng_(i,j)*enr_funct_(j),j));

            

            // get velocity norm
            const double vel_norm = sqrt(blitz::sum(velint*velint));

            // normed velocity at element centre
            static BlitzVec velino(3);
            if (vel_norm>=1e-6)
            {
                velino = velint/vel_norm;
            }
            else
            {
                velino = 0.;
                velino(0) = 1;
            }

            // get streamlength
            const double val = blitz::sum(blitz::abs(blitz::sum(velino(j)*derxy(j,i),j)));
            const double strle = 2.0/val;

            double tau_stab_M;
            double tau_stab_Mp;
            double tau_stab_C;
            if (instationary)
            {
                // calculate tau: stabilization parameters for stationary case
                
                const double visceff = visc;
                /* viscous : reactive forces */
                const double re1 = 4.0 * timefac * visceff / (mk * DSQR(strle));

                /* convective : viscous forces */
                const double re2 = mk * vel_norm * strle / (2.0 * visceff);

                const double xi1 = max(re1,1.0);
                const double xi2 = max(re2,1.0);

                tau_stab_M = DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

                // compute tau_Mp
                //    stability parameter definition according to Franca and Valentin (2000)
                //                                       and Barrenechea and Valentin (2002)

                 /* viscous : reactive forces */
                const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk));
                /* convective : viscous forces */
                const double re_convect = mk * vel_norm * hk / (2.0 * visceff);

                const double xi_viscous = max(re_viscous,1.0);
                const double xi_convect = max(re_convect,1.0);

                /*
                                xi1,xi2 ^
                                        |      /
                                        |     /
                                        |    /
                                      1 +---+
                                        |
                                        |
                                        |
                                        +--------------> re1,re2
                                            1
                */
                tau_stab_Mp = DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visceff/mk) * xi_convect);

                /*------------------------------------------------------ compute tau_C ---*/
                /*-- stability parameter definition according to Codina (2002), CMAME 191
                 *
                 * Analysis of a stabilized finite element approximation of the transient
                 * convection-diffusion-reaction equation using orthogonal subscales.
                 * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
                 *
                 * */
                //tau[2] = sqrt(DSQR(visc)+DSQR(0.5*vel_norm*hk));

                // Wall Diss. 99
                /*
                                    xi2 ^
                                        |
                                      1 |   +-----------
                                        |  /
                                        | /
                                        |/
                                        +--------------> Re2
                                            1
                */
                const double xi_tau_c = min(re2,1.0);
                tau_stab_C = vel_norm * hk * 0.5 * xi_tau_c /timefac;
            }
            else
            {
                // calculate tau: stabilization parameters for stationary case
                
                // compute tau_Mu
                const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
                const double xi_tau_mu = max(re_tau_mu, 1.0);
                tau_stab_M = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);
    
                // compute tau_Mp
                const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
                const double xi_tau_mp = max(re_tau_mp,1.0);
                tau_stab_Mp = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);
    
                // compute tau_C
                const double xi_tau_c = min(re_tau_mp, 1.0);
                tau_stab_C = 0.5*vel_norm*hk*xi_tau_c;
            }
            
            // stabilisation parameter
            const double tau_M  = fac*tau_stab_M;
            const double tau_Mp = fac*tau_stab_Mp;
            const double tau_C  = fac*tau_stab_C;
            
            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            static BlitzVec rhsint(3);
            rhsint = histvec + bodyforce*timefac;

            /*----------------- get numerical representation of single operators ---*/

            /* Convective term  u_old * grad u_old: */
            static BlitzVec conv_old(3);
            conv_old = blitz::sum(vderxy(i, j)*velint(j), j);
            
            /* Viscous term  div epsilon(u_old) */
            static BlitzVec visc_old(3);
            visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
            visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
            visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));
            
            // evaluate residual once for all stabilisation right hand sides
            static BlitzVec res_old(3);
            res_old = -rhsint+timefac*(conv_old+gradp-2.0*visc*visc_old);
            if (instationary)
                res_old += velint;
      
            /* Reactive term  u:  funct */
            /* linearise convective term */

            /*--- convective part u_old * grad (funct) --------------------------*/
            /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
             with  N .. form function matrix                                   */
            const BlitzVec enr_conv_c_(blitz::sum(enr_derxy(j,i)*velint(j), j));

//              /*--- convective grid part u_G * grad (funct) -----------------------*/
//              /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
//              enr_conv_g_ = 0.0;


          /*--- viscous term  - grad * epsilon(u): ----------------------------*/
          /*   /                                                \
               |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
             1 |                                                |
             - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
             2 |                                                |
               |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
               \                                                /
    
               with N_x .. x-line of N
               N_y .. y-line of N                                             */
            blitz::Array<double,3> enr_viscs2_(3,3,numparamvelx);
            enr_viscs2_(0,0,_) = 0.5 * (2.0 * shp_dxdx + shp_dydy + shp_dzdz);
            enr_viscs2_(0,1,_) = 0.5 *  shp_dxdy;
            enr_viscs2_(0,2,_) = 0.5 *  shp_dxdz;
            enr_viscs2_(1,0,_) = 0.5 *  shp_dydx;
            enr_viscs2_(1,1,_) = 0.5 * (shp_dxdx + 2.0 * shp_dydy + shp_dzdz);
            enr_viscs2_(1,2,_) = 0.5 *  shp_dydz;
            enr_viscs2_(2,0,_) = 0.5 *  shp_dzdx;
            enr_viscs2_(2,1,_) = 0.5 *  shp_dzdy;
            enr_viscs2_(2,2,_) = 0.5 * (shp_dxdx + shp_dydy + 2.0 * shp_dzdz);


            /*--------------------------------- now build single stiffness terms ---*/

            //----------------------------------------------------------------------
            //                            GALERKIN PART

            if (instationary)
            {
                /* inertia (contribution to mass matrix) */
                /*
                                     /        \
                                    |          |
                                    |  Du , v  |
                                    |          |
                                     \        /
                */
                assembler.template Matrix<Velx,Velx>(shp, fac, shp);
                assembler.template Matrix<Vely,Vely>(shp, fac, shp);
                assembler.template Matrix<Velz,Velz>(shp, fac, shp);
                
                assembler.template Vector<Velx>(shp, -fac*velint(0));
                assembler.template Vector<Vely>(shp, -fac*velint(1));
                assembler.template Vector<Velz>(shp, -fac*velint(2));
            }
            
            /* convection, convective part */
            /*
                         /                       \
                        |      / n+1       \      |
                        | v , | u   o nabla | Du  |
                        |      \ (i)       /      |
                         \                       /
            */
            assembler.template Matrix<Velx,Velx>(shp, timefacfac, enr_conv_c_);
            assembler.template Matrix<Vely,Vely>(shp, timefacfac, enr_conv_c_);
            assembler.template Matrix<Velz,Velz>(shp, timefacfac, enr_conv_c_);
            
            assembler.template Vector<Velx>(shp, -timefacfac*(velint(0)*vderxy(0,0) // check order
                                                             +velint(1)*vderxy(0,1)
                                                             +velint(2)*vderxy(0,2)));
            assembler.template Vector<Vely>(shp, -timefacfac*(velint(0)*vderxy(1,0)
                                                             +velint(1)*vderxy(1,1)
                                                             +velint(2)*vderxy(1,2)));
            assembler.template Vector<Velz>(shp, -timefacfac*(velint(0)*vderxy(2,0)
                                                             +velint(1)*vderxy(2,1)
                                                             +velint(2)*vderxy(2,2)));
            
            if (newton)
            {
                /*  convection, reactive part */
                /*
                       /                         \
                      |      /          \   n+1   |
                      | v , | Du o nabla | u      |
                      |      \          /   (i)   |
                       \                         /
                */
                assembler.template Matrix<Velx,Velx>(shp, timefacfac*vderxy(0,0), shp);
                assembler.template Matrix<Velx,Vely>(shp, timefacfac*vderxy(0,1), shp);
                assembler.template Matrix<Velx,Velz>(shp, timefacfac*vderxy(0,2), shp);
                assembler.template Matrix<Vely,Velx>(shp, timefacfac*vderxy(1,0), shp);
                assembler.template Matrix<Vely,Vely>(shp, timefacfac*vderxy(1,1), shp);
                assembler.template Matrix<Vely,Velz>(shp, timefacfac*vderxy(1,2), shp);
                assembler.template Matrix<Velz,Velx>(shp, timefacfac*vderxy(2,0), shp);
                assembler.template Matrix<Velz,Vely>(shp, timefacfac*vderxy(2,1), shp);
                assembler.template Matrix<Velz,Velz>(shp, timefacfac*vderxy(2,2), shp);
            }
            
            /* Viskositaetsterm */
            /*
                          /                        \
                         |       / \         /  \   |
                         |  eps | v | , tau | Du |  |
                         |       \ /         \  /   |
                          \                        /
            */
            assembler.template Matrix<Velx,Velx>(shp_dx, 2.0*visc*timefacfac, shp_dx);
            assembler.template Matrix<Velx,Velx>(shp_dy,     visc*timefacfac, shp_dy);
            assembler.template Matrix<Velx,Velx>(shp_dz,     visc*timefacfac, shp_dz);
            assembler.template Matrix<Velx,Vely>(shp_dy,     visc*timefacfac, shp_dx);
            assembler.template Matrix<Velx,Velz>(shp_dz,     visc*timefacfac, shp_dx);
            
            assembler.template Matrix<Vely,Velx>(shp_dx,     visc*timefacfac, shp_dy);
            assembler.template Matrix<Vely,Vely>(shp_dx,     visc*timefacfac, shp_dx);
            assembler.template Matrix<Vely,Vely>(shp_dy, 2.0*visc*timefacfac, shp_dy);
            assembler.template Matrix<Vely,Vely>(shp_dz,     visc*timefacfac, shp_dz);
            assembler.template Matrix<Vely,Velz>(shp_dz,     visc*timefacfac, shp_dy);
            
            assembler.template Matrix<Velz,Velx>(shp_dx,     visc*timefacfac, shp_dz);
            assembler.template Matrix<Velz,Vely>(shp_dy,     visc*timefacfac, shp_dz);
            assembler.template Matrix<Velz,Velz>(shp_dx,     visc*timefacfac, shp_dx);
            assembler.template Matrix<Velz,Velz>(shp_dy,     visc*timefacfac, shp_dy);
            assembler.template Matrix<Velz,Velz>(shp_dz, 2.0*visc*timefacfac, shp_dz);
            
            assembler.template Vector<Velx>(shp_dx,     -visc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
            assembler.template Vector<Velx>(shp_dy,     -visc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
            assembler.template Vector<Velx>(shp_dz,     -visc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
            
            assembler.template Vector<Vely>(shp_dx,     -visc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
            assembler.template Vector<Vely>(shp_dy,     -visc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
            assembler.template Vector<Vely>(shp_dz,     -visc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));
            
            assembler.template Vector<Velz>(shp_dx,     -visc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
            assembler.template Vector<Velz>(shp_dy,     -visc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
            assembler.template Vector<Velz>(shp_dz,     -visc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));
            
            /* Druckterm */
            /*
                            /                \
                           |                  |
                         - |  nabla o v , Dp  |
                           |                  |
                            \                /
            */
            assembler.template Matrix<Velx,Pres>(shp_dx, -timefacfac, shp);
            assembler.template Matrix<Vely,Pres>(shp_dy, -timefacfac, shp);
            assembler.template Matrix<Velz,Pres>(shp_dz, -timefacfac, shp);
            
            assembler.template Vector<Velx>(shp_dx, timefacfac*press);
            assembler.template Vector<Vely>(shp_dy, timefacfac*press);
            assembler.template Vector<Velz>(shp_dz, timefacfac*press);
            
            /* Divergenzfreiheit - continuity equation*/
            /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
            */
            assembler.template Matrix<Pres,Velx>(shp, timefacfac, shp_dx);
            assembler.template Matrix<Pres,Vely>(shp, timefacfac, shp_dy);
            assembler.template Matrix<Pres,Velz>(shp, timefacfac, shp_dz);
            
            const double trace_gamma = (vderxy(0, 0) + vderxy(1, 1) + vderxy(2, 2));
            assembler.template Vector<Pres>(shp, -timefacfac*trace_gamma);
            
            
            
            // source term of the right hand side
            assembler.template Vector<Velx>(shp, fac*rhsint(0));
            assembler.template Vector<Vely>(shp, fac*rhsint(1));
            assembler.template Vector<Velz>(shp, fac*rhsint(2));
            
            
            // Hellinger-Reissner terms
            if (stress_unknowns_present)
            {
                const double reciproke_visc = 1.0 / visc;
                /* stress-stress coupling, quadratic part */
                /*
                               /                    \
                              |                      |
                            - | virt tau , eps(Dtau) |
                              |                      |
                               \                    /
                */
                const double factor = 1.0;//reciproke_visc;
                assembler.template Matrix<Sigmaxx,Sigmaxx>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmaxy,Sigmaxy>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmaxz,Sigmaxz>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmayx,Sigmayx>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmayy,Sigmayy>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmayz,Sigmayz>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmazx,Sigmazx>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmazy,Sigmazy>(shp_sigma, -factor*timefacfac, shp_sigma);
                assembler.template Matrix<Sigmazz,Sigmazz>(shp_sigma, -factor*timefacfac, shp_sigma);
                
                assembler.template Vector<Sigmaxx>(shp_sigma, factor*timefacfac*sigma(0,0));
                assembler.template Vector<Sigmaxy>(shp_sigma, factor*timefacfac*sigma(0,1));
                assembler.template Vector<Sigmaxz>(shp_sigma, factor*timefacfac*sigma(0,2));
                assembler.template Vector<Sigmayx>(shp_sigma, factor*timefacfac*sigma(1,0));
                assembler.template Vector<Sigmayy>(shp_sigma, factor*timefacfac*sigma(1,1));
                assembler.template Vector<Sigmayz>(shp_sigma, factor*timefacfac*sigma(1,2));
                assembler.template Vector<Sigmazx>(shp_sigma, factor*timefacfac*sigma(2,0));
                assembler.template Vector<Sigmazy>(shp_sigma, factor*timefacfac*sigma(2,1));
                assembler.template Vector<Sigmazz>(shp_sigma, factor*timefacfac*sigma(2,2));
                
                /* stress-stress coupling, rectangular part */
                /*
                               /                  \
                              |                    |
                              | virt tau , eps(Du) |
                              |                    |
                               \                  /
                */
                const double rect_factor = visc;
                assembler.template Matrix<Sigmaxx,Velx>(shp_sigma, 2.0*timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmaxy,Velx>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dy);
                assembler.template Matrix<Sigmaxy,Vely>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmaxz,Velx>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dz);
                assembler.template Matrix<Sigmaxz,Velz>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmayx,Velx>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dy);
                assembler.template Matrix<Sigmayx,Vely>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmayy,Vely>(shp_sigma, 2.0*timefacfac*0.5*rect_factor, shp_dy);
                assembler.template Matrix<Sigmayz,Velz>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dz);
                assembler.template Matrix<Sigmayz,Vely>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmazx,Velx>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dz);
                assembler.template Matrix<Sigmazx,Velz>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dx);
                assembler.template Matrix<Sigmazy,Vely>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dz);
                assembler.template Matrix<Sigmazy,Velz>(shp_sigma,     timefacfac*0.5*rect_factor, shp_dy);
                assembler.template Matrix<Sigmazz,Velz>(shp_sigma, 2.0*timefacfac*0.5*rect_factor, shp_dz);
                
                assembler.template Vector<Sigmaxx>(shp_sigma, -2.0*timefacfac*0.5*rect_factor*vderxy(0, 0));
                assembler.template Vector<Sigmaxy>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(0, 1));
                assembler.template Vector<Sigmaxy>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(1, 0));
                assembler.template Vector<Sigmaxz>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(0, 2));
                assembler.template Vector<Sigmaxz>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(2, 0));
                assembler.template Vector<Sigmayx>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(1, 0));
                assembler.template Vector<Sigmayx>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(0, 1));
                assembler.template Vector<Sigmayy>(shp_sigma, -2.0*timefacfac*0.5*rect_factor*vderxy(1, 1));
                assembler.template Vector<Sigmayz>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(1, 2));
                assembler.template Vector<Sigmayz>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(2, 1));
                assembler.template Vector<Sigmazx>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(2, 0));
                assembler.template Vector<Sigmazx>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(0, 2));
                assembler.template Vector<Sigmazy>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(2, 1));
                assembler.template Vector<Sigmazy>(shp_sigma, -    timefacfac*0.5*rect_factor*vderxy(1, 2));
                assembler.template Vector<Sigmazz>(shp_sigma, -2.0*timefacfac*0.5*rect_factor*vderxy(2, 2));
                
                /*
                               /                \
                              |                  |
                            - | virt sigma , p I |
                              |                  |
                               \                /
                */
                const double taue_p_fac = 1.0;//reciproke_visc;
                assembler.template Matrix<Sigmaxx,Pres>(shp_sigma, -taue_p_fac*timefacfac, shp);
                assembler.template Matrix<Sigmayy,Pres>(shp_sigma, -taue_p_fac*timefacfac, shp);
                assembler.template Matrix<Sigmazz,Pres>(shp_sigma, -taue_p_fac*timefacfac, shp);
                
                assembler.template Vector<Sigmaxx>(shp_sigma, taue_p_fac*timefacfac*press);
                assembler.template Vector<Sigmayy>(shp_sigma, taue_p_fac*timefacfac*press);
                assembler.template Vector<Sigmazz>(shp_sigma, taue_p_fac*timefacfac*press);
            }
            
            //----------------------------------------------------------------------
            //                 PRESSURE STABILISATION PART
            if(pstab)
            {
                const double timetauMp  = timefac * tau_Mp;
                if (instationary)
                {
                    /* pressure stabilisation: inertia */
                    /*
                                /              \
                               |                |
                               |  Du , nabla q  |
                               |                |
                                \              /
                    */
                    assembler.template Matrix<Pres,Velx>(shp_dx, timetauMp, shp);
                    assembler.template Matrix<Pres,Vely>(shp_dy, timetauMp, shp);
                    assembler.template Matrix<Pres,Velz>(shp_dz, timetauMp, shp);
                }
                const double ttimetauMp = timefac * timefac * tau_Mp;
                /* pressure stabilisation: convection, convective part */
                /*
                          /                             \
                         |             / n+1       \     |
                         | nabla q ,  | u   o nabla | Du |
                         |             \ i         /     |
                          \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp_dx, ttimetauMp, enr_conv_c_);
                assembler.template Matrix<Pres,Vely>(shp_dy, ttimetauMp, enr_conv_c_);
                assembler.template Matrix<Pres,Velz>(shp_dz, ttimetauMp, enr_conv_c_);
                
                if (newton)
                {
                    /*  pressure stabilisation: convection, reactive part
                          /                             \
                         |           /          \   n+1  |
                         | grad q , | Du o nabla | u     |
                         |           \          /   (i)  |
                          \                             /
                    */
                    assembler.template Matrix<Pres,Velx>(shp_dx, ttimetauMp*vderxy(0,0), shp);
                    assembler.template Matrix<Pres,Velx>(shp_dy, ttimetauMp*vderxy(1,0), shp);
                    assembler.template Matrix<Pres,Velx>(shp_dz, ttimetauMp*vderxy(2,0), shp);
                    
                    assembler.template Matrix<Pres,Vely>(shp_dx, ttimetauMp*vderxy(0,1), shp);
                    assembler.template Matrix<Pres,Vely>(shp_dy, ttimetauMp*vderxy(1,1), shp);
                    assembler.template Matrix<Pres,Vely>(shp_dz, ttimetauMp*vderxy(2,1), shp);
                    
                    assembler.template Matrix<Pres,Velz>(shp_dx, ttimetauMp*vderxy(0,2), shp);
                    assembler.template Matrix<Pres,Velz>(shp_dy, ttimetauMp*vderxy(1,2), shp);
                    assembler.template Matrix<Pres,Velz>(shp_dz, ttimetauMp*vderxy(2,2), shp);
                }
                
                /* pressure stabilisation: viscosity (-L_visc_u) */
                /*
                           /                             \
                          |                         /  \  |
                        - |  nabla q , nabla o tau | Du | |
                          |                         \  /  |
                           \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp_dx, -2.0*visc*ttimetauMp, enr_viscs2_(0, 0, _));
                assembler.template Matrix<Pres,Vely>(shp_dx, -2.0*visc*ttimetauMp, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Pres,Velz>(shp_dx, -2.0*visc*ttimetauMp, enr_viscs2_(0, 2, _));
                
                assembler.template Matrix<Pres,Velx>(shp_dy, -2.0*visc*ttimetauMp, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Pres,Vely>(shp_dy, -2.0*visc*ttimetauMp, enr_viscs2_(1, 1, _));
                assembler.template Matrix<Pres,Velz>(shp_dy, -2.0*visc*ttimetauMp, enr_viscs2_(1, 2, _));
                
                assembler.template Matrix<Pres,Velx>(shp_dz, -2.0*visc*ttimetauMp, enr_viscs2_(0, 2, _));
                assembler.template Matrix<Pres,Vely>(shp_dz, -2.0*visc*ttimetauMp, enr_viscs2_(1, 2, _));
                assembler.template Matrix<Pres,Velz>(shp_dz, -2.0*visc*ttimetauMp, enr_viscs2_(2, 2, _));
                      
                /* pressure stabilisation: pressure( L_pres_p) */
                /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
                */
                assembler.template Matrix<Pres,Pres>(shp_dx, ttimetauMp, shp_dx);
                assembler.template Matrix<Pres,Pres>(shp_dy, ttimetauMp, shp_dy);
                assembler.template Matrix<Pres,Pres>(shp_dz, ttimetauMp, shp_dz);
                
                // pressure stabilization
                assembler.template Vector<Pres>(shp_dx, -timetauMp*res_old(0));
                assembler.template Vector<Pres>(shp_dy, -timetauMp*res_old(1));
                assembler.template Vector<Pres>(shp_dz, -timetauMp*res_old(2));
            }
            
            //----------------------------------------------------------------------
            //                     SUPG STABILISATION PART
            if(supg)
            {
                const double timetauM   = timefac * tau_M;
                if (instationary)
                {
                    /* supg stabilisation: inertia  */
                    /*
                              /                        \
                             |        / n+1       \     |
                             |  Du , | u   o nabla | v  |
                             |        \ (i)       /     |
                              \                        /
                    */
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, timetauM, shp);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, timetauM, shp);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, timetauM, shp);
                    
                    if (newton)
                    {
                        /* supg stabilisation: inertia, linearisation of testfunction  */
                        /*
                                   /                           \
                                  |   n+1      /          \     |
                                  |  u      , | Du o nabla | v  |
                                  |   (i)      \          /     |
                                   \                           /
    
                        */
                        assembler.template Matrix<Velx,Velx>(shp_dx, timetauM*velint(0), shp);
                        assembler.template Matrix<Velx,Vely>(shp_dy, timetauM*velint(0), shp);
                        assembler.template Matrix<Velx,Velz>(shp_dz, timetauM*velint(0), shp);
                        assembler.template Matrix<Vely,Velx>(shp_dx, timetauM*velint(1), shp);
                        assembler.template Matrix<Vely,Vely>(shp_dy, timetauM*velint(1), shp);
                        assembler.template Matrix<Vely,Velz>(shp_dz, timetauM*velint(1), shp);
                        assembler.template Matrix<Velz,Velx>(shp_dx, timetauM*velint(2), shp);
                        assembler.template Matrix<Velz,Vely>(shp_dy, timetauM*velint(2), shp);
                        assembler.template Matrix<Velz,Velz>(shp_dz, timetauM*velint(2), shp);
                    }
                }
                const double ttimetauM  = timefac * timefac * tau_M;
                /* supg stabilisation: convective part ( L_conv_u) */
                /*
                     /                                          \
                    |  / n+1        \        / n+1        \      |
                    | | u    o nabla | v ,  | u    o nabla | Du  |
                    |  \ (i)        /        \ (i)        /      |
                     \                                          /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM, enr_conv_c_);
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM, enr_conv_c_);
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM, enr_conv_c_);
                /* supg stabilisation: pressure part  ( L_pres_p) */
                /*
                          /                             \
                         |   / n+1       \               |
                         |  | u   o nabla | v , nabla Dp |
                         |   \ (i)       /               |
                          \                             /
                */
                assembler.template Matrix<Velx,Pres>(enr_conv_c_, ttimetauM, shp_dx);
                assembler.template Matrix<Vely,Pres>(enr_conv_c_, ttimetauM, shp_dy);
                assembler.template Matrix<Velz,Pres>(enr_conv_c_, ttimetauM, shp_dz);

                /* supg stabilisation: viscous part  (-L_visc_u) */
                /*
                      /                                        \
                     |               /  \    / n+1        \     |
                   - |  nabla o eps | Du |, | u    o nabla | v  |
                     |               \  /    \ (i)        /     |
                      \                                        /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(0, 0, _));
                assembler.template Matrix<Velx,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Velx,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(0, 2, _));
    
                assembler.template Matrix<Vely,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(1, 1, _));
                assembler.template Matrix<Vely,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(1, 2, _));
    
                assembler.template Matrix<Velz,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(0, 2, _));
                assembler.template Matrix<Velz,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(1, 2, _));
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_(2, 2, _));
                
                if (newton)
                {
                    /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                    /*
                               /                                           \
                              |    /          \   n+1    / n+1        \     |
                              |   | Du o nabla | u    , | u    o nabla | v  |
                              |    \          /   (i)    \ (i)        /     |
                               \                                           /
                    */
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*vderxy(0,0), shp);
                    assembler.template Matrix<Velx,Vely>(enr_conv_c_, ttimetauM*vderxy(0,1), shp);
                    assembler.template Matrix<Velx,Velz>(enr_conv_c_, ttimetauM*vderxy(0,2), shp);
                    
                    assembler.template Matrix<Vely,Velx>(enr_conv_c_, ttimetauM*vderxy(1,0), shp);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*vderxy(1,1), shp);                    
                    assembler.template Matrix<Vely,Velz>(enr_conv_c_, ttimetauM*vderxy(1,2), shp);                    
                    
                    assembler.template Matrix<Velz,Velx>(enr_conv_c_, ttimetauM*vderxy(2,0), shp);
                    assembler.template Matrix<Velz,Vely>(enr_conv_c_, ttimetauM*vderxy(2,1), shp);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*vderxy(2,2), shp);
                    
                    /*
                             /                                           \
                            |    / n+1        \   n+1    /          \     |
                            |   | u    o nabla | u    , | Du o nabla | v  |
                            |    \ (i)        /   (i)    \          /     |
                             \                                           /
                    */
                    const double con0 = ttimetauM*(velint(0)*vderxy(0,0) + velint(1)*vderxy(0,1) + velint(2)*vderxy(0,2));
                    assembler.template Matrix<Velx,Velx>(shp_dx, con0, shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, con0, shp);                        
                    assembler.template Matrix<Velx,Velz>(shp_dz, con0, shp); 
                    
                    const double con1 = ttimetauM*(velint(0)*vderxy(1,0) + velint(1)*vderxy(1,1) + velint(2)*vderxy(1,2));
                    assembler.template Matrix<Vely,Velx>(shp_dx, con1, shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, con1, shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, con1, shp);
                    
                    const double con2 = ttimetauM*(velint(0)*vderxy(2,0) + velint(1)*vderxy(2,1) + velint(2)*vderxy(2,2));
                    assembler.template Matrix<Velz,Velx>(shp_dx, con2, shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, con2, shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, con2, shp);
                    
                    /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
                    /*
                                    /                               \
                                   |         n+1    /          \     |
                                   |  nabla p    , | Du o nabla | v  |
                                   |         (i)    \          /     |
                                    \                               /
                    */
                    assembler.template Matrix<Velx,Velx>(shp_dx, ttimetauM*gradp(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, ttimetauM*gradp(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, ttimetauM*gradp(0), shp);

                    assembler.template Matrix<Vely,Velx>(shp_dx, ttimetauM*gradp(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, ttimetauM*gradp(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, ttimetauM*gradp(1), shp);
                        
                    assembler.template Matrix<Velz,Velx>(shp_dx, ttimetauM*gradp(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, ttimetauM*gradp(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, ttimetauM*gradp(2), shp);

                      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                      /*
                              /                                         \
                             |               / n+1 \    /          \     |
                           - |  nabla o eps | u     |, | Du o nabla | v  |
                             |               \ (i) /    \          /     |
                              \                                         /
                      */
                    assembler.template Matrix<Velx,Velx>(shp_dx, -2.0*visc*ttimetauM*visc_old(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, -2.0*visc*ttimetauM*visc_old(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, -2.0*visc*ttimetauM*visc_old(0), shp);

                    assembler.template Matrix<Vely,Velx>(shp_dx, -2.0*visc*ttimetauM*visc_old(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, -2.0*visc*ttimetauM*visc_old(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, -2.0*visc*ttimetauM*visc_old(1), shp);
                        
                    assembler.template Matrix<Velz,Velx>(shp_dx, -2.0*visc*ttimetauM*visc_old(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, -2.0*visc*ttimetauM*visc_old(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, -2.0*visc*ttimetauM*visc_old(2), shp);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                                  /                             \
                                 |              /          \     |
                               - |  rhsint   , | Du o nabla | v  |
                                 |              \          /     |
                                  \                             /

                    */
                    assembler.template Matrix<Velx,Velx>(shp_dx, -timetauM*rhsint(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, -timetauM*rhsint(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, -timetauM*rhsint(0), shp);
                    
                    assembler.template Matrix<Vely,Velx>(shp_dx, -timetauM*rhsint(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, -timetauM*rhsint(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, -timetauM*rhsint(1), shp);
                    
                    assembler.template Matrix<Velz,Velx>(shp_dx, -timetauM*rhsint(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, -timetauM*rhsint(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, -timetauM*rhsint(2), shp);
                } // if newton
                
                // supg stabilisation
                assembler.template Vector<Velx>(enr_conv_c_, -timetauM*res_old(0));
                assembler.template Vector<Vely>(enr_conv_c_, -timetauM*res_old(1));
                assembler.template Vector<Velz>(enr_conv_c_, -timetauM*res_old(2));
            }

            
            //----------------------------------------------------------------------
            //                     STABILISATION, CONTINUITY PART
            if(cstab)
            {
                const double timefac_timefac_tau_C=timefac*timefac*tau_C;
                const double timefac_timefac_tau_C_divunp=timefac_timefac_tau_C*(vderxy(0, 0)+vderxy(1, 1)+vderxy(2, 2));
                /* continuity stabilisation on left hand side */
                /*
                         /                        \
                        |                          |
                        | nabla o Du  , nabla o v  |
                        |                          |
                         \                        /
                */
                assembler.template Matrix<Velx,Velx>(shp_dx, timefac_timefac_tau_C, shp_dx);
                assembler.template Matrix<Velx,Vely>(shp_dx, timefac_timefac_tau_C, shp_dy);
                assembler.template Matrix<Velx,Velz>(shp_dx, timefac_timefac_tau_C, shp_dz);
                
                assembler.template Matrix<Vely,Velx>(shp_dy, timefac_timefac_tau_C, shp_dx);
                assembler.template Matrix<Vely,Vely>(shp_dy, timefac_timefac_tau_C, shp_dy);
                assembler.template Matrix<Vely,Velz>(shp_dy, timefac_timefac_tau_C, shp_dz);
                
                assembler.template Matrix<Velz,Velx>(shp_dz, timefac_timefac_tau_C, shp_dx);
                assembler.template Matrix<Velz,Vely>(shp_dz, timefac_timefac_tau_C, shp_dy);
                assembler.template Matrix<Velz,Velz>(shp_dz, timefac_timefac_tau_C, shp_dz);
                
                assembler.template Vector<Velx>(shp_dx, -timefac_timefac_tau_C_divunp);
                assembler.template Vector<Vely>(shp_dy, -timefac_timefac_tau_C_divunp);
                assembler.template Vector<Velz>(shp_dz, -timefac_timefac_tau_C_divunp);
            } // endif cstab
        } // end loop over gauss points
    } // end loop over integration cells
    

    //if (false)
    if (ASSTYPE == XFEM::xfem_assembly)
    {
        // for now, I don't try to compare to elements wiout stress unknowns, since they lock anyway
        if (stress_unknowns_present){
    
    // information about boundary integration cells
    const XFEM::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());
    
    
    
    bool printcell = true;
    
    if (printcell)
    {
        printcell = false;
        std::ofstream f_system("example.pos");
        
        
        {
            // draw elements with associated gid
            std::stringstream gmshfilecontent;
            gmshfilecontent << "View \" " << "DomainIntegrationCells \" {" << endl;
            XFEM::DomainIntCells::const_iterator cell;
            for (cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
            {
                //gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXYZ(*ele), 0.0, cell->Shape()) << endl;
                gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXiDomainBlitz(), 0.0, cell->Shape()) << endl;
            };
            gmshfilecontent << "};" << endl;
            f_system << gmshfilecontent.str();
        }
        {
            // draw elements with associated gid
            std::stringstream gmshfilecontent;
            gmshfilecontent << "View \" " << "BoundaryIntegrationCells \" {" << endl;
            XFEM::BoundaryIntCells::const_iterator cell;
            for (cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
            {
                //gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXYZ(*ele), 0.0, cell->Shape()) << endl;
                gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXiDomainBlitz(), 0.0, cell->Shape()) << endl;
            };
            gmshfilecontent << "};" << endl;
            f_system << gmshfilecontent.str();
        }
        f_system.close();

    }
    
    
    
    
    
    // loop over boundary integration cells
    for (XFEM::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
        
        //TODO: get the right element
        const DRT::Element* boundaryele = ih->GetBoundaryEle(cell->GetSurfaceEleGid());
        //const DRT::Element* boundaryele = ih->GetBoundaryEle(0);
        //cout << (*boundaryele) << endl;
        
        // get node coordinates
        const BlitzMat xyze_boundary(DRT::UTILS::PositionArrayBlitz(boundaryele));

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta^\boundary
            XFEM::PosEtaBoundary pos_eta_boundary(2);
            pos_eta_boundary(0) = intpoints.qxg[iquad][0];
            pos_eta_boundary(1) = intpoints.qxg[iquad][1];
      
            // coordinates of the current integration point in element coordinates \xi
            const XFEM::PosXiBoundary posXiBoundary(mapEtaBToXiB(*cell, pos_eta_boundary));
            const XFEM::PosXiDomain posXiDomain(mapEtaBToXiD(*cell, pos_eta_boundary));
            const double detcell = abs(detEtaBToXiB(*cell, pos_eta_boundary)); //TODO: check normals

            // shape functions and their first derivatives
            const BlitzVec funct_boundary(DRT::UTILS::shape_function_2D(posXiBoundary(0),posXiBoundary(1),boundaryele->Shape()));
            const BlitzMat deriv_boundary(DRT::UTILS::shape_function_2D_deriv1(posXiBoundary(0),posXiBoundary(1),boundaryele->Shape()));
            
            // shape functions and their first derivatives
            const BlitzVec funct(DRT::UTILS::shape_function<DISTYPE>(posXiDomain));
      
            const BlitzVec funct_stress(DRT::UTILS::shape_function_3D(posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype));

            // position of the gausspoint in physical coordinates
            const BlitzVec gauss_pos_xyz(blitz::sum(funct_boundary(j)*xyze_boundary(i,j),j));
      
            // get jacobian matrix d x / d \xi  (3x2)
            const BlitzMat dxyzdrs(blitz::sum(xyze_boundary(i,k)*deriv_boundary(j,k),k));
            
            // compute covariant metric tensor G for surface element (2x2)
            const BlitzMat metric(blitz::sum(dxyzdrs(k,i)*dxyzdrs(k,j),k));
            //const BlitzMat metric = computeMetricTensor(xyze_boundary,deriv_boundary);
            
            // compute global derivates
            //const BlitzMat derxy(blitz::sum(dxyzdrs(i,k)*deriv_boundary(k,j),k));
            //const BlitzMat derxy_stress(blitz::sum(xji(i,k)*deriv_stress(k,j),k));
            
            const double detmetric = sqrt(metric(0,0)*metric(1,1) - metric(0,1)*metric(1,0));
            const double fac = intpoints.qwgt[iquad]*detmetric*detcell;

            if (detmetric < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detmetric);
            }

            // compute second global derivative
            static BlitzMat derxy2(6,numnode,blitz::ColumnMajorArray<2>());
            derxy2 = 0.;

            // after this call, one should only use the enriched shape functions and derivatives!
            BlitzVec enr_funct(numparamvelx);
            BlitzMat enr_derxy(3,numparamvelx,blitz::ColumnMajorArray<2>());
            BlitzMat enr_derxy2(6,numparamvelx,blitz::ColumnMajorArray<2>());
            BlitzVec enr_funct_stress(numparamsigmaxx);
            
            // shape function for nodal dofs
            XFEM::ComputeEnrichedNodalShapefunction(
                    *ele,
                    ih,
                    dofman,
                    Velx,
                    gauss_pos_xyz,
                    XFEM::Enrichment::approachFromPlus,
                    funct,
                    enr_funct);
            
//            // shape function for nodal dofs
//            XFEM::ComputeEnrichedNodalShapefunction(
//                    *ele,
//                    ih,
//                    dofman,
//                    Velx,
//                    gauss_pos_xyz,
//                    XFEM::Enrichment::approachFromPlus,
//                    funct,
//                    derxy,
//                    derxy2,
//                    enr_funct,
//                    enr_derxy,
//                    enr_derxy2);

            // shape functions for element dofs
            XFEM::ComputeEnrichedElementShapefunction(
                    *ele,
                    ih,
                    dofman,
                    Sigmaxx,
                    gauss_pos_xyz,
                    XFEM::Enrichment::approachFromPlus,
                    funct_stress,
                    enr_funct_stress);
                
            // perform integration for entire matrix and rhs
            // create vievs on the shape function arrays for easy handling in the assembly process
            const BlitzVec shp(enr_funct);
//            const BlitzVec shp_dx(enr_derxy(0,_));
//            const BlitzVec shp_dy(enr_derxy(1,_));
//            const BlitzVec shp_dz(enr_derxy(2,_));
            const BlitzVec shp_sigma(enr_funct_stress);
            
            // get normal vector (in x coordinates) to surface element at integration point
            static BlitzVec3 normalvec;
            computeNormalToBoundaryElement(boundaryele, deriv_boundary, xyze_boundary, normalvec);
      
            // get velocities (n+g,i) at integration point
            const BlitzVec velint(blitz::sum(shp(j)*evelnp(i,j),j));
            
            // get pressure
            //const double press(blitz::sum(shp*eprenp));

            // get viscous stress unknowns
            static BlitzMat sigma(3,3);
            sigma(0,0) = blitz::sum(shp_sigma(_)*esigma(0,_)); sigma(0,1) = blitz::sum(shp_sigma(_)*esigma(3,_)); sigma(0,2) = blitz::sum(shp_sigma(_)*esigma(4,_));
            sigma(1,0) = sigma(0,1);                           sigma(1,1) = blitz::sum(shp_sigma(_)*esigma(1,_)); sigma(1,2) = blitz::sum(shp_sigma(_)*esigma(5,_));
            sigma(2,0) = sigma(0,2);                           sigma(2,1) = sigma(1,2);                           sigma(2,2) = blitz::sum(shp_sigma(_)*esigma(2,_));


//            // get velocity norm
//            const double vel_norm = sqrt(blitz::sum(velint*velint));
//
//            // normed velocity at element centre
//            static BlitzVec velino(3);
//            if (vel_norm>=1e-6)
//            {
//                velino = velint/vel_norm;
//            }
//            else
//            {
//                velino = 0.;
//                velino(0) = 1;
//            }
//
//            // get streamlength
//            const double val = blitz::sum(blitz::abs(blitz::sum(velino(j)*derxy(j,i),j)));
//            const double strle = 2.0/val;
//
//            // calculate tau: stabilization parameters for stationary case
//
//            // compute tau_Mu
//            const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
//            const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
//            const double tau_stab_0 = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);
//
//            // compute tau_Mp
//            const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
//            const double xi_tau_mp = DMAX(re_tau_mp,1.0);
//            const double tau_stab_1 = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);
//
//            // compute tau_C
//            const double xi_tau_c = DMIN(re_tau_mp, 1.0);
//            const double tau_stab_2 = 0.5*vel_norm*hk*xi_tau_c;
//      
//            // stabilisation parameter
//            const double tau_M  = tau_stab_0*fac;
//            const double tau_Mp = tau_stab_1*fac;
//            const double tau_C  = tau_stab_2*fac;
            
            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;
            

            /*--------------------------------- now build single stiffness terms ---*/

            /*
                           /                   \
                          |                     |
                        - | (virt tau) * n , Du |
                          |                     |
                           \                   /
            */
            const double taue_u_factor = 1.0;
            assembler.template Matrix<Sigmaxx,Velx>(shp_sigma, -taue_u_factor*timefacfac*normalvec(0), shp);
            assembler.template Matrix<Sigmaxy,Velx>(shp_sigma, -taue_u_factor*timefacfac*normalvec(1), shp);
            assembler.template Matrix<Sigmaxz,Velx>(shp_sigma, -taue_u_factor*timefacfac*normalvec(2), shp);
            assembler.template Matrix<Sigmayx,Vely>(shp_sigma, -taue_u_factor*timefacfac*normalvec(0), shp);
            assembler.template Matrix<Sigmayy,Vely>(shp_sigma, -taue_u_factor*timefacfac*normalvec(1), shp);
            assembler.template Matrix<Sigmayz,Vely>(shp_sigma, -taue_u_factor*timefacfac*normalvec(2), shp);
            assembler.template Matrix<Sigmazx,Velz>(shp_sigma, -taue_u_factor*timefacfac*normalvec(0), shp);
            assembler.template Matrix<Sigmazy,Velz>(shp_sigma, -taue_u_factor*timefacfac*normalvec(1), shp);
            assembler.template Matrix<Sigmazz,Velz>(shp_sigma, -taue_u_factor*timefacfac*normalvec(2), shp);
            
            assembler.template Vector<Sigmaxx>(shp_sigma, taue_u_factor*timefacfac*normalvec(0)*velint(0));
            assembler.template Vector<Sigmaxy>(shp_sigma, taue_u_factor*timefacfac*normalvec(1)*velint(0));
            assembler.template Vector<Sigmaxz>(shp_sigma, taue_u_factor*timefacfac*normalvec(2)*velint(0));
            assembler.template Vector<Sigmayx>(shp_sigma, taue_u_factor*timefacfac*normalvec(0)*velint(1));
            assembler.template Vector<Sigmayy>(shp_sigma, taue_u_factor*timefacfac*normalvec(1)*velint(1));
            assembler.template Vector<Sigmayz>(shp_sigma, taue_u_factor*timefacfac*normalvec(2)*velint(1));
            assembler.template Vector<Sigmazx>(shp_sigma, taue_u_factor*timefacfac*normalvec(0)*velint(2));
            assembler.template Vector<Sigmazy>(shp_sigma, taue_u_factor*timefacfac*normalvec(1)*velint(2));
            assembler.template Vector<Sigmazz>(shp_sigma, taue_u_factor*timefacfac*normalvec(2)*velint(2));
            
            // symmetric part from partial integration of stress term
            /*
                           /            \
                          |              |
                        - | v , Dtau * n |
                          |              |
                           \            /
            */
            assembler.template Matrix<Velx,Sigmaxx>(shp, -timefacfac*normalvec(0), shp_sigma);
            assembler.template Matrix<Velx,Sigmaxy>(shp, -timefacfac*normalvec(1), shp_sigma);
            assembler.template Matrix<Velx,Sigmaxz>(shp, -timefacfac*normalvec(2), shp_sigma);
            assembler.template Matrix<Vely,Sigmayx>(shp, -timefacfac*normalvec(0), shp_sigma);
            assembler.template Matrix<Vely,Sigmayy>(shp, -timefacfac*normalvec(1), shp_sigma);
            assembler.template Matrix<Vely,Sigmayz>(shp, -timefacfac*normalvec(2), shp_sigma);
            assembler.template Matrix<Velz,Sigmazx>(shp, -timefacfac*normalvec(0), shp_sigma);
            assembler.template Matrix<Velz,Sigmazy>(shp, -timefacfac*normalvec(1), shp_sigma);
            assembler.template Matrix<Velz,Sigmazz>(shp, -timefacfac*normalvec(2), shp_sigma);
            
            const BlitzVec tau_times_n(blitz::sum(sigma(i,j)*normalvec(j),j));
            assembler.template Vector<Velx>(shp, timefacfac*tau_times_n(0));
            assembler.template Vector<Vely>(shp, timefacfac*tau_times_n(1));
            assembler.template Vector<Velz>(shp, timefacfac*tau_times_n(2));
        
            
            
            
            if (pstab)
            {
                // symmetric part from partial integration of stress term
                /*
                               /                 \
                              |                   |
                              | nabla q , tau * n |
                              |                   |
                               \                 /
                */
//                assembler.template Matrix<Pres,Sigmaxx>(shp_dx, -tau_Mp*normalvec(0), shp_sigma);
//                assembler.template Matrix<Pres,Sigmaxy>(shp_dx, -tau_Mp*normalvec(1), shp_sigma);
//                assembler.template Matrix<Pres,Sigmaxz>(shp_dx, -tau_Mp*normalvec(2), shp_sigma);
//                
//                assembler.template Matrix<Pres,Sigmayx>(shp_dy, -tau_Mp*normalvec(0), shp_sigma);
//                assembler.template Matrix<Pres,Sigmayy>(shp_dy, -tau_Mp*normalvec(1), shp_sigma);
//                assembler.template Matrix<Pres,Sigmayz>(shp_dy, -tau_Mp*normalvec(2), shp_sigma);
//                
//                assembler.template Matrix<Pres,Sigmazx>(shp_dz, -tau_Mp*normalvec(0), shp_sigma);
//                assembler.template Matrix<Pres,Sigmazy>(shp_dz, -tau_Mp*normalvec(1), shp_sigma);
//                assembler.template Matrix<Pres,Sigmazz>(shp_dz, -tau_Mp*normalvec(2), shp_sigma);
//                
//                const BlitzVec tau_times_n(blitz::sum(sigma(i,j)*normalvec(j),j));
//                assembler.template Vector<Pres>(shp_dx, tau_Mp*tau_times_n(0));
//                assembler.template Vector<Pres>(shp_dy, tau_Mp*tau_times_n(1));
//                assembler.template Vector<Pres>(shp_dz, tau_Mp*tau_times_n(2));
            }
        } // end loop over gauss points
    } // end loop over boundary integration cells
    }
    } // if (ASSTYPE == XFEM::xfem_assembly)
    return;
}

/*!
 * \brief entry point for Sysmat call
 * at one point, one has to call specific template intantiations of Sysmat using the current Shape() of the element.
 * This is the point.
 */
void callSysmat(
        const XFEM::AssemblyType          assembly_type,
        const DRT::ELEMENTS::XFluid3*     ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const std::vector<double>&        locval,
        const std::vector<double>&        locval_hist,
        BlitzMat&                         estif,
        BlitzVec&                         eforce,
        const struct _MATERIAL*           material,
        const double                      time,          ///< current time (pseudotime for stationary formulation)
        const double                      timefac,       ///< One-step-Theta: theta*dt, BDF2: 2/3 * dt
        const bool                        newton ,
        const bool                        pstab  ,
        const bool                        supg   ,
        const bool                        cstab  ,
        const bool                        instationary
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            case DRT::Element::tet4:
                Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            default:
                dserror("Sysmat not templated yet");
        };
    }
    else
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, locval_hist, estif, eforce,
                        material, time, timefac, newton, pstab, supg, cstab, instationary);
                break;
            default:
                dserror("Sysmat not templated yet");
        };
    }
}
} // end namespace XFLUID

#endif

#endif
#endif
