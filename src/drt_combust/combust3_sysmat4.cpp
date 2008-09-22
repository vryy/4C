/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat4.cpp

\brief element formulations for 3d Combustfluid element

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat.H"
#include "combust3_utils.H"
#include "../drt_f3/fluid3_stabilization.H"
#include "combust3_local_assembler.H"
#include "combust3_interpolation.H"
#include "../drt_geometry/coordinate_transformation.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_xfem/spacetime_boundary.H"
#include "../drt_lib/drt_utils.H"

template <int shpVecSize>
struct Shp
{
  blitz::TinyVector<double,shpVecSize>  d0;
  blitz::TinyVector<double,shpVecSize>  dx;
  blitz::TinyVector<double,shpVecSize>  dy;
  blitz::TinyVector<double,shpVecSize>  dz;
  blitz::TinyVector<double,shpVecSize>  dxdx;
  blitz::TinyVector<double,shpVecSize>  dxdy;
  blitz::TinyVector<double,shpVecSize>  dxdz;
  blitz::TinyVector<double,shpVecSize>  dydx;
  blitz::TinyVector<double,shpVecSize>  dydy;
  blitz::TinyVector<double,shpVecSize>  dydz;
  blitz::TinyVector<double,shpVecSize>  dzdx;
  blitz::TinyVector<double,shpVecSize>  dzdy;
  blitz::TinyVector<double,shpVecSize>  dzdz;
};

using namespace XFEM::PHYSICS;

//namespace COMBUST
//{

  //! size factor to allow static arrays
  ///
  /// to allow static arrays for a unknown number of unknowns, we make them bigger than necessary
  /// this factor is multiplied times numnode(distype) to get the size of many arrays
  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFac {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFac<XFEM::standard_assembly> {static const int fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFac<XFEM::xfem_assembly>     {static const int fac = 3;};

  //! interpolate from nodal vector array to integration point vector using the shape function
  template <class M, class MS>
  static BlitzVec3 interpolateVectorFieldToIntPoint(
      const M&  eleVectorField,       ///< array with nodal vector values
      const MS& shp,                  ///< array with nodal shape function
      const int numparam              ///< number of parameters
      )
  {
    BlitzVec3 v;
    const int nsd = 3;
    for (int isd = 0; isd < nsd; ++isd)
    {
        v(isd) = 0.0;
        for (int iparam = 0; iparam < numparam; ++iparam)
            v(isd) += eleVectorField(isd,iparam)*shp(iparam);
    }
    return v;
  }
  
  template<class M>
  static bool modifyOldTimeStepsValues(
      const DRT::Element*                        ele,           ///< the element those matrix is calculated
      const Teuchos::RCP<XFEM::InterfaceHandle>  ih,   ///< connection to the interface handler
      const M&                                   xyze,
      const BlitzVec3&                           posXiDomain,
      const int                                  labelnp,
      const Epetra_Vector&                       ivelcoln,
      const Epetra_Vector&                       iacccoln,
      BlitzVec3&                                 gpveln,
      BlitzVec3&                                 gpaccn
      )
  {
    GEO::PosX posx_gp;
    GEO::elementToCurrentCoordinates(ele, xyze, posXiDomain, posx_gp);
    
    const bool is_in_fluid = (labelnp == 0);
    
    if (not is_in_fluid)
    {
      //std::cout << "should I arrive here?" << std::endl;
      return false;
    }

    const bool was_in_fluid = (ih->PositionWithinConditionN(posx_gp) == 0);
    
    const bool in_space_time_slab_area = (is_in_fluid and (not was_in_fluid));
    
    if (in_space_time_slab_area)
    {
      XFEM::SpaceTimeBoundaryCell slab;
      BlitzVec3 rst;
      
      const bool found_cell = ih->FindSpaceTimeLayerCell(posx_gp,slab,rst);
      
      if (found_cell)
      {
        
        const double delta_slab = -(rst(2)-1.0)*0.5;
  
        if (delta_slab > (1.0+1.0e-7) or -1.0e-7 > delta_slab)
        {
          cout << rst(2) <<  "  " << delta_slab << endl;
          cout << slab.toString() << endl << endl;
          dserror("wrong value of delta_slab");
        }
//              theta_dt = theta_dt_pure;// * delta_slab;
      
        DRT::Element* boundaryele = ih->cutterdis()->gElement(slab.getBeleId());
        const int numnode_boundary = boundaryele->NumNode();
        
        BlitzVec3 iveln;
        BlitzVec3 iaccn;
        
        if (ivelcoln.GlobalLength() > 3)
        {
          // get interface velocities at the boundary element nodes
          BlitzMat veln_boundary(3,numnode_boundary*2);
          BlitzMat accn_boundary(3,numnode_boundary*2);
          if (numnode_boundary != 4)
            dserror("needs more generalizashun!");
          const DRT::Node*const* nodes = boundaryele->Nodes();
          
          for (int inode = 0; inode < numnode_boundary; ++inode)
          {
            const DRT::Node* node = nodes[inode];
            const std::vector<int> lm = ih->cutterdis()->Dof(node);
            std::vector<double> myvel(3);
            DRT::UTILS::ExtractMyValues(ivelcoln,myvel,lm);
            veln_boundary(0,inode) = myvel[0];
            veln_boundary(1,inode) = myvel[1];
            veln_boundary(2,inode) = myvel[2];
            DRT::UTILS::ExtractMyValues(ivelcoln,myvel,lm);
            veln_boundary(0,inode+4) = myvel[0];
            veln_boundary(1,inode+4) = myvel[1];
            veln_boundary(2,inode+4) = myvel[2];
            
            std::vector<double> myacc(3);
            DRT::UTILS::ExtractMyValues(iacccoln,myacc,lm);
            accn_boundary(0,inode) = myacc[0];
            accn_boundary(1,inode) = myacc[1];
            accn_boundary(2,inode) = myacc[2];
            DRT::UTILS::ExtractMyValues(iacccoln,myacc,lm);
            accn_boundary(0,inode+4) = myacc[0];
            accn_boundary(1,inode+4) = myacc[1];
            accn_boundary(2,inode+4) = myacc[2];
          }
  //                cout << "veln_boundary: " << veln_boundary << endl;
  //                cout << "accn_boundary: " << accn_boundary << endl;
          
          BlitzVec funct_ST(numnode_boundary*2);
          DRT::UTILS::shape_function_3D(funct_ST,rst(0),rst(1),rst(2),DRT::Element::hex8);
          iveln  = interpolateVectorFieldToIntPoint(veln_boundary , funct_ST, 8);
          iaccn  = interpolateVectorFieldToIntPoint(accn_boundary , funct_ST, 8);
  //                cout << "iveln " << iveln << endl;
  //                cout << "iaccn " << iaccn << endl;
        }
        else
        {
          iveln = 0.0;
          iaccn = 0.0;
        }
      
        //cout << "using space time boundary values " << ele->Id() << endl; 
        gpveln(0) = iveln(0);
        gpveln(1) = iveln(1);
        gpveln(2) = iveln(2);
        
        gpaccn(0) = iaccn(0);
        gpaccn(1) = iaccn(1);
        gpaccn(2) = iaccn(2);
      }
    }
    return true;
  }
  
  //! fill a number of arrays with unknown values from the unknown vector given by the discretization
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            class M1, class V1, class M2>
  void fillElementUnknownsArrays4(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::Combust3::MyState mystate,
          M1& evelnp,
          M1& eveln,
          M1& evelnm,
          M1& eaccn,
          V1& eprenp,
          M2& etau
          )
  {
      
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      
      // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
      //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
      const int numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
      const int numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
      const int numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
      const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
      // put one here to create arrays of size 1, since they are not needed anyway
      // in the xfem assembly, the numparam is determined by the dofmanager
      //const int numparamtauxx = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 1);
      const int numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
           
      const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
      const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
      const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
      const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());
      
      for (int iparam=0; iparam<numparamvelx; ++iparam)
      {
          evelnp(0,iparam) = mystate.velnp[velxdof[iparam]];
          eveln( 0,iparam) = mystate.veln[ velxdof[iparam]];
          evelnm(0,iparam) = mystate.velnm[velxdof[iparam]];
          eaccn( 0,iparam) = mystate.accn[ velxdof[iparam]];
      }
      for (int iparam=0; iparam<numparamvely; ++iparam)
      {
          evelnp(1,iparam) = mystate.velnp[velydof[iparam]];
          eveln( 1,iparam) = mystate.veln[ velydof[iparam]];
          evelnm(1,iparam) = mystate.velnm[velydof[iparam]];
          eaccn( 1,iparam) = mystate.accn[ velydof[iparam]];
      }
      for (int iparam=0; iparam<numparamvelz; ++iparam)
      {
          evelnp(2,iparam) = mystate.velnp[velzdof[iparam]];
          eveln( 2,iparam) = mystate.veln[ velzdof[iparam]];
          evelnm(2,iparam) = mystate.velnm[velzdof[iparam]];
          eaccn( 2,iparam) = mystate.accn[ velzdof[iparam]];
      }
      for (int iparam=0; iparam<numparampres; ++iparam)
          eprenp(iparam) = mystate.velnp[presdof[iparam]];
      const bool tauele_unknowns_present = (XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 0) > 0);
      if (tauele_unknowns_present)
      {
          const int numparamtauyy = XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayy, 1);
          const int numparamtauzz = XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmazz, 1);
          const int numparamtauxy = XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxy, 1);
          const int numparamtauxz = XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxz, 1);
          const int numparamtauyz = XFEM::comgetNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayz, 1);
          const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxx>());
          const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayy>());
          const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmazz>());
          const std::vector<int>& tauxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxy>());
          const std::vector<int>& tauxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxz>());
          const std::vector<int>& tauyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayz>());
          for (int iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = mystate.velnp[tauxxdof[iparam]];
          for (int iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = mystate.velnp[tauyydof[iparam]];
          for (int iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = mystate.velnp[tauzzdof[iparam]];
          for (int iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = mystate.velnp[tauxydof[iparam]];
          for (int iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = mystate.velnp[tauxzdof[iparam]];
          for (int iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = mystate.velnp[tauyzdof[iparam]];
      }
  }
  
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M1, class V1, class M2>
static void SysmatDomain4(
    const DRT::Element*                 ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandle>  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const V1&                           eprenp,
    const M2&                           etau,
    const Teuchos::RCP<const Epetra_Vector>   ivelcol,       ///< velocity for interface nodes
    const Teuchos::RCP<Epetra_Vector>   iforcecol,     ///< reaction force due to given interface velocity
    const struct _MATERIAL*             material,      ///< fluid material
    const FLUID_TIMEINTTYPE             timealgo,      ///< time discretization type
    const double                        dt,            ///< delta t (time step size)
    const double                        theta,         ///< factor for one step theta scheme
    const bool                          newton,        ///< full Newton or fixed-point-like
    const bool                          pstab,         ///< flag for stabilization
    const bool                          supg,          ///< flag for stabilization
    const bool                          cstab,         ///< flag for stabilization
    const bool                          instationary,  ///< switch between stationary and instationary formulation
    const LocalAssembler<DISTYPE, ASSTYPE>& assembler
)
{
    TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat4 - domain");
  
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);
    
    // get node coordinates of the current element
    static blitz::TinyMatrix<double,nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // rigid body hack - assume structure is rigid and has uniform acceleration and velocity
//    const Epetra_Vector& ivelcolnp = *ih->cutterdis()->GetState("ivelcolnp");
    const Epetra_Vector& ivelcoln  = *ih->cutterdis()->GetState("ivelcoln");
    const Epetra_Vector& iacccoln  = *ih->cutterdis()->GetState("iacccoln");
    
    // dead load in element nodes
    /////////////////////////////////////////////////// , BlitzMat edeadng_(BodyForce(ele->Nodes(),time));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double visc = material->m.fluid->viscosity;

    // flag for higher order elements
    const bool higher_order_ele = XFEM::isHigherOrderElement<DISTYPE>();
    //const bool higher_order_ele = secondDerivativesAvailable<DISTYPE>();
    
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    
    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::comgetNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
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
    //const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
    const int numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
    const int numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
    
    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,xyze);
    const double mk = XFLUID::MK<DISTYPE>();
    
    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele->Id(),DISTYPE));
    //cout << "Element "<< ele->Id() << ": ";
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {
        const BlitzVec3 cellcenter(cell->GetPhysicalCenterPosition(*ele));
        
        // shortcut for intersected elements: if cell is only in solid domains for all influencing enrichments, skip it
        const int labelnp = ih->PositionWithinConditionNP(cellcenter);
        const std::set<int> xlabelset(dofman.getUniqueEnrichmentLabels());
        bool compute = false;
        if (labelnp == 0) // fluid
        {
          compute = true;
        }
        else if (xlabelset.size() > 1) // multiple interface labels
        {
          compute = true;
        }
        else if (xlabelset.find(labelnp) == xlabelset.end()) // ???
        {
          compute = true;
        }
        if (not compute)
        {
          continue;
        }
            
        const std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
              ih,
              dofman.getUniqueEnrichments(),
              cellcenter,
              XFEM::Enrichment::approachUnknown));
        
        const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule(ih->ElementIntersected(ele->Id()),cell->Shape(),ele->Shape());
        
        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            GEO::PosEtaDomain pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            GEO::PosXiDomain posXiDomain;
            GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static BlitzVec funct(numnode);
            static BlitzMat deriv(3, numnode, blitz::ColumnMajorArray<2>());
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
            DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            // discontinuous stress shape functions
            static BlitzVec funct_stress(DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement);
            static BlitzMat deriv_stress(3, DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,blitz::ColumnMajorArray<2>());
            if (ASSTYPE == XFEM::xfem_assembly)
            {
              if (tauele_unknowns_present)
              {
                DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
                DRT::UTILS::shape_function_3D_deriv1(deriv_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
              }
              else
              {
                funct_stress = 0.0;
                deriv_stress = 0.0;
              }
            }

            // position of the gausspoint in physical coordinates
            BlitzVec3 gauss_pos_xyz;
            BLITZTINY::MV_product<3,numnode>(xyze,funct,gauss_pos_xyz);
      
            // get transposed of the jacobian matrix d x / d \xi
            //xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);
            static BlitzMat3x3 xjm;
            BLITZTINY::MMt_product<3,3,numnode>(deriv,xyze,xjm);

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
            static BlitzMat3x3 xji;
            GEO::Inverse3x3(xjm, det, xji);

            // compute global derivates
            static BlitzMat derxy(3, numnode);
            //static BlitzMat derxy_stress(3, DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,blitz::ColumnMajorArray<2>());
            //static BlitzMat derxy_discpres(3, DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement,blitz::ColumnMajorArray<2>());
            //derxy          = blitz::sum(xji(i,k)*deriv(k,j),k);
            for (int isd = 0; isd < nsd; ++isd)
            {
              for (int inode = 0; inode < numnode; ++inode)
              {
                derxy(isd,inode) = 0.0;
                for (int jsd = 0; jsd < nsd; ++jsd)
                {
                   derxy(isd,inode) += xji(isd,jsd)*deriv(jsd,inode);
                }
              }
            }
//            for (int isd = 0; isd < nsd; ++isd)
//            {
//              for (int inode = 0; inode < DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement; ++inode)
//              {
//                derxy_discpres(isd,inode) = 0.0;
//                for (int jsd = 0; jsd < nsd; ++jsd)
//                {
//                  derxy_discpres(isd,inode) += xji(isd,jsd)*deriv_discpres(jsd,inode);
//                }
//              }
//            }

            // compute second global derivative
            static BlitzMat derxy2(6,numnode,blitz::ColumnMajorArray<2>());
            if (higher_order_ele)
            {
                static BlitzMat deriv2(6,numnode,blitz::ColumnMajorArray<2>());
                DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
                XFEM::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
            }
            else
            {
                derxy2 = 0.;
            }

            const int shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
            const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
            
            static Shp<shpVecSize> shp;
            
            typedef blitz::TinyVector<double,shpVecSize> ShpVec;
/*
            static ShpVec shp;
            static ShpVec shp_dx;
            static ShpVec shp_dy;
            static ShpVec shp_dz;
            static ShpVec shp_dxdx;
            static ShpVec shp_dxdy;
            static ShpVec shp_dxdz;
            static ShpVec shp_dydx;
            static ShpVec shp_dydy;
            static ShpVec shp_dydz;
            static ShpVec shp_dzdx;
            static ShpVec shp_dzdy;
            static ShpVec shp_dzdz;
*/
            
            static blitz::TinyVector<double,shpVecSizeStress> shp_tau;
            
            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // temporary arrays
                static BlitzVec enr_funct(3*numnode);
                static BlitzMat enr_derxy(3,3*numnode);
                static BlitzMat enr_derxy2(6,3*numnode);
                
//                const std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
//                      ih,
//                      dofman.getUniqueEnrichments(),
//                      gauss_pos_xyz,
//                      XFEM::Enrichment::approachUnknown));
              
                // shape function for nodal dofs
                XFEM::ComputeEnrichedNodalShapefunction(
                        *ele,
                        ih,
                        dofman,
                        Velx,
                        enrvals,
                        funct,
                        derxy,
                        derxy2, 
                        enr_funct,
                        enr_derxy,
                        enr_derxy2);
          
                for (int iparam = 0; iparam < numparamvelx; ++iparam)
                {
                  shp.d0(iparam) = enr_funct(iparam);
                  shp.dx(iparam) = enr_derxy(0,iparam);
                  shp.dy(iparam) = enr_derxy(1,iparam);
                  shp.dz(iparam) = enr_derxy(2,iparam);
                  shp.dxdx(iparam) = enr_derxy2(0,iparam);
                  shp.dxdy(iparam) = enr_derxy2(3,iparam);
                  shp.dxdz(iparam) = enr_derxy2(4,iparam);
                  shp.dydx(iparam) = shp.dxdy(iparam);
                  shp.dydy(iparam) = enr_derxy2(1,iparam);
                  shp.dydz(iparam) = enr_derxy2(5,iparam);
                  shp.dzdx(iparam) = shp.dxdz(iparam);
                  shp.dzdy(iparam) = shp.dydz(iparam);
                  shp.dzdz(iparam) = enr_derxy2(2,iparam);
                }
                

                
                if (tauele_unknowns_present)
                {
                    static BlitzVec enr_funct_stress(3*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement);
                  
                    // shape functions for element dofs
                    XFEM::ComputeEnrichedElementShapefunction(
                            *ele,
                            ih,
                            dofman,
                            Sigmaxx,
                            enrvals,
                            funct_stress,
                            enr_funct_stress);
                    
                    for (int iparam = 0; iparam < numparamtauxx; ++iparam)
                    {
                      shp_tau(iparam) = enr_funct_stress(iparam);
                    }
                }
                else
                {
                  shp_tau = 0.0;
                }
            }
            else
            {
              for (int iparam = 0; iparam < numnode; ++iparam)
              {
                shp.d0(iparam) = funct(iparam);
                shp.dx(iparam) = derxy(0,iparam);
                shp.dy(iparam) = derxy(1,iparam);
                shp.dz(iparam) = derxy(2,iparam);
                shp.dxdx(iparam) = derxy2(0,iparam);
                shp.dxdy(iparam) = derxy2(3,iparam);
                shp.dxdz(iparam) = derxy2(4,iparam);
                shp.dydx(iparam) = shp.dxdy(iparam);
                shp.dydy(iparam) = derxy2(1,iparam);
                shp.dydz(iparam) = derxy2(5,iparam);
                shp.dzdx(iparam) = shp.dxdz(iparam);
                shp.dzdy(iparam) = shp.dydz(iparam);
                shp.dzdz(iparam) = derxy2(2,iparam);
              }
              
              if (tauele_unknowns_present)
              {
                dserror("no stress enrichments without xfem assembly");
              }
            }
            
            // get velocities and accelerations at integration point
            const BlitzVec3 gpvelnp = interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
            BlitzVec3 gpveln  = interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
            const BlitzVec3 gpvelnm = interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
            BlitzVec3 gpaccn  = interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);
            
            
            if (ASSTYPE == XFEM::xfem_assembly)
            {
              const bool valid_spacetime_cell_found = modifyOldTimeStepsValues(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, iacccoln, gpveln, gpaccn);
              if (not valid_spacetime_cell_found)
                continue;
            }
//            cout << gpvelnp << endl;
//            cout << evelnp << endl;
//            cout << shp << endl;

            // get history data (n) at integration point
//            BlitzVec3 histvec;
//            //histvec = blitz::sum(enr_funct(j)*evelnp_hist(i,j),j);
//            for (int isd = 0; isd < nsd; ++isd)
//            {
//                histvec(isd) = 0.0;
//                for (int iparam = 0; iparam < numparamvelx; ++iparam)
//                    histvec(isd) += evelnp_hist(isd,iparam)*shp.d0(iparam);
//            }
            const BlitzVec3 histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
                gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
            
            // get velocity (np,i) derivatives at integration point
            BlitzMat3x3 vderxy;
            //vderxy = blitz::sum(enr_derxy(j,k)*evelnp(i,k),k);
            vderxy = 0.0;
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
              for (int isd = 0; isd < nsd; ++isd)
              {
                vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
                vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
                vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
              }
            }
            
            //cout << "eps_xy" << (0.5*(vderxy(0,1)+vderxy(1,0))) << ", "<< endl;
      
            // calculate 2nd velocity derivatives at integration point
            static blitz::TinyMatrix<double,3,6> vderxy2;
            if (higher_order_ele)
            {
              //vderxy2 = blitz::sum(enr_derxy2(j,k)*evelnp(i,k),k);
              vderxy2 = 0.0;
              for (int iparam = 0; iparam < numparamvelx; ++iparam)
              {
                for (int isd = 0; isd < nsd; ++isd)
                {
                  vderxy2(isd,0) += evelnp(isd,iparam)*shp.dxdx(iparam);
                  vderxy2(isd,1) += evelnp(isd,iparam)*shp.dydy(iparam);
                  vderxy2(isd,2) += evelnp(isd,iparam)*shp.dzdz(iparam);
                  vderxy2(isd,3) += evelnp(isd,iparam)*shp.dxdy(iparam);
                  vderxy2(isd,4) += evelnp(isd,iparam)*shp.dxdz(iparam);
                  vderxy2(isd,5) += evelnp(isd,iparam)*shp.dydz(iparam);
                }
              }
            }
            else
            {
                vderxy2 = 0.;
            }

            // get pressure gradients
            BlitzVec3 gradp;
            //gradp = blitz::sum(enr_derxy(i,j)*eprenp(j),j);
            gradp = 0.0;
            for (int iparam = 0; iparam < numparampres; ++iparam)
            {
              gradp(0) += shp.dx(iparam)*eprenp(iparam);
              gradp(1) += shp.dy(iparam)*eprenp(iparam);
              gradp(2) += shp.dz(iparam)*eprenp(iparam);
            }
    
//            // get discont. pressure gradients
//            static BlitzVec3 graddiscp;
//            //gradp = blitz::sum(enr_derxy(i,j)*eprenp(j),j);
//            for (int isd = 0; isd < nsd; ++isd)
//            {
//                graddiscp(isd) = 0.0;
//                for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
//                    graddiscp(isd) += enr_derxy_discpres(isd,iparam)*ediscprenp(iparam);
//            }
            
            // get pressure
            //const double pres(blitz::sum(enr_funct*eprenp));
            double pres = 0.0;
            for (int iparam = 0; iparam < numparampres; ++iparam)
              pres += shp.d0(iparam)*eprenp(iparam);
            
            // get viscous stress unknowns
            BlitzMat3x3 tau;
            if (tauele_unknowns_present)
            {
              XFEM::fill_tau(numparamtauxx, shp_tau, etau, tau);
            }
            else
            {
              tau = 0.0;
            }
          
//            BlitzVec3 nabla_dot_tau;
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
//            BlitzVec3 bodyforce;
//            bodyforce = 0.0;
//            cout << bodyforce << endl;
            //////////////////////////////////////////BlitzVec bodyforce_(blitz::sum(enr_edeadng_(i,j)*enr_funct_(j),j));

            // compute stabilization parameters (3 taus)
            double tau_stab_M  = 0.0;
            double tau_stab_Mp = 0.0;
            double tau_stab_C  = 0.0;
            XFLUID::computeStabilization(derxy, gpvelnp, numparamvelx, instationary, visc, hk, mk, timefac,
                tau_stab_M, tau_stab_Mp, tau_stab_C);

            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            BlitzVec3 rhsint;
            BlitzVec3 bodyforce;
            bodyforce = 0.0;
            //bodyforce(0) = 1.0;
            for (int isd = 0; isd < nsd; ++isd)
                rhsint(isd) = histvec(isd) + bodyforce(isd)*timefac;

            /*----------------- get numerical representation of single operators ---*/

            /* Convective term  u_old * grad u_old: */
            BlitzVec3 conv_old;
            //conv_old = blitz::sum(vderxy(i, j)*gpvelnp(j), j);
            BLITZTINY::MV_product<3,3>(vderxy,gpvelnp,conv_old);
            
            /* Viscous term  div epsilon(u_old) */
            BlitzVec3 visc_old;
            visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
            visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
            visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));
            
            // evaluate residual once for all stabilisation right hand sides
            BlitzVec3 res_old;
            //res_old = -rhsint+timefac*(conv_old+gradp-2.0*visc*visc_old);
            for (int isd = 0; isd < nsd; ++isd)
                res_old(isd) = -rhsint(isd)+timefac*(conv_old(isd)+gradp(isd)-2.0*visc*visc_old(isd));  
            
            if (instationary)
                res_old += gpvelnp;
      
            /* Reactive term  u:  funct */
            /* linearise convective term */

            /*--- convective part u_old * grad (funct) --------------------------*/
            /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
             with  N .. form function matrix                                   */
            //const BlitzVec enr_conv_c_(blitz::sum(enr_derxy(j,i)*gpvelnp(j), j));
            static ShpVec enr_conv_c_;
            enr_conv_c_ = 0.0;
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
                enr_conv_c_(iparam) += shp.dx(iparam)*gpvelnp(0);
                enr_conv_c_(iparam) += shp.dy(iparam)*gpvelnp(1);
                enr_conv_c_(iparam) += shp.dz(iparam)*gpvelnp(2);
            }
            

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
            static ShpVec enr_viscs2_xx;
            static ShpVec enr_viscs2_xy;
            static ShpVec enr_viscs2_xz;
            static ShpVec enr_viscs2_yx;
            static ShpVec enr_viscs2_yy;
            static ShpVec enr_viscs2_yz;
            static ShpVec enr_viscs2_zx;
            static ShpVec enr_viscs2_zy;
            static ShpVec enr_viscs2_zz;

            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
              enr_viscs2_xx(iparam) = 0.5 * (2.0 * shp.dxdx(iparam) + shp.dydy(iparam) + shp.dzdz(iparam));
              enr_viscs2_xy(iparam) = 0.5 *  shp.dxdy(iparam);
              enr_viscs2_xz(iparam) = 0.5 *  shp.dxdz(iparam);
              enr_viscs2_yx(iparam) = 0.5 *  shp.dydx(iparam);
              enr_viscs2_yy(iparam) = 0.5 * (shp.dxdx(iparam) + 2.0 * shp.dydy(iparam) + shp.dzdz(iparam));
              enr_viscs2_yz(iparam) = 0.5 *  shp.dydz(iparam);
              enr_viscs2_zx(iparam) = 0.5 *  shp.dzdx(iparam);
              enr_viscs2_zy(iparam) = 0.5 *  shp.dzdy(iparam);
              enr_viscs2_zz(iparam) = 0.5 * (shp.dxdx(iparam) + shp.dydy(iparam) + 2.0 * shp.dzdz(iparam));
            }



            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

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
                assembler.template Matrix<Velx,Velx>(shp.d0, fac, shp.d0);
                assembler.template Matrix<Vely,Vely>(shp.d0, fac, shp.d0);
                assembler.template Matrix<Velz,Velz>(shp.d0, fac, shp.d0);
                
                assembler.template Vector<Velx>(shp.d0, -fac*gpvelnp(0));
                assembler.template Vector<Vely>(shp.d0, -fac*gpvelnp(1));
                assembler.template Vector<Velz>(shp.d0, -fac*gpvelnp(2));
            }
            
            /* convection, convective part */
            /*
                         /                       \
                        |      / n+1       \      |
                        | v , | u   o nabla | Du  |
                        |      \ (i)       /      |
                         \                       /
            */
            assembler.template Matrix<Velx,Velx>(shp.d0, timefacfac, enr_conv_c_);
            assembler.template Matrix<Vely,Vely>(shp.d0, timefacfac, enr_conv_c_);
            assembler.template Matrix<Velz,Velz>(shp.d0, timefacfac, enr_conv_c_);
            
            assembler.template Vector<Velx>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(0,0) // check order
                                                             +gpvelnp(1)*vderxy(0,1)
                                                             +gpvelnp(2)*vderxy(0,2)));
            assembler.template Vector<Vely>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(1,0)
                                                             +gpvelnp(1)*vderxy(1,1)
                                                             +gpvelnp(2)*vderxy(1,2)));
            assembler.template Vector<Velz>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(2,0)
                                                             +gpvelnp(1)*vderxy(2,1)
                                                             +gpvelnp(2)*vderxy(2,2)));
            
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
                assembler.template Matrix<Velx,Velx>(shp.d0, timefacfac*vderxy(0,0), shp.d0);
                assembler.template Matrix<Velx,Vely>(shp.d0, timefacfac*vderxy(0,1), shp.d0);
                assembler.template Matrix<Velx,Velz>(shp.d0, timefacfac*vderxy(0,2), shp.d0);
                assembler.template Matrix<Vely,Velx>(shp.d0, timefacfac*vderxy(1,0), shp.d0);
                assembler.template Matrix<Vely,Vely>(shp.d0, timefacfac*vderxy(1,1), shp.d0);
                assembler.template Matrix<Vely,Velz>(shp.d0, timefacfac*vderxy(1,2), shp.d0);
                assembler.template Matrix<Velz,Velx>(shp.d0, timefacfac*vderxy(2,0), shp.d0);
                assembler.template Matrix<Velz,Vely>(shp.d0, timefacfac*vderxy(2,1), shp.d0);
                assembler.template Matrix<Velz,Velz>(shp.d0, timefacfac*vderxy(2,2), shp.d0);
            }
            
            /* Viskositaetsterm */
            /*
                          /                        \
                         |       / \         /  \   |
                         |  eps | v | , tau | Du |  |
                         |       \ /         \  /   |
                          \                        /
            */
            assembler.template Matrix<Velx,Velx>(shp.dx,   2*visc*timefacfac, shp.dx);
            assembler.template Matrix<Velx,Velx>(shp.dy,     visc*timefacfac, shp.dy);
            assembler.template Matrix<Velx,Vely>(shp.dy,     visc*timefacfac, shp.dx);
            assembler.template Matrix<Velx,Velx>(shp.dz,     visc*timefacfac, shp.dz);
            assembler.template Matrix<Velx,Velz>(shp.dz,     visc*timefacfac, shp.dx);
            
            assembler.template Matrix<Vely,Vely>(shp.dx,     visc*timefacfac, shp.dx);
            assembler.template Matrix<Vely,Velx>(shp.dx,     visc*timefacfac, shp.dy);
            assembler.template Matrix<Vely,Vely>(shp.dy,   2*visc*timefacfac, shp.dy);
            assembler.template Matrix<Vely,Vely>(shp.dz,     visc*timefacfac, shp.dz);
            assembler.template Matrix<Vely,Velz>(shp.dz,     visc*timefacfac, shp.dy);
            
            assembler.template Matrix<Velz,Velz>(shp.dx,     visc*timefacfac, shp.dx);
            assembler.template Matrix<Velz,Velx>(shp.dx,     visc*timefacfac, shp.dz);
            assembler.template Matrix<Velz,Velz>(shp.dy,     visc*timefacfac, shp.dy);
            assembler.template Matrix<Velz,Vely>(shp.dy,     visc*timefacfac, shp.dz);
            assembler.template Matrix<Velz,Velz>(shp.dz,   2*visc*timefacfac, shp.dz);
            
            assembler.template Vector<Velx>(shp.dx,     -visc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
            assembler.template Vector<Velx>(shp.dy,     -visc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
            assembler.template Vector<Velx>(shp.dz,     -visc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
            
            assembler.template Vector<Vely>(shp.dx,     -visc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
            assembler.template Vector<Vely>(shp.dy,     -visc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
            assembler.template Vector<Vely>(shp.dz,     -visc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));
            
            assembler.template Vector<Velz>(shp.dx,     -visc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
            assembler.template Vector<Velz>(shp.dy,     -visc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
            assembler.template Vector<Velz>(shp.dz,     -visc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));
            
            /* Druckterm */
            /*
                            /                \
                           |                  |
                         - |  nabla o v , Dp  |
                           |                  |
                            \                /
            */
            assembler.template Matrix<Velx,Pres>(shp.dx, -timefacfac, shp.d0);
            assembler.template Matrix<Vely,Pres>(shp.dy, -timefacfac, shp.d0);
            assembler.template Matrix<Velz,Pres>(shp.dz, -timefacfac, shp.d0);
            
            assembler.template Vector<Velx>(shp.dx, timefacfac*pres);
            assembler.template Vector<Vely>(shp.dy, timefacfac*pres);
            assembler.template Vector<Velz>(shp.dz, timefacfac*pres);
            
            /* Divergenzfreiheit - continuity equation*/
            /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
            */
            assembler.template Matrix<Pres,Velx>(shp.d0, timefacfac, shp.dx);
            assembler.template Matrix<Pres,Vely>(shp.d0, timefacfac, shp.dy);
            assembler.template Matrix<Pres,Velz>(shp.d0, timefacfac, shp.dz);
            
            const double trace_gamma = (vderxy(0, 0) + vderxy(1, 1) + vderxy(2, 2));
            assembler.template Vector<Pres>(shp.d0, -timefacfac*trace_gamma);
            
            // source term of the right hand side
            assembler.template Vector<Velx>(shp.d0, fac*rhsint(0));
            assembler.template Vector<Vely>(shp.d0, fac*rhsint(1));
            assembler.template Vector<Velz>(shp.d0, fac*rhsint(2));
            
            
            // Hellinger-Reissner terms
            if (tauele_unknowns_present)
            {
                
                                 /*                     \
                              - |  virt tau , eps(Dtau)  |
                                 \                     */
              
                const double reciproke_viscfac = 1.0/(2.0*visc);
                assembler.template Matrix<Sigmaxx,Sigmaxx>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
                assembler.template Matrix<Sigmaxy,Sigmaxy>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
                assembler.template Matrix<Sigmaxz,Sigmaxz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
                assembler.template Matrix<Sigmayy,Sigmayy>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
                assembler.template Matrix<Sigmayz,Sigmayz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
                assembler.template Matrix<Sigmazz,Sigmazz>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
                
                assembler.template Vector<Sigmaxx>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,0));
                assembler.template Vector<Sigmaxy>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,1)*2.0);
                assembler.template Vector<Sigmaxz>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,2)*2.0);
                assembler.template Vector<Sigmayy>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,1));
                assembler.template Vector<Sigmayz>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,2)*2.0);
                assembler.template Vector<Sigmazz>(shp_tau,  reciproke_viscfac*timefacfac*tau(2,2));
                
                //             /                  |
                //            | virt tau , eps(Du) |
                //             |                  |
                
                assembler.template Matrix<Sigmaxx,Velx>(shp_tau,     timefacfac    , shp.dx);
                assembler.template Matrix<Sigmaxy,Velx>(shp_tau,     timefacfac    , shp.dy);
                assembler.template Matrix<Sigmaxy,Vely>(shp_tau,     timefacfac    , shp.dx);
                assembler.template Matrix<Sigmaxz,Velx>(shp_tau,     timefacfac    , shp.dz);
                assembler.template Matrix<Sigmaxz,Velz>(shp_tau,     timefacfac    , shp.dx);
                assembler.template Matrix<Sigmayy,Vely>(shp_tau,     timefacfac    , shp.dy);
                assembler.template Matrix<Sigmayz,Vely>(shp_tau,     timefacfac    , shp.dz);
                assembler.template Matrix<Sigmayz,Velz>(shp_tau,     timefacfac    , shp.dy);
                assembler.template Matrix<Sigmazz,Velz>(shp_tau,     timefacfac    , shp.dz);
                
                assembler.template Vector<Sigmaxx>(shp_tau,    - timefacfac*vderxy(0, 0));
                assembler.template Vector<Sigmaxy>(shp_tau,    - timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
                assembler.template Vector<Sigmaxz>(shp_tau,    - timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
                assembler.template Vector<Sigmayy>(shp_tau,    - timefacfac*vderxy(1, 1));
                assembler.template Vector<Sigmayz>(shp_tau,    - timefacfac*(vderxy(1, 2) + vderxy(1, 2)));
                assembler.template Vector<Sigmazz>(shp_tau,    - timefacfac*vderxy(2, 2));
                
                
                /* pressure-pressure coupling, rectangular part */
                /*
                               /                    \
                              |                      |
                            - | tr(virt tau^e) , p I |
                              |                      |
                               \                    /
                */
                assembler.template Matrix<Sigmaxx,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
                assembler.template Matrix<Sigmayy,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
                assembler.template Matrix<Sigmazz,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
                
                assembler.template Vector<Sigmaxx>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
                assembler.template Vector<Sigmayy>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
                assembler.template Vector<Sigmazz>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
                
            }
            
            //----------------------------------------------------------------------
            //                 PRESSURE STABILISATION PART
            if(pstab)
            {
                const double timetauMp  = timefac * tau_stab_Mp * fac;
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
                    assembler.template Matrix<Pres,Velx>(shp.dx, timetauMp, shp.d0);
                    assembler.template Matrix<Pres,Vely>(shp.dy, timetauMp, shp.d0);
                    assembler.template Matrix<Pres,Velz>(shp.dz, timetauMp, shp.d0);
                }
                const double ttimetauMp = timefac * timefac * tau_stab_Mp * fac;
                /* pressure stabilisation: convection, convective part */
                /*
                          /                             \
                         |             / n+1       \     |
                         | nabla q ,  | u   o nabla | Du |
                         |             \ i         /     |
                          \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp.dx, ttimetauMp, enr_conv_c_);
                assembler.template Matrix<Pres,Vely>(shp.dy, ttimetauMp, enr_conv_c_);
                assembler.template Matrix<Pres,Velz>(shp.dz, ttimetauMp, enr_conv_c_);
                
                if (newton)
                {
                    /*  pressure stabilisation: convection, reactive part
                          /                             \
                         |           /          \   n+1  |
                         | grad q , | Du o nabla | u     |
                         |           \          /   (i)  |
                          \                             /
                    */
                    assembler.template Matrix<Pres,Velx>(shp.dx, ttimetauMp*vderxy(0,0), shp.d0);
                    assembler.template Matrix<Pres,Velx>(shp.dy, ttimetauMp*vderxy(1,0), shp.d0);
                    assembler.template Matrix<Pres,Velx>(shp.dz, ttimetauMp*vderxy(2,0), shp.d0);
                    
                    assembler.template Matrix<Pres,Vely>(shp.dx, ttimetauMp*vderxy(0,1), shp.d0);
                    assembler.template Matrix<Pres,Vely>(shp.dy, ttimetauMp*vderxy(1,1), shp.d0);
                    assembler.template Matrix<Pres,Vely>(shp.dz, ttimetauMp*vderxy(2,1), shp.d0);
                    
                    assembler.template Matrix<Pres,Velz>(shp.dx, ttimetauMp*vderxy(0,2), shp.d0);
                    assembler.template Matrix<Pres,Velz>(shp.dy, ttimetauMp*vderxy(1,2), shp.d0);
                    assembler.template Matrix<Pres,Velz>(shp.dz, ttimetauMp*vderxy(2,2), shp.d0);
                }
                
                /* pressure stabilisation: viscosity (-L_visc_u) */
                /*
                           /                             \
                          |                         /  \  |
                        - |  nabla q , nabla o tau | Du | |
                          |                         \  /  |
                           \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2_xx);
                assembler.template Matrix<Pres,Vely>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2_xy);
                assembler.template Matrix<Pres,Velz>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2_xz);
                
                assembler.template Matrix<Pres,Velx>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2_xy);
                assembler.template Matrix<Pres,Vely>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2_yy);
                assembler.template Matrix<Pres,Velz>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2_yz);
                
                assembler.template Matrix<Pres,Velx>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2_xz);
                assembler.template Matrix<Pres,Vely>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2_yz);
                assembler.template Matrix<Pres,Velz>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2_zz);
                      
                /* pressure stabilisation: pressure( L_pres_p) */
                /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
                */
                assembler.template Matrix<Pres,Pres>(shp.dx, ttimetauMp, shp.dx);
                assembler.template Matrix<Pres,Pres>(shp.dy, ttimetauMp, shp.dy);
                assembler.template Matrix<Pres,Pres>(shp.dz, ttimetauMp, shp.dz);
                
                // pressure stabilization
                assembler.template Vector<Pres>(shp.dx, -timetauMp*res_old(0));
                assembler.template Vector<Pres>(shp.dy, -timetauMp*res_old(1));
                assembler.template Vector<Pres>(shp.dz, -timetauMp*res_old(2));
                
            }
            
            //----------------------------------------------------------------------
            //                     SUPG STABILISATION PART
            if(supg)
            {
                const double timetauM   = timefac * tau_stab_M * fac;
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
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, timetauM, shp.d0);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, timetauM, shp.d0);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, timetauM, shp.d0);
                    
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
                        assembler.template Matrix<Velx,Velx>(shp.dx, timetauM*gpvelnp(0), shp.d0);
                        assembler.template Matrix<Velx,Vely>(shp.dy, timetauM*gpvelnp(0), shp.d0);
                        assembler.template Matrix<Velx,Velz>(shp.dz, timetauM*gpvelnp(0), shp.d0);
                        assembler.template Matrix<Vely,Velx>(shp.dx, timetauM*gpvelnp(1), shp.d0);
                        assembler.template Matrix<Vely,Vely>(shp.dy, timetauM*gpvelnp(1), shp.d0);
                        assembler.template Matrix<Vely,Velz>(shp.dz, timetauM*gpvelnp(1), shp.d0);
                        assembler.template Matrix<Velz,Velx>(shp.dx, timetauM*gpvelnp(2), shp.d0);
                        assembler.template Matrix<Velz,Vely>(shp.dy, timetauM*gpvelnp(2), shp.d0);
                        assembler.template Matrix<Velz,Velz>(shp.dz, timetauM*gpvelnp(2), shp.d0);
                    }
                }
                const double ttimetauM  = timefac * timefac * tau_stab_M * fac;
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
                assembler.template Matrix<Velx,Pres>(enr_conv_c_, ttimetauM, shp.dx);
                assembler.template Matrix<Vely,Pres>(enr_conv_c_, ttimetauM, shp.dy);
                assembler.template Matrix<Velz,Pres>(enr_conv_c_, ttimetauM, shp.dz);

                /* supg stabilisation: viscous part  (-L_visc_u) */
                /*
                      /                                        \
                     |               /  \    / n+1        \     |
                   - |  nabla o eps | Du |, | u    o nabla | v  |
                     |               \  /    \ (i)        /     |
                      \                                        /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_xx);
                assembler.template Matrix<Velx,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_xy);
                assembler.template Matrix<Velx,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_xz);
    
                assembler.template Matrix<Vely,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_yx);
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_yy);
                assembler.template Matrix<Vely,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_yz);
    
                assembler.template Matrix<Velz,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_zx);
                assembler.template Matrix<Velz,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_zy);
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2_zz);
                
                if (newton)
                {
                    /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                    /*
                               /                                           |
                              |   |            |   n+1  | n+1          |   |
                              |   | Du o nabla | u    , | u    o nabla | v  |
                              |   |            |   (i)  | (i)          |   |
                               \                                           |
                    */
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*vderxy(0,0), shp.d0);
                    assembler.template Matrix<Velx,Vely>(enr_conv_c_, ttimetauM*vderxy(0,1), shp.d0);
                    assembler.template Matrix<Velx,Velz>(enr_conv_c_, ttimetauM*vderxy(0,2), shp.d0);
                    
                    assembler.template Matrix<Vely,Velx>(enr_conv_c_, ttimetauM*vderxy(1,0), shp.d0);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*vderxy(1,1), shp.d0);                    
                    assembler.template Matrix<Vely,Velz>(enr_conv_c_, ttimetauM*vderxy(1,2), shp.d0);                    
                    
                    assembler.template Matrix<Velz,Velx>(enr_conv_c_, ttimetauM*vderxy(2,0), shp.d0);
                    assembler.template Matrix<Velz,Vely>(enr_conv_c_, ttimetauM*vderxy(2,1), shp.d0);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*vderxy(2,2), shp.d0);
                    
                    /*
                             /                                            |
                            |    / n+1         |   n+1    /          |    |
                            |   | u    o nabla | u    , | Du o nabla | v  |
                            |   | (i)          |   (i)   |           |    |
                             \                                            |
                    */
                    const double con0 = ttimetauM*(gpvelnp(0)*vderxy(0,0) + gpvelnp(1)*vderxy(0,1) + gpvelnp(2)*vderxy(0,2));
                    assembler.template Matrix<Velx,Velx>(shp.dx, con0, shp.d0);
                    assembler.template Matrix<Velx,Vely>(shp.dy, con0, shp.d0);                        
                    assembler.template Matrix<Velx,Velz>(shp.dz, con0, shp.d0); 
                    
                    const double con1 = ttimetauM*(gpvelnp(0)*vderxy(1,0) + gpvelnp(1)*vderxy(1,1) + gpvelnp(2)*vderxy(1,2));
                    assembler.template Matrix<Vely,Velx>(shp.dx, con1, shp.d0);
                    assembler.template Matrix<Vely,Vely>(shp.dy, con1, shp.d0);
                    assembler.template Matrix<Vely,Velz>(shp.dz, con1, shp.d0);
                    
                    const double con2 = ttimetauM*(gpvelnp(0)*vderxy(2,0) + gpvelnp(1)*vderxy(2,1) + gpvelnp(2)*vderxy(2,2));
                    assembler.template Matrix<Velz,Velx>(shp.dx, con2, shp.d0);
                    assembler.template Matrix<Velz,Vely>(shp.dy, con2, shp.d0);
                    assembler.template Matrix<Velz,Velz>(shp.dz, con2, shp.d0);
                    
                    /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
                    /*
                                    /                               \
                                   |         n+1    /          \     |
                                   |  nabla p    , | Du o nabla | v  |
                                   |         (i)    \          /     |
                                    \                               /
                    */
                    assembler.template Matrix<Velx,Velx>(shp.dx, ttimetauM*gradp(0), shp.d0);
                    assembler.template Matrix<Velx,Vely>(shp.dy, ttimetauM*gradp(0), shp.d0);
                    assembler.template Matrix<Velx,Velz>(shp.dz, ttimetauM*gradp(0), shp.d0);

                    assembler.template Matrix<Vely,Velx>(shp.dx, ttimetauM*gradp(1), shp.d0);
                    assembler.template Matrix<Vely,Vely>(shp.dy, ttimetauM*gradp(1), shp.d0);
                    assembler.template Matrix<Vely,Velz>(shp.dz, ttimetauM*gradp(1), shp.d0);
                        
                    assembler.template Matrix<Velz,Velx>(shp.dx, ttimetauM*gradp(2), shp.d0);
                    assembler.template Matrix<Velz,Vely>(shp.dy, ttimetauM*gradp(2), shp.d0);
                    assembler.template Matrix<Velz,Velz>(shp.dz, ttimetauM*gradp(2), shp.d0);

                      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                      /*
                              /                                         \
                             |               / n+1 \    /          \     |
                           - |  nabla o eps | u     |, | Du o nabla | v  |
                             |               \ (i) /    \          /     |
                              \                                         /
                      */
                    assembler.template Matrix<Velx,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(0), shp.d0);
                    assembler.template Matrix<Velx,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(0), shp.d0);
                    assembler.template Matrix<Velx,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(0), shp.d0);

                    assembler.template Matrix<Vely,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(1), shp.d0);
                    assembler.template Matrix<Vely,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(1), shp.d0);
                    assembler.template Matrix<Vely,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(1), shp.d0);
                        
                    assembler.template Matrix<Velz,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(2), shp.d0);
                    assembler.template Matrix<Velz,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(2), shp.d0);
                    assembler.template Matrix<Velz,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(2), shp.d0);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                                  /                             \
                                 |              /          \     |
                               - |  rhsint   , | Du o nabla | v  |
                                 |              \          /     |
                                  \                             /

                    */
                    assembler.template Matrix<Velx,Velx>(shp.dx, -timetauM*rhsint(0), shp.d0);
                    assembler.template Matrix<Velx,Vely>(shp.dy, -timetauM*rhsint(0), shp.d0);
                    assembler.template Matrix<Velx,Velz>(shp.dz, -timetauM*rhsint(0), shp.d0);
                    
                    assembler.template Matrix<Vely,Velx>(shp.dx, -timetauM*rhsint(1), shp.d0);
                    assembler.template Matrix<Vely,Vely>(shp.dy, -timetauM*rhsint(1), shp.d0);
                    assembler.template Matrix<Vely,Velz>(shp.dz, -timetauM*rhsint(1), shp.d0);
                    
                    assembler.template Matrix<Velz,Velx>(shp.dx, -timetauM*rhsint(2), shp.d0);
                    assembler.template Matrix<Velz,Vely>(shp.dy, -timetauM*rhsint(2), shp.d0);
                    assembler.template Matrix<Velz,Velz>(shp.dz, -timetauM*rhsint(2), shp.d0);
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
                const double timefac_timefac_tau_C=timefac*timefac*tau_stab_C * fac;
                const double timefac_timefac_tau_C_divunp=timefac_timefac_tau_C*(vderxy(0, 0)+vderxy(1, 1)+vderxy(2, 2));
                /* continuity stabilisation on left hand side */
                /*
                         /                        \
                        |                          |
                        | nabla o Du  , nabla o v  |
                        |                          |
                         \                        /
                */
                assembler.template Matrix<Velx,Velx>(shp.dx, timefac_timefac_tau_C, shp.dx);
                assembler.template Matrix<Velx,Vely>(shp.dx, timefac_timefac_tau_C, shp.dy);
                assembler.template Matrix<Velx,Velz>(shp.dx, timefac_timefac_tau_C, shp.dz);
                
                assembler.template Matrix<Vely,Velx>(shp.dy, timefac_timefac_tau_C, shp.dx);
                assembler.template Matrix<Vely,Vely>(shp.dy, timefac_timefac_tau_C, shp.dy);
                assembler.template Matrix<Vely,Velz>(shp.dy, timefac_timefac_tau_C, shp.dz);
                
                assembler.template Matrix<Velz,Velx>(shp.dz, timefac_timefac_tau_C, shp.dx);
                assembler.template Matrix<Velz,Vely>(shp.dz, timefac_timefac_tau_C, shp.dy);
                assembler.template Matrix<Velz,Velz>(shp.dz, timefac_timefac_tau_C, shp.dz);
                
                assembler.template Vector<Velx>(shp.dx, -timefac_timefac_tau_C_divunp);
                assembler.template Vector<Vely>(shp.dy, -timefac_timefac_tau_C_divunp);
                assembler.template Vector<Velz>(shp.dz, -timefac_timefac_tau_C_divunp);
            } // endif cstab
        } // end loop over gauss points
    } // end loop over integration cells

    return;
}
  
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M1, class V1, class M2>
static void SysmatBoundary4(
    const DRT::Element*               ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandle>  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const M1&                         eveln,
    const M1&                         evelnm,
    const M1&                         eaccn,
    const V1&                         eprenp,
    const M2&                         etau,
    const Teuchos::RCP<const Epetra_Vector> ivelcol,       ///< velocity for interface nodes
    const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
    const struct _MATERIAL*           material,      ///< fluid material
    const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
    const double                      dt,            ///< delta t (time step size)
    const double                      theta,         ///< factor for one step theta scheme
    const bool                        newton,        ///< full Newton or fixed-point-like
    const bool                        pstab,         ///< flag for stabilization
    const bool                        supg,          ///< flag for stabilization
    const bool                        cstab,         ///< flag for stabilization
    const bool                        instationary,  ///< switch between stationary and instationary formulation
    const LocalAssembler<DISTYPE, ASSTYPE>& assembler,
    const bool                        ifaceForceContribution
)
{
    TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat4 - boundary");
  
    //if (false)
    if (ASSTYPE == XFEM::xfem_assembly)
    {
      const int nsd = 3;
      
      // time integration constant
      const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);
      
      // number of nodes for element
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      
      // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
      const int numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
      //const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
      // put one here to create arrays of size 1, since they are not needed anyway
      // in the xfem assembly, the numparam is determined by the dofmanager
      //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
      const int numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
      
      
      const bool tauele_unknowns_present = (XFEM::comgetNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
        // for now, I don't try to compare to elements without stress unknowns, since they lock anyway
        if (tauele_unknowns_present){
    
    // information about boundary integration cells
    const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());
    
    
    // loop over boundary integration cells
    for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
    {
      
        // gaussian points
        const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_37point);
        
        // get the right boundary element
        const DRT::Element* boundaryele = ih->GetBoundaryEle(cell->GetSurfaceEleGid());
        //cout << (*boundaryele) << endl;
        const int numnode_boundary = boundaryele->NumNode();
//        cout << "numnode_boundary: " << numnode_boundary << endl;
        
        // get current node coordinates
        const std::map<int,blitz::TinyVector<double,3> >* positions = ih->cutterposnp();
        const BlitzMat xyze_boundary(GEO::getCurrentNodalPositions(boundaryele, *positions));
        
        // get interface velocities at the boundary element nodes
        BlitzMat vel_boundary(3,numnode_boundary);
        const DRT::Node*const* nodes = boundaryele->Nodes();
        {
          for (int inode = 0; inode < numnode_boundary; ++inode)
          {
            const DRT::Node* node = nodes[inode];
            std::vector<int> lm = ih->cutterdis()->Dof(node);
            std::vector<double> myvel(3);
            DRT::UTILS::ExtractMyValues(*ivelcol,myvel,lm);
            vel_boundary(0,inode) = myvel[0];
            vel_boundary(1,inode) = myvel[1];
            vel_boundary(2,inode) = myvel[2];
          }
        }
        
        BlitzMat force_boundary(3,numnode_boundary);
        force_boundary = 0.0;

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta^\boundary
            GEO::PosEtaBoundary pos_eta_boundary;
            pos_eta_boundary(0) = intpoints.qxg[iquad][0];
            pos_eta_boundary(1) = intpoints.qxg[iquad][1];
//            cout << pos_eta_boundary << endl;
            
            // coordinates of the current integration point in element coordinates \xi
            GEO::PosXiBoundary posXiBoundary;
            mapEtaBToXiB(*cell, pos_eta_boundary, posXiBoundary);
            //cout << posXiBoundary << endl;
            
            GEO::PosXiDomain posXiDomain;
            mapEtaBToXiD(*cell, pos_eta_boundary, posXiDomain);
//            cout << cell->toString() << endl;
//            cout << posXiDomain << endl;
            const double detcell = fabs(detEtaBToXiB(*cell, pos_eta_boundary)); //TODO: check normals
            if (detcell <= 0.0)
            {
              cout << "detcel :  " << detcell << endl;
              dserror("negative detcell! should be a bug!");
            }

            // shape functions and their first derivatives
            BlitzVec funct_boundary(DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            DRT::UTILS::shape_function_2D(funct_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());
            BlitzMat deriv_boundary(3, DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());
            
            // shape functions and their first derivatives
            static BlitzVec funct(DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement);
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            // stress shape function
            const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
            static BlitzVec funct_stress(DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement);
            DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
            
            // position of the gausspoint in physical coordinates
//            gauss_pos_xyz = blitz::sum(funct_boundary(j)*xyze_boundary(i,j),j);
            BlitzVec3 gauss_pos_xyz;
            for (int isd = 0; isd < 3; ++isd)
            {
                gauss_pos_xyz(isd) = 0.0;
                for (int inode = 0; inode < numnode_boundary; ++inode)
                {
                    gauss_pos_xyz(isd) += funct_boundary(inode)*xyze_boundary(isd,inode);
                }
            }
      
            // get jacobian matrix d x / d \xi  (3x2)
            static BlitzMat3x2 dxyzdrs;
            //dxyzdrs = blitz::sum(xyze_boundary(i,k)*deriv_boundary(j,k),k);
            for (int isd = 0; isd < 3; ++isd)
            {
                for (int j = 0; j < 2; ++j)
                {
                    dxyzdrs(isd,j) = 0.0;
                    for (int k = 0; k < numnode_boundary; ++k)
                    {
                        dxyzdrs(isd,j) += xyze_boundary(isd,k)*deriv_boundary(j,k);
                    }
                }
            }
            
            // compute covariant metric tensor G for surface element (2x2)
            static BlitzMat2x2 metric;
            //metric = blitz::sum(dxyzdrs(k,i)*dxyzdrs(k,j),k);
            BLITZTINY::MtM_product<2,2,3>(dxyzdrs,dxyzdrs,metric);
            //const BlitzMat metric = computeMetricTensor(xyze_boundary,deriv_boundary);
            
            // compute global derivates
            //const BlitzMat derxy(blitz::sum(dxyzdrs(i,k)*deriv_boundary(k,j),k));
            //const BlitzMat derxy_stress(blitz::sum(xji(i,k)*deriv_stress(k,j),k));
            
            const double detmetric = sqrt(metric(0,0)*metric(1,1) - metric(0,1)*metric(1,0));
            if (detmetric <= 0.0)
            {
              dserror("negative detmetric! should be a bug!");
            }
            
            const double fac = intpoints.qwgt[iquad]*detmetric*detcell;
            if (fac <= 0.0)
            {
              dserror("negative fac! should be a bug!");
            }
            
            // temporary arrays
            BlitzVec enr_funct(numparamvelx);

            BlitzVec enr_funct_stress(numparamtauxx);
            
            if (dofman.getUniqueEnrichments().size() > 1)
              dserror("for an intersected element, we assume only 1 enrichment for now!");
            const std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
                  ih,
                  dofman.getUniqueEnrichments(),
                  gauss_pos_xyz,
                  XFEM::Enrichment::approachFromPlus));
            
            // shape function for nodal dofs
            XFEM::ComputeEnrichedNodalShapefunction(
                    *ele,
                    ih,
                    dofman,
                    Velx,
                    enrvals,
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
                    enrvals,
                    funct_stress,
                    enr_funct_stress);
                
            // perform integration for entire matrix and rhs
            const int shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
            const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
            typedef blitz::TinyVector<double,shpVecSize> ShpVec;
            static ShpVec shp;
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
              shp(iparam) = enr_funct(iparam);
            }
            static blitz::TinyVector<double,shpVecSizeStress> shp_tau;
            for (int iparam = 0; iparam < numparamtauxx; ++iparam)
            {
              shp_tau(iparam) = enr_funct_stress(iparam);
            }
            
            // get normal vector (in x coordinates) to surface element at integration point
            BlitzVec3 normalvec_solid;
            GEO::computeNormalToSurfaceElement(boundaryele, xyze_boundary, posXiBoundary, normalvec_solid);
//            cout << "normalvec " << normalvec << ", " << endl;
            BlitzVec3 normalvec_fluid;
            normalvec_fluid = -normalvec_solid;
//            cout << "normalvec : ";
//            cout << normalvec_fluid << endl;
      
            // get velocities (n+g,i) at integration point
//            gpvelnp = blitz::sum(evelnp(i,j)*shp(j),j);
            BlitzVec3 gpvelnp;
            //gpvelnp = blitz::sum(enr_funct(j)*evelnp(i,j),j);
            for (int isd = 0; isd < nsd; ++isd)
            {
                gpvelnp(isd) = 0.0;
                for (int iparam = 0; iparam < numparamvelx; ++iparam)
                    gpvelnp(isd) += evelnp(isd,iparam)*shp(iparam);
            }
            
            // get interface velocity
            BlitzVec3 interface_gpvelnp;
            if (timealgo == timeint_stationary)
            {
              interface_gpvelnp = 0.0;
            }
            else
            {
              for (int isd = 0; isd < 3; ++isd)
              {
                  interface_gpvelnp(isd) = 0.0;
                  for (int inode = 0; inode < numnode_boundary; ++inode)
                  {
                      interface_gpvelnp(isd) += vel_boundary(isd,inode)*funct_boundary(inode);
                  }
              }
            }

            // get viscous stress unknowns
            static BlitzMat3x3 tau;
            XFEM::fill_tau(numparamtauxx, shp_tau, etau, tau);
            
            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;
            
            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////
            
               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */
            
            const double taue_u_factor = 1.0;
            assembler.template Matrix<Sigmaxx,Velx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmaxy,Velx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmaxz,Velx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Sigmayx,Vely>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmayy,Vely>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmayz,Vely>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Sigmazx,Velz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmazy,Velz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmazz,Velz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2), shp);
            
            assembler.template Vector<Sigmaxx>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(0)*gpvelnp(0));
            assembler.template Vector<Sigmaxy>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(1)*gpvelnp(0));
            assembler.template Vector<Sigmaxz>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(2)*gpvelnp(0));
            assembler.template Vector<Sigmayx>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(0)*gpvelnp(1));
            assembler.template Vector<Sigmayy>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(1)*gpvelnp(1));
            assembler.template Vector<Sigmayz>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(2)*gpvelnp(1));
            assembler.template Vector<Sigmazx>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(0)*gpvelnp(2));
            assembler.template Vector<Sigmazy>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(1)*gpvelnp(2));
            assembler.template Vector<Sigmazz>(shp_tau, taue_u_factor*timefacfac*normalvec_fluid(2)*gpvelnp(2));
            
            
               /*                            \
              |  (virt tau) * n^f , u^\iface  |
               \                            */
            
            assembler.template Vector<Sigmaxx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0)*interface_gpvelnp(0));
            assembler.template Vector<Sigmaxy>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1)*interface_gpvelnp(0));
            assembler.template Vector<Sigmaxz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2)*interface_gpvelnp(0));
            assembler.template Vector<Sigmayx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0)*interface_gpvelnp(1));
            assembler.template Vector<Sigmayy>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1)*interface_gpvelnp(1));
            assembler.template Vector<Sigmayz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2)*interface_gpvelnp(1));
            assembler.template Vector<Sigmazx>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(0)*interface_gpvelnp(2));
            assembler.template Vector<Sigmazy>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(1)*interface_gpvelnp(2));
            assembler.template Vector<Sigmazz>(shp_tau, -taue_u_factor*timefacfac*normalvec_fluid(2)*interface_gpvelnp(2));
            
               /*               \
            - |  v , Dtau * n^f  |
               \               */

            const double vtaun_fac = 1.0;
            assembler.template Matrix<Velx,Sigmaxx>(shp, -vtaun_fac*timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velx,Sigmaxy>(shp, -vtaun_fac*timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velx,Sigmaxz>(shp, -vtaun_fac*timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Vely,Sigmayx>(shp, -vtaun_fac*timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Vely,Sigmayy>(shp, -vtaun_fac*timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Vely,Sigmayz>(shp, -vtaun_fac*timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Velz,Sigmazx>(shp, -vtaun_fac*timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velz,Sigmazy>(shp, -vtaun_fac*timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velz,Sigmazz>(shp, -vtaun_fac*timefacfac*normalvec_fluid(2), shp_tau);
            
            BlitzVec3 disctau_times_n;
            BLITZTINY::MV_product<3,3>(tau,normalvec_fluid,disctau_times_n);
            //cout << "sigmaijnj : " << disctau_times_n << endl;
            assembler.template Vector<Velx>(shp, vtaun_fac*timefacfac*disctau_times_n(0));
            assembler.template Vector<Vely>(shp, vtaun_fac*timefacfac*disctau_times_n(1));
            assembler.template Vector<Velz>(shp, vtaun_fac*timefacfac*disctau_times_n(2));
            
            // here the interface force is integrated
            // this is done using test shape functions of the boundary mesh
            // hence, we can't use the local assembler here
            for (int inode = 0; inode < numnode_boundary; ++inode)
            {
              force_boundary(0,inode) += funct_boundary(inode) * (disctau_times_n(0) * timefacfac);
              force_boundary(1,inode) += funct_boundary(inode) * (disctau_times_n(1) * timefacfac);
              force_boundary(2,inode) += funct_boundary(inode) * (disctau_times_n(2) * timefacfac);
            }
            
        } // end loop over gauss points
        
        // here we need to assemble into the global force vector of the boundary discretization
        // note that we assemble into a overlapping vector, hence we later have to figure out,
        // how the values get into the right places of the force vector in unique distribution 
        if (ifaceForceContribution)
        {
          const Epetra_Map* dofcolmap = ih->cutterdis()->DofColMap();
          for (int inode = 0; inode < numnode_boundary; ++inode)
          {
            const DRT::Node* node = nodes[inode];
            const std::vector<int> gdofs(ih->cutterdis()->Dof(node));
            (*iforcecol)[dofcolmap->LID(gdofs[0])] += force_boundary(0,inode);
            (*iforcecol)[dofcolmap->LID(gdofs[1])] += force_boundary(1,inode);
            (*iforcecol)[dofcolmap->LID(gdofs[2])] += force_boundary(2,inode);
          }
        }
        
        
      } // end loop over boundary integration cells
    }
    } // if (ASSTYPE == XFEM::xfem_assembly)
    return;
}
  
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
static void Sysmat4(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandle>  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const DRT::ELEMENTS::Combust3::MyState  mystate,  ///< element state variables
        const Teuchos::RCP<const Epetra_Vector> ivelcol,       ///< velocity for interface nodes
        const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,         ///< element matrix to calculate
        Epetra_SerialDenseVector&         eforce,        ///< element rhs to calculate
        const struct _MATERIAL*           material,      ///< fluid material
        const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
        const double                      dt,            ///< delta t (time step size)
        const double                      theta,         ///< factor for one step theta scheme
        const bool                        newton,        ///< full Newton or fixed-point-like
        const bool                        pstab,         ///< flag for stabilisation
        const bool                        supg,          ///< flag for stabilisation
        const bool                        cstab,         ///< flag for stabilisation
        const bool                        instationary,  ///< switch between stationary and instationary formulation
        const bool                        ifaceForceContribution
        )
{
    // initialize arrays
    estif.Scale(0.0);
    eforce.Scale(0.0);
    
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    
    // dead load in element nodes
    //////////////////////////////////////////////////// , BlitzMat edeadng_(BodyForce(ele->Nodes(),time));

    const LocalAssembler<DISTYPE, ASSTYPE> assembler(dofman, estif, eforce);
    
    // split velocity and pressure (and stress)
    const int shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    blitz::TinyVector<double,shpVecSize> eprenp;
    blitz::TinyMatrix<double,3,shpVecSize> evelnp;
    blitz::TinyMatrix<double,3,shpVecSize> eveln;
    blitz::TinyMatrix<double,3,shpVecSize> evelnm;
    blitz::TinyMatrix<double,3,shpVecSize> eaccn;
    blitz::TinyMatrix<double,6,shpVecSizeStress> etau;
    
    fillElementUnknownsArrays4<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp, etau);
    
    SysmatDomain4<DISTYPE,ASSTYPE>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, etau, ivelcol, iforcecol,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler);
    
    if (ASSTYPE == XFEM::xfem_assembly)
    {
      SysmatBoundary4<DISTYPE,ASSTYPE>(
          ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, etau, ivelcol, iforcecol,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler, ifaceForceContribution);
    }
}




/*!
 * \brief entry point for Sysmat call

 */
void COMBUST::callSysmat4(
        const XFEM::AssemblyType          assembly_type,
        const DRT::ELEMENTS::Combust3*     ele,
        const Teuchos::RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const DRT::ELEMENTS::Combust3::MyState  mystate,   ///< element state variables
        const Teuchos::RCP<const Epetra_Vector> ivelcol,
        const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,
        Epetra_SerialDenseVector&         eforce,
        const struct _MATERIAL*           material,
        const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
        const double                      dt,            ///< delta t (time step size)
        const double                      theta,         ///< factor for one step theta scheme
        const bool                        newton ,
        const bool                        pstab  ,
        const bool                        supg   ,
        const bool                        cstab  ,
        const bool                        instationary,
        const bool                        ifaceForceContribution
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat4<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::hex27:
//                Sysmat4<DRT::Element::hex27,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary);
//                break;
            case DRT::Element::tet4:
                Sysmat4<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::tet10:
//                Sysmat4<DRT::Element::tet4,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary);
//                break;
            default:
                dserror("Sysmat not templated yet");
        };
    }
    else
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat4<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::hex27:
//                Sysmat4<DRT::Element::hex27,XFEM::xfem_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary);
//                break;
//            case DRT::Element::tet10:
//                Sysmat4<DRT::Element::tet4,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary);
//                break;
            default:
                dserror("Sysmat not templated yet");
        };
    }
}
//} // end namespace COMBUST

#endif

#endif
