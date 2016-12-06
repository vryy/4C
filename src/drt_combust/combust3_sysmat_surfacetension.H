/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_surfacetension.H

\brief system matrix formulation for surface tension boundary integrals / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_SURFACETENSION_H
#define COMBUST3_SYSMAT_SURFACETENSION_H


#include "combust3_sysmat.H"
#include "combust3_utils.H"
#include "combust_time_integration_element.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_normals.H"

#ifdef DEBUG
#include "../drt_mat/material.H"
#include "../drt_inpar/inpar_material.H"
#endif

namespace COMBUST
{


/*!
  Calculate surface tension 1D-boundary integrals for a line segment of a 2D-boundary integration cell
  for kink erichments in velocity
 */
template <DRT::Element::DiscretizationType DISTYPE,
DRT::Element::DiscretizationType CELLDISTYPE,
XFEM::AssemblyType ASSTYPE,
int NUMDOF,
class M1, class M2, class M3, class V1>
void RHS_LaplaceBeltrami_1DLineIntegrals_TPF(
    const DRT::ELEMENTS::Combust3*             ele,           ///< the element those matrix is calculated
    const XFEM::ElementDofManager&             dofman,
    const V1&                                  ephi,
    M1&                                        cellXiDomain,
    M2&                                        xyze,          ///< xyz coordinates of element ele
    const M3&                                  grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const size_t&                              numparamvelx,
    const double&                              timefac,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
)
{
  // additional 1D surface tension integrals get assembled

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  //--------------------------------------------------------------------------------------------
  // get all line segments of current 2D boundary integration cell
  //--------------------------------------------------------------------------------------------

  // shape functions in 1D => only lines are assumed
  const size_t numvertices_line = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::line2>::numNodePerElement;

  // define vector of line segments
  // each line is represented by a (nsd x 2)-matrix
  std::vector< LINALG::Matrix<nsd,numvertices_line> > lineVec;

  //--------------------------------------------------------------------------------------------
  // set elements of lineVec
  //--------------------------------------------------------------------------------------------
  // loop over nodes of boundary cell
  for (size_t inode=0; inode<numvertices; inode++)
  {
    LINALG::Matrix<nsd,numvertices_line> linetmp;
    linetmp.Clear();

    for (size_t isd = 0; isd<nsd; isd++){
      // second point
      linetmp(isd,1) = cellXiDomain(isd,(inode+1)%numvertices);
      // first point
      linetmp(isd,0) = cellXiDomain(isd,(inode)%numvertices);
    }
    lineVec.push_back(linetmp);
  }



  //============================   loop over boundary lines ==================================

  // define vector iterator
  typedef std::vector<LINALG::Matrix<nsd,numvertices_line> > LinesVec;
  typedef LinesVec::iterator LinesVec_Iterator;

  for (LinesVec_Iterator line_it = lineVec.begin(); line_it != lineVec.end(); line_it++ )
  {
    // check if line lies in boundary of the domain (bound(Omega))
    // REMARK:
    // we integrate all boundary lines because isbound-check is not so easy
    // interior lines are integrated twice but with contrarious directions -> cancel each other out

    bool isbound =true;

    if (isbound == true) // line is a boundary line => integrate boundary integral for surface tension at Interface-boundaries
    {
      //============================ INTEGRATE BOUNDARY TERM ==============================
      //        std::cout << "boundary line is found:\n" << *line_it;

      // get integration points
      // TODO get the right intrule
      const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);


      LINALG::Matrix<nsd,numvertices_line> currLine_Xi(*line_it);
      //----------------------------------------------------------------------------------------------
      // integration loop over Gaussian points
      //----------------------------------------------------------------------------------------------
      for (int iquad=0; iquad<intpoint_line.nquad; ++iquad)
      {

        // 1D-integration points
        const LINALG::Matrix<1,1> posGammaBound(intpoint_line.qxg[iquad]);
        const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);

        // jacobian for coupled transformation
        // get derivatives dxi_3D/dline_1D
        static LINALG::Matrix<1,numvertices_line> deriv_gamma;
        DRT::UTILS::shape_function_1D_deriv1(deriv_gamma,posGammaBound(0,0),DRT::Element::line2);

        // calculate dxi3Ddline_1D
        static LINALG::Matrix<3,1> dXi3Ddline1D;
        dXi3Ddline1D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 1; j++) // derivatives
            for (int k = 0; k < (int)numvertices_line; k++)
              dXi3Ddline1D(i,j) += currLine_Xi(i,k)*deriv_gamma(j,k);

        // transform Gauss point to xi3D space (element parameter space)
        static LINALG::Matrix<nsd,1> gpinXi3D;
        gpinXi3D.Clear();

        static LINALG::Matrix<numvertices_line,1> shape_gamma;
        DRT::UTILS::shape_function_1D(shape_gamma,posGammaBound(0,0),DRT::Element::line2);

        gpinXi3D.Multiply(currLine_Xi,shape_gamma);


        static LINALG::Matrix<3,numnode> deriv_xi3D;
        DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);


        // calculate dx3Ddxi3D
        static LINALG::Matrix<3,3> dX3DdXi3D;
        dX3DdXi3D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 3; j++) // derivatives
            for (int k = 0; k < (int)numnode; k++)
              dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

        // get the coupled Jacobian dx3Ddeta2D
        static LINALG::Matrix<3,1> dx3Ddline1D;
        dx3Ddline1D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 1; j++) // derivatives
            for (int k = 0; k < 3; k++)
              dx3Ddline1D(i,j) += dX3DdXi3D(i,k) * dXi3Ddline1D(k,j);

        // get deformation factor
        static LINALG::Matrix<1,1> Jac_tmp; // J^T*J
        Jac_tmp.Clear();
        Jac_tmp.MultiplyTN(dx3Ddline1D,dx3Ddline1D);

        if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

        const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

        const double fac = intpoint_line.qwgt[iquad]*deform_factor;


        // compute total (time and spatial) integration factor (and coefficients of single terms)?
        const double timefacfac_gamma = timefac * fac;



        //--------------------------------------------------------------------------------------------
        // evaluate element shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<numnode,1> funct_gp;
        funct_gp.Clear();
        DRT::UTILS::shape_function_3D(funct_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);

        static LINALG::Matrix<nsd,numnode> deriv_gp;
        deriv_gp.Clear();
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);

        //--------------------------------------------------------------------------------------------
        // procedures involving Jacobian matrix for domain mapping
        //--------------------------------------------------------------------------------------------
        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.Clear();
        xjm_gp.MultiplyNT(deriv_gp,xyze);
        // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
        //        const double detXtoXi = xjm_gp.Determinant();
        //        //#ifdef DEBUG
        //        if (detXtoXi < 0.0)
        //        {
        //          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
        //        }
        //#endif
        // inverse of jacobian
        static LINALG::Matrix<nsd,nsd> xji_gp;
        xji_gp.Clear();
        xji_gp.Invert(xjm_gp);

        //--------------------------------------------------------------------------------------------
        // compute global derivates of shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        // compute first global derivative
        static LINALG::Matrix<3,numnode> derxy_gp;
        derxy_gp.Clear();
        // derxy(i,j) = xji(i,k) * deriv(k,j)
        derxy_gp.Multiply(xji_gp,deriv_gp);

        // compute second global derivative
        static LINALG::Matrix<6,numnode> derxy2_gp;
        if (false) // if (higher_order_ele)
        {
          static LINALG::Matrix<6,numnode> deriv2_gp;
          DRT::UTILS::shape_function_3D_deriv2(deriv2_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);
          DRT::UTILS::gder2<DISTYPE,numnode>(xjm_gp, derxy_gp, deriv2_gp, xyze, derxy2_gp);
        }
        else
        {
          derxy2_gp.Clear();
        }

        //--------------------------------------------------------------------------------------------
        // rearrange (enriched) shape functions as approximation functions
        //--------------------------------------------------------------------------------------------
        const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

        // see remark: SysmatBoundarySurfaceTensionIntCell
        const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi,funct_gp,derxy_gp,derxy2_gp);


        // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
        static LINALG::Matrix<shpVecSize,1>       enrfunct;
        enrfunct.Clear();


        // shape functions for nodal parameters (dofs) on plus and minus side
        enrvals.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct_gp, enrfunct);

        // (enriched) shape functions = approximation functions (P = N * \Psi)
        static XFEM::ApproxFunc<0,shpVecSize> shp_mean_inv;


        // fill approximation functions
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          shp_mean_inv.d0(iparam) = enrfunct(iparam);
        }


        //================= get normal to gamma-line  =======================================
        // normal vector has to lie in the tangential space to the interface
        // normal has to lie orthogonal to gamma-line
        // normal has to be normalized
        // normal must have the right orientation

        // get tangential space to the interface at current gaussian point

        //------------------------------------------------------------------------------------------
        // get position of line in XYZ system

        LINALG::Matrix<numnode,1> shape_tmp;
        shape_tmp.Clear();

        // get position of line in XYZ system

        // first point of line
        LINALG::Matrix<nsd,1> firstpoint;
        firstpoint.Clear();

        DRT::UTILS::shape_function_3D(shape_tmp,(*line_it)(0,0), (*line_it)(1,0), (*line_it)(2,0), DISTYPE);
        firstpoint.Multiply(xyze,shape_tmp);

        // second point of line
        LINALG::Matrix<nsd,1> secondpoint;
        secondpoint.Clear();

        DRT::UTILS::shape_function_3D(shape_tmp,(*line_it)(0,1), (*line_it)(1,1), (*line_it)(2,1), DISTYPE);
        secondpoint.Multiply(xyze,shape_tmp);

        // define line in XYZ-system
        LINALG::Matrix<nsd,numvertices_line> currLine_XYZ;
        currLine_XYZ.Clear();
        for(int isd=0; isd<3; isd++)
        {
          currLine_XYZ(isd,0) = firstpoint(isd);
          currLine_XYZ(isd,1) = secondpoint(isd);
        }
        //------------------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------------------
        // 1. define vector for line-direction line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> line;
        line.Clear();
        line(0) = currLine_XYZ(0,1) - currLine_XYZ(0,0);
        line(1) = currLine_XYZ(1,1) - currLine_XYZ(1,0);
        line(2) = currLine_XYZ(2,1) - currLine_XYZ(2,0);


        //--------------------------------------------------------------------------------------------
        // 2. get normal_Phi = grad(phi)/| grad_phi | at gp
        //--------------------------------------------------------------------------------------------
        LINALG::Matrix<nsd,1> normal_Phi;
        normal_Phi.Clear();

        for(size_t i=0; i<numnode; i++)
        {

          // interpolate nodal phi gradients
          normal_Phi(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
          normal_Phi(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
          normal_Phi(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);

        }

        // normalize final gradient
        normal_Phi.Scale(-1.0/normal_Phi.Norm2());

        //--------------------------------------------------------------------------------------------
        // 3. get normal_gamma: <n_gamma,line > = 0, <n_gamma,normal_Phi> = 0, orientation!!!
        // set n_gamma pointing outward from domain n_gamma = normal_Phi x line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> normal_gamma;
        normal_gamma.Clear();
        normal_gamma(0) = normal_Phi(1)*line(2) - normal_Phi(2)*line(1);
        normal_gamma(1) = normal_Phi(2)*line(0) - normal_Phi(0)*line(2);
        normal_gamma(2) = normal_Phi(0)*line(1) - normal_Phi(1)*line(0);

        normal_gamma.Scale(1.0/normal_gamma.Norm2());


        // =============================  assemble vectors ==================
        //----------------------------------------     |                                 |
        // surface tension term boundary of interface  |  < v > , surftens*n_gamma >      |
        //----------------------------------------     |                                 |

        assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(0));
        assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(1));
        assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(2));


      } // ==================================== end loop gaussian points ==========================


    } // ================================= end line is boundary / end of assembly ============


  } //=========================== end loop over lines==================================
  return;
}


/*!
 Calculate surface tension 1D-boundary integrals for a line segment of a 2D-boundary integration cell
 for jump enrichments in velocity!!!
 */
template <DRT::Element::DiscretizationType DISTYPE,
DRT::Element::DiscretizationType CELLDISTYPE,
XFEM::AssemblyType ASSTYPE,
int NUMDOF,
class M1, class M2, class M3>
void RHS_LaplaceBeltrami_1DLineIntegrals_Nitsche(
    M1&                                        cellXiDomain,
    M2&                                        xyze,          ///< xyz coordinates of element ele
    const M3&                                  grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const XFEM::ElementEnrichmentValues&       enrvals_plus,
    const XFEM::ElementEnrichmentValues&       enrvals_minus,
    const double&                              kappa_plus,
    const double&                              kappa_minus,
    const size_t&                              numparamvelx,
    const double&                              timefac,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
)
{
  // additional 1D surface tension integrals get assembled

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  //--------------------------------------------------------------------------------------------
  // get all line segments of current 2D boundary integration cell
  //--------------------------------------------------------------------------------------------

  // shape functions in 1D => only lines are assumed
  const size_t numvertices_line = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::line2>::numNodePerElement;

  // define vector of line segments
  // each line is represented by a (nsd x 2)-matrix
  std::vector< LINALG::Matrix<nsd,numvertices_line> > lineVec;

  //--------------------------------------------------------------------------------------------
  // set elements of lineVec
  //--------------------------------------------------------------------------------------------
  // loop over nodes of boundary cell
  for (size_t inode=0; inode<numvertices; inode++)
  {
    LINALG::Matrix<nsd,numvertices_line> linetmp;
    linetmp.Clear();

    for (size_t isd = 0; isd<nsd; isd++){
      // second point
      linetmp(isd,1) = cellXiDomain(isd,(inode+1)%numvertices);
      // first point
      linetmp(isd,0) = cellXiDomain(isd,(inode)%numvertices);
    }
    lineVec.push_back(linetmp);
  }



  //============================   loop over boundary lines ==================================

  // define vector iterator
  typedef std::vector<LINALG::Matrix<nsd,numvertices_line> > LinesVec;
  typedef LinesVec::iterator LinesVec_Iterator;

  for (LinesVec_Iterator line_it = lineVec.begin(); line_it != lineVec.end(); line_it++ )
  {
    // check if line lies in boundary of the domain (bound(Omega))
    // REMARK:
    // we integrate all boundary lines because isbound-check is not so easy
    // interior lines are integrated twice but with contrarious directions -> cancel each other out

    bool isbound =true;

    if (isbound == true) // line is a boundary line => integrate boundary integral for surface tension at Interface-boundaries
    {
      //============================ INTEGRATE BOUNDARY TERM ==============================
      //        std::cout << "boundary line is found:\n" << *line_it;

      // get integration points
      // TODO get the right intrule
      const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);


      LINALG::Matrix<nsd,numvertices_line> currLine_Xi(*line_it);
      //----------------------------------------------------------------------------------------------
      // integration loop over Gaussian points
      //----------------------------------------------------------------------------------------------
      for (int iquad=0; iquad<intpoint_line.nquad; ++iquad)
      {

        // 1D-integration points
        const LINALG::Matrix<1,1> posGammaBound(intpoint_line.qxg[iquad]);
        const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);

        // jacobian for coupled transformation
        // get derivatives dxi_3D/dline_1D
        static LINALG::Matrix<1,numvertices_line> deriv_gamma;
        DRT::UTILS::shape_function_1D_deriv1(deriv_gamma,posGammaBound(0,0),DRT::Element::line2);

        // calculate dxi3Ddline_1D
        static LINALG::Matrix<3,1> dXi3Ddline1D;
        dXi3Ddline1D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 1; j++) // derivatives
            for (int k = 0; k < (int)numvertices_line; k++)
              dXi3Ddline1D(i,j) += currLine_Xi(i,k)*deriv_gamma(j,k);

        // transform Gauss point to xi3D space (element parameter space)
        static LINALG::Matrix<nsd,1> gpinXi3D;
        gpinXi3D.Clear();

        static LINALG::Matrix<numvertices_line,1> shape_gamma;
        DRT::UTILS::shape_function_1D(shape_gamma,posGammaBound(0,0),DRT::Element::line2);

        gpinXi3D.Multiply(currLine_Xi,shape_gamma);


        static LINALG::Matrix<3,numnode> deriv_xi3D;
        DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);


        // calculate dx3Ddxi3D
        static LINALG::Matrix<3,3> dX3DdXi3D;
        dX3DdXi3D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 3; j++) // derivatives
            for (int k = 0; k < (int)numnode; k++)
              dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

        // get the coupled Jacobian dx3Ddeta2D
        static LINALG::Matrix<3,1> dx3Ddline1D;
        dx3Ddline1D.Clear();
        for (int i = 0; i < 3; i++)   // dimensions
          for (int j = 0; j < 1; j++) // derivatives
            for (int k = 0; k < 3; k++)
              dx3Ddline1D(i,j) += dX3DdXi3D(i,k) * dXi3Ddline1D(k,j);

        // get deformation factor
        static LINALG::Matrix<1,1> Jac_tmp; // J^T*J
        Jac_tmp.Clear();
        Jac_tmp.MultiplyTN(dx3Ddline1D,dx3Ddline1D);

        if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

        const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

        const double fac = intpoint_line.qwgt[iquad]*deform_factor;


        // compute total (time and spatial) integration factor (and coefficients of single terms)?
        const double timefacfac_gamma = timefac * fac;



        //--------------------------------------------------------------------------------------------
        // evaluate element shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<numnode,1> funct_gp;
        funct_gp.Clear();
        DRT::UTILS::shape_function_3D(funct_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);

        static LINALG::Matrix<nsd,numnode> deriv_gp;
        deriv_gp.Clear();
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);

        //--------------------------------------------------------------------------------------------
        // procedures involving Jacobian matrix for domain mapping
        //--------------------------------------------------------------------------------------------
        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.Clear();
        xjm_gp.MultiplyNT(deriv_gp,xyze);
        // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
        //        const double detXtoXi = xjm_gp.Determinant();
        //        //#ifdef DEBUG
        //        if (detXtoXi < 0.0)
        //        {
        //          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
        //        }
        //#endif
        // inverse of jacobian
        static LINALG::Matrix<nsd,nsd> xji_gp;
        xji_gp.Clear();
        xji_gp.Invert(xjm_gp);

        //--------------------------------------------------------------------------------------------
        // compute global derivates of shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        // compute first global derivative
        static LINALG::Matrix<3,numnode> derxy_gp;
        derxy_gp.Clear();
        // derxy(i,j) = xji(i,k) * deriv(k,j)
        derxy_gp.Multiply(xji_gp,deriv_gp);

        // compute second global derivative
        static LINALG::Matrix<6,numnode> derxy2_gp;
        if (false) // if (higher_order_ele)
        {
          static LINALG::Matrix<6,numnode> deriv2_gp;
          DRT::UTILS::shape_function_3D_deriv2(deriv2_gp,gpinXi3D(0),gpinXi3D(1),gpinXi3D(2),DISTYPE);
          DRT::UTILS::gder2<DISTYPE,numnode>(xjm_gp, derxy_gp, deriv2_gp, xyze, derxy2_gp);
        }
        else
        {
          derxy2_gp.Clear();
        }



        //--------------------------------------------------------------------------------------------
        // rearrange (enriched) shape functions as approximation functions
        //--------------------------------------------------------------------------------------------
        const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

        // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
        static LINALG::Matrix<shpVecSize,1>       enrfunct_plus;
        enrfunct_plus.Clear();
        static LINALG::Matrix<shpVecSize,1>       enrfunct_minus;
        enrfunct_minus.Clear();


        // shape functions for nodal parameters (dofs) on plus and minus side
        enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct_gp, enrfunct_plus);
        enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct_gp, enrfunct_minus);

        // (enriched) shape functions = approximation functions (P = N * \Psi)
        static XFEM::ApproxFunc<0,shpVecSize> shp_mean_inv;


        // fill approximation functions
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          // <>-average
          shp_mean_inv.d0(iparam) = kappa_minus * enrfunct_plus(iparam) + kappa_plus * enrfunct_minus(iparam);
        }

        //================= get normal to gamma-line  =======================================
        // normal vector has to lie in the tangential space to the interface
        // normal has to lie orthogonal to gamma-line
        // normal has to be normalized
        // normal must have the right orientation

        // get tangential space to the interface at current gaussian point

        //------------------------------------------------------------------------------------------
        // get position of line in XYZ system

        LINALG::Matrix<numnode,1> shape_tmp;
        shape_tmp.Clear();

        // get position of line in XYZ system

        // first point of line
        LINALG::Matrix<nsd,1> firstpoint;
        firstpoint.Clear();

        DRT::UTILS::shape_function_3D(shape_tmp,(*line_it)(0,0), (*line_it)(1,0), (*line_it)(2,0), DISTYPE);
        firstpoint.Multiply(xyze,shape_tmp);

        // second point of line
        LINALG::Matrix<nsd,1> secondpoint;
        secondpoint.Clear();

        DRT::UTILS::shape_function_3D(shape_tmp,(*line_it)(0,1), (*line_it)(1,1), (*line_it)(2,1), DISTYPE);
        secondpoint.Multiply(xyze,shape_tmp);

        // define line in XYZ-system
        LINALG::Matrix<nsd,numvertices_line> currLine_XYZ;
        currLine_XYZ.Clear();
        for(int isd=0; isd<3; isd++)
        {
          currLine_XYZ(isd,0) = firstpoint(isd);
          currLine_XYZ(isd,1) = secondpoint(isd);
        }
        //------------------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------------------
        // 1. define vector for line-direction line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> line;
        line.Clear();
        line(0) = currLine_XYZ(0,1) - currLine_XYZ(0,0);
        line(1) = currLine_XYZ(1,1) - currLine_XYZ(1,0);
        line(2) = currLine_XYZ(2,1) - currLine_XYZ(2,0);


        //--------------------------------------------------------------------------------------------
        // 2. get normal_Phi = grad(phi)/| grad_phi | at gp
        //--------------------------------------------------------------------------------------------
        LINALG::Matrix<nsd,1> normal_Phi;
        normal_Phi.Clear();

        for(size_t i=0; i<numnode; i++)
        {

          // interpolate nodal phi gradients
          normal_Phi(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
          normal_Phi(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
          normal_Phi(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);

        }

        // normalize final gradient
        // remark: n = -grad(\phi)
        normal_Phi.Scale(-1.0/normal_Phi.Norm2());


        //--------------------------------------------------------------------------------------------
        // 3. get normal_gamma: <n_gamma,line > = 0, <n_gamma,normal_Phi> = 0, orientation!!!
        // set n_gamma pointing outward from domain n_gamma = normal_Phi x line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> normal_gamma;
        normal_gamma.Clear();
        normal_gamma(0) = normal_Phi(1)*line(2) - normal_Phi(2)*line(1);
        normal_gamma(1) = normal_Phi(2)*line(0) - normal_Phi(0)*line(2);
        normal_gamma(2) = normal_Phi(0)*line(1) - normal_Phi(1)*line(0);

        normal_gamma.Scale(1.0/normal_gamma.Norm2());

        // =============================  assemble vectors ==================
        //----------------------------------------     |                                 |
        // surface tension term boundary of interface  |  < v > , surftens*n_gamma >     |
        //----------------------------------------     |                                 |

        assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(0));
        assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(1));
        assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, - timefacfac_gamma*surftenscoeff*normal_gamma(2));


      } // ==================================== end loop gaussian points ==========================


    } // ================================= end line is boundary / end of assembly ============


  } //=========================== end loop over lines==================================
  return;
}


/*!
 Calculate surface tension 2D-boundary integrals
 */
template <DRT::Element::DiscretizationType DISTYPE,
XFEM::AssemblyType ASSTYPE,
int NUMDOF,
class M1, class M2, class M3, class M4, class M5, class M6, class V1>
void RHS_SurfaceTension2DIntegrals(
    M1&                                        posXiDomain, // position of gaussian point for evaluating of surface tension terms
    M5&                                        shp_mean_inv,
    const DRT::ELEMENTS::Combust3*             ele,           ///< the element those matrix is calculated
    M2&                                        xyze,          ///< xyz coordinates of element ele
    const M3&                                  grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const M6&                                  grad_ephi2_smoothed, ///<smoothed nodal G-function gradient values for this element
    const V1&                                  ecurv,              ///<curvature for this element?
    M4&                                        normal,
    const double&                              timefacfac,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               variablesurftens,
    const bool                                 second_deriv,
    const double                               surftenscoeff,   ///<,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
)
{
  // get derivative of surface tension coefficient
  // for the time being only linear in x available
  double surftenscoeffx=variablesurftens;
  double surftenscoeffy=0.0;
  double surftenscoeffz=0.0;
  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  if(surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami)
  {

    //-----------------------------------------------------------     |                                           |
    // surface tension (Laplace-Beltrami surface tension) term      - |  < nabla_BEL v > ° nabla_BEL id_gamma     |
    //-----------------------------------------------------------     |                                           |
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(1.0 - normal(0)*normal(0)));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(0)));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(0)));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(1)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(1.0 - normal(1)*normal(1)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(1)));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(2)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(2)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(1.0 - normal(2)*normal(2)));

  }

  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed)
  {

    //-----------------------------------------------------------     |                                           |
    // surface tension (Laplace-Beltrami surface tension) term      - |  < nabla_BEL v > ° nabla_BEL id_gamma     |
    //-----------------------------------------------------------     |                                           |

    //-----------------------------------------------------------------------------------------------------
    //additional terms to non-smoothed Laplace-Beltrami-version
    //better approximation of curvature
    //-----------------------------------------------------------------------------------------------------

    //    1. needs interpolated normal vectors ( normal_Phi = nabla(Phi)/ |nabla(Phi)|)

    //here are used unenriched shape functions,
    // because scatra-field and fluid-field have the same discretization we can use the variables
    // derxy which is multiplied with the jacobian
    // ephi contains the nodal values of phi in element ele
    // numnode is the number of nodes in the 3D element


    //exact normals for static bubble 3D test example
    static LINALG::Matrix<nsd,1> normal_Phi;
    normal_Phi.Clear();

    // get shape functions at gaussian point
    static LINALG::Matrix<numnode,1> funct_gp;
    funct_gp.Clear();
    DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    // loop over nodes
    for(size_t i = 0; i< numnode; i++)
    {
      normal_Phi(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
      normal_Phi(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
      normal_Phi(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);
    }


    // normalize final gradient
    if (normal_Phi.Norm2() == 0.0 ) dserror("normal norm is zero !!! ");
    normal_Phi.Scale(-1.0/normal_Phi.Norm2());

    // new variant of implementation
    LINALG::Matrix<nsd,nsd> Proj;
    Proj.Clear();
    LINALG::Matrix<nsd,nsd> Proj_smooth;
    Proj_smooth.Clear();
    LINALG::Matrix<nsd,nsd> Proj_combine;
    Proj_combine.Clear();

    for (size_t isd = 0; isd < nsd; isd++) //rows
    {
      for(size_t jsd= 0; jsd < nsd; jsd++) // cols
      {
        if(isd==jsd){
          Proj(isd,jsd) = 1.0 - normal(isd)*normal(jsd);
          Proj_smooth(isd,jsd) = 1.0 - normal_Phi(isd)*normal_Phi(jsd);
        }
        else{
          Proj(isd,jsd) = -normal(isd)*normal(jsd);
          Proj_smooth(isd,jsd) = -normal_Phi(isd)*normal_Phi(jsd);
        }
      }
    }

    Proj_combine.MultiplyNN(Proj,Proj_smooth);

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,0));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,1));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,2));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,2));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,2));


  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_fixed_curvature)
  {

    //----------------------------------------   |                                 |
    // surface tension (constant jump) term      |  < v > , surftens*Curv*n >      |
    //----------------------------------------   |                                 |
    double rad_1 = 0.25;
    double rad_2 = 0.25;

    // 2D
    //rad_2 = 0.0; // not needed; in fact, rad_2 is infinity
    //double curvature = 1.0/rad_1;
    // 3D
    double curvature = 1.0/rad_1 + 1.0/rad_2;

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(0) );
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(1) );
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(2) );

  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_divgrad
          or surftensapprox == INPAR::COMBUST::surface_tension_approx_divgrad_normal
          or surftensapprox == INPAR::COMBUST::surface_tension_approx_nodal_curvature)
  {

    //----------------------------------------     |                                 |
    // surface tension (constant jump) term      - |  < v > , surftens*Curv*n >      |
    //----------------------------------------     |                                 |



    double curvature = 0.0;

    if (surftensapprox == INPAR::COMBUST::surface_tension_approx_nodal_curvature)
    {
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
       DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      for(size_t i=0; i<numnode; i++)
        curvature -= ecurv(i) * funct_gp(i); // Do not ask me why ecurv contains the negative surface tension...
                                             // ... answer curvature is multiplied by -1 in flamefront
    }
    else
    {
      // compute curvature at gausspoint
      COMBUST::CalcCurvature<DISTYPE>(curvature,posXiDomain,xyze,grad_ephi_smoothed,grad_ephi2_smoothed,second_deriv);

      if (fabs(curvature)<1.0E-9) // spurious velocities for almost planar surfaces observed
      {                           // -> set curvature to zero
        curvature=0.0;
//        std::cout << "small curvature value < 1.0e-9 set to zero" << std::endl;
      }
    }


    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(0) );
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(1) );
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, timefacfac*curvature*surftenscoeff*normal(2) );
  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_none){;}
  else dserror("warning: wrong type of surface tension approximation");

  //add term for variable surface tension coefficient
 if (variablesurftens != 0.0)
 {
   assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, - timefacfac*surftenscoeffx*(1.0 - normal(0)*normal(0)));
   assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, - timefacfac*surftenscoeffy*(    - normal(1)*normal(0)));
   assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, - timefacfac*surftenscoeffz*(    - normal(2)*normal(0)));

   assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, - timefacfac*surftenscoeffx*(    - normal(0)*normal(1)));
   assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, - timefacfac*surftenscoeffy*(1.0 - normal(1)*normal(1)));
   assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, - timefacfac*surftenscoeffz*(    - normal(2)*normal(1)));

   assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, - timefacfac*surftenscoeffx*(    - normal(0)*normal(2)));
   assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, - timefacfac*surftenscoeffy*(    - normal(1)*normal(2)));
   assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, - timefacfac*surftenscoeffz*(1.0 - normal(2)*normal(2)));
 }

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
DRT::Element::DiscretizationType CELLDISTYPE,
XFEM::AssemblyType ASSTYPE,
size_t NUMDOF,
class M1, class M2, class M3, class M4, class V1, class V2, class V3>
void SysmatBoundarySurfaceTensionIntCell(
    const XFEM::ElementDofManager&    dofman,
    GEO::BoundaryIntCell              cell,
    const DRT::UTILS::GaussRule2D&    intrule2D,
    M2&                               xyze,
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const bool&                       higher_order_ele,
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const M4&                         grad_ephi2_smoothed, ///<smoothed nodal G-function gradient values for this element
    const V3&                         ecurv,              ///<curvature for this element?
    const double&                       timefac,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const size_t&                              numparamvelx,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const double                               variablesurftens, ///<
    const bool                                 second_deriv,
    const bool                                 connected_interface,
    const bool                                  smoothed_boundary_integration
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

  // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
  LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
  // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
  LINALG::SerialDenseMatrix cellXYZDomaintmp = cell.CellNodalPosXYZ();
  // transform to fixed size format
  const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
  const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

  // assembly line integrals for surface tension if interface is not connected
  // and has intersection with domain boundary
  // no line integrals for not laplace-beltrami surface tension approximation types
  if(!connected_interface
      && ( surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami ||
           surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed )
  ) // only for laplace_beltrami approximation
  {
    COMBUST::RHS_LaplaceBeltrami_1DLineIntegrals_TPF<DISTYPE,CELLDISTYPE,ASSTYPE,NUMDOF>(
        ele, dofman, ephi, cellXiDomain, xyze, grad_ephi_smoothed,
        numparamvelx,timefac,surftensapprox,surftenscoeff,assembler
    );
  }

  //----------------------------------------------------------------------------------------------
  // integration loop over Gaussian points
  //----------------------------------------------------------------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // new transformation for boundary integrals
    // 1. define a coupled transformation x_3D(xi_3D(eta_2D)): transformation from 2D->3D
    // 2. compute the corresponding Jacobian J_eta2D->x_3D and
    // 3. the corresponding surface integral factor sqrt(det(J_eta2D->x_3D^T * J_eta2D->x_3D))
    // 4. approximate integral with Gauss rule in eta coordinate system
    // 5. evaluate the transformed integrand f(x(xi(eta)))

    const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

    //  std::cout << "numnode" << numnode << std::endl;

    // jacobian for coupled transformation
    // get derivatives dxi_3D/deta_2D
    static LINALG::Matrix<2,numvertices> deriv_eta2D;
    DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),CELLDISTYPE);

    // calculate dxi3Ddeta2D
    static LINALG::Matrix<3,2> dXi3Ddeta2D;
    dXi3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < (int)numvertices; k++)
          dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

    //    std::cout << "dXi3Ddeta2D: " << dXi3Ddeta2D << std::endl;

    // transform Gauss point to xi3D space (element parameter space)
    static LINALG::Matrix<nsd,1> gpinXi3D;
    gpinXi3D.Clear();
    // coordinates of this integration point in element coordinates \xi^domain
    GEO::mapEtaBToXiD(cell, gpinEta2D, gpinXi3D);

    //        std::cout << "gpinXi3D" << gpinXi3D << std::endl;

    static LINALG::Matrix<3,numnode> deriv_xi3D;
    DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);

    //    std::cout << "deriv_xi3D" << deriv_xi3D << std::endl;

    // calculate dx3Ddxi3D
    static LINALG::Matrix<3,3> dX3DdXi3D;
    dX3DdXi3D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 3; j++) // derivatives
        for (int k = 0; k < (int)numnode; k++)
          dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);
    //    std::cout << "dX3DdXi3D: " << dX3DdXi3D << std::endl;

    // get the couped Jacobian dx3Ddeta2D
    static LINALG::Matrix<3,2> dx3Ddeta2D;
    dx3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < 3; k++)
          dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);
    //    std::cout << "dx3Ddeta2D: " << dx3Ddeta2D << std::endl;
    // get deformation factor
    static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
    Jac_tmp.Clear();
    Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

    if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");
    //    std::cout << "jac_det" << Jac_tmp.Determinant() << std::endl;
    const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

    const double fac = intpoints.qwgt[iquad]*deform_factor;

    LINALG::Matrix<2,1> posEtaBoundary;
    posEtaBoundary.Clear();
    for (int i= 0; i< 2; i++){
      posEtaBoundary(i,0) = gpinEta2D(i,0);
    }

    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    for (int i= 0; i< 3; i++){
      posXiDomain(i,0) = gpinXi3D(i,0);
    }

    // compute total (time and spatial) integration factor (and coefficients of single terms)?
    const double timefacfac = timefac * fac;
    //std::cout << "boundary time and spatial factor " << timefacfac << std::endl;

    //--------------------------------------------------------------------------------------------
    // compute normal vector (in physical coordinates)
    // remark: for linear boundary integration cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();

    // normal points from Omega^+ to Omega^-
    if ((smoothed_boundary_integration == true)
        and (surftensapprox != INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed
        and surftensapprox != INPAR::COMBUST::surface_tension_approx_divgrad_normal
        and surftensapprox != INPAR::COMBUST::surface_tension_approx_nodal_curvature))
    {
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      normal.Clear();
      LINALG::Matrix<nsd,1> grad_xi;
      grad_xi.Clear();

      // loop over nodes
      for(size_t i = 0; i< numnode; i++)
      {
        normal(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
        normal(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
        normal(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);
      }

      // 4. normalize final gradient
      if (normal.Norm2() == 0.0 ) dserror("normal norm is zero !!! ");
//      normal.Scale(1.0/normal.Norm2());
      normal.Scale(-1.0/normal.Norm2());
    }
    else //(smoothed_boundary_integration == false)
    {
      GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);
    }

    //--------------------------------------------------------------------------------------------
    // evaluate shape functions and their first derivatives at this Gaussian point
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<numnode,1> funct;
    funct.Clear();
    static LINALG::Matrix<nsd,numnode> deriv;
    deriv.Clear();
    DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    //--------------------------------------------------------------------------------------------
    // procedures involving Jacobian matrix for domain mapping
    //--------------------------------------------------------------------------------------------
    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd> xjm;
    xjm.Clear();
    xjm.MultiplyNT(deriv,xyze);
    // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
    const double detXtoXi = xjm.Determinant();

    if (detXtoXi < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
    }
    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> xji;
    xji.Clear();
    xji.Invert(xjm);

    //--------------------------------------------------------------------------------------------
    // compute global derivates of shape functions at this Gaussian point
    //--------------------------------------------------------------------------------------------
    // compute first global derivative
    static LINALG::Matrix<3,numnode> derxy;
    derxy.Clear();
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    derxy.Multiply(xji,deriv);

    // compute second global derivative
    static LINALG::Matrix<6,numnode> derxy2;
    if (higher_order_ele) // if (higher_order_ele)
    {
      static LINALG::Matrix<6,numnode> deriv2;
      DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }

    // get global coordinates to compute location dependent surface tension coefficient
    // for the time being only linear in x available
    LINALG::Matrix<3,1> xgp(true);
    for (size_t i = 0; i < nsd; ++i)
      for (size_t j = 0; j < numnode; ++j)
        xgp(i,0)+=funct(j,0)*xyze(i,j);
    // gamma(x) = gamma - dgamma/dx*x
    double complete_surftenscoeff = surftenscoeff - variablesurftens*xgp(0,0);

    //==========================================================================================
    //------------------------------------------------------------------------------------------
    // evaluate the enrichment function at the interface (boundary integration cells)
    // combined evaluation for kink enrichments in velocity and jump enrichments in pressure
    //------------------------------------------------------------------------------------------

    const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // remark: we need the velocity shape function here
    //         and we have kink enrichments for velocity
    //         hence, XFEM::Enrichment::approachFromPlus is not used
    //         we transfer this parameter because this class also computes pressure shape functions and otherwise
    //         we get a dserror
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi,funct,derxy,derxy2);

    //--------------------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------------------

    // temporary variables holding enriched shape functions (N * \Psi) on either side of the interface
    // for velocity (kink enriched)
    static LINALG::Matrix<shpVecSize,1>       enrfunct_vel;
    enrfunct_vel.Clear();

    static LINALG::Matrix<3,shpVecSize> enrderxy_vel;
    enrderxy_vel.Clear();

    static LINALG::Matrix<6,shpVecSize> enrderxy2_vel;
    enrderxy2_vel.Clear();

    // shape functions for nodal parameters (dofs) on plus and minus side
    // Velx okay, all shapefunctions for Velx, Vely, Velz are the same (kink enriched)
    enrvals.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
        enrfunct_vel, enrderxy_vel, enrderxy2_vel);

    // (enriched) shape functions = approximation functions (P = N * \Psi)

    static XFEM::ApproxFunc<2,shpVecSize> shp_mean_inv;

#if 0
    // version for jumps, unused
    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_mean_inv.d0(iparam) = kappa_minus * enrfunct_plus_vel(iparam) + kappa_plus * enrfunct_minus_vel(iparam);
      //std::cout << "jump   " << shp_mean_inv.d0(iparam) << std::endl;

      shp_mean_inv.dx(iparam) = kappa_minus * enrderxy_plus_vel(0,iparam) + kappa_plus * enrderxy_minus_vel(0,iparam);
      shp_mean_inv.dy(iparam) = kappa_minus * enrderxy_plus_vel(1,iparam) + kappa_plus * enrderxy_minus_vel(1,iparam);
      shp_mean_inv.dz(iparam) = kappa_minus * enrderxy_plus_vel(2,iparam) + kappa_plus * enrderxy_minus_vel(2,iparam);
      //std::cout << "jump   " << shp_mean_inv.dx(iparam) << std::endl;
    }
#endif
    // new reduced version for kink enrichments
    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_mean_inv.d0(iparam) = enrfunct_vel(iparam);

      shp_mean_inv.dx(iparam) = enrderxy_vel(0,iparam);
      shp_mean_inv.dy(iparam) = enrderxy_vel(1,iparam);
      shp_mean_inv.dz(iparam) = enrderxy_vel(2,iparam);
    }

    //================================ SURFACE TENSION ASSEMBLY =================================
    //---------------------assemble 2D surface integrals -------------------------
    COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(
    posXiDomain,shp_mean_inv,ele,xyze,
    grad_ephi_smoothed,grad_ephi2_smoothed,ecurv,normal,timefacfac,
    surftensapprox,variablesurftens, second_deriv, complete_surftenscoeff,
    assembler);

    //===========================================================================================
  }// loop Gaussian points

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
XFEM::AssemblyType ASSTYPE,
int NUMDOF,
class M1, class M2, class M3, class M4, class V1, class V2, class V3>
void SysmatBoundarySurfaceTension(
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const M4&                         grad_ephi2_smoothed, ///<smoothed nodal G-function gradient values for this element
    const V3&                         ecurv,              ///<curvature for this element?
    const M2&                         etau,
    Teuchos::RCP<const MAT::Material> material,      ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    const double                      ga_alphaF,
    const double                      ga_alphaM,
    const double                      ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const double                      flamespeed,
    const double                      nitschevel,
    const double                      nitschepres,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               variablesurftens,
    const bool                                 second_deriv,
    const bool                                 connected_interface,
    const bool                                 smoothed_boundary_integration
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - boundary");

#ifdef DEBUG
  if (ASSTYPE != XFEM::xfem_assembly)
    dserror("boundary integration integrals must only be added to intersected XFEM elements");
#endif

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration constant
  double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  if (timealgo == INPAR::FLUID::timeint_afgenalpha) // here the rhs timefac is needed
    timefac /= ga_alphaF;

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  //------------------------------------------------------------------------------------------------
  // get material parameters for all boundary integration cell of this element
  //------------------------------------------------------------------------------------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  double surftenscoeff = 0.0;
  // set surface tension
  COMBUST::GetMaterialParams(material, surftenscoeff);

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();


  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->ElementBoundaryIntCells(ele->Id());

  //------------------------------------------------------------------------------------------------
  // loop over boundary integration cells
  //------------------------------------------------------------------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {

    case DRT::Element::tri3:
      COMBUST::SysmatBoundarySurfaceTensionIntCell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          dofman,
          *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed, grad_ephi2_smoothed, ecurv,
          timefac, assembler,
          numparamvelx, surftensapprox, surftenscoeff, variablesurftens, second_deriv,
          connected_interface, smoothed_boundary_integration
      );
      break;
    case DRT::Element::quad4:
      COMBUST::SysmatBoundarySurfaceTensionIntCell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          dofman,
          *cell, DRT::UTILS::intrule_quad_36point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed, grad_ephi2_smoothed, ecurv,
          timefac, assembler,
          numparamvelx, surftensapprox, surftenscoeff, variablesurftens, second_deriv,
          connected_interface, smoothed_boundary_integration
      );
      break;
    default:
      dserror("cell distype not implemented yet");
    }

  } // loop boundary integration cells

  return;
}


} // namespace COMBUST

#endif
