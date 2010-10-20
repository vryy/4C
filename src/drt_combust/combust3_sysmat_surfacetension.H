/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_TPF_surfacetension.H

\brief system matrix formulation for surface tension boundary integrals / 3D combustion element

<pre>
Maintainer: Benedikt Schott
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat.H"
#include "combust3_sysmat_benedikt.H"
#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_f3/xfluid3_utils.H"
#include "../drt_f3/fluid3_stabilization.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_xfem/enrichment.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"




using namespace XFEM::PHYSICS;

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
    M1&                                        cellXYZDomain,
    M2&                                        xyze,          ///< xyz coordinates of element ele
    const M3&                                  grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const double&                              kappa_plus,
    const double&                              kappa_minus,
    const size_t&                              numparamvelx,
    const double&                              timefac,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
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

  // assume tri3 boundary integration cells
  // if (cellShape != DRT::Element::tri3) dserror("invalid type of boundary integration cell");
  //const size_t numlines = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

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
      linetmp(isd,1) = cellXYZDomain(isd,(inode+1)%numvertices);
      // first point
      linetmp(isd,0) = cellXYZDomain(isd,(inode)%numvertices);
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

    if(isbound==false);
    else // line is a boundary line => integrate boundary integral for surface tension at Interface-boundaries
    {
      //============================ INTEGRATE BOUNDARY TERM ==============================
      //    		cout << "boundary line is found:\n" << *line_it;

      // get integration points
      // TODO get the right intrule
      const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);


      LINALG::Matrix<nsd,numvertices_line> currLine(*line_it);
      //----------------------------------------------------------------------------------------------
      // integration loop over Gaussian points
      //----------------------------------------------------------------------------------------------
      for (int iquad=0; iquad<intpoint_line.nquad; ++iquad)
      {
        // 1D-integration points
        const LINALG::Matrix<1,1> posGammaBound(intpoint_line.qxg[iquad]);

        //--------------------------------------------------------------------------------------------
        // transformation from GammaBoundary to XYZ
        //--------------------------------------------------------------------------------------------
        double detXYZtoGamma = 0.0;


        static LINALG::Matrix<1,numvertices_line> deriv_gamma;
        DRT::UTILS::shape_function_1D_deriv1(deriv_gamma,posGammaBound(0,0),DRT::Element::line2);


        // det for line surface integrals  sqrt(det(D^T*D)) = norm2(D) in 1D case
        static LINALG::Matrix<3,1> dXYZdGamma;
        dXYZdGamma.Clear();
        for (size_t i = 0; i < nsd; i++)		//dimensions
          for (size_t k = 0; k < numvertices_line; k++)		// number of vertices at line
            dXYZdGamma(i,0) += currLine(i,k)*deriv_gamma(0,k);

        detXYZtoGamma = dXYZdGamma.Norm2();

        //--------------------------------------------------------------------------------------------
        // now map posGammaBound to 3D-Xi-coordinates
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> posGammainXi;
        posGammainXi.Clear();

        static LINALG::Matrix<nsd,1> posGammainXYZ;
        posGammainXYZ.Clear();
        static LINALG::Matrix<numvertices_line,1> shape_gamma;
        DRT::UTILS::shape_function_1D(shape_gamma,posGammaBound(0,0),DRT::Element::line2);

        posGammainXYZ.Multiply(currLine,shape_gamma);


        //get Xi-coordinates of posEtaBoundary 3D -> 3D transformation in element, not in cell
        COMBUST::MapXYZToXiDomain<DISTYPE>(xyze,posGammainXYZ,posGammainXi);


        //--------------------------------------------------------------------------------------------
        // compute integration factors
        //--------------------------------------------------------------------------------------------
        // compute spatial integration factor
        const double fac_gamma = intpoint_line.qwgt[iquad]*detXYZtoGamma;

        // compute total (time and spatial) integration factor (and coefficients of single terms)?
        const double timefacfac_gamma = timefac * fac_gamma;


        //--------------------------------------------------------------------------------------------
        // evaluate element shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<numnode,1> funct_gp;
        funct_gp.Clear();
        DRT::UTILS::shape_function_3D(funct_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);

        static LINALG::Matrix<nsd,numnode> deriv_gp;
        deriv_gp.Clear();
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);

        //--------------------------------------------------------------------------------------------
        // procedures involving Jacobian matrix for domain mapping
        //--------------------------------------------------------------------------------------------
        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.Clear();
        xjm_gp.MultiplyNT(deriv_gp,xyze);
        // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
        const double detXtoXi = xjm_gp.Determinant();
        //#ifdef DEBUG
        if (detXtoXi < 0.0)
        {
          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
        }
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
          DRT::UTILS::shape_function_3D_deriv2(deriv2_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);
          DRT::UTILS::gder2<DISTYPE>(xjm_gp, derxy_gp, deriv2_gp, xyze, derxy2_gp);
        }
        else
        {
          derxy2_gp.Clear();
        }

        //--------------------------------------------------------------------------------------------
        // rearrange (enriched) shape functions as approximation functions
        //--------------------------------------------------------------------------------------------
        const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

        const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi,funct_gp,derxy_gp,derxy2_gp);
        const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi,funct_gp,derxy_gp,derxy2_gp);


        // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
        static LINALG::Matrix<shpVecSize,1>       enrfunct_plus;
        enrfunct_plus.Clear();
        static LINALG::Matrix<shpVecSize,1>       enrfunct_minus;
        enrfunct_minus.Clear();


        // shape functions for nodal parameters (dofs) on plus and minus side
        enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct_gp, enrfunct_plus);
        enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct_gp, enrfunct_minus);

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

        //--------------------------------------------------------------------------------------------
        // 1. define vector for line-direction line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> line;
        line.Clear();
        line(0) = currLine(0,1) - currLine(0,0);
        line(1) = currLine(1,1) - currLine(1,0);
        line(2) = currLine(2,1) - currLine(2,0);

        //--------------------------------------------------------------------------------------------
        // 2. get normal_Phi = grad(phi)/| grad_phi | at gp
        //--------------------------------------------------------------------------------------------
        LINALG::Matrix<nsd,1> normal_Phi;
        normal_Phi.Clear();

        // xi-gradient at gaussian point
        LINALG::Matrix<nsd,1> grad_xi;
        grad_xi.Clear();

        for(size_t i=0; i<numnode; i++)
        {
#if(1)
          LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
          LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
          nodal_grad_xyz_tmp.Clear();
          nodal_grad_xi_tmp.Clear();


          // get smoothed xyz-gradient at node i
          // funct = shape functions at Gaussian point
          //nodal_grad_xyz_tmp.Multiply(grad_ephi_smoothed,funct);

          // transform xyz-gradient
          // 1. transformt nodal xyz-gradient (nodal_grad_xyz_tmp) to nodal xi-gradient => nodal_grad_xi_tmp
          // 2. interpolate xi-gradient at Gaussian point          => nodal_grad_xi
          // 3. transform xi-gradient at Gaussian point to xyz-gradient
          // 4. normalize final gradient

          // 1. transform nodal xyz-gradient to nodal xi-gradient
          // get smoothed xyz-gradient at node i
          nodal_grad_xyz_tmp(0) =  grad_ephi_smoothed(0,i);
          nodal_grad_xyz_tmp(1) =  grad_ephi_smoothed(1,i);
          nodal_grad_xyz_tmp(2) =  grad_ephi_smoothed(2,i);

          // interpolate nodal phi gradients
          normal_Phi(0)+=nodal_grad_xyz_tmp(0)*funct_gp(i);
          normal_Phi(1)+=nodal_grad_xyz_tmp(1)*funct_gp(i);
          normal_Phi(2)+=nodal_grad_xyz_tmp(2)*funct_gp(i);


          //     	          	// transform to xi-coordinates
          //     	          	// jacobi at node
          //
          //     	          	// get Xi-coordinates of current node in current element
          //     	          	static LINALG::Matrix<nsd,1> node_Xicoordinates;
          //     	          	node_Xicoordinates.Clear();
          //     	          	node_Xicoordinates = DRT::UTILS::getNodeCoordinates(i, DISTYPE);
          //
          //     	          	//cout << node_Xicoordinates << endl;
          //
          //     	          	// get derivatives of shapefunctions at node
          //     	          	static LINALG::Matrix<nsd,numnode> deriv_node;
          //     	          	deriv_node.Clear();
          //     	          	DRT::UTILS::shape_function_3D_deriv1(deriv_node,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);
          //
          //     	          	// get jacobian at node
          //
          //     	              // get transposed of the jacobian matrix d x / d \xi
          //     	              // xjm(i,j) = deriv(i,k)*xyze(j,k)
          //     	              static LINALG::Matrix<nsd,nsd> xjm_node;
          //     	              xjm_node.Clear();
          //     	              xjm_node.MultiplyNT(deriv_node,xyze);
          //
          //     	    //            static LINALG::Matrix<nsd,nsd> xji_node;
          //     	    //            xji_node.Clear();
          //     	    //            xji_node.Invert(xjm_node);
          //
          //
          //     	    //            nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);
          //     	              nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);
          //
          //     	          	// funct = shape functions at Gaussian point
          //
          //     	          	grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
          //     	          	grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
          //     	          	grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);
          //
          //

#else        	
          LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
          nodal_grad_xi_tmp.Clear();

          // get smoothed xi-gradient at node i
          nodal_grad_xi_tmp(0) =  grad_ephi_smoothed(0,i);
          nodal_grad_xi_tmp(1) =  grad_ephi_smoothed(1,i);
          nodal_grad_xi_tmp(2) =  grad_ephi_smoothed(2,i);


          grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
          grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
          grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);
#endif
        }

        // transform xi-gradient at Gaussian point to xyz-gradient
        // 	            static LINALG::Matrix<nsd,numnode> deriv_gp;
        // 	        	deriv_gp.Clear();
        // 	        	DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);
        //
        //
        // 	            static LINALG::Matrix<nsd,nsd> xjm_gp;
        // 	            xjm_gp.Clear();
        // 	            xjm_gp.MultiplyNT(deriv_gp,xyze);
        //
        // 	            static LINALG::Matrix<nsd,nsd> xji_gp;
        // 	            xji_gp.Clear();
        // 	            xji_gp.Invert(xjm_gp);

        // 	            normal_Phi.Multiply(xji_gp,grad_xi);

        // normalize final gradient
        normal_Phi.Scale(1.0/normal_Phi.Norm2());

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

        assembler.template Vector<Velx>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(0));
        assembler.template Vector<Vely>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(1));
        assembler.template Vector<Velz>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(2));


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
    M1&                                        cellXYZDomain,
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
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
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

  // assume tri3 boundary integration cells
  //if (cellShape != DRT::Element::tri3) dserror("invalid type of boundary integration cell");
  //const size_t numlines = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

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
      linetmp(isd,1) = cellXYZDomain(isd,(inode+1)%numvertices);
      // first point
      linetmp(isd,0) = cellXYZDomain(isd,(inode)%numvertices);
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

    if(isbound==false);
    else // line is a boundary line => integrate boundary integral for surface tension at Interface-boundaries
    {
      //============================ INTEGRATE BOUNDARY TERM ==============================
      //    		cout << "boundary line is found:\n" << *line_it;

      // get integration points
      // TODO get the right intrule
      const DRT::UTILS::IntegrationPoints1D intpoint_line(DRT::UTILS::intrule_line_2point);


      LINALG::Matrix<nsd,numvertices_line> currLine(*line_it);
      //----------------------------------------------------------------------------------------------
      // integration loop over Gaussian points
      //----------------------------------------------------------------------------------------------
      for (int iquad=0; iquad<intpoint_line.nquad; ++iquad)
      {
        // 1D-integration points
        const LINALG::Matrix<1,1> posGammaBound(intpoint_line.qxg[iquad]);

        //--------------------------------------------------------------------------------------------
        // transformation from GammaBoundary to XYZ
        //--------------------------------------------------------------------------------------------
        double detXYZtoGamma = 0.0;


        static LINALG::Matrix<1,numvertices_line> deriv_gamma;
        DRT::UTILS::shape_function_1D_deriv1(deriv_gamma,posGammaBound(0,0),DRT::Element::line2);


        // det for line surface integrals  sqrt(det(D^T*D)) = norm2(D) in 1D case
        static LINALG::Matrix<3,1> dXYZdGamma;
        dXYZdGamma.Clear();
        for (size_t i = 0; i < nsd; i++)		//dimensions
          for (size_t k = 0; k < numvertices_line; k++)		// number of vertices at line
            dXYZdGamma(i,0) += currLine(i,k)*deriv_gamma(0,k);

        detXYZtoGamma = dXYZdGamma.Norm2();

        //--------------------------------------------------------------------------------------------
        // now map posGammaBound to 3D-Xi-coordinates
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> posGammainXi;
        posGammainXi.Clear();

        static LINALG::Matrix<nsd,1> posGammainXYZ;
        posGammainXYZ.Clear();
        static LINALG::Matrix<numvertices_line,1> shape_gamma;
        DRT::UTILS::shape_function_1D(shape_gamma,posGammaBound(0,0),DRT::Element::line2);

        posGammainXYZ.Multiply(currLine,shape_gamma);


        //get Xi-coordinates of posEtaBoundary 3D -> 3D transformation in element, not in cell
        COMBUST::MapXYZToXiDomain<DISTYPE>(xyze,posGammainXYZ,posGammainXi);


        //--------------------------------------------------------------------------------------------
        // compute integration factors
        //--------------------------------------------------------------------------------------------
        // compute spatial integration factor
        const double fac_gamma = intpoint_line.qwgt[iquad]*detXYZtoGamma;

        // compute total (time and spatial) integration factor (and coefficients of single terms)?
        const double timefacfac_gamma = timefac * fac_gamma;


        //--------------------------------------------------------------------------------------------
        // evaluate element shape functions at this Gaussian point
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<numnode,1> funct_gp;
        funct_gp.Clear();
        DRT::UTILS::shape_function_3D(funct_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);


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
        enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct_gp, enrfunct_plus);
        enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct_gp, enrfunct_minus);

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

        //--------------------------------------------------------------------------------------------
        // 1. define vector for line-direction line
        //--------------------------------------------------------------------------------------------
        static LINALG::Matrix<nsd,1> line;
        line.Clear();
        line(0) = currLine(0,1) - currLine(0,0);
        line(1) = currLine(1,1) - currLine(1,0);
        line(2) = currLine(2,1) - currLine(2,0);

        //--------------------------------------------------------------------------------------------
        // 2. get normal_Phi = grad(phi)/| grad_phi | at gp
        //--------------------------------------------------------------------------------------------
        LINALG::Matrix<nsd,1> normal_Phi;
        normal_Phi.Clear();

        // xi-gradient at gaussian point
        LINALG::Matrix<nsd,1> grad_xi;
        grad_xi.Clear();

        for(size_t i=0; i<numnode; i++)
        {
#if(1)
          LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
          LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
          nodal_grad_xyz_tmp.Clear();
          nodal_grad_xi_tmp.Clear();


          // get smoothed xyz-gradient at node i
          // funct = shape functions at Gaussian point
          //nodal_grad_xyz_tmp.Multiply(grad_ephi_smoothed,funct);

          // transform xyz-gradient
          // 1. transformt nodal xyz-gradient (nodal_grad_xyz_tmp) to nodal xi-gradient => nodal_grad_xi_tmp
          // 2. interpolate xi-gradient at Gaussian point          => nodal_grad_xi
          // 3. transform xi-gradient at Gaussian point to xyz-gradient
          // 4. normalize final gradient

          // 1. transform nodal xyz-gradient to nodal xi-gradient
          // get smoothed xyz-gradient at node i
          nodal_grad_xyz_tmp(0) =  grad_ephi_smoothed(0,i);
          nodal_grad_xyz_tmp(1) =  grad_ephi_smoothed(1,i);
          nodal_grad_xyz_tmp(2) =  grad_ephi_smoothed(2,i);

          // transform to xi-coordinates
          // jacobi at node

          // get Xi-coordinates of current node in current element
          static LINALG::Matrix<nsd,1> node_Xicoordinates;
          node_Xicoordinates.Clear();
          node_Xicoordinates = DRT::UTILS::getNodeCoordinates(i, DISTYPE);

          //cout << node_Xicoordinates << endl;

          // get derivatives of shapefunctions at node
          static LINALG::Matrix<nsd,numnode> deriv_node;
          deriv_node.Clear();
          DRT::UTILS::shape_function_3D_deriv1(deriv_node,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);

          // get jacobian at node

          // get transposed of the jacobian matrix d x / d \xi
          // xjm(i,j) = deriv(i,k)*xyze(j,k)
          static LINALG::Matrix<nsd,nsd> xjm_node;
          xjm_node.Clear();
          xjm_node.MultiplyNT(deriv_node,xyze);

          //            static LINALG::Matrix<nsd,nsd> xji_node;
          //            xji_node.Clear();
          //            xji_node.Invert(xjm_node);


          //            nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);
          nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);

          // funct = shape functions at Gaussian point

          grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
          grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
          grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);



#else        	
          LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
          nodal_grad_xi_tmp.Clear();

          // get smoothed xi-gradient at node i
          nodal_grad_xi_tmp(0) =  grad_ephi_smoothed(0,i);
          nodal_grad_xi_tmp(1) =  grad_ephi_smoothed(1,i);
          nodal_grad_xi_tmp(2) =  grad_ephi_smoothed(2,i);


          grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
          grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
          grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);
#endif
        }

        // transform xi-gradient at Gaussian point to xyz-gradient
        static LINALG::Matrix<nsd,numnode> deriv_gp;
        deriv_gp.Clear();
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posGammainXi(0),posGammainXi(1),posGammainXi(2),DISTYPE);


        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.Clear();
        xjm_gp.MultiplyNT(deriv_gp,xyze);

        static LINALG::Matrix<nsd,nsd> xji_gp;
        xji_gp.Clear();
        xji_gp.Invert(xjm_gp);

        normal_Phi.Multiply(xji_gp,grad_xi);

        // normalize final gradient
        normal_Phi.Scale(1.0/normal_Phi.Norm2());

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

        assembler.template Vector<Velx>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(0));
        assembler.template Vector<Vely>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(1));
        assembler.template Vector<Velz>(shp_mean_inv.d0, + timefacfac_gamma*surftenscoeff*normal_gamma(2));


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
class M1, class M2, class M3, class M4, class M5>
void RHS_SurfaceTension2DIntegrals(
    M1&                                        posXiDomain, // position of gaussian point for evaluating of surface tension terms
    M5&                                        shp_mean_inv,
    const DRT::ELEMENTS::Combust3*             ele,           ///< the element those matrix is calculated
    M2&                                        xyze,          ///< xyz coordinates of element ele
    const M3&                                  grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    M4&                                        normal,
    const double&                              timefacfac,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler
)
{
  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  if(surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami)
  {

    //-----------------------------------------------------------     |                                           |
    // surface tension (Laplace-Beltrami surface tension) term      - |  < nabla_BEL v > ° nabla_BEL id_gamma     |
    //-----------------------------------------------------------     |                                           |

    assembler.template Vector<Velx>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(1.0 - normal(0)*normal(0)));
    assembler.template Vector<Velx>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(0)));
    assembler.template Vector<Velx>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(0)));

    assembler.template Vector<Vely>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(1)));
    assembler.template Vector<Vely>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(1.0 - normal(1)*normal(1)));
    assembler.template Vector<Vely>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(1)));

    assembler.template Vector<Velz>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(2)));
    assembler.template Vector<Velz>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(2)));
    assembler.template Vector<Velz>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(1.0 - normal(2)*normal(2)));

  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed)
  {

    //-----------------------------------------------------------     |                                           |
    // surface tension (Laplace-Beltrami surface tension) term      - |  < nabla_BEL v > ° nabla_BEL id_gamma     |
    //-----------------------------------------------------------     |                                           |

    assembler.template Vector<Velx>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(1.0 - normal(0)*normal(0)));
    assembler.template Vector<Velx>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(0)));
    assembler.template Vector<Velx>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(0)));

    assembler.template Vector<Vely>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(1)));
    assembler.template Vector<Vely>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(1.0 - normal(1)*normal(1)));
    assembler.template Vector<Vely>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(    - normal(2)*normal(1)));

    assembler.template Vector<Velz>(shp_mean_inv.dx, - timefacfac*surftenscoeff*(    - normal(0)*normal(2)));
    assembler.template Vector<Velz>(shp_mean_inv.dy, - timefacfac*surftenscoeff*(    - normal(1)*normal(2)));
    assembler.template Vector<Velz>(shp_mean_inv.dz, - timefacfac*surftenscoeff*(1.0 - normal(2)*normal(2)));

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


    LINALG::Matrix<nsd,1> grad_xi;
    grad_xi.Clear();

    // loop over nodes
    for(size_t i = 0; i< numnode; i++)
    {
#if(1)
      LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
      LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
      nodal_grad_xyz_tmp.Clear();
      nodal_grad_xi_tmp.Clear();

      // transform xyz-gradient
      // 1. transformt nodal xyz-gradient (nodal_grad_xyz_tmp) to nodal xi-gradient => nodal_grad_xi_tmp
      // 2. interpolate xi-gradient at Gaussian point          => nodal_grad_xi
      // 3. transform xi-gradient at Gaussian point to xyz-gradient
      // 4. normalize final gradient

      // 1. transform nodal xyz-gradient to nodal xi-gradient
      // get smoothed xyz-gradient at node i
      nodal_grad_xyz_tmp(0) =  grad_ephi_smoothed(0,i);
      nodal_grad_xyz_tmp(1) =  grad_ephi_smoothed(1,i);
      nodal_grad_xyz_tmp(2) =  grad_ephi_smoothed(2,i);


      // new
      normal_Phi(0)+=nodal_grad_xyz_tmp(0)*funct_gp(i);
      normal_Phi(1)+=nodal_grad_xyz_tmp(1)*funct_gp(i);
      normal_Phi(2)+=nodal_grad_xyz_tmp(2)*funct_gp(i);


      //      	// 2. transform to xi-coordinates
      //      	// jacobi at node
      //
      //      	// get Xi-coordinates of current node in current element
      //      	static LINALG::Matrix<nsd,1> node_Xicoordinates;
      //      	node_Xicoordinates.Clear();
      //      	node_Xicoordinates = DRT::UTILS::getNodeCoordinates(i, DISTYPE);
      //
      //
      //      	// get derivatives of shapefunctions at node
      //      	static LINALG::Matrix<nsd,numnode> deriv_node;
      //      	deriv_node.Clear();
      //      	DRT::UTILS::shape_function_3D_deriv1(deriv_node,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);
      //
      //      	// get jacobian at node
      //
      //      	// get transposed of the jacobian matrix d x / d \xi
      //      	// xjm(i,j) = deriv(i,k)*xyze(j,k)
      //      	static LINALG::Matrix<nsd,nsd> xjm_node;
      //      	xjm_node.Clear();
      //      	xjm_node.MultiplyNT(deriv_node,xyze);
      //
      //      	nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);
      //
      //      	// funct_gp = shape functions at Gaussian point
      //
      //      	grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
      //      	grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
      //      	grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);
#else // trafo with xi-grad stored grad_ephi_smoothed
      LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
      nodal_grad_xi_tmp.Clear();

      // get smoothed xi-gradient at node i
      nodal_grad_xi_tmp(0) =  grad_ephi_smoothed(0,i);
      nodal_grad_xi_tmp(1) =  grad_ephi_smoothed(1,i);
      nodal_grad_xi_tmp(2) =  grad_ephi_smoothed(2,i);

      //        	cout << "nodal_grad_xi_tmp:\n: " << nodal_grad_xi_tmp << endl;

      grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
      grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
      grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);

#endif


    }

    //      // 3. transform xi-gradient at Gaussian point to xyz-gradient
    //
    //      // TODO: already defined above
    //      static LINALG::Matrix<nsd,numnode> deriv_gp;
    //      deriv_gp.Clear();
    //      DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    //
    //      static LINALG::Matrix<nsd,nsd> xjm_gp;
    //      xjm_gp.Clear();
    //      xjm_gp.MultiplyNT(deriv_gp,xyze);
    //
    //      static LINALG::Matrix<nsd,nsd> xji_gp;
    //      xji_gp.Clear();
    //      xji_gp.Invert(xjm_gp);


    //      normal_Phi.Multiply(xji_gp,grad_xi);

    // 4. normalize final gradient
    normal_Phi.Scale(1.0/normal_Phi.Norm2());

#if(0)
    // new variant of implementation
    static LINALG::Matrix<nsd,nsd> Proj;
    Proj.Clear();
    static LINALG::Matrix<nsd,nsd> Proj_smooth;
    Proj_smooth.Clear();
    static LINALG::Matrix<nsd,nsd> Proj_combine;
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


    assembler.template Vector<Velx>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,0));
    assembler.template Vector<Velx>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,0));
    assembler.template Vector<Velx>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,0));

    assembler.template Vector<Vely>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,1));
    assembler.template Vector<Vely>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,1));
    assembler.template Vector<Vely>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,1));

    assembler.template Vector<Velz>(shp_mean_inv.dx, - timefacfac * surftenscoeff*Proj_combine(0,2));
    assembler.template Vector<Velz>(shp_mean_inv.dy, - timefacfac * surftenscoeff*Proj_combine(1,2));
    assembler.template Vector<Velz>(shp_mean_inv.dz, - timefacfac * surftenscoeff*Proj_combine(2,2));




#endif
    //------------------------assemble smoothed Laplace-Beltrami-terms---------------------------

    // TODO: implement normal if scatra-field and fluid-field have different discretizations
    // a) for each gp load the phi-element in which x_gp lies
    // b) calculate xi_domain of x_gp in phi-element
    // c) calculate shape function derivates in phi-element
    // d) interpolate nabla(Phi) at gp
    // e) scale normal=nabla(Phi)/ |nabla(Phi)|


    //    2. define n°normal_Phi =: normal_dotp
    double normal_dotp = normal(0)*normal_Phi(0) + normal(1)*normal_Phi(1) + normal(2)*normal_Phi(2);

    //    3. assemble additional terms
    assembler.template Vector<Velx>(shp_mean_inv.dx, - timefacfac * surftenscoeff*( normal_dotp*normal(0) -normal_Phi(0))
        * normal_Phi(0) );
    assembler.template Vector<Velx>(shp_mean_inv.dy, - timefacfac * surftenscoeff*( normal_dotp*normal(1) -normal_Phi(1))
        * normal_Phi(0) );
    assembler.template Vector<Velx>(shp_mean_inv.dz, - timefacfac * surftenscoeff*( normal_dotp*normal(2) -normal_Phi(2))
        * normal_Phi(0) );

    assembler.template Vector<Vely>(shp_mean_inv.dx, - timefacfac * surftenscoeff*( normal_dotp*normal(0) -normal_Phi(0))
        * normal_Phi(1) );
    assembler.template Vector<Vely>(shp_mean_inv.dy, - timefacfac * surftenscoeff*( normal_dotp*normal(1) -normal_Phi(1))
        * normal_Phi(1) );
    assembler.template Vector<Vely>(shp_mean_inv.dz, - timefacfac * surftenscoeff*( normal_dotp*normal(2) -normal_Phi(2))
        * normal_Phi(1) );

    assembler.template Vector<Velz>(shp_mean_inv.dx, - timefacfac * surftenscoeff*( normal_dotp*normal(0) -normal_Phi(0))
        * normal_Phi(2) );
    assembler.template Vector<Velz>(shp_mean_inv.dy, - timefacfac * surftenscoeff*( normal_dotp*normal(1) -normal_Phi(1))
        * normal_Phi(2) );
    assembler.template Vector<Velz>(shp_mean_inv.dz, - timefacfac * surftenscoeff*( normal_dotp*normal(2) -normal_Phi(2))
        * normal_Phi(2) );

  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_fixed_curvature)
  {

    //----------------------------------------     |                                 |
    // surface tension (constant jump) term      - |  < v > , surftens*Curv*n >      |
    //----------------------------------------     |                                 |
    double rad_1 = 0.3;
    double rad_2 = 0.3;
    double curvature = 1.0/rad_1 + 1.0/rad_2;


    assembler.template Vector<Velx>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(0) );
    assembler.template Vector<Vely>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(1) );
    assembler.template Vector<Velz>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(2) );

  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_divgrad)
  {

    //----------------------------------------     |                                 |
    // surface tension (constant jump) term      - |  < v > , surftens*Curv*n >      |
    //----------------------------------------     |                                 |

    //==========================reconstruct normal_phi = grad(phi) / |grad(phi)| ===========

    static LINALG::Matrix<nsd,1> normal_Phi;
    normal_Phi.Clear();

    // get shape functions at gaussian point
    static LINALG::Matrix<numnode,1> funct_gp;
    funct_gp.Clear();
    DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    // TODO: already defined above
    static LINALG::Matrix<nsd,numnode> deriv_gp;
    deriv_gp.Clear();
    DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    static LINALG::Matrix<nsd,nsd> xjm_gp;
    xjm_gp.Clear();
    xjm_gp.MultiplyNT(deriv_gp,xyze);

    static LINALG::Matrix<nsd,nsd> xji_gp;
    xji_gp.Clear();
    xji_gp.Invert(xjm_gp);

    // derivatives in xyz-direction
    static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
    deriv_gp_xyz.Clear();
    deriv_gp_xyz.Multiply(xji_gp,deriv_gp);


    LINALG::Matrix<nsd,1> grad_xi;
    grad_xi.Clear();

    // get second derivatives of phi
    LINALG::Matrix<9,1> grad_phi2;
    grad_phi2.Clear();

    // loop over nodes
    for(size_t i = 0; i< numnode; i++)
    {
#if(1)
      LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
      LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
      nodal_grad_xyz_tmp.Clear();
      nodal_grad_xi_tmp.Clear();


      // transform xyz-gradient
      // 1. transformt nodal xyz-gradient (nodal_grad_xyz_tmp) to nodal xi-gradient => nodal_grad_xi_tmp
      // 2. interpolate xi-gradient at Gaussian point          => nodal_grad_xi
      // 3. transform xi-gradient at Gaussian point to xyz-gradient
      // 4. normalize final gradient

      // 1. transform nodal xyz-gradient to nodal xi-gradient
      // get smoothed xyz-gradient at node i
      nodal_grad_xyz_tmp(0) =  grad_ephi_smoothed(0,i);
      nodal_grad_xyz_tmp(1) =  grad_ephi_smoothed(1,i);
      nodal_grad_xyz_tmp(2) =  grad_ephi_smoothed(2,i);

      normal_Phi(0) += nodal_grad_xyz_tmp(0) *funct_gp(i);
      normal_Phi(1) += nodal_grad_xyz_tmp(1) *funct_gp(i);
      normal_Phi(2) += nodal_grad_xyz_tmp(2) *funct_gp(i);

      // transform to xi-coordinates
      // jacobi at node

      // get Xi-coordinates of current node in current element
      static LINALG::Matrix<nsd,1> node_Xicoordinates;
      node_Xicoordinates.Clear();
      node_Xicoordinates = DRT::UTILS::getNodeCoordinates(i, DISTYPE);


      // get derivatives of shapefunctions at node
      static LINALG::Matrix<nsd,numnode> deriv_node;
      deriv_node.Clear();
      DRT::UTILS::shape_function_3D_deriv1(deriv_node,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);

      // get jacobian at node

      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm_node;
      xjm_node.Clear();
      xjm_node.MultiplyNT(deriv_node,xyze);

      nodal_grad_xi_tmp.Multiply(xjm_node,nodal_grad_xyz_tmp);

      // funct_gp = shape functions at Gaussian point

      grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
      grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
      grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);


      grad_phi2(0) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(0); // ,xx
      grad_phi2(1) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(1); // ,yy
      grad_phi2(2) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(2); // ,zz
      grad_phi2(3) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(0); // ,xy
      grad_phi2(4) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(0); // ,xz
      grad_phi2(5) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(1); // ,yz
      grad_phi2(6) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(1); // ,yx
      grad_phi2(7) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(2); // ,zx
      grad_phi2(8) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(2); // ,zy

#else // trafo with xi-grad stored grad_ephi_smoothed
      LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
      nodal_grad_xi_tmp.Clear();

      // get smoothed xi-gradient at node i
      nodal_grad_xi_tmp(0) =  grad_ephi_smoothed(0,i);
      nodal_grad_xi_tmp(1) =  grad_ephi_smoothed(1,i);
      nodal_grad_xi_tmp(2) =  grad_ephi_smoothed(2,i);

      //        	cout << "nodal_grad_xi_tmp:\n: " << nodal_grad_xi_tmp << endl;

      grad_xi(0) += nodal_grad_xi_tmp(0)*funct_gp(i);
      grad_xi(1) += nodal_grad_xi_tmp(1)*funct_gp(i);
      grad_xi(2) += nodal_grad_xi_tmp(2)*funct_gp(i);

#endif


    }

    // 3. transform xi-gradient at Gaussian point to xyz-gradient
    //      normal_Phi.Multiply(xji_gp,grad_xi);

    // 4. normalize final gradient
    //normal_Phi.Scale(1.0/normal_Phi.Norm2());


    //=========================================== calculate curvature =========================
    double curvature = 0.0;

    double grad_phi_norm = normal_Phi.Norm2();

    curvature = -1.0/pow(grad_phi_norm,3)*(  normal_Phi(0)*normal_Phi(0)*grad_phi2(0)
        + normal_Phi(1)*normal_Phi(1)*grad_phi2(1)
        + normal_Phi(2)*normal_Phi(2)*grad_phi2(2)  )
        -1.0/pow(grad_phi_norm,3)*(  normal_Phi(0)*normal_Phi(1)*(grad_phi2(3)+grad_phi2(6) )
            + normal_Phi(0)*normal_Phi(2)*(grad_phi2(4)+grad_phi2(7) )
            + normal_Phi(1)*normal_Phi(2)*(grad_phi2(5)+grad_phi2(8)) )
            +1.0/grad_phi_norm * ( grad_phi2(0)+grad_phi2(1)+grad_phi2(2) );

    assembler.template Vector<Velx>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(0) );
    assembler.template Vector<Vely>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(1) );
    assembler.template Vector<Velz>(shp_mean_inv.d0, - timefacfac*curvature*surftenscoeff*normal(2) );
  }
  else if(surftensapprox == INPAR::COMBUST::surface_tension_approx_none){;}
  else dserror("warning: wrong type of surface tension approximation");

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
DRT::Element::DiscretizationType CELLDISTYPE,
XFEM::AssemblyType ASSTYPE,
size_t NUMDOF,
class M1, class M2, class M3, class V1, class V2>
void Sysmat_TPF_SurfaceTension_Boundary_Cell(
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
    const double&                     jt,
    const double&                       timefac,
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const double&                        kappa_plus,    /// weights of average in plus domain
    const double&                        kappa_minus,    /// weights of average in minus domain
    const size_t&                              numparamvelx,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::NormalTensionJumpType normaltensionjumptype
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
      && ( normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_SurfaceTension)
      && (        surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami
          || surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed )
  ) // only for laplace_beltrami approximation
  {
    COMBUST::RHS_LaplaceBeltrami_1DLineIntegrals_TPF<DISTYPE,CELLDISTYPE,ASSTYPE,NUMDOF>(
        ele, dofman, ephi, cellXYZDomain, xyze, grad_ephi_smoothed,
        kappa_plus,kappa_minus,numparamvelx,timefac,surftensapprox,surftenscoeff,assembler
    );
  }

  //----------------------------------------------------------------------------------------------
  // integration loop over Gaussian points
  //----------------------------------------------------------------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {



    // 1. cell coordinates are stored in cellXYZDomain <nsd,numvertices>

    // 2. 2D-shape functions -> evaluate at each integration point

    // coordinates of this integration point in boundary cell coordinates \eta^boundary
    const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);
    double detXYZtoEta1Eta2 = 0.0;

    static LINALG::Matrix<2,numvertices> deriv_boundary;
    DRT::UTILS::shape_function_2D_deriv1(deriv_boundary,posEtaBoundary(0,0),posEtaBoundary(1,0),CELLDISTYPE);


    // det for surface integrals  ||  dxyz/deta1 x dxyz/deta2  ||
    static LINALG::Matrix<3,2> dXYZdEta1Eta2;
    dXYZdEta1Eta2.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < (int)numvertices; k++)
          dXYZdEta1Eta2(i,j) += cellXYZDomain(i,k)*deriv_boundary(j,k);

    static LINALG::Matrix<3,1> tmp_vec;     // vector for dxyz/deta1 x dxyz/deta2
    tmp_vec.Clear();
    tmp_vec(0) = dXYZdEta1Eta2(1,0) * dXYZdEta1Eta2(2,1) - dXYZdEta1Eta2(2,0) * dXYZdEta1Eta2(1,1);
    tmp_vec(1) = dXYZdEta1Eta2(2,0) * dXYZdEta1Eta2(0,1) - dXYZdEta1Eta2(0,0) * dXYZdEta1Eta2(2,1);
    tmp_vec(2) = dXYZdEta1Eta2(0,0) * dXYZdEta1Eta2(1,1) - dXYZdEta1Eta2(1,0) * dXYZdEta1Eta2(0,1);

    detXYZtoEta1Eta2 = tmp_vec.Norm2();


    // now map posEtaBoundary in 3D-Xi-coordinates
    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    // coordinates of this integration point in element coordinates \xi^domain

    // the same results like GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain

    //TODO: check this transformation
    // schott Apr 30, 2010
    //GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

    // schott Apr 30, 2010
    static LINALG::Matrix<nsd,1> posEtaBoundaryinXYZ;
    posEtaBoundaryinXYZ.Clear();
    static LINALG::Matrix<numvertices,1> shape_boundary;
    DRT::UTILS::shape_function_2D(shape_boundary,posEtaBoundary(0,0),posEtaBoundary(1,0),cell.Shape());

    posEtaBoundaryinXYZ.Multiply(cellXYZDomain,shape_boundary);
    //get Xi-coordinates of posEtaBoundary 3D -> 3D transformation in element, not in cell
    COMBUST::MapXYZToXiDomain<DISTYPE>(xyze,posEtaBoundaryinXYZ,posXiDomain);

    //        GEO::currentToVolumeElementCoordinates(DISTYPE, xyze, posEtaBoundaryinXYZ, posXiDomain);


    //--------------------------------------------------------------------------------------------
    // compute integration factors
    //--------------------------------------------------------------------------------------------
    // compute spatial integration factor
    const double fac = intpoints.qwgt[iquad]*detXYZtoEta1Eta2;
    //cout << "boundary spatial factor " << fac << endl;
    // compute total (time and spatial) integration factor (and coefficients of single terms)?
    const double timefacfac = timefac * fac;
    //cout << "boundary time and spatial factor " << timefacfac << endl;


    //--------------------------------------------------------------------------------------------
    // compute normal vector (in physical coordinates)
    // remark: for linear boundary integrastion cells this could be done before the loop over all
    //         Gaussian points
    // TODO: is this normal normed?
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();

    // here normal points from Omega^+ to Omega^-
    GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);
    // but we need a normal vector pointing from Omega^- to Omega^+ -> Scale(-1.0) => n = n^- = -n^+
    normal.Scale(-1.0);

    static LINALG::Matrix<nsd,1> jump_tension;
    jump_tension.Clear();

    if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_Const)
    {
      // normal tension jump is set to jt
      jump_tension(0,0) = jt * normal(0);
      jump_tension(1,0) = jt * normal(1);
      jump_tension(2,0) = jt * normal(2);
    }
    // no jump value necessary for surface tension  => another implementation (jump value implicit given)

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
      DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }
    //==========================================================================================
    //------------------------------------------------------------------------------------------
    // evaluate the enrichment function at the interface (boundary integration cells)
    // combined evaluation for kink enrichments in velocity and jump enrichments in pressure
    //------------------------------------------------------------------------------------------

    const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi,funct,derxy,derxy2);
    const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi,funct,derxy,derxy2);


    //--------------------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------------------

    // temporary variables holding enriched shape functions (N * \Psi) on either side of the interface
    // for velocity (kink enriched)
    static LINALG::Matrix<shpVecSize,1>       enrfunct_plus_vel;
    enrfunct_plus_vel.Clear();
    static LINALG::Matrix<shpVecSize,1>       enrfunct_minus_vel;
    enrfunct_minus_vel.Clear();

    static LINALG::Matrix<3,shpVecSize> enrderxy_plus_vel;
    enrderxy_plus_vel.Clear();
    static LINALG::Matrix<3,shpVecSize> enrderxy_minus_vel;
    enrderxy_minus_vel.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_plus_vel;
    enrderxy2_plus_vel.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_minus_vel;
    enrderxy2_minus_vel.Clear();

    // shape functions for nodal parameters (dofs) on plus and minus side
    // Velx okay, all shapefunctions for Velx, Vely, Velz are the same (kink enriched)
    enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_plus_vel, enrderxy_plus_vel, enrderxy2_plus_vel);
    enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_minus_vel, enrderxy_minus_vel, enrderxy2_minus_vel);

    // (enriched) shape functions = approximation functions (P = N * \Psi)

    static XFEM::ApproxFunc<2,shpVecSize> shp_mean_inv;

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_mean_inv.d0(iparam) = kappa_minus * enrfunct_plus_vel(iparam) + kappa_plus * enrfunct_minus_vel(iparam);

      shp_mean_inv.dx(iparam) = kappa_minus * enrderxy_plus_vel(0,iparam) + kappa_plus * enrderxy_minus_vel(0,iparam);
      shp_mean_inv.dy(iparam) = kappa_minus * enrderxy_plus_vel(1,iparam) + kappa_plus * enrderxy_minus_vel(1,iparam);
      shp_mean_inv.dz(iparam) = kappa_minus * enrderxy_plus_vel(2,iparam) + kappa_plus * enrderxy_minus_vel(2,iparam);
    }

    if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_SurfaceTension)
    {
      //================================ SURFACE TENSION ASSEMBLY =================================
      //---------------------assemble 2D surface integrals -------------------------
      COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(posXiDomain,shp_mean_inv,ele,xyze,grad_ephi_smoothed,normal,timefacfac,surftensapprox,surftenscoeff,assembler);

      //===========================================================================================
    }
    else
    { // all other normal tension jumps

      //----------------------------------------     |                     |
      // surface tension (constant jump) term      - |  < v > , J_t >      |
      //----------------------------------------     |                     |


      assembler.template Vector<Velx>(shp_mean_inv.d0, - timefacfac*jump_tension(0,0) );
      assembler.template Vector<Vely>(shp_mean_inv.d0, - timefacfac*jump_tension(1,0) );
      assembler.template Vector<Velz>(shp_mean_inv.d0, - timefacfac*jump_tension(2,0) );
    }

  }// loop Gaussian points

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
XFEM::AssemblyType ASSTYPE,
int NUMDOF,
class M1, class M2, class M3, class V1, class V2>
void Sysmat_TPF_SurfaceTension_Boundary(
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const M2&                         eepsilon,
    Teuchos::RCP<const MAT::Material> material,      ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const double                      flamespeed,
    const double                      nitschevel,
    const double                      nitschepres,
    const double                      ele_meas_plus,    /// measure of element in plus domain
    const double                      ele_meas_minus,    /// measure of element in minus domain
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType veljumptype,
    const INPAR::COMBUST::NormalTensionJumpType normaltensionjumptype
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
  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

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


#ifdef TPF_NITSCHEPRES
  const double alphap = nitschepres;
  //cout << "alpha_p: " << alphap << endl;
#endif

  //=====================================================================================================
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus

  double ele_meas = ele_meas_plus + ele_meas_minus;
  double kappa_plus  = ele_meas_plus/ele_meas;
  double kappa_minus = ele_meas_minus/ele_meas;

  //=====================================================================================================

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();


#if 0
  const bool epsilonele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
  // for now, I don't try to compare to elements without stress unknowns, since they lock anyway
#endif

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
#if 0
  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, the numparam is determined by the dofmanager
  //const int numparamepsilonxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
#endif

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

  //------------------------------------------------------------------------------------------------
  // compute jump values
  //------------------------------------------------------------------------------------------------
  double jt = 0.0;

  // jump in normal tension with value jt
  if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_Const) jt = 3.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);



  //------------------------------------------------------------------------------------------------
  // loop over boundary integration cells
  //------------------------------------------------------------------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {

    case DRT::Element::tri3:
      COMBUST::Sysmat_TPF_SurfaceTension_Boundary_Cell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          dofman,
          *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed,
          jt,
          timefac, assembler,
          kappa_plus, kappa_minus,
          numparamvelx, surftensapprox, surftenscoeff, connected_interface, veljumptype, normaltensionjumptype
      );
      break;
    case DRT::Element::quad4:
      // TODO check this intrule !!!
      COMBUST::Sysmat_TPF_SurfaceTension_Boundary_Cell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          dofman,
          *cell, DRT::UTILS::intrule_quad_9point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed,
          jt,
          timefac, assembler,
          kappa_plus, kappa_minus,
          numparamvelx, surftensapprox, surftenscoeff, connected_interface, veljumptype, normaltensionjumptype
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
#endif
