/*----------------------------------------------------------------------*/
/*!
\file combust3_transformations_benedikt.H

\brief system matrix formulation for surface tension boundary integrals / 3D combustion element

<pre>
Maintainer: Benedikt Schott
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "combust3_sysmat.H"
#include "combust_defines.H"
#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust3_interpolation.H"
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
  Transforms Gaussian point of cell from 3D-Eta-coordinates in 3D-XYZ-coordinates
 */
template <DRT::Element::DiscretizationType CELL_DISTYPE,
          class M1, class M2>
void GetGPinXYZ(
    GEO::DomainIntCell          cell,
    M1&                         pos_eta_domain,
    M2&                         posGPinXYZ)
{
  const int nsd = 3;

  // number of nodes of this cell
  const size_t numnode3Dcell = DRT::UTILS::DisTypeToNumNodePerEle<CELL_DISTYPE>::numNodePerElement;

  // cel 3D-shape functions
  LINALG::Matrix<numnode3Dcell,1> funct3Dcell;
  DRT::UTILS::shape_function_3D(funct3Dcell,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),CELL_DISTYPE);

  LINALG::SerialDenseMatrix xyze_cell_tmp = cell.CellNodalPosXYZ();
  const LINALG::Matrix<nsd,numnode3Dcell> xyze_cell(xyze_cell_tmp);

  posGPinXYZ.Multiply(xyze_cell,funct3Dcell);
  return;
}



/*!
  Transforms XYZ-coordinates in Xi-coordinates (3D->3D) via newton iteration
 */
template <DRT::Element::DiscretizationType DISTYPE,
          class M1, class M2>
void MapXYZToXiDomain(
    M2&                  xyze,           // element coordinates
    M1&                  posXYZDomain,   // position of point in xyz-coordinates to transform
    M1&                  posXiDomain     // position in xi-coordinates, to set
)
{
  const size_t nsd = 3;
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  //=====================================================================================================
  // implements a nonlinear solution for Xi_gp in 3D element coordinates via Newton iteration
  //=====================================================================================================
  bool newton_converged = false;

  static LINALG::Matrix<nsd,1> xi_gp;
  xi_gp.Clear();

  static LINALG::Matrix<nsd,1> incr_tmp;
  incr_tmp(0)=1.0;
  incr_tmp(1)=0.0;
  incr_tmp(2)=0.0;	//incr_tmp gets norm 1.0

  static LINALG::Matrix<nsd,1> residuum;
  residuum.Clear();

  //initialize Newton
  static LINALG::Matrix<numnode,1> shape_ele_Xi;
  static LINALG::Matrix<nsd,numnode> deriv3D_ele_Xi;

  static LINALG::Matrix<nsd,nsd> xjm_Element_XiToXYZ;		//jacobian at start vector
  static LINALG::Matrix<nsd,nsd> xji_Element_XiToXYZ;		// invers of jacobian

  int iteration_count = 0;

  const double lambda_min = 1.0/1024.0;
  const int max_iter = 15;
  // absolut incremental TOL
  const double TOL_incr = 1e-013;

  //prepare first step
  if (DISTYPE == DRT::Element::tet4){
    // start vector for newton in tet4 element (1/4*sqrt(2))
    xi_gp(0)= 0.25*1.41;
    xi_gp(1)= 0.25*1.41;
    xi_gp(2)= 0.25*1.41;
    cout << "tet4 element" << endl;
  }
  else if(DISTYPE == DRT::Element::hex8){
    // xi_gp is cleared
    //cout << "hex8 element" << endl;
  }
  // ==================== get jacobian and residuum of startvector for first iteration ====================
  DRT::UTILS::shape_function_3D(shape_ele_Xi, xi_gp(0),xi_gp(1),xi_gp(2),DISTYPE);
  DRT::UTILS::shape_function_3D_deriv1(deriv3D_ele_Xi, xi_gp(0),xi_gp(1),xi_gp(2),DISTYPE);

  xjm_Element_XiToXYZ.MultiplyNT(deriv3D_ele_Xi,xyze); // jacobian
  xji_Element_XiToXYZ.Invert(xjm_Element_XiToXYZ);	 // jacobian inverted

  residuum.Clear();
  residuum.Multiply(xyze, shape_ele_Xi);  // x_k at start position
  residuum -= posXYZDomain;        // x_k - x_gp
  residuum.Scale(-1.0);		            // negative residuum -> RHS for Newton iteration

  //======================= iteration loop for newton method with damping =================================
  while(iteration_count < max_iter)
  {
    iteration_count += 1; // count iterations
    double lambda = 1.0;  // initialize damping factor

    //================== solve Newton iteration ============================
    incr_tmp.Clear();
    incr_tmp.MultiplyTN(xji_Element_XiToXYZ,residuum); // J^(-T)*residuum
    //=====================================================================

    // abort criterion
    if (incr_tmp.Norm2()< TOL_incr) newton_converged = true;

    //==================== // natural monotonicity test =========================================
    // test || d(xi_k+1) || := || J^-T(x_k) * (x_k+1 - x_gp)|| <= (1.0 - lambda_k / 2.0)*|| J^-T * (x_k - x_gp) || = || incr_tmp ||

    static LINALG::Matrix<nsd,1> xi_gp_tmp; // needed for damping, only tmp variable
    static LINALG::Matrix<nsd,1> x_gp_tmp;
    static LINALG::Matrix<nsd,1> dxi_monot;

    // do-while-loop, natural monotonicity test is calculated at least once
    do
    {

      if(lambda < lambda_min) cout << "natural monotonicity test failed for newton iteration with standard start vector, incr.Norm2()= " << incr_tmp.Norm2() << endl;

      // update iteration
      xi_gp_tmp(0,0) = xi_gp(0,0) + lambda*incr_tmp(0,0);
      xi_gp_tmp(1,0) = xi_gp(1,0) + lambda*incr_tmp(1,0);
      xi_gp_tmp(2,0) = xi_gp(2,0) + lambda*incr_tmp(2,0);
      //cout << "incr_tmp:\n" << incr_tmp << endl;

      if(newton_converged == true) break;

      // clear tmp variables
      x_gp_tmp.Clear();
      dxi_monot.Clear();
      residuum.Clear();

      // get x_gp_tmp corresponing to xi_gp_tmp

      //new shape functions
      DRT::UTILS::shape_function_3D(shape_ele_Xi, xi_gp_tmp(0),xi_gp_tmp(1),xi_gp_tmp(2),DISTYPE);
      // current x-position
      x_gp_tmp.Multiply(xyze, shape_ele_Xi);

      //cout << "xi_gp_tmp" << xi_gp_tmp << endl;

      // get (x_k+1 - x_gp)
      residuum(0,0) = x_gp_tmp(0,0) - posXYZDomain(0,0);
      residuum(1,0) = x_gp_tmp(1,0) - posXYZDomain(1,0);
      residuum(2,0) = x_gp_tmp(2,0) - posXYZDomain(2,0);

      // get d(xi_k+1) = J^-T(x_k) * (x_k+1 - x_gp)
      dxi_monot.MultiplyTN(xji_Element_XiToXYZ, residuum);

      // reduce lambda, needed if monotonicity test failed
      lambda /= 2.0;
    } while(dxi_monot.Norm2() > ((1.0-lambda)*incr_tmp.Norm2()));

    // update xi_gp with tmp variable
    xi_gp(0,0) = xi_gp_tmp(0,0);
    xi_gp(1,0) = xi_gp_tmp(1,0);
    xi_gp(2,0) = xi_gp_tmp(2,0);

    // break newton iterations if relative increment is smaller then TOL_rel_incr
    if (newton_converged == true) break;

    // else prepare next iteration
    // get residuum of current iteration with current lambda
    residuum.Scale(-1.0);		// negative residuum -> RHS for Newton iteration (new step)


    // get derivatives of shape functions to update jacobian
    DRT::UTILS::shape_function_3D_deriv1(deriv3D_ele_Xi, xi_gp(0),xi_gp(1),xi_gp(2),DISTYPE);

    // update jacobian and inverted jacobian
    xjm_Element_XiToXYZ.MultiplyNT(deriv3D_ele_Xi,xyze);
    xji_Element_XiToXYZ.Invert(xjm_Element_XiToXYZ);
  }
  if(iteration_count == max_iter){dserror("newton not converged in max iterations");}
  posXiDomain(0) = xi_gp(0);
  posXiDomain(1) = xi_gp(1);
  posXiDomain(2) = xi_gp(2);

  //===============================================================================================
  return;
}




/*!
  Get domain integral features like jacobian, fac, shape functions evaluated for each Gaussian point
 */
template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType CELL_DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M1, class M2,class M3, class M4, class M5, class M6>
void GetDomainIntegrationFeatures(
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    M1&                            pos_eta_domain,
    M2&                            xyze,
    GEO::DomainIntCell             cell,
    double&                        det_cell_EtaToXYZ,
    M3&                            funct3Dele,
    M4&                            deriv3Dele_xyz,
    M5&                            deriv3Dele_xyz_2nd,
    M6&                            xji_ele_XiToXYZ
)
{
  const size_t NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  const int nsd = 3;

  //============================================================
  //      1. Gaussian points for cell in eta coordinates (cell)
  //============================================================
  // Gaussian points of 3D-cell in eta coordinates
  // pos_eta_domain

  //=========================================================
  //      2. Gaussian points for cell in XYZ coordinates
  //=========================================================
  // coordinates of Gaussian points in xyz-domain
  static LINALG::Matrix<nsd,1> posGPinXYZ;
  posGPinXYZ.Clear();

  // number of nodes of this cell
  const size_t numnode3Dcell = DRT::UTILS::DisTypeToNumNodePerEle<CELL_DISTYPE>::numNodePerElement;
  //cout << numnode3Dcell << endl;

  LINALG::Matrix<numnode3Dcell,1> funct3Dcell;
  DRT::UTILS::shape_function_3D(funct3Dcell,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),CELL_DISTYPE);



  LINALG::SerialDenseMatrix xyze_cell_tmp = cell.CellNodalPosXYZ();
  const LINALG::Matrix<nsd,numnode3Dcell> xyze_cell(xyze_cell_tmp);

  posGPinXYZ.Multiply(xyze_cell,funct3Dcell);

  //================================================================
  //      3. Gaussian points for cell in XI coordinates (element)
  //================================================================
  // coordinates of Gaussian points in xi_1, xi_2, xi_3 domain
  static LINALG::Matrix<nsd,1> posGPinXi;
  posGPinXi.Clear();

  // transform coordinates from xyz to xi-elementcoordinates of element -> nonlinear system
  COMBUST::MapXYZToXiDomain<DISTYPE>(xyze, posGPinXYZ, posGPinXi);


  //======================================================================================
  //      4. Jacobi-Matrix for cell in eta coordinates (cell) (Transformation eta->XYZ)
  //======================================================================================
  static LINALG::Matrix<nsd,numnode3Dcell> deriv3Dcell;
  DRT::UTILS::shape_function_3D_deriv1(deriv3Dcell,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),CELL_DISTYPE);

  // get transposed of the jacobian matrix dXYZ/dETA
  // xjm(i,j) = deriv(i,k)*xyze(j,k)
  static LINALG::Matrix<nsd,nsd> xjm_cell_EtaToXYZ;
  xjm_cell_EtaToXYZ.MultiplyNT(deriv3Dcell,xyze_cell);

  //const double det_cell_EtaToXYZ = xjm_cell_EtaToXYZ.Determinant();
  det_cell_EtaToXYZ = xjm_cell_EtaToXYZ.Determinant();

  //=================================================================================
  //      5. not enriched shape functions in XI coordinates (element) at Gaussian point
  //=================================================================================

  //--------------------------------------------------------------------------------------
  // evaluate shape functions (not enriched) and their first derivatives at this Gaussian point (XI-coordinates)
  //--------------------------------------------------------------------------------------
  // numnode is number of element nodes (DISTYPE)
  static LINALG::Matrix<nsd,NUMNODE> deriv3Dele;
  static LINALG::Matrix<6,NUMNODE> deriv3Dele_2nd;

  DRT::UTILS::shape_function_3D(funct3Dele,posGPinXi(0),posGPinXi(1),posGPinXi(2),DISTYPE);
  DRT::UTILS::shape_function_3D_deriv1(deriv3Dele,posGPinXi(0),posGPinXi(1),posGPinXi(2),DISTYPE);

  //=================================================================================
  //      6. Jacobi-Matrix for element in XI coordinates (element) (Transformation XI->XYZ)
  //=================================================================================

  // no additional entry in fac, only for transformation of derivatives of enriched shape functions
  static LINALG::Matrix<nsd,nsd> xjm_ele_XiToXYZ;
  xjm_ele_XiToXYZ.MultiplyNT(deriv3Dele,xyze);

  // inverse of jacobian
  xji_ele_XiToXYZ.Invert(xjm_ele_XiToXYZ);

  //=================================================================================
  //      7. Global derivatives of shape functions (Transformation XI->XYZ)
  //=================================================================================

  //--------------------------------------------------------------------------------------
  // compute global derivates of shape functions at this Gaussian point
  //--------------------------------------------------------------------------------------
  // derxy(i,j) = xji(i,k) * deriv(k,j)
  deriv3Dele_xyz.Multiply(xji_ele_XiToXYZ,deriv3Dele);

  // compute second global derivative
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();
  if (higher_order_ele)
  {
    static LINALG::Matrix<6,NUMNODE> deriv3Dele_2nd;
    DRT::UTILS::shape_function_3D_deriv2(deriv3Dele_2nd,posGPinXi(0),posGPinXi(1),posGPinXi(2),DISTYPE);
    DRT::UTILS::gder2<DISTYPE>(xjm_ele_XiToXYZ, deriv3Dele_xyz, deriv3Dele_2nd, xyze, deriv3Dele_xyz_2nd);
  }
  else
  {
    deriv3Dele_xyz_2nd.Clear();
  }

  return;
}


} // namespace COMBUST


#endif
#endif
