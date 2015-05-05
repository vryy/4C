/*----------------------------------------------------------------------*/
/*!
\file j_integral.H

\brief Calculates J-integral around crack tip

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef J_INTEGRAL_H_
#define J_INTEGRAL_H_

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_integration.H"

namespace DRT {
namespace CRACK {

class J_Integral {
public:

  /*!
  \brief Constructor
   */
  J_Integral( Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<std::vector<char> >& stressdata,
      Teuchos::RCP<std::vector<char> >& straindata,
      Teuchos::RCP<Epetra_Vector>& disp_col,
      Teuchos::RCP<Epetra_Vector>& iforce_col,
      std::set<int> cracknodes,
      std::set<int> boundary_nodes,
      std::vector<int>& tipnodes,
      LINALG::Matrix<3,1> normal,
      LINALG::Matrix<3,1> tangent,
      int segment_id )
  :discret_( dis ),
   gausspts_stress_( stressdata ),
   gausspts_strain_( straindata ),
   disp_col_( disp_col ),
   iforce_col_( iforce_col ),
   cracknodes_( cracknodes ),
   boundary_nodes_( boundary_nodes ),
   tipnodes_( tipnodes ),
   normal_( normal ),
   tangent_( tangent ),
   segment_id_( segment_id ),
   comm_( discret_->Comm() ),
   myrank_( comm_.MyPID() )
  {
  }

  /*!
  \brief Destructor
   */
  ~J_Integral(){};

  /*!
  \brief Performs all operations for computing J-integral
   */
  std::vector<double> compute_j_integral();

private:
  /*!
  \brief Identify all surfaces over which integration should be performed
   */
  void build_contour_surfaces( const std::set<int>& all_ele, const std::set<int>& outer_nodes );

  void clear_tip_ele( std::set<int>&allele );

  /*!
  \brief Build the data structure of structural surfaces used in the integration
   */
  void store_element_surfaces_in_outer_or_boundary( DRT::Element * ele, const std::set<int>& outer_nodes );

  /*!
  \brief Returns true if all the nodes given in "surnodes" is positioned in the boundary or in outer nodes
   */
  bool is_all_nodes_in_outer_or_boundary( const std::vector<int>& surnodes, const std::set<int>& outer_nodes );

  /*!
  \brief Return true if  all nodes are on the crack surface
   */
  bool is_all_nodes_crack_surfaces( const std::vector<int>& surnodes );

  /*!
  \brief Write GMSH output of domain used for computation of J-integral
   */
  void Write_GMSH_output_surface();

  void Write_GMSH_output_domain( std::set<int>& allele );

  void Initialize_support_function( std::set<int>& allnodes, std::set<int>& outer_nodes );

  double getMinOuterDist( std::set<int>& outer_nodes, LINALG::Matrix<2,1>& tipcoo_2d );

  double getMaxTipDist( LINALG::Matrix<2,1>& tipcoo_2d, std::set<int>& attached );

  void perform_integration_domain( std::set<int>& allele, std::set<int>&outer_nodes );

  void perform_integration_surfaces( std::set<int>& allele, std::set<int>&outer_nodes );

  DRT::UTILS::GaussRule3D getGaussRuleElement( DRT::Element * ele );

  Teuchos::RCP<DRT::Element> getCrackSurfaceEle( DRT::Element* parele );

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! Cauchy stress at Gauss points
  Teuchos::RCP<std::vector<char> >& gausspts_stress_;

  //! Green-Lagrange strain at Gauss points
  Teuchos::RCP<std::vector<char> >& gausspts_strain_;

  //! Holds the displacement in column vector
  Teuchos::RCP<Epetra_Vector> disp_col_;

  //! Holds the interface forces in column vector
  Teuchos::RCP<Epetra_Vector> iforce_col_;

  //! nodes that are on crack surface corresponding to this tip
  std::set<int> cracknodes_;

  std::set<int> boundary_nodes_;

  //! all the nodes that defines the crack tip
  std::vector<int> tipnodes_;

  //! unit normal vector at crack tip
  LINALG::Matrix<3,1> normal_;

  //! unit tangent vector at crack tip
  LINALG::Matrix<3,1> tangent_;

  //! segment id of this crack tip
  int segment_id_;

  //! communicator
  const Epetra_Comm& comm_;

  //! processor id
  int myrank_;

  //! holds the J-integral values
  std::vector<double> Jvector_;

  //! data structure that holds the surfaces that lie on J-integral contour
  // key -- sorted node ids of surface ele
  // std::pair<>.first -- parent FE ele id that has this surface
  // std::pair<>.second -- correct ordering surface nodes to get outward pointing normal
  std::map<std::set<int>,std::pair<int, std::vector<int> > > surfaces_;

  //! store surface number within given element
  std::map<std::set<int>, int> surf_no_in_ele_;

  std::map<int,double> supp_func_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* J_INTEGRAL_H_ */
