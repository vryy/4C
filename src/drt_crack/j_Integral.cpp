/*----------------------------------------------------------------------*/
/*!
\file j_integral.cpp

\brief Calculates J-integral around crack tip

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "j_Integral.H"
#include "crackUtils.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/so3_material.H"
#include "../drt_so3/so_surface.H"
#include "../drt_io/io_gmsh.H"

/*----------------------------------------------------------------------------------------*
 * Main routine that perform all operations of computing J-integral               sudhakar 12/14
 * The contour J-integral is transformed into domain form after applying the divergence
 * theorem. Domain form is evaluated
 *----------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::J_Integral::compute_j_integral()
{
  surfaces_.clear();
  Jvector_.clear();

  std::set<int> allnodes, outer_nodes, all_ele;
  for( std::vector<int>::iterator it = tipnodes_.begin(); it != tipnodes_.end(); it++ )
  {
    outer_nodes.insert( *it );
    allnodes.insert( *it );
  }

  // Get all the elements that are located around a few layers around the crack tip
  // Integration will be carried out on these elements
  //TODO: read number of layers  from input file
  int num_layers = 5;
  for( int i=0; i<num_layers; i++ )
  {
    DRT::CRACK::UTILS::get_layer_ele_nodes( allnodes, outer_nodes, all_ele, discret_, myrank_, boundary_nodes_, cracknodes_ );
    LINALG::GatherAll( all_ele, comm_ );
    LINALG::GatherAll( outer_nodes, comm_ );
  }

  //build_contour_surfaces( all_ele, outer_nodes );

  // Set support function on all elements
  Initialize_support_function( allnodes, outer_nodes );

  // Erase all elements connected to the crack tip
  clear_tip_ele( all_ele );


  // Write output of domain used for integration together with the support function
  //Write_GMSH_output_domain( all_ele );

  // Evaluate domain form of J-integral
  perform_integration_domain( all_ele, outer_nodes );

  // For crack surfaces subjected to external traction (as in crack-FSI examples)
  // it is necessary to evaluate additional contour integral
  perform_integration_surfaces( all_ele, outer_nodes );

  return Jvector_;
}

/*--------------------------------------------------------------------------------------------------------------------------*
 * Initialize the scalar support function used in the calculation of domain integral                                sudhakar 12/14
 * This function gets value of 1.0 on in outer layer nodes and 0 in inner layer
 *--------------------------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::Initialize_support_function( std::set<int>& allnodes, std::set<int>& outer_nodes )
{
  if( allnodes.size() == 0 )
    dserror("No nodes found\n");

  supp_func_.clear();

  LINALG::Matrix<2,1> tipcoo_2d;
  double rmax = getMinOuterDist( outer_nodes, tipcoo_2d );
  std::set<int> attached;
  double rmin = getMaxTipDist( tipcoo_2d, attached );

  if(0)//if( fabs(rmin-rmax) < 1e-12 )
  {
    std::cout<<"WARNING!!! Only one layer of elements for J-integral?\n";
    for( std::set<int>::iterator it = allnodes.begin(); it != allnodes.end(); it++ )
    {
      int nodid = *it;
      supp_func_[nodid] = 1.0;

      // If the node falls either on boundary nodes or outer nodes, then support function is zero
      if( outer_nodes.find( nodid ) != outer_nodes.end() )
        supp_func_[nodid] = 0.0;
    }
  }

  else
  {
    for( std::set<int>::iterator it = allnodes.begin(); it != allnodes.end(); it++ )
    {
      int nodid = *it;
      if(  attached.find( nodid ) != attached.end() )
      {
        supp_func_[nodid] = 1.0;
        continue;
      }

      // If the node falls either on boundary nodes or outer nodes, then support function is zero
      if( outer_nodes.find( nodid ) != outer_nodes.end() )
      {
        supp_func_[nodid] = 0.0;
        continue;
      }

      if( discret_->HaveGlobalNode( nodid ) )
      {
        DRT::Node * nod = discret_->gNode(nodid);
        if( nod->Owner() == myrank_ )
        {
          const double * outx = nod->X();
          double tipdist = pow( (outx[0]-tipcoo_2d(0,0)), 2 ) + pow( (outx[1]-tipcoo_2d(1,0)), 2 );
          tipdist = sqrt( tipdist );

          if( tipdist < rmin )
            supp_func_[nodid] = 1.0;
          else if( tipdist > rmax or fabs(tipdist-rmax) < 1e-10 )
            supp_func_[nodid] = 0.0;
          else
          {
            double val = (rmax-tipdist)/(rmax-rmin);
            supp_func_[nodid] = val;
          }
        }
      }
    }
  }
  LINALG::GatherAll( supp_func_, comm_ );
}

/*------------------------------------------------------------------------------------------*
 * From the given set of elements, remove all the elements that are                 sudhakar 12/14
 * connected to the crack tip
 *------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::clear_tip_ele( std::set<int>&allele )
{
  std::set<int> tipele;

  // get all elements connected to tip nodes and store them in tipele
  for( std::vector<int>::iterator it = tipnodes_.begin(); it != tipnodes_.end(); it++ )
  {
    int tipid = *it;
    if( discret_->HaveGlobalNode( tipid ) )
    {
      DRT::Node * tipnode = discret_->gNode( tipid );
      if( tipnode->Owner() == myrank_ )
      {
        DRT::Element ** eleptr = tipnode->Elements();
        int numele = tipnode->NumElement();
        for( int eleno = 0; eleno < numele; eleno++ )
        {
          int eleid = eleptr[eleno]->Id();
          tipele.insert( eleid );
        }
      }
    }
  }

  LINALG::GatherAll( tipele, comm_ );

  // remove all tip elements from the given set
  std::vector<std::set<int>::iterator > deliter;
  for( std::set<int>::iterator it = allele.begin(); it != allele.end(); it++ )
  {
    if( tipele.find( *it ) != tipele.end() )
      deliter.push_back( it );
  }

  for( unsigned ii = 0; ii < deliter.size(); ii++ )
    allele.erase( deliter[ii] );

}

/*-----------------------------------------------------------------------------------------------------------------*
 * Evaluate the all the domain form of J-integral                                             sudhakar 01/15
 *-----------------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::perform_integration_domain( std::set<int>& allele, std::set<int>&outer_nodes )
{
  LINALG::Matrix<2,1> tipcoo_2d;
  double rmax = getMinOuterDist( outer_nodes, tipcoo_2d );
  std::set<int> attached;
  double rmin = getMaxTipDist( tipcoo_2d, attached );

  double j1_local = 0.0, j2_local = 0.0;

  for( std::set<int>::iterator it = allele.begin(); it != allele.end(); it++ )
  {
    int eleid = *it;
    if( discret_->HaveGlobalElement( eleid ) )
    {
      DRT::Element * ele = discret_->gElement( eleid );
      if( ele->Owner() == myrank_ )
      {
        //---
        // Get the reference coordinates, displacement and current configuration at each node of element
        //---
        int numnode = ele->NumNode();
        Epetra_SerialDenseMatrix xcurr(3,numnode);   // current  coord.
        Epetra_SerialDenseMatrix xdisp(3,numnode);   // displacements
        Epetra_SerialDenseMatrix xref(3,numnode);   // reference coordinates
        const int * nodeids = ele->NodeIds();
        DRT::Node** nodes = ele->Nodes();
        for( int pn = 0; pn < numnode; pn++ )
        {
          std::vector<double> dn = DRT::CRACK::UTILS::getDisplacementNode( discret_, nodeids[pn], disp_col_ );
          const double* x = nodes[pn]->X();

          for( int dim = 0; dim < 3; dim++ )
          {
            xdisp(dim,pn) = dn[dim];
            xcurr(dim,pn) = x[dim] + xdisp(dim,pn);
            xref(dim,pn) = x[dim];
          }
        }

        // to get shape functions and its derivatives in natural coordinate system
        LINALG::SerialDenseVector xsi(numnode);
        LINALG::SerialDenseMatrix deriv(3,numnode);

#if 1
        // Get Gauss integration rule
        DRT::UTILS::GaussRule3D gaussrule = getGaussRuleElement( ele );
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        const int ngp = intpoints.nquad;

        for (int gp=0; gp<ngp; gp++)
        {
          // Gauss points in  natural coordinate system
          const double e0 = intpoints.qxg[gp][0];
          const double e1 = intpoints.qxg[gp][1];
          const double e2 = intpoints.qxg[gp][2];
          const double wei = intpoints.qwgt[gp];
#else
        int gp_no = -1;
        DRT::UTILS::GaussIntegration gi_temp = DRT::UTILS::GaussIntegration( ele->Shape(), 5 );
        for ( DRT::UTILS::GaussIntegration::iterator iquad=gi_temp.begin(); iquad!=gi_temp.end(); ++iquad )
        {
          gp_no++;
          const LINALG::Matrix<3,1> eta( iquad.Point() );
          const double e0 = eta(0,0);
          const double e1 = eta(1,0);
          const double e2 = eta(2,0);
          const double wei = iquad.Weight();
#endif

          // get shape functions and derivatives of surface
          DRT::UTILS::shape_function_3D(xsi,e0,e1,e2,ele->Shape());
          DRT::UTILS::shape_function_3D_deriv1(deriv,e0,e1,e2,ele->Shape());

          // Find Jacobian and its inverse
          LINALG::SerialDenseMatrix xjm(3,3);
          LINALG::SerialDenseMatrix xji(3,3);
          xjm.Multiply('N','T',1.0,deriv,xref,0.0);

          double det = xjm(0,0)*(xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2)) + xjm(0,1)*(- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2)) + xjm(0,2)*(xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1));
          xji(0,0) = (xjm(1,1)*xjm(2,2)-xjm(2,1)*xjm(1,2))/det;
          xji(1,1) = (xjm(0,0)*xjm(2,2)-xjm(2,0)*xjm(0,2))/det;
          xji(2,2) = (xjm(0,0)*xjm(1,1)-xjm(1,0)*xjm(0,1))/det;
          xji(0,1) = (xjm(2,1)*xjm(0,2)-xjm(0,1)*xjm(2,2))/det;
          xji(0,2) = (xjm(0,1)*xjm(1,2)-xjm(1,1)*xjm(0,2))/det;
          xji(1,0) = (xjm(2,0)*xjm(1,2)-xjm(1,0)*xjm(2,2))/det;
          xji(1,2) = (xjm(1,0)*xjm(0,2)-xjm(0,0)*xjm(1,2))/det;
          xji(2,0) = (xjm(1,0)*xjm(2,1)-xjm(2,0)*xjm(1,1))/det;
          xji(2,1) = (xjm(2,0)*xjm(0,1)-xjm(0,0)*xjm(2,1))/det;

          LINALG::SerialDenseMatrix pderxy(3,numnode);
          LINALG::SerialDenseMatrix defgrad(3,3);      // deformation gradient
          LINALG::SerialDenseMatrix dispgrad(3,3);      // usual displacement gradient
          pderxy.Multiply('N','N',1.0,xji,deriv,0.0);
          dispgrad.Multiply('N','T',1.0,xdisp,pderxy,0.0);

          //---
          // compute gradiens of support functions
          //---
          // Get support functions at nodal points
          LINALG::SerialDenseVector supp(numnode);
          for( int pn = 0; pn < numnode; pn++ )
          {
            int nid = nodeids[pn];
            supp(pn) = supp_func_[nid];
          }
          LINALG::SerialDenseVector suppgrad(3);        // gradient of support functions
          suppgrad.Multiply('N','N',1.0,pderxy,supp,0.0);
#if 1
          double glo_gaus[3] = {0.0, 0.0, 0.0};
          for( int nn = 0; nn < numnode; nn++)
          {
            for( int dim = 0; dim < 3; dim++ )
              glo_gaus[dim] += xref(dim,nn) *  xsi(nn);
          }

          double tipdist = pow( (glo_gaus[0]-tipcoo_2d(0,0)), 2 ) + pow( (glo_gaus[1]-tipcoo_2d(1,0)), 2 );
          tipdist = sqrt( tipdist );

          suppgrad(0) = 0.0; suppgrad(1) = 0.0; suppgrad(2) = 0.0;
          if( tipdist < rmax )
          {
            suppgrad(0) = -(glo_gaus[0]-tipcoo_2d(0,0)) /( tipdist * (rmax-rmin) );
            suppgrad(1) = -(glo_gaus[1]-tipcoo_2d(1,0)) /( tipdist * (rmax-rmin) );
          }
#endif

          for( int ii=0; ii<3; ii++ )
          {
            for( int jj=0; jj<3; jj++ )
              defgrad(ii,jj) = dispgrad(ii,jj);
          }

          for( int dim = 0; dim < 3; dim++ )
            defgrad(dim,dim) = dispgrad(dim,dim) + 1.0;

          //---
          // Get strain energy
          //---
          LINALG::SerialDenseMatrix cauchygreen(3,3);
          cauchygreen.Multiply('T','N',1.0,defgrad,defgrad,0.0);
          LINALG::Matrix<6,1> glstrain;
          glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
          glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
          glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
          glstrain(3) = cauchygreen(0,1);
          glstrain(4) = cauchygreen(1,2);
          glstrain(5) = cauchygreen(2,0);

          double psi = 0.0;
          Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
          so3mat->StrainEnergy(glstrain,psi,ele->Id());

          //---
          // Get second Piola Kirchhoff stress tensor
          //---
          Teuchos::ParameterList params;
          params.set<int>("gp",gp);
          //params.set<int>("gp",gp_no);

          LINALG::Matrix<6,1> stress2pk;
          LINALG::Matrix<6,6> cmat;
          params.set<int>("iostress", INPAR::STR::stress_2pk); // needed for activefiber material; if output is requested only active stresses are written
          LINALG::Matrix<3,3> defgrad_temp(defgrad.A(),true);
          so3mat->Evaluate(&defgrad_temp,&glstrain,params,&stress2pk,&cmat,ele->Id());

          //---
          // Find first Piola Kirchhoff stress tensor
          // P1 = F(P2)
          //---
          LINALG::SerialDenseMatrix first_piola(3,3);
          for( int dim=0; dim<3; dim++ )
          {
            first_piola(dim,0) = defgrad(dim,0)*stress2pk(0,0) + defgrad(dim,1)*stress2pk(3,0) + defgrad(dim,2)*stress2pk(5,0);
            first_piola(dim,1) = defgrad(dim,0)*stress2pk(3,0) + defgrad(dim,1)*stress2pk(1,0) + defgrad(dim,2)*stress2pk(4,0);
            first_piola(dim,2) = defgrad(dim,0)*stress2pk(5,0) + defgrad(dim,1)*stress2pk(4,0) + defgrad(dim,2)*stress2pk(2,0);
          }

          //---
          // Find first part of J-integral
          //---
          LINALG::SerialDenseMatrix jfirst(3,3);
          jfirst.Multiply('T','N',1.0,dispgrad,first_piola,0.0);

          //---
          // Find second part of J-integral
          //---
          for( int dim = 0; dim < 3; dim++ )
            jfirst(dim,dim) -= psi;

          LINALG::SerialDenseVector j_xyz(3);
          j_xyz.Multiply('N','N',1.0,jfirst,suppgrad,0.0);

          for( int dim=0; dim<3; dim++ )
            j_xyz(dim) = j_xyz(dim) * wei * det;

          for( int dim=0; dim<3; dim++ )
          {
            j1_local += j_xyz(dim) * normal_(dim,0);
            j2_local += j_xyz(dim) * tangent_(dim,0);
          }
        }
      }
    }
  }

  double i1=0.0, i2=0.0;
  comm_.SumAll( &j1_local, &i1, 1 );
  comm_.SumAll( &j2_local, &i2, 1 );

  Jvector_.push_back( i1 );
  Jvector_.push_back( i2 );
}

/*----------------------------------------------------------------------------------------------------------*
 * Identify all surfaces over which integration should be performed                                 sudhakar 10/14
 * (Currently not used. This was the initial attempt when trying to compute the J-integral in
 * its original contour form without transforming to domain form.)
 *----------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::build_contour_surfaces( const std::set<int>& all_ele, const std::set<int>& outer_nodes )
{
  if( all_ele.size() == 0 )
    std::cout<<"No elements found for this crack tip\n";

  for( std::set<int>::iterator it = all_ele.begin(); it != all_ele.end(); it++ )
  {
    int eleid = *it;
    if( discret_->HaveGlobalElement( eleid ) )
    {
      DRT::Element * ele = discret_->gElement( eleid );

      if( ele->Owner() == myrank_ )
      {
        store_element_surfaces_in_outer_or_boundary( ele, outer_nodes );
      }
    }
  }

  // check if surfaces_ are stored at least in one processor
  // otherwise something is wrong
  int sur_local = surfaces_.size();
  int sur_global = 0;
  comm_.SumAll( &sur_local, &sur_global, 1 );
  if( sur_global == 0 )
    dserror( "not a single surface contour element found \n" );
}

/*----------------------------------------------------------------------------------------------------------*
 * Build the data structure of structural surfaces used in the integration                          sudhakar 10/14
 * (Currently not used. This was the initial attempt when trying to compute the J-integral in
 * its original contour form without transforming to domain form.)
 *----------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::store_element_surfaces_in_outer_or_boundary( DRT::Element * ele,
                                                                          const std::set<int>& outer_nodes )
{
  const DRT::Element::DiscretizationType distype = ele->Shape();
  std::vector< std::vector<int> > allsurf_nodes = DRT::UTILS::getEleNodeNumberingSurfaces(distype);

  unsigned int numsurf = ele->NumSurface();
  if( numsurf != allsurf_nodes.size() )
    dserror( "Some surface nodes are missing?\n" );

  const int * elenodes = ele->NodeIds();

  for( unsigned surno = 0; surno < allsurf_nodes.size(); surno++ )
  {
    std::vector<int> surnodes_index = allsurf_nodes[surno];

    std::vector<int> surnodes;
    for( unsigned ind = 0; ind < surnodes_index.size(); ind++ )
      surnodes.push_back( elenodes[surnodes_index[ind]] );


    if( is_all_nodes_in_outer_or_boundary( surnodes, outer_nodes ) )
    {
      std::set<int> sur_sorted;
      sur_sorted.insert( surnodes.begin(), surnodes.end() );
      surfaces_[sur_sorted] = std::make_pair( ele->Id(), surnodes );
      surf_no_in_ele_[sur_sorted] = surno;
    }
  }
}

/*------------------------------------------------------------------------------------------------------*
 * Returns true if all the nodes passed as input are located in the crack surfaces              sudhakar 12/14
 *------------------------------------------------------------------------------------------------------*/
bool DRT::CRACK::J_Integral::is_all_nodes_crack_surfaces( const std::vector<int>& surnodes )
{
  for( unsigned nodeno = 0; nodeno < surnodes.size(); nodeno++ )
  {
    if( cracknodes_.find( surnodes[nodeno] ) == cracknodes_.end() )
      return false;
  }
  return true;
}

/*----------------------------------------------------------------------------------------------------------*
 * Returns true if all the nodes given in "surnodes" are in outer surface                           sudhakar 10/14
 * or fall over the boundary
 *----------------------------------------------------------------------------------------------------------*/
bool DRT::CRACK::J_Integral::is_all_nodes_in_outer_or_boundary( const std::vector<int>& surnodes,
                                                                const std::set<int>& outer_nodes )
{
  bool is_contour_surface = true;

  //---------------
  // Case 1: All nodes are lying on crack surfaces
  //---------------
  // If only few nodes are on crack surface, then this does not satisfy other two conditions as well
  /*int oncrack = 0;
  for( unsigned nodeno = 0; nodeno < surnodes.size(); nodeno++ )
  {
    if( cracknodes_.find( surnodes[nodeno] ) == cracknodes_.end() )
      is_contour_surface = false;
    else
      oncrack++;
  }
  if( is_contour_surface )
    return true;
  else if( oncrack > 0 )
    return false;
  else
    is_contour_surface = true;*/

  //---------------
  // Case 2: All nodes are lying on outer nodes
  //---------------
  for( unsigned nodeno = 0; nodeno < surnodes.size(); nodeno++ )
  {
    if( outer_nodes.find( surnodes[nodeno] ) == outer_nodes.end() )
    {
      is_contour_surface = false;
      break;
    }
  }
  if( is_contour_surface )
    return true;
  else
    is_contour_surface = true;

  return false; //???

  //---------------
  // Case 3: All nodes are lying on boundary nodes
  //---------------
  for( unsigned nodeno = 0; nodeno < surnodes.size(); nodeno++ )
  {
    if( boundary_nodes_.find( surnodes[nodeno] ) == boundary_nodes_.end() )
    {
      is_contour_surface = false;
      break;
    }
  }

  return is_contour_surface;
}

/*------------------------------------------------------------------------------------------------------------*
 * Of all the nodes located on the outer surface of the J-integral domain, get the node               sudhakar 01/15
 * which is at the minimum distance from the crack tip node
 *------------------------------------------------------------------------------------------------------------*/
double DRT::CRACK::J_Integral::getMinOuterDist( std::set<int>& outer_nodes, LINALG::Matrix<2,1>& tipcoo_2d )
{
  //---
  // Get tip coordinates in 2D on all processors
  //---
  int tipid = tipnodes_[0];
  int lmaster = 0, gmaster = 0;
  if( discret_->HaveGlobalNode( tipid ) )
  {
    DRT::Node * node = discret_->gNode( tipid );
    if( node->Owner() == myrank_ )
    {
      lmaster = myrank_;
      const double *tipx = node->X();
      for( int dim=0; dim<2; dim++ )
        tipcoo_2d(dim,0) = tipx[dim];
    }
  }
  comm_.SumAll( &lmaster, &gmaster, 1 );
  comm_.Broadcast( &tipcoo_2d(0,0), 3, gmaster );

  //---
  // Get minimum outer node distance from tip node
  //---
  std::map<int,double> out_dist; // maximum outer distance at each processor
  double temp_min_dist = 1000000.0;
  // Get the farthest outer node to form an analytical function
  for( std::set<int>::iterator it = outer_nodes.begin(); it != outer_nodes.end(); it++ )
  {
    int nodid = *it;
    if( discret_->HaveGlobalNode( nodid ) )
    {
      DRT::Node * node = discret_->gNode( nodid );
      if( node->Owner() == myrank_ )
      {
        const double * outx = node->X();
        double tipdist = pow( (outx[0]-tipcoo_2d(0,0)), 2 ) + pow( (outx[1]-tipcoo_2d(1,0)), 2 );
        tipdist = sqrt( tipdist );
        if( tipdist < temp_min_dist )
          temp_min_dist = tipdist;
      }
    }
  }
  out_dist[myrank_] = temp_min_dist;
  LINALG::GatherAll( out_dist, comm_ );
  // Go through all processor values to get the minimum;
  double min_outer_dist = 1000000000000.0;
  for( std::map<int,double>::iterator it = out_dist.begin(); it != out_dist.end(); it++)
  {
    if( (it->second) < min_outer_dist )
      min_outer_dist = it->second;
  }
  return min_outer_dist;
}

/*------------------------------------------------------------------------------------------------------------*
 * Of all the nodes located on the elements attached to the crack tip node, get the node              sudhakar 01/15
 * which is at the maximum distance from the crack tip node
 *------------------------------------------------------------------------------------------------------------*/
double DRT::CRACK::J_Integral::getMaxTipDist( LINALG::Matrix<2,1>& tipcoo_2d, std::set<int>& attached )
{
  attached.clear();
  int tipid = tipnodes_[0];
  if( discret_->HaveGlobalNode( tipid ) )
  {
    DRT::Node * node = discret_->gNode( tipid );
    if( node->Owner() == myrank_ )
    {
      DRT::Element** eles = node->Elements();
      const int numeles = node->NumElement();
      for( int j=0; j<numeles; j++ )
      {
        const DRT::Element* e = eles[j];
        const int* nids = e->NodeIds();
        const int numnodes = e->NumNode();
        for( int k=0; k<numnodes; k++ )
          attached.insert( nids[k] );
      }
    }
  }
  LINALG::GatherAll( attached, comm_ );

  //---
  // Get maximum distance from tip node
  //---
  std::map<int,double> out_dist; // maximum outer distance at each processor
  double temp_max_dist = 0.0;
  // Get the farthest outer node to form an analytical function
  for( std::set<int>::iterator it = attached.begin(); it != attached.end(); it++ )
  {
    int nodid = *it;
    if( discret_->HaveGlobalNode( nodid ) )
    {
      DRT::Node * node = discret_->gNode( nodid );
      if( node->Owner() == myrank_ )
      {
        const double * outx = node->X();
        double tipdist = pow( (outx[0]-tipcoo_2d(0,0)), 2 ) + pow( (outx[1]-tipcoo_2d(1,0)), 2 );
        tipdist = sqrt( tipdist );
        if( tipdist > temp_max_dist )
          temp_max_dist = tipdist;
      }
    }
  }
  out_dist[myrank_] = temp_max_dist;
  LINALG::GatherAll( out_dist, comm_ );
  // Go through all processor values to get the maximum;
  double min_outer_dist = 0.0;
  for( std::map<int,double>::iterator it = out_dist.begin(); it != out_dist.end(); it++)
  {
    if( (it->second) > min_outer_dist )
      min_outer_dist = it->second;
  }
  return min_outer_dist;
}


/*---------------------------------------------------------------------------------------------------------------*
 * For simulation that involve application of external traction forces on the crack surfaces,             sudhakar 01/15
 * an additional surface integral is evaluated over the crack surfaces.
 *---------------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::perform_integration_surfaces( std::set<int>& allele, std::set<int>&outer_nodes )
{

  if ( DRT::Problem::Instance()->ProblemType() != prb_fsi_crack )
    return;

  double j1_local = 0.0, j2_local = 0.0;
  for( std::set<int>::iterator it = allele.begin(); it != allele.end(); it++ )
  {
    int eleid = *it;
    if( discret_->HaveGlobalElement( eleid ) )
    {
      DRT::Element * ele = discret_->gElement( eleid );
      if( ele->Owner() == myrank_ )
      {
        Teuchos::RCP<DRT::Element> surfele = getCrackSurfaceEle( ele );

        if( surfele != Teuchos::null )
        {
          Teuchos::RCP<DRT::ELEMENTS::StructuralSurface> surface_ele = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::StructuralSurface>(surfele);
          std::vector<double> J = surface_ele->Evaluate_J_integral_crack( discret_, disp_col_, supp_func_, normal_, tangent_ );

          j1_local += J[0];
          j2_local += J[1];
        }
      }
    }
  }

  double i1=0.0, i2=0.0;
  comm_.SumAll( &j1_local, &i1, 1 );
  comm_.SumAll( &j2_local, &i2, 1 );

  if( myrank_ == 0 )
  {
    std::cout<<"-------------Jvector before = "<<Jvector_[0]<<" "<<Jvector_[1]<<"\n";
    std::cout<<"adding to j-vector = "<<i1<<" "<<i2<<"\n";
  }

  Jvector_[0] -= i1;
  Jvector_[1] -= i2;

  if( myrank_ == 0 )
  {
    std::cout<<"-------------Jvector after = "<<Jvector_[0]<<" "<<Jvector_[1]<<"\n";
  }
}

/*----------------------------------------------------------------------------------------------------------*
 * For the given  element, get the surface which is located on the crack surface                    sudhakar 01/15
 * It is assumed in this work that only one surface per element can fall on the crack surface
 *----------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::CRACK::J_Integral::getCrackSurfaceEle( DRT::Element* parele )
{
  const DRT::Element::DiscretizationType distype = parele->Shape();
  std::vector< std::vector<int> > allsurf_nodes = DRT::UTILS::getEleNodeNumberingSurfaces(distype);

  unsigned int numsurf = parele->NumSurface();
  if( numsurf != allsurf_nodes.size() )
    dserror( "Some surface nodes are missing?\n" );

  const int * elenodes = parele->NodeIds();

  for( unsigned surno = 0; surno < allsurf_nodes.size(); surno++ )
  {
    std::vector<int> surnodes_index = allsurf_nodes[surno];

    std::vector<int> surnodeids;
    for( unsigned ind = 0; ind < surnodes_index.size(); ind++ )
      surnodeids.push_back( elenodes[surnodes_index[ind]] );

    //TODO: It is assumed here that there exists only one crack surface element per parent element
    // Needs to be checked thoroughly
    if ( is_all_nodes_crack_surfaces( surnodeids ) )
    {
      std::vector<DRT::Node*> surnodes;
      for( unsigned nodno = 0; nodno < surnodeids.size(); nodno++ )
      {
        int gid = surnodeids[nodno];
        surnodes.push_back( discret_->gNode( gid ) );
      }

      Teuchos::RCP<DRT::Element> surfele = Teuchos::rcp(new DRT::ELEMENTS::StructuralSurface( 0, myrank_, surnodes.size(),
                                                                                              &surnodeids[0], &surnodes[0], parele, surno ));

      if( surfele == Teuchos::null )
        dserror("can't create a surface element for the considered parent element\n");

      return surfele;
    }
  }
  return Teuchos::null;
}



/*--------------------------------------------------------------------------------------*
 * Write domain used for J-integral computation  into GMSH output               sudhakar 10/14
 * For debugging resons
 *--------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::Write_GMSH_output_surface()
{
  static int step = 0;

  std::string str1 = "Jsurface";
  std::ostringstream ostr;
  ostr<<str1<<segment_id_;

  std::string jfile = ostr.str();
  const std::string filename = IO::GMSH::GetFileName(jfile, step, false, myrank_);

  std::ofstream fj(filename.c_str());
  // add 'View' to Gmsh postprocessing file
  fj << "View \" " << "J_Surface \" {" << std::endl;

  for( std::map<std::set<int>,std::pair<int, std::vector<int> > >::iterator it = surfaces_.begin();
                                                                            it != surfaces_.end(); it++ )
  {
    int par_ele_id = it->second.first;
    DRT::Element * parele = discret_->gElement( par_ele_id );
    if( discret_->HaveGlobalElement( par_ele_id ) )
    {
      if( parele->Owner() == myrank_ )
      {
        int sur_no_ele = surf_no_in_ele_[it->first];

        std::vector<int> surfnodeids = it->second.second;
        int numsurfnodes = surfnodeids.size();

        if( numsurfnodes == 3 )
          fj<<"ST(";
        else if( numsurfnodes == 4 )
          fj<<"SQ(";
        else
          dserror("Surface of an element should be either Tri3 or Quad4\n");

        for( int i=0; i<numsurfnodes; i++ )
        {
          int nodeid = surfnodeids[i];
          DRT::Node* node = discret_->gNode( nodeid );
          const double * coo = node->X();

          std::vector<double> disp_node = DRT::CRACK::UTILS::getDisplacementNode( discret_, nodeid, disp_col_ );
          for( int dim=0; dim<3; dim++ )
            disp_node[dim] += coo[dim];
          fj<<disp_node[0]<<","<<disp_node[1]<<","<<disp_node[2];
          if( i != numsurfnodes-1 )
            fj<<",";
        }
        if( numsurfnodes == 3 )
          fj<<"){"<<sur_no_ele<<","<<sur_no_ele<<","<<sur_no_ele<<"};\n";
        else if( numsurfnodes == 4 )
          fj<<"){"<<sur_no_ele<<","<<sur_no_ele<<","<<sur_no_ele<<","<<sur_no_ele<<"};\n";

      }
    }
  }

  fj << "};" << std::endl;
  fj.close();
  step++;
}

/*------------------------------------------------------------------------------------------------*
 * Write domain used for interaction integral computation  into GMSH output               sudhakar 10/14
 * For debugging resons
 *------------------------------------------------------------------------------------------------*/
void DRT::CRACK::J_Integral::Write_GMSH_output_domain( std::set<int>& allele )
{
  static int step = 0;

  std::string str1 = "Intdomain";
  std::ostringstream ostr;
  ostr<<str1<<segment_id_;

  std::string jfile = ostr.str();
  const std::string filename = IO::GMSH::GetFileName(jfile, step, false, myrank_);

  std::ofstream fj(filename.c_str());
  // add 'View' to Gmsh postprocessing file
  fj << "View \" " << "Jdomain \" {" << std::endl;

  for( std::set<int>::iterator it = allele.begin(); it != allele.end(); it++ )
  {
    if( discret_->HaveGlobalElement( *it ) )
    {
      DRT::Element* ele = discret_->gElement( *it );
      if( ele->Owner() == myrank_ )
      {
        const int * nodeids = ele->NodeIds();
        const int numnodes = ele->NumNode();
        DRT::Node ** nodes = ele->Nodes();

        if( numnodes == 8 )
          fj<<"SH(";
        else if ( numnodes == 6 )
          fj<<"SI(";
        else
          dserror( "Elements must be a Hex8 or Wedge6\n" );

        for( int nodno = 0; nodno < numnodes; nodno++ )
        {
          DRT::Node * node = nodes[nodno];
          const double * coo = node->X();

          std::vector<double> disp_node = DRT::CRACK::UTILS::getDisplacementNode( discret_, node, disp_col_ );
          for( int dim=0; dim<3; dim++ )
            disp_node[dim] += coo[dim];
          fj<<disp_node[0]<<","<<disp_node[1]<<","<<disp_node[2];
          if( nodno != numnodes-1 )
            fj<<",";
        }

        fj<<"){";
        for( int nodno = 0; nodno < numnodes; nodno++ )
        {
          int nodeid = nodeids[nodno];
          fj<<supp_func_[nodeid];
          if( nodno != numnodes-1 )
            fj<<",";
        }
        fj<<"};\n";
      }
    }
  }

  fj << "};" << std::endl;
  fj.close();
  step++;
}

/*--------------------------------------------------------------------------------------------------*
 * Obtain Gaussian quadrature of the given element                                          sudhakar 02/15
 *--------------------------------------------------------------------------------------------------*/
DRT::UTILS::GaussRule3D DRT::CRACK::J_Integral::getGaussRuleElement( DRT::Element * ele )
{
  DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;

  switch( ele->Shape() )
  {
  case DRT::Element::hex8:
  {
    gaussrule = DRT::UTILS::intrule_hex_8point;
    //gaussrule = DRT::UTILS::intrule_hex_18point;
    break;
  }
  case DRT::Element::wedge6:
  {
    gaussrule = DRT::UTILS::intrule_wedge_6point;
    //gaussrule = DRT::UTILS::intrule_wedge_9point;
    break;
  }
  default:
  {
    dserror( "Crack modules only only for Hex8 and Wedge6 elements\n" );
    break;
  }
  }

  return gaussrule;
}

/*--------------------------------------------------------------------------------------*
 * Perform actual integration over the surfaces                                 sudhakar 10/14
 *--------------------------------------------------------------------------------------*/
/*void DRT::CRACK::J_Integral::perform_integration_surfaces()
{
  //----
  // Interpolate stresses and strains at FE nodes from their Gauss point values
  //----
  //const Epetra_Map* noderowmap = discret_->NodeRowMap();

  double Jvec_loc[2];
  Jvec_loc[0] = 0.0; Jvec_loc[1] = 0.0;
  for( std::map<std::set<int>,std::pair<int, std::vector<int> > >::iterator it = surfaces_.begin();
                                                                            it != surfaces_.end(); it++ )
  {
    int par_ele_id = it->second.first;
    if( discret_->HaveGlobalElement( par_ele_id ) )
    {
      DRT::Element * parele = discret_->gElement( par_ele_id );
      if( parele->Owner() == myrank_ )
      {
        std::vector<int> surfnodeids = it->second.second;
        std::vector<DRT::Node*> surfnodes;

        for( unsigned nodno = 0; nodno < surfnodeids.size(); nodno++ )
        {
          int gid = surfnodeids[nodno];
          surfnodes.push_back( discret_->gNode( gid ) );
        }

        int surid_in_ele = surf_no_in_ele_[it->first];
        // now we have all quantities required to compute J-integral
        // we create temporary surface elements here
        Teuchos::RCP<DRT::Element> surfele = Teuchos::rcp(new DRT::ELEMENTS::StructuralSurface( 0, myrank_, surfnodeids.size(),
                                                                                                &surfnodeids[0], &surfnodes[0], parele, surid_in_ele ));
        Teuchos::RCP<DRT::ELEMENTS::StructuralSurface> surface_ele = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::StructuralSurface>(surfele);
        std::vector<double> J = surface_ele->Evaluate_J_integral_crack( discret_, disp_col_, normal_, tangent_ );

        Jvec_loc[0] += J[0];
        Jvec_loc[1] += J[1];
      }
    }
  }

  Jvector_.push_back(0.0);
  Jvector_.push_back(0.0);
  comm_.SumAll( &Jvec_loc[0], &Jvector_[0], 2 );
}*/
