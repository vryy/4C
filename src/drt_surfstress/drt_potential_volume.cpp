/*!-------------------------------------------------------------------
\file drt_potential_volume.cpp

\brief  Class controlling volume stresses due to potential forces
        between interfaces of mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_potential_volume.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/element_normals.H"
#include <cstdlib>



/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 05/09|
 *-------------------------------------------------------------------*/
POTENTIAL::VolumePotential::VolumePotential(
    Teuchos::RCP<DRT::Discretization>     discretRCP,
    DRT::Discretization&                  discret,
    const GEO::TreeType&                  treetype):
    Potential(discretRCP, discret)
{
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*discretRCP_);
  DRT::UTILS::CollectElementsByConditionLabel(*discretRCP_, elementsByLabel_,"Potential" );
  
  //if(prob_dim_ == 3)
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype);
  //else if (prob_dim_ == 2)
  //	searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::QUADTREE));
  //else
  //    	    dserror("problem dimension not correct");

  // std::cout << "Potential manager constructor done" << endl;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::EvaluatePotential(  ParameterList& p,
                                                  RefCountPtr<Epetra_Vector> disp,
                                                  RefCountPtr<Epetra_Vector> fint,
                                                  RefCountPtr<LINALG::SparseMatrix> stiff)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  
  // update displacement for volume discretization discretRCP
  UpdateDisplacementsOfPotentialDiscretization(disp);

  EvaluatePotentialCondition(p,stiff,null,fint,null,null,"Potential");

  return;
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| for volume elements                                                |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule3D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);

  // find nodal ids influencing a given element
  const int     label     = cond->GetInt("label");
  const double  cutOff    = cond->GetDouble("cutOff");
  std::map<int,std::set<int> > potentialElementIds;
  
  /*
  for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
  {
    // compute AABB von jedem Element plus cutoff radius
    const DRT::Node* node = element->Nodes()[i];
    const LINALG::Matrix<3,1> x_node = currentpositions_.find(node->Id())->second;
    // octtree search
    treeSearchElementsInCutOffRadius(discretRCP_, currentpositions_, x_node, potentialElementIds, cutOff, label);
    // serial search
    //searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
    // searchElementsInCutOffRadius(discretRCP_, currentpositions_, x_node, potentialElementIds, cutOff);
  }
  */
  treeSearchElementsInCutOffRadius(discretRCP_, elemXAABBList_, element, potentialElementIds, cutOff, label);
  
  // initialize time variables
  const int    curvenum = cond->GetInt("curve");
  const double time     = params.get<double>("total time",-1.0);
  //const double dt     = params.get<double>("delta time",0.0);
  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
  double curvefac       = 1.0;
  // apply potential forces gradually
  if (time <= t_end)
    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // compute internal force and stiffness matrix
  // TODO if poteles empty don t do assembly
  computeFandK(element, gaussrule, potentialElementIds, lm, K_surf, F_int, cond, label, curvefac);
  // cout << "stiffness stop" << endl;
  return;
}




/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
	// initialize Lennard Jones potential constant variables
	  RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);

	  // find nodal ids influencing a given element
	  const int     label     = cond->GetInt("label");
	  const double  cutOff    = cond->GetDouble("cutOff");
	  std::map<int,std::set<int> > potentialElementIds;
	  

	  for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
	  {
	    // compute AABB von jedem Element plus cutoff radius
	    const DRT::Node* node = element->Nodes()[i];
	    const LINALG::Matrix<3,1> x_node = currentpositions_.find(node->Id())->second;
	    // octtree search
	    treeSearchElementsInCutOffRadius(discretRCP_, currentpositions_, x_node, potentialElementIds, cutOff, label);
	    // serial search
	    //searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
	    // searchElementsInCutOffRadius(discretRCP_, currentpositions_, x_node, potentialElementIds, cutOff);
	  }

	//  treeSearchElementsInCutOffRadius(discretRCP_, elemXAABBList_, element, potentialElementIds, cutOff, label);
	  
	  // initialize time variables
	  const int    curvenum = cond->GetInt("curve");
	  const double time     = params.get<double>("total time",-1.0);
	  //const double dt     = params.get<double>("delta time",0.0);
	  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
	  double curvefac       = 1.0;
	  // apply potential forces gradually
	  if (time <= t_end)
	    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);

	  // compute internal force and stiffness matrix
	  // TODO if poteles empty don t do assembly
	  computeFandK(element, gaussrule, potentialElementIds, lm, K_surf, F_int, cond, label, curvefac);
	  // cout << "stiffness stop" << endl;
	  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::UpdateDisplacementsOfPotentialDiscretization(
    Teuchos::RCP<Epetra_Vector>     idisp_solid
    )
{
  // extract displacements associated with potential elements from solid discretization
  //idisp_onproc_ = potboundary_.ExtractCondVector(idisp_solid);
  //idisp_total_->Scale(0.0);

  // import
  // int err = idisp_total_->Import((*idisp_onproc_), (*importer_),Insert);
  // if(err) dserror("Import using importer returned err=%d",err);

  currentpositions_.clear();
  
  // run over volume discretization
  {
    for (int lid = 0; lid < discretRCP_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = discretRCP_->lColNode(lid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      discretRCP_->Dof(node, lm);
   
      vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;
      DRT::UTILS::ExtractMyValues(*idisp_solid,mydisp,lm);
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentpositions_[node->Id()] = currpos;
    }
  }

  //hier muss noch eine Abfrage rein welcher Suchbaum verwendet wird, bei dem alten Suchb. wird die Abfrage neg.
  
  // reinitialize search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*discretRCP_, currentpositions_);
  if(prob_dim_ == 2)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::QUADTREE));
  else if(prob_dim_ == 3)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  else
    dserror("problem dimension not correct");
  
  //build boxes around every element
  // if abfrage
  for (int lid = 0; lid < discretRCP_->NumMyColElements(); ++lid)
  {
  	const DRT::Element*  element = discretRCP_->lColElement(lid);
  	const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element,currentpositions_));
	  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
	  GEO::checkRoughGeoType(element, xyze, eleGeoType);
	  elemXAABBList_[lid]=GEO::computeFastXAABB(element->Shape(), xyze, eleGeoType);	
  }
  return;
}


/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
/*
void POTENTIAL::Potential::treeSearchVolumParallel(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    std::map<int,LINALG::Matrix<3,2> >&         elemXAABBList,
    const DRT::Element*                         element,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius,
    const int                                   label)
{

  LINALG::Matrix<3,2> eleXAABB = elemXAABBList[element->LID()];

  // enlarge box by cut off radius
  for(int dim = 0; dim < 3; dim++)
  {
    eleXAABB(dim,0) = eleXAABB(dim,0) - radius;
    eleXAABB(dim,1) = eleXAABB(dim,1) + radius;
  }
  
  1. local search
  searchTree_->queryPotentialElements(elemXAABBList, eleXAABB, potentialElementIds, label);
  
  create empty element list for each proc
  
  for(int i = 0; i < num_procs-1; i++)
  {
  
    MPI_Barrier();
  
    2. Receive incoming eleXAABB from previous processor
    MPI_Receive();
  
    3. send eleXAABB to the next processor
    MPI_Send();
    
    
    // send elementlist to next proc
    // receive elelist from previous proc
    
    searchTree_->queryPotentialElements(elemXAABBList, eleXAABB, potentialElementIds, label);
    
    // get elements with potentialElementIds 
    // pack elements
    // store in element list
  }
  
  2. Receive incoming eleXAABB from previous processor
  MPI_Receive();

  3. send eleXAABB to the next processor
  MPI_Send();
  
   // send elementlist to next proc
   // receive elelist from previous proc

  return;
}
*/


/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (volume)       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule3D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   vector<int>&                     lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RefCountPtr<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{

  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, potElements, lm); //MINE potentialdis_
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");	
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(3,numnode);
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule3D rule_pot = DRT::UTILS::intrule3D_undefined;
         GetGaussRule3D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints3D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(3,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
         
           const int numdof = 3;
	
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 3; dim++)
			     F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
	
	         // computation of stiffness matrix (possibly with non-local values)
	         for (int inode = 0;inode < numnode; ++inode)
	           for(int dim = 0; dim < 3; dim++)
	           {
	             // k,ii
	             for (int jnode = 0; jnode < numnode; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
	                   funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
	
	             // k,ij
	             for (int jnode = 0;jnode < numnode_pot; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
	                    funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);
	
	            }
	         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule2D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   vector<int>&                     lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RefCountPtr<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(2,numnode);
    LINALG::Matrix<2,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
         GetGaussRule2D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(2,numnode_pot);
           LINALG::Matrix<2,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<2,1>  potderiv1;
           LINALG::Matrix<2,2>  potderiv2;
           
           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
            
           const int numdof = 2;
  
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 2; dim++)
               F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
  
           // computation of stiffness matrix (possibly with non-local values)
           for (int inode = 0;inode < numnode; ++inode)
             for(int dim = 0; dim < 2; dim++)
             {
               // k,ii
               for (int jnode = 0; jnode < numnode; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                     funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
  
               // k,ij
               for (int jnode = 0;jnode < numnode_pot; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

            }
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for volume elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints3D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];	
  const double e1 = intpoints.qxg[gp][1];
  const double e2 = intpoints.qxg[gp][2];
  
  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_3D(funct,e0,e1,e2,element->Shape());
  DRT::UTILS::shape_function_3D_deriv1(deriv,e0,e1,e2,element->Shape());

  
  LINALG::SerialDenseMatrix Jacobi(3,3);
  LINALG::SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);		
  LINALG::SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_,element,x,3);	
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);
  //LINALG::SerialDenseMatrix  metrictensor(3,3);
  //metrictensor.Multiply('N','T',1.0,dXYZdrs,dXYZdrs,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*(
		  				Jacobi(1,1)*Jacobi(2,2) - Jacobi(2,1)*Jacobi(1,2)  ) -
		  			  Jacobi(0,1)*(
		  				Jacobi(1,0)*Jacobi(2,2) - Jacobi(2,0)*Jacobi(1,2)  ) +
		  			  Jacobi(0,2)*(
						Jacobi(1,0)*Jacobi(2,1) - Jacobi(2,0)*Jacobi(1,1)  );                 		
                            		
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for surface elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<2,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  LINALG::SerialDenseMatrix Jacobi(2,2);
  LINALG::SerialDenseMatrix X(numnode,2);
  ReferenceConfiguration(element,X,2);
  LINALG::SerialDenseMatrix x(numnode,2);
  SpatialConfiguration(currentpositions_,element,x,2);
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*Jacobi(1,1) - Jacobi(1,0)*Jacobi(0,1);
  		  				               		
                              		
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 2; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}



/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 09/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::GetGaussRule3D(
    const DRT::Element::DiscretizationType& distype,
    DRT::UTILS::GaussRule3D&                rule_pot)
{
  // obtain Gaussrule and integration points
  switch (distype)
  {
    case DRT::Element::hex8:
      rule_pot = DRT::UTILS::intrule_hex_8point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return;
}



/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 04/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::GetGaussRule2D(
    const DRT::Element::DiscretizationType& distype,
    DRT::UTILS::GaussRule2D&                rule_pot)
{
  // obtain Gaussrule and integration points
  switch (distype)
  {
    case DRT::Element::quad4:
      rule_pot = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      rule_pot = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3:
      rule_pot = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      rule_pot = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return;
}

#endif

