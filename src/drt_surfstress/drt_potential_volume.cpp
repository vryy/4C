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
    Teuchos::RCP<DRT::Discretization> discretRCP,
    DRT::Discretization& discret):
    Potential(discretRCP, discret)
{
  // create discretization from potential boundary elements and distribute them
  // on every processor
/*  
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("Potential");

  if(prob_dim_ == 2)
    potentialdis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE3", conditions_to_copy);
  else if(prob_dim_ == 3)
    potentialdis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE3", conditions_to_copy);
  else
    dserror("problem dimension not correct");
  dsassert(potentialdis_->NumGlobalNodes() > 0, "empty discretization detected. Potential conditions applied?");

  // set new dof set
  RCP<POTENTIAL::PotentialDofSet> pdofset = rcp(new POTENTIAL::PotentialDofSet(discretRCP_));
  (*potentialdis_).ReplaceDofSet(pdofset);
  (*potentialdis_).FillComplete(false, false, false);

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *(potentialdis_->NodeRowMap());
  const Epetra_Map elemrowmap = *(potentialdis_->ElementRowMap());

  // put all boundary nodes and elements onto all processors
  const Teuchos::RCP<const Epetra_Map> nodecolmap = LINALG::AllreduceEMap(*(potentialdis_->NodeRowMap()));
  const Teuchos::RCP<const Epetra_Map> elemcolmap = LINALG::AllreduceEMap(*(potentialdis_->ElementRowMap()));

  // redistribute nodes and elements to column (ghost) map
  potentialdis_->ExportColumnNodes(*nodecolmap);
  potentialdis_->ExportColumnElements(*elemcolmap);

  // Now we are done. :)
  // fill complete calls the assgin degrees of fredom method
  // provided by PotentialDofSet class
  const int err = potentialdis_->FillComplete();
  if (err) dserror("FillComplete() returned err=%d",err);

  // split dof vector of soliddiscretization into a dof vector of potential boundary condition and
  // remaining dofs
  DRT::UTILS::SetupNDimExtractor(*discretRCP_ ,"Potential", potboundary_);

  // create potential surface dof row map using the solid parallel distribution
  // for all potential surface elements belonging to the solid discretization on a single
  // processor
  const Teuchos::RCP< const Epetra_Map > potsurface_condmap_onOneProc = potboundary_.CondMap();
  if(! potsurface_condmap_onOneProc->UniqueGIDs() ) dserror("cond_map not unique");
  if(! (potsurface_condmap_onOneProc->PointSameAs(*(*potentialdis_).DofRowMap())))
    dserror("maps are not point equal");
  if(! (potsurface_condmap_onOneProc->SameAs(*(*potentialdis_).DofRowMap())))
    dserror("maps are not equal");
  // create disp vector for boundary discretization corresponding to the
  // potential condition elements stored on one proc
  idisp_onproc_    = LINALG::CreateVector(*potsurface_condmap_onOneProc,true);

  // create potential surface dof col map based on the boundayr dis
  // that is distributed redundantly on all processors
  const Epetra_Map* potsurface_dofcolmap_total = potentialdis_->DofColMap();
  // create disp vector for boundary discretization redundant on all procs
  idisp_total_    = LINALG::CreateVector(*potsurface_dofcolmap_total,true);

  // create importer
  importer_ = rcp(new Epetra_Import(idisp_total_->Map(),idisp_onproc_->Map()));

  // set up tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*potentialdis_);
  DRT::UTILS::CollectElementsByConditionLabel(*potentialdis_, elementsByLabel_,"Potential" );
  if(prob_dim_ == 2)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::QUADTREE));
  else if(prob_dim_ == 3)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  else
    dserror("problem dimension not correct");


  // std::cout << "Potential manager constructor done" << endl;
   */
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
    const DRT::Node* node = element->Nodes()[i];
    //const LINALG::Matrix<3,1> x_node = currentpositions_.find(node->Id())->second;
    // octtree search
    //treeSearchElementsInCutOffRadius(potentialdis_, currentpositions_, x_node, potentialElementIds, cutOff, label);
    
    // serial search neu für volume enlemente
    //searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
  }

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
| for line elements                  (nicht beachten)                |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule1D&  gaussrule,
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
    const DRT::Node* node = element->Nodes()[i];
    const LINALG::Matrix<3,1> x_node = currentpositions_.find(node->Id())->second;
    // octtree search
    treeSearchElementsInCutOffRadius(potentialdis_, currentpositions_, x_node, potentialElementIds, cutOff, label);
    // serial search
    // searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
  }

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
    for (int lid = 0; lid < potentialdis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = potentialdis_->lColNode(lid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      potentialdis_->Dof(node, lm);
      vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;
      DRT::UTILS::ExtractMyValues(*idisp_total_,mydisp,lm);
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentpositions_[node->Id()] = currpos;
    }
  }

  // reinitialize search tree
  /*const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*potentialdis_, currentpositions_);
  if(prob_dim_ == 2)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::QUADTREE));
  else if(prob_dim_ == 3)
    searchTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  else
    dserror("problem dimension not correct");
    */
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
  CollectLmcol(potentialdis_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = 25000000;
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
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);

    // compute normal n_gp to act ele in x_gp
    LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = potentialdis_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner);

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
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // compute normal of surface element in x_pot_gp
           LINALG::Matrix<3,1> n_pot_gp = ComputeNormalInGP(element_pot, x_pot_gp);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           bool contact = true;
           bool valid = false;
           if(contact)
           {
             double welldepth = 0.02;
             LINALG::Matrix<3,1> x_gp_con(true);
             x_gp_con = x_gp;

             //Add: this = scalarThis * this + scalarOther * other.
             x_gp_con.Update( (-1.0)*welldepth, n_gp, 1.0);
             x_pot_gp.Update( (-1.0)*welldepth, n_pot_gp, 1.0);

             EvaluatePotentialfromCondition(cond, x_gp_con, x_pot_gp, potderiv1, potderiv2);
             valid = DetermineValidContribution(x_gp_con, x_pot_gp, n_gp, n_pot_gp);
           }
           else
           {
             EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
             valid = DetermineValidContribution(x_gp, x_pot_gp, n_gp, n_pot_gp);
           }

          // if valid contribution
          if(valid)
          {
             //cout << "potderiv1 = " << potderiv1 << endl;
             //cout << "potderiv2 = " << potderiv2 << endl;

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
             }// valid contribution
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix for line        |
| elements (2D problems)                                             |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule1D&   gaussrule,
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
  CollectLmcol(potentialdis_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = 25000000;
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(1,numnode);
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = potentialdis_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule1D rule_pot = DRT::UTILS::intrule1D_undefined;
         switch (element_pot->Shape())
         {
            case DRT::Element::line2:
               rule_pot = DRT::UTILS::intrule_line_2point;
               break;
            case DRT::Element::line3:
              rule_pot = DRT::UTILS::intrule_line_3point;
              break;
            default:
               dserror("unknown number of nodes for gaussrule initialization");
         }
         const DRT::UTILS::IntegrationPoints1D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(1,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

            // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;
           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
            //cout << "potderiv1 = " << potderiv1 << endl;
            //cout << "potderiv2 = " << potderiv2 << endl;

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



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for surface elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  LINALG::SerialDenseMatrix dXYZdrs(2,3);
  LINALG::SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  LINALG::SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_,element,x,3);
  dXYZdrs.Multiply('N','N',1.0,deriv,X,0.0);
  LINALG::SerialDenseMatrix  metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dXYZdrs,dXYZdrs,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = sqrt(  metrictensor(0,0)*metrictensor(1,1)
                            -metrictensor(0,1)*metrictensor(1,0));
  double factor = intpoints.qwgt[gp] * detA * curve_fac;

  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp,  for line elements    u.may 02/09|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints1D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_1D(funct,e0,element->Shape());
  DRT::UTILS::shape_function_1D_deriv1(deriv,e0,element->Shape());

  LINALG::SerialDenseMatrix dXYZdr(1,3);
  LINALG::SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  LINALG::SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_, element,x,3);
  dXYZdr.Multiply('N','N',1.0,deriv,X,0.0);
  LINALG::SerialDenseMatrix  metrictensor(1,1);
  metrictensor.Multiply('N','T',1.0,dXYZdr,dXYZdr,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = metrictensor(0,0);
  double factor = intpoints.qwgt[gp] * detA * curve_fac;

  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
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



/*----------------------------------------------------------------------*
 |  compute surface normal in gauss point                    u.may 04/09|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> POTENTIAL::VolumePotential::ComputeNormalInGP(
    const DRT::Element*                     element,
    const LINALG::Matrix<3,1>&              x_gp)
{

  const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element, currentpositions_));
  LINALG::Matrix<2,1> elecoord_gp(true);
  GEO::CurrentToSurfaceElementCoordinates(element->Shape(), xyze, x_gp, elecoord_gp);
  LINALG::Matrix<3,1> normal_gp(true);
  GEO::computeNormalToSurfaceElement(element->Shape(), xyze, elecoord_gp, normal_gp);

  return normal_gp;
}



/*----------------------------------------------------------------------*
 |  compute surface normal in gauss point                    u.may 04/09|
 *----------------------------------------------------------------------*/
bool POTENTIAL::VolumePotential::DetermineValidContribution(
    const LINALG::Matrix<3,1>&              x_gp,
    const LINALG::Matrix<3,1>&              x_pot_gp,
    const LINALG::Matrix<3,1>&              normal_gp,
    const LINALG::Matrix<3,1>&              normal_pot_gp)
{
  // compare normals with radius in order to delete senseless contributions passing through the element
  // compute distance vector
  LINALG::Matrix<3,1> radius(true);
  radius.Update(1.0, x_pot_gp, -1.0, x_gp);

  LINALG::Matrix<3,1> radius_pot(true);
  radius_pot.Update(1.0, x_gp, -1.0, x_pot_gp);

  // compare normals
  const double scalarproduct = radius(0)*normal_gp(0) + radius(1)*normal_gp(1) + radius(2)*normal_gp(2);
  const double scalarproduct_pot =  radius_pot(0)*normal_pot_gp(0) +
                                    radius_pot(1)*normal_pot_gp(1) +
                                    radius_pot(2)*normal_pot_gp(2);

  // if valid contribution
  if(scalarproduct > (-1)*GEO::TOL13 && scalarproduct_pot > (-1)*GEO::TOL13 )
    return true;

  return false;
}



#endif

