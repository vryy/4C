/*!-------------------------------------------------------------------
\file drt_potential_manager.cpp

\brief  Class controlling surface stresses due to potential forces
        between interfaces of mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#include "drt_potential_surface.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/element_normals.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_bele3/bele2.H"
#include "../drt_lib/drt_condition_selector.H"
#include <cstdlib>

using namespace GEO;

/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
POTENTIAL::SurfacePotential::SurfacePotential(
    Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                discret,
    const GEO::TreeType&                treetype):
    Potential (discretRCP, discret)

{
  //std::cout << "pot man"  << std::endl;

  // create discretization from potential boundary elements and distribute them
  // on every processor
  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("Potential");

  if(prob_dim_ == 2)
    potentialdis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE2", conditions_to_copy);
  else if(prob_dim_ == 3)
    potentialdis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE3", conditions_to_copy);
  else
    dserror("problem dimension not correct");
  dsassert(potentialdis_->NumGlobalNodes() > 0, "empty discretization detected. Potential conditions applied?");

  // set new dof set
  RCP<POTENTIAL::PotentialDofSet> pdofset = Teuchos::rcp(new POTENTIAL::PotentialDofSet(discretRCP_));
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

  // split dof vector of soliddiscretization into a dof vector of potential boundary
  // condition and remaining dofs
  const int ndim = DRT::Problem::Instance()->NDim();
  {
    DRT::UTILS::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(*discretRCP_,"Potential",0,ndim)));
    mcs.SetupExtractor(*discretRCP_,*(discretRCP_->DofRowMap()),potboundary_);
  }

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
  importer_ = Teuchos::rcp(new Epetra_Import(idisp_total_->Map(),idisp_onproc_->Map()));

  // set up tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*potentialdis_);
  DRT::UTILS::CollectElementsByConditionLabel(*potentialdis_, elementsByLabel_,"Potential" );
  InvertElementsPerLabel(elementsByLabel_,labelByElement_);

  treetype_ = treetype;
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype_);



  // std::cout << "Potential manager constructor done" << std::endl;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::EvaluatePotential(  Teuchos::ParameterList& p,
                                                  RCP<Epetra_Vector> disp,
                                                  RCP<Epetra_Vector> fint,
                                                  RCP<LINALG::SparseMatrix> stiff)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  UpdateDisplacementsOfPotentialDiscretization(disp);

  EvaluateSurfacePotentialCondition(p,stiff,Teuchos::null,fint,Teuchos::null,Teuchos::null,"Potential");

  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| evaluate potential conditions based on a Epetra_FecrsMatrix        |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::EvaluateSurfacePotentialCondition(
    Teuchos::ParameterList&                          params,
    RCP<LINALG::SparseMatrix>       systemmatrix1,
    RCP<LINALG::SparseMatrix>       systemmatrix2,
    RCP<Epetra_Vector>              systemvector1,
    Teuchos::RCP<Epetra_Vector>             systemvector2,
    Teuchos::RCP<Epetra_Vector>             systemvector3,
    const std::string&                           condstring)
{
  if (!discret_.Filled()) dserror("FillComplete() was not called");
  if (!discret_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time < 0.0) usetime = false;

  if(time < 0.0)
    std::cout <<  "no time curve set " << std::endl;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through potential conditions and evaluate them
  //----------------------------------------------------------------------
  std::vector<DRT::Condition*> potentialcond;
  discret_.GetCondition(condstring, potentialcond);
  for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++ condIter)
  {
    std::map<int,RCP<DRT::Element> >& geom = (*condIter)->Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,RCP<DRT::Element> >::iterator curr;

    // Evaluate Loadcurve if defined. Put current load factor in parameterlist
    const std::vector<int>*    curve  = (*condIter)->Get<std::vector<int> >("curve");
    int                   curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double                curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

    params.set("LoadCurveFactor",curvefac);

    params.set<RCP<DRT::Condition> >("condition", Teuchos::rcp(*condIter,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmrowowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(discret_,lm,lmrowowner,lmstride);
      const int rowsize = lm.size();

      // Reshape element matrices and vectors and init to zero in element->evaluate
      // call the element specific evaluate method
      int err = curr->second->Evaluate( params,discret_,lm, elematrix1,elematrix2,
                                        elevector1,elevector2,elevector3);


      if (err) dserror("error while evaluating elements");

      // specify lm row and lm col
      std::vector<int> lmrow;
      std::vector<int> lmcol;
      // only local values appeared
      if((int) lm.size() == rowsize)
      {
        lmrow.resize(lm.size());
        lmcol.resize(lm.size());
        lmrow = lm;
        lmcol = lm;
      }
      // non-local values appeared
      else if((int) lm.size() > rowsize)
      {
        for(int i = 0; i < rowsize; i++)
          lmrow.push_back(lm[i]);

        lmcol.resize(lm.size());
        lmcol = lm;
      }
      else
        dserror("lm is not properly filled");

      // assembly
      int eid = curr->second->Id();
      if (assemblemat1) systemmatrix1->FEAssemble(eid,elematrix1,lmrow,lmrowowner,lmcol);
      if (assemblemat2) systemmatrix2->FEAssemble(eid,elematrix2,lmrow,lmrowowner,lmcol);

      if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lmrow,lmrowowner);
      if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lmrow,lmrowowner);
      if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lmrow,lmrowowner);
    }
  } // for(std::vector<DRT::Condition*>::iterator condIter = potentialcond->begin() ; condIter != potentialcond->end(); ++ condIter)
  return;
} // end of EvaluatePotentialCondition




/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  RCP<DRT::Condition> cond = params.get<RCP<DRT::Condition> >("condition",Teuchos::null);

  // find nodal ids influencing a given element
  const int     label     = cond->GetInt("label");		//jeder Körper besitzt ein label
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
  // std::cout << "stiffness stop" << std::endl;
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  11/09|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level (volume approximation 1)                          |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::StiffnessAndInternalForcesPotentialApprox1(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int
    )
{

  // get potential paramters
  RCP<DRT::Condition> cond = params.get<RCP<DRT::Condition> >("condition",Teuchos::null);
  const int     label     = cond->GetInt("label");		//jeder Körper besitzt ein label
  const double  cutOff    = cond->GetDouble("cutOff");

  // initialize time variables
  const int    curvenum = cond->GetInt("curve");
  const double time     = params.get<double>("total time",-1.0);
  //const double dt     = params.get<double>("delta time",0.0);

  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
  double curvefac       = 1.0;
  if (time <= t_end)
    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);


  // get element ids of influencing structures
  std::map<int,std::set<int> > potentialElementIds;

  for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
  {
    const DRT::Node* node = element->Nodes()[i];
    const LINALG::Matrix<3,1> x_node = currentpositions_.find(node->Id())->second;
    treeSearchElementsInCutOffRadius(potentialdis_, currentpositions_, x_node, potentialElementIds, cutOff, label);
  }

  // compute internal force and stiffness matrix
  computeFandK_Approx1_new(element, gaussrule, potentialElementIds, lm, K_surf, F_int, cond, label, curvefac);

  //K_surf.Print(std::cout);
  //F_int.Print(std::cout);


  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  11/09|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level (volume approximation 2)                          |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::StiffnessAndInternalForcesPotentialApprox2(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // get potential condition
  RCP<DRT::Condition> cond = params.get<RCP<DRT::Condition> >("condition",Teuchos::null);
  const int     label     = cond->GetInt("label");    //jeder Körper besitzt ein label
  const double  cutOff    = cond->GetDouble("cutOff");

  // initialize time variables
  const int    curvenum = cond->GetInt("curve");
  const double time     = params.get<double>("total time",-1.0);
  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
  double curvefac       = 1.0;
  if (time <= t_end)
    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get for each gauss point , each structure the nearest object
  std::map< int, std::map<int, GEO::NearestObject> > potentialObjectsAtGP = GetPotentialObjectsAtGP(element, gaussrule,label, cutOff);

  // compute internal force and stiffness matrix
  computeFandK_Approx2_new(element, gaussrule, potentialObjectsAtGP, lm, K_surf, F_int, cond, label, curvefac);

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| for line elements                                                  |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule1D&  gaussrule,
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  RCP<DRT::Condition> cond = params.get<RCP<DRT::Condition> >("condition",Teuchos::null);

  // find nodal ids influencing a given element
  const int     label     = cond->GetInt("label");
  const double  cutOff    = cond->GetDouble("cutOff");

  // initialize time variables
  const int    curvenum = cond->GetInt("curve");
  const double time     = params.get<double>("total time",-1.0);
  const double t_end    = DRT::Problem::Instance()->Curve(curvenum).end();
  double curvefac       = 1.0;

  if (time <= t_end)
    curvefac      = DRT::Problem::Instance()->Curve(curvenum).f(time);

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

  // compute internal force and stiffness matrix
  // TODO if poteles empty don t do assembly
  computeFandK(element, gaussrule, potentialElementIds, lm, K_surf, F_int, cond, label, curvefac);

  return;
}


/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::UpdateDisplacementsOfPotentialDiscretization(
    Teuchos::RCP<Epetra_Vector>     idisp_solid
    )
{
  // extract displacements associated with potential elements from solid discretization
  idisp_onproc_ = potboundary_.ExtractCondVector(idisp_solid);
  idisp_total_->Scale(0.0);

 // LINALG::Export(*idisp_onproc_,*idisp_total_);

  // import
  int err = idisp_total_->Import((*idisp_onproc_), (*importer_),Insert);
  if(err) dserror("Import using importer returned err=%d",err);

  currentpositions_.clear();
  {
    for (int lid = 0; lid < potentialdis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = potentialdis_->lColNode(lid);
      std::vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      potentialdis_->Dof(node, lm);
      std::vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;
      DRT::UTILS::ExtractMyValues(*idisp_total_,mydisp,lm);
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentpositions_[node->Id()] = currpos;
    }
  }

  // reinitialize search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*potentialdis_, currentpositions_);
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype_);
  /*if(prob_dim_ == 2)
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
void POTENTIAL::SurfacePotential::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule2D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   std::vector<int>&                lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RCP<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(potentialdis_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta"); // = 25000000;
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    Epetra_SerialDenseVector   funct(numnode);
    Epetra_SerialDenseMatrix   deriv(2,numnode);
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
         const int numnode_pot = element_pot->NumNode();
         const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

         // obtain current potential dofs
         std::vector<int> lmpot;
         std::vector<int> lmowner;
         std::vector<int> lmstride;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);

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
           Epetra_SerialDenseVector  funct_pot(numnode_pot);
           Epetra_SerialDenseMatrix  deriv_pot(2,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // compute normal of surface element in x_pot_gp
           LINALG::Matrix<3,1> n_pot_gp = ComputeNormalInGP(element_pot, x_pot_gp);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           bool validContribution = false;
           // contact
           //fabs(cond->GetDouble("exvollength")) > 1e-7
           if(fabs(cond->GetDouble("exvollength")) > 1e-7)
           {
             double exvollength = cond->GetDouble("exvollength");
             LINALG::Matrix<3,1> x_gp_con(true);
             LINALG::Matrix<3,1> x_pot_gp_con(true);
             x_gp_con = x_gp;
             x_pot_gp_con = x_pot_gp;

             //Add: this = scalarThis * this + scalarOther * other.
             x_gp_con.Update( (-1.0)*exvollength, n_gp, 1.0);
             x_pot_gp_con.Update( (-1.0)*exvollength, n_pot_gp, 1.0);

             EvaluatePotentialfromCondition(cond, x_gp_con, x_pot_gp_con, potderiv1, potderiv2);
             validContribution = DetermineValidContribution(x_gp_con, x_pot_gp_con, n_gp, n_pot_gp);
           }
           else
           {
             EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
             validContribution = DetermineValidContribution(x_gp, x_pot_gp, n_gp, n_pot_gp);
           }

          // if valid contribution
          if(validContribution)
          {
             //std::cout << "potderiv1 = " << potderiv1 << std::endl;
             //std::cout << "potderiv2 = " << potderiv2 << std::endl;
             const int numdof = 3;
             // computation of internal forces (possibly with non-local values)
             for (int inode = 0; inode < numnode; inode++)
               for(int dim = 0; dim < 3; dim++)
                 F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

             // computation of stiffness matrix (possibly with non-local values)
             for (int inode = 0;inode < numnode; ++inode)
               for(int dim = 0; dim < 3; dim++)
               {
                 // k,ii
                 for (int jnode = 0; jnode < numnode; ++jnode)
                   for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                     K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                       funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

                 // k,ij
                 for (int jnode = 0;jnode < numnode_pot; ++jnode)
                   for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                     K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                        funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

                }
             }// valid contribution
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}

/*-------------------------------------------------------------------*
| (private)                                               umay  10/09|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
| approximation                                                      |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeFandK_Approx1(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule2D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   std::vector<int>&                lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RCP<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac
   )
{

  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(potentialdis_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta"); // = 50 e21;

  //gaussrule sagt, wieviel gausspunkte verwendet werden sollen
  //intpoints enthält die Koordinaten und die Gewichtungsfaktoren
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------

  //nquad ist Attribut, bei vier Gausspunkten also 4
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    // zusätzlich Berechnung von FInvers und X_gp
    const int numnode = actEle->NumNode();
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(2,numnode);
    Epetra_SerialDenseMatrix 	  FInvers(3,3);
    LINALG::Matrix<3,1>         x_gp(true);
    LINALG::Matrix<3,1>         X_gp(true);

    const double fac = ComputeFactorApprox(actEle, funct, deriv, FInvers, intpoints, gp, x_gp, X_gp, curvefac);

    // compute normal n_gp to act ele in x_gp
    LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);
    //Normale auf das aktive Element in gp_X berechnen
    LINALG::Matrix<3,1> N_gp = ComputeNormalInGP_Initialconf(actEle, X_gp);


    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------

    //Schleife über die verschiedenen Körper, die nicht gleich dem eigenen label sind?
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)

      //Schleife über die potentiellen Elemente
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = potentialdis_->gElement(*eleIter);

         const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

         // obtain current potential dofs
         std::vector<int> lmpot;
         std::vector<int> lmowner;
         std::vector<int> lmstride;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);

         // obtain Gaussrule and integration points
         //warum muss Gaussregel so kompliziert hergeholt werden?
         DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
         GetGaussRule2D(element_pot->Shape(), rule_pot);
         //hier werden wieder die Gausspunkte + Gewichtungsfaktoren gespeichert
         const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           //Zusätzlich FInvers_pot und X_pot_gp
           const int numnode_pot = element_pot->NumNode();
           Epetra_SerialDenseVector  funct_pot(numnode_pot);
           Epetra_SerialDenseMatrix  deriv_pot(2,numnode_pot);
           Epetra_SerialDenseMatrix  FInvers_pot(3,3);
           LINALG::Matrix<3,1>        X_pot_gp(true);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactorApprox(element_pot, funct_pot, deriv_pot, FInvers_pot, intpoints_pot,
        		   									  gp_pot, x_pot_gp, X_pot_gp, curvefac);

           // compute normal of surface element in x_pot_gp
           LINALG::Matrix<3,1> n_pot_gp = ComputeNormalInGP(element_pot, x_pot_gp);
           //Normale in Gausspunkt für potentielles Element berechnen
           LINALG::Matrix<3,1> N_pot_gp = ComputeNormalInGP_Initialconf(element_pot, X_pot_gp);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           // contact
           //fabs(cond->GetDouble("exvollength")) > 1e-7
           //was soll das sein?
           if(fabs(cond->GetDouble("exvollength")) > 1e-7)
           {
             double exvollength = cond->GetDouble("exvollength");
             LINALG::Matrix<3,1> x_gp_con(true);
             LINALG::Matrix<3,1> x_pot_gp_con(true);
             x_gp_con = x_gp;
             x_pot_gp_con = x_pot_gp;

             //Add: this = scalarThis * this + scalarOther * other.
             x_gp_con.Update( (-1.0)*exvollength, n_gp, 1.0);
             x_pot_gp_con.Update( (-1.0)*exvollength, n_pot_gp, 1.0);

             EvaluatePotentialfromCondition_Approx1(cond, x_gp_con, x_pot_gp_con, potderiv1, potderiv2);
           }
           else
           {
        	 //Testet welches Potential verwendet wird Zeta oder LJ
        	 //potderiv1 speichert dF/dr*Einheitsdistanz potderiv2 entsprechend
             EvaluatePotentialfromCondition_Approx1(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
             //Testen ob die Elemente auf sich gegenüber liegenden Flächen liegen
             //Noch prüfen, ob das bei mir nötig ist!!!
             //validContribution = DetermineValidContribution(x_gp, x_pot_gp, n_gp, n_pot_gp);
           }

          // if valid contribution
          //if(validContribution)
          //{
             //std::cout << "potderiv1 = " << potderiv1 << std::endl;
             //std::cout << "potderiv2 = " << potderiv2 << std::endl;
             const int numdof = 3;
             //Zwischen speichert F_Invers*N
             Epetra_SerialDenseVector 	Zwischen (3);
             Epetra_SerialDenseVector 	Zwischen_pot (3);
             Epetra_SerialDenseVector 	N_gp_SD (3);
             Epetra_SerialDenseVector 	N_pot_gp_SD (3);
             LINALG::Matrix<3,1> 	      Einheitsabstand(true);
             double Teta=0;
             double Teta_pot=0;

             //Überschreiben von Matix<> in SerialDenseMatrix
             for(int i=0;i<3;i++)
             {
            	 N_gp_SD(i)=N_gp(i);
            	 N_pot_gp_SD(i)=N_pot_gp(i);
             }
             //Berechnung der zusätzlichen Teta Terme
             Zwischen.Multiply('T','N',1.0,FInvers,N_gp_SD,0.0);
             Zwischen_pot.Multiply('T','N',1.0,FInvers_pot,N_pot_gp_SD,0.0);
             computeDistanceVector(x_gp, x_pot_gp,Einheitsabstand);

             for(int i=0;i<3;i++)
             {
            	 Teta+=(-1)*Einheitsabstand(i)*Zwischen(i);
            	 Teta_pot+=Einheitsabstand(i)*Zwischen_pot(i);
             }

             // computation of internal forces (possibly with non-local values)
             for (int inode = 0; inode < numnode; inode++)
               for(int dim = 0; dim < 3; dim++)
                 F_int[inode*numdof+dim] += funct(inode)*beta*fac*Teta*(beta_pot*potderiv1(dim)*fac_pot*Teta_pot);

             // computation of stiffness matrix (possibly with non-local values)
             for (int inode = 0;inode < numnode; ++inode)
               for(int dim = 0; dim < 3; dim++)
               {
                 // k,ii
                 for (int jnode = 0; jnode < numnode; ++jnode)
                   {for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                     K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                       funct(inode)*beta*fac*Teta*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot*Teta_pot);}

                 // k,ij
                 for (int jnode = 0;jnode < numnode_pot; ++jnode)
                   for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                     {K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                        funct(inode)*beta*fac*Teta*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot*Teta_pot);}

                }
             //K_surf.Print(std::cout);
             //}// valid contribution
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element
  //F_int.Print(std::cout);

  return;
}


/*-------------------------------------------------------------------*
| (private)                                               umay  10/09|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
| approximation                                                      |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeFandK_Approx1_new(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule2D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   std::vector<int>&                lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RCP<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(potentialdis_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  const double beta = cond->GetDouble("beta"); // = 50 e21;
  //intpoints enthält die Koordinaten und die Gewichtungsfaktoren
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(2,numnode);
    LINALG::Matrix<3,1>         x_gp(true);
    LINALG::Matrix<3,1>         n_xsi(true);

    //std::cout << " act element" << std::endl;
    //actEle->Print(std::cout);
    const double fac = ComputeFactorApprox_new(actEle, funct, deriv, intpoints, gp, x_gp, n_xsi, curvefac);

    //std::cout << "act fac = "<< fac << std::endl;
    // compute normal n_gp to act ele in x_gp
    LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = potentialdis_->gElement(*eleIter);
         const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

         //std::cout << " influencing element" << std::endl;
         //element_pot->Print(std::cout);
         // obtain current potential dofs
         std::vector<int> lmpot;
         std::vector<int> lmowner;
         std::vector<int> lmstride;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);

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
           Epetra_SerialDenseVector  funct_pot(numnode_pot);
           Epetra_SerialDenseMatrix  deriv_pot(2,numnode_pot);
           LINALG::Matrix<3,1>       x_pot_gp(true);
           LINALG::Matrix<3,1>       n_xsi_pot(true);

           //std::cout << "fac_pot = "<< fac << std::endl;
           const double fac_pot = ComputeFactorApprox_new(element_pot,funct_pot, deriv_pot,
                                                      intpoints_pot, gp_pot, x_pot_gp,
                                                      n_xsi_pot, curvefac);

           // compute normal of surface element in x_pot_gp
           LINALG::Matrix<3,1> n_pot_gp = ComputeNormalInGP(element_pot, x_pot_gp);

           //n_pot_gp.Print(std::cout);
           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;
           LINALG::Matrix<3,1>  radius_unit(true);

           //bool validContribution = false;
           // if contact
/*           if(fabs(cond->GetDouble("exvollength")) > 1e-7)
           {
             double exvollength = cond->GetDouble("exvollength");
             LINALG::Matrix<3,1> x_gp_con(true);
             LINALG::Matrix<3,1> x_pot_gp_con(true);
             x_gp_con = x_gp;
             x_pot_gp_con = x_pot_gp;

             //Add: this = scalarThis * this + scalarOther * other.
             x_gp_con.Update( (-1.0)*exvollength, n_gp, 1.0);
             x_pot_gp_con.Update( (-1.0)*exvollength, n_pot_gp, 1.0);

             EvaluatePotentialfromCondition_Approx1(cond, x_gp_con, x_pot_gp_con,  potderiv1, potderiv2);
             validContribution = DetermineValidContribution(x_gp_con, x_pot_gp_con, n_gp, n_pot_gp);
           }
           else
           {
           */
             EvaluatePotentialfromCondition_Approx1(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
             //Testen ob die Elemente auf sich gegenüber liegenden Flächen liegen
             //Noch prüfen, ob das bei mir nötig ist!!!
             //validContribution = DetermineValidContribution(x_gp, x_pot_gp, n_gp, n_pot_gp);
        //   }

          // if valid contribution
          //if(validContribution)
          //{
             const int numdof = 3;

             // compute theta and theta_pot
             // TODO use LINALG MULTIPLY
             computeDistanceVector(x_gp, x_pot_gp,radius_unit);

             double theta = 0.0;
             double theta_pot = 0.0;
             for(int i=0;i<3;i++)
             {
               theta     +=  (-1.0)*radius_unit(i)*n_xsi(i);
               theta_pot +=  radius_unit(i)*n_xsi_pot(i);
             }

             theta = 1.0;
             theta_pot = 1.0;

             //std::cout << "theta  = "<< theta  <<  std::endl;
             //std::cout << "theta_pot  = "<< theta_pot  <<  std::endl;
             //std::cout << "beta = "<<  beta << std::endl;
             //std::cout << "beta_pot = "<<  beta_pot << std::endl;
             //std::cout << "fac = "<<  beta << std::endl;
             //std::cout << "fac_pot = "<<  beta_pot << std::endl;

             // computation of internal forces
             for (int inode = 0; inode < numnode; inode++)
               for(int dim = 0; dim < 3; dim++)
                 F_int[inode*numdof+dim] += funct(inode)*beta*fac*theta*(beta_pot*potderiv1(dim)*fac_pot*theta_pot);
             // computation of stiffness matrix (possibly with non-local values)
             for (int inode = 0;inode < numnode; ++inode)
               for(int dim = 0; dim < 3; dim++)
               {
                 // k,ii
                 for (int jnode = 0; jnode < numnode; ++jnode)
                   for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                     K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                       funct(inode)*beta*fac*theta*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot*theta_pot);

                 // k,ij
                 for (int jnode = 0;jnode < numnode_pot; ++jnode)
                   for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                      K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                        funct(inode)*beta*fac*theta*(beta_pot*(-1.0)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot*theta_pot);

                }
             //}// valid contribution
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  10/09|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
| approximation                                                      |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeFandK_Approx2(
   const DRT::Element*                                        actEle,
   const DRT::UTILS::GaussRule2D&                             gaussrule,
   const std::map< int, std::map<int, GEO::NearestObject> >&  potentialObjects,
   std::vector<int>&                                          lm,
   Epetra_SerialDenseMatrix&                                  K_surf,
   Epetra_SerialDenseVector&                                  F_int,
   RCP<DRT::Condition>                                cond,
   const int                                                  label,
   const double                                               curvefac)
{
  if((int)potentialObjects.size() == 0)
  {
    F_int.Size(lm.size());
    K_surf.Shape(lm.size(), lm.size());
    return;
  }

  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(potentialdis_, potentialObjects, lm);

  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  //Atomdichte
  const double beta = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    const int numnode = actEle->NumNode();
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(2,numnode);
    Epetra_SerialDenseMatrix    FInvers(3,3);
    LINALG::Matrix<3,1>         x_gp(true);
    LINALG::Matrix<3,1>         X_gp(true);
    const double fac = ComputeFactorApprox(actEle, funct, deriv, FInvers, intpoints, gp, x_gp, X_gp, curvefac);
    // compute normal n_gp to act ele in x_gp
    LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);
    //Normale auf das aktive Element in gp_X berechnen
    LINALG::Matrix<3,1> N_gp = ComputeNormalInGP_Initialconf(actEle, X_gp);

    //----------------------------------------------------------------------
    // loop over all structures
    //----------------------------------------------------------------------
    if(potentialObjects.find(gp) != potentialObjects.end())
      for(std::map<int, GEO::NearestObject >::const_iterator labelIter = (potentialObjects.find(gp)->second).begin();
      labelIter != (potentialObjects.find(gp)->second).end(); labelIter++)
      {
        GEO::NearestObject potObject = labelIter->second;
        const LINALG::Matrix<3,1> xp = potObject.getPhysCoord();
        const DRT::Element* element_pot = potentialdis_->gElement(GetElementId(potentialdis_, potObject));
        const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

        // obtain current potential dofs
        std::vector<int> lmpot;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);

        //Compute normal and detF in xp
        LINALG::Matrix<3,1> np;
        double detF = ComputeNormalAndDetFinXp(element_pot,xp,np);

        LINALG::Matrix<3,1>   Fs;
        LINALG::Matrix<3,3>   Fsderiv;
        LINALG::Matrix<3,1>   x_gp_con(true);
        LINALG::Matrix<3,1>   xp_con(true);

        bool validContribution = false;
        // contact
        //fabs(cond->GetDouble("exvollength")) > 1e-7
        if(fabs(cond->GetDouble("exvollength")) > 1e-7)
        {
          double exvollength = cond->GetDouble("exvollength");
          x_gp_con = x_gp;
          xp_con    = xp;
          //Add: this = scalarThis * this + scalarOther * other.
          x_gp_con.Update( (-1.0)*exvollength, n_gp, 1.0);
          xp_con.Update( (-1.0)*exvollength, np, 1.0);

          EvaluatePotentialfromCondition_Approx2(cond, x_gp_con, xp_con, Fs, Fsderiv);
          validContribution = DetermineValidContribution(x_gp_con, xp_con, n_gp, np);
        }
        else
        {
          EvaluatePotentialfromCondition_Approx2(cond, x_gp, xp, Fs, Fsderiv);
          //Testen ob die Elemente auf sich gegenüber liegenden Flächen liegen
          //Noch prüfen, ob das bei mir nötig ist!!!
          validContribution = DetermineValidContribution(x_gp, xp, n_gp, np);
        }

        // if valid contribution
        if(validContribution)
        {
          //std::cout << "potderiv1 = " << potderiv1 << std::endl;
          //std::cout << "potderiv2 = " << potderiv2 << std::endl;
          const int numdof = 3;
          //Zwischen speichert F_Invers(T)*N
          Epetra_SerialDenseVector  Zwischen(3);
          Epetra_SerialDenseVector  N_gp_SD(3);
          LINALG::Matrix<3,1>      Einheitsabstand(true);
          double Teta=0;

          //Überschreiben von Matix<> in SerialDenseMatrix
          for(int i=0;i<3;i++)
            N_gp_SD(i)=N_gp(i);

          //Berechnung des zusätzlichen Teta Terms
          Zwischen.Multiply('T','N',1.0,FInvers,N_gp_SD,0.0);

          for(int i=0;i<3;i++)
            Teta+=np(i)*Zwischen(i);

          // computation of internal forces (possibly with non-local values)
          for (int inode = 0; inode < numnode; inode++)
            for(int dim = 0; dim < 3; dim++)
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*Teta*(beta_pot*Fs(dim)*(1/detF));

          DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
          GetGaussRule2D(element_pot->Shape(), rule_pot);
          //hier werden wieder die Gausspunkte + Gewichtungsfaktoren gespeichert
          const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
          const int numnode_pot = element_pot->NumNode();
          //----------------------------------------------------------------------
          // run over all gauss points of a influencing element
          //----------------------------------------------------------------------
          for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
          {
            Epetra_SerialDenseVector  funct_pot(numnode_pot);
            const double e0 = intpoints_pot.qxg[gp_pot][0];
            const double e1 = intpoints_pot.qxg[gp_pot][1];
            // get shape functions of the element
            DRT::UTILS::shape_function_2D(funct_pot,e0,e1,element_pot->Shape());

            // computation of stiffness matrix (possibly with non-local values)
            for (int inode = 0;inode < numnode; ++inode)
              for(int dim = 0; dim < 3; dim++)
              {
                //was ist funct(jnode) funct_pot(jnode)
                // k,ii
                for (int jnode = 0; jnode < numnode; ++jnode)
                  for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                    K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                      funct(inode)*beta*fac*Teta*(beta_pot*(1.0/detF)*Fsderiv(dim,dim_pot)*funct(jnode));

                // k,ij
                for (int jnode = 0;jnode < numnode_pot; ++jnode)
                  for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                    K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*Teta*(beta_pot*(1.0/detF)*(-1)*Fsderiv(dim,dim_pot)*funct_pot(jnode));
              }
          }
          //K_surf.Print(std::cout);
        }// valid contribution
      }//einzelne Körper
    } // loop over all gauss points of the actual element
  //F_int.Print(std::cout);
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  10/09|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
| approximation                                                      |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeFandK_Approx2_new(
   const DRT::Element*                                        actEle,
   const DRT::UTILS::GaussRule2D&                             gaussrule,
   const std::map< int, std::map<int, GEO::NearestObject> >&  potentialObjects,
   std::vector<int>&                                          lm,
   Epetra_SerialDenseMatrix&                                  K_surf,
   Epetra_SerialDenseVector&                                  F_int,
   RCP<DRT::Condition>                                cond,
   const int                                                  label,
   const double                                               curvefac)
{
  if((int)potentialObjects.size() == 0)
  {
    F_int.Size(lm.size());
    K_surf.Shape(lm.size(), lm.size());
    return;
  }

  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(potentialdis_, potentialObjects, lm);

  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  //Atomdichte
  const double beta = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //std::cout << "ELEMENT = "<< actEle->Id()<<  std::endl;

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    const int numnode = actEle->NumNode();
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(2,numnode);
    //Epetra_SerialDenseMatrix    FInvers(3,3);
    LINALG::Matrix<3,1>         x_gp(true);
    //LINALG::Matrix<3,1>         X_gp(true);
    // const double fac = ComputeFactorApprox(actEle, funct, deriv, FInvers, intpoints, gp, x_gp, X_gp, curvefac);

    LINALG::Matrix<3,1>         n_xsi(true);
    const double fac = ComputeFactorApprox_new(actEle, funct, deriv, intpoints, gp, x_gp, n_xsi, curvefac);
    //std::cout << "fac = "  << fac << std::endl;
    // compute normal n_gp to act ele in x_gp
    //LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);
    //Normale auf das aktive Element in gp_X berechnen
    // LINALG::Matrix<3,1> N_gp = ComputeNormalInGP_Initialconf(actEle, X_gp);

    //----------------------------------------------------------------------
    // loop over all structures
    //----------------------------------------------------------------------
    if(potentialObjects.find(gp) != potentialObjects.end())
      for(std::map<int, GEO::NearestObject >::const_iterator labelIter = (potentialObjects.find(gp)->second).begin();
      labelIter != (potentialObjects.find(gp)->second).end(); labelIter++)
      {
        GEO::NearestObject potObject = labelIter->second;
        const LINALG::Matrix<3,1> xp = potObject.getPhysCoord();
        const DRT::Element* element_pot = potentialdis_->gElement(GetElementId(potentialdis_, potObject));
        const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

        // obtain current potential dofs
        std::vector<int> lmpot;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);

        //Compute normal on potential elment in xp
        LINALG::Matrix<3,1> np;
        ComputeNormalOnPotentialElement(element_pot,xp,np);

        //std::cout << "Normal on potential element" << std::endl;
        //np.Print(std::cout);

        LINALG::Matrix<3,1>   Fs;
        LINALG::Matrix<3,3>   Fsderiv;
        LINALG::Matrix<3,1>   x_gp_con(true);
        LINALG::Matrix<3,1>   xp_con(true);

        bool validContribution = false;
        // contact
        //fabs(cond->GetDouble("exvollength")) > 1e-7
        const LINALG::Matrix<3,1> n_gp = ComputeNormalInGP(actEle, x_gp);

        //std::cout << "Normal in Gauss Point" << std::endl;
        //n_gp.Print(std::cout);

        if(fabs(cond->GetDouble("exvollength")) > 1e-7)
        {
          double exvollength = cond->GetDouble("exvollength");
          x_gp_con = x_gp;
          xp_con    = xp;
          //Add: this = scalarThis * this + scalarOther * other.
          x_gp_con.Update( (-1.0)*exvollength, n_gp, 1.0);
          xp_con.Update( (-1.0)*exvollength, np, 1.0);

          EvaluatePotentialfromCondition_Approx2(cond, x_gp_con, xp_con, Fs, Fsderiv);
          validContribution = DetermineValidContribution(x_gp_con, xp_con, n_gp, np);
        }
        else
        {
          EvaluatePotentialfromCondition_Approx2(cond, x_gp, xp, Fs, Fsderiv);
          //Testen ob die Elemente auf sich gegenüber liegenden Flächen liegen
          //Noch prüfen, ob das bei mir nötig ist!!!
          validContribution = DetermineValidContribution(x_gp, xp, n_gp, np);
        }

        // if valid contribution
        if(validContribution)
        {
          //std::cout << "potderiv1 = " << potderiv1 << std::endl;
          //std::cout << "potderiv2 = " << potderiv2 << std::endl;
          const int numdof = 3;

          // computation of theta
          double theta = 0.0;
          for(int i=0;i<3;i++)
            theta += np(i)*n_xsi(i);

          //std::cout << "theta = " << theta << std::endl;
          // computation of internal forces (possibly with non-local values)
          for (int inode = 0; inode < numnode; inode++)
            for(int dim = 0; dim < 3; dim++)
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*theta*beta_pot*Fs(dim);

          DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
          GetGaussRule2D(element_pot->Shape(), rule_pot);
          const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
          const int numnode_pot = element_pot->NumNode();
          //----------------------------------------------------------------------
          // run over all gauss points of a influencing element
          //----------------------------------------------------------------------
          for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
          {
            Epetra_SerialDenseVector  funct_pot(numnode_pot);
            const double e0 = intpoints_pot.qxg[gp_pot][0];
            const double e1 = intpoints_pot.qxg[gp_pot][1];
            // get shape functions of the element
            DRT::UTILS::shape_function_2D(funct_pot,e0,e1,element_pot->Shape());

            // computation of stiffness matrix (possibly with non-local values)
            for (int inode = 0;inode < numnode; ++inode)
              for(int dim = 0; dim < 3; dim++)
              {
                //was ist funct(jnode) funct_pot(jnode)
                // k,ii
                for (int jnode = 0; jnode < numnode; ++jnode)
                  for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                    K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                      funct(inode)*beta*fac*theta*(beta_pot*Fsderiv(dim,dim_pot)*funct(jnode));

                // k,ij
                for (int jnode = 0;jnode < numnode_pot; ++jnode)
                  for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                    K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*theta*(beta_pot*(-1)*Fsderiv(dim,dim_pot)*funct_pot(jnode));
              }
          }
          //K_surf.Print(std::cout);
        }// valid contribution
      }//einzelne Körper
    } // loop over all gauss points of the actual element

  //F_int.Print(std::cout);
  return;
}





/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix for line        |
| elements (2D problems)                                             |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule1D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   std::vector<int>&                lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RCP<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{

  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
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
    Epetra_SerialDenseVector   funct(numnode);
    Epetra_SerialDenseMatrix   deriv(1,numnode);
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
         std::vector<int> lmpot;
         std::vector<int> lmowner;
         std::vector<int> lmstride;
         element_pot->LocationVector(*potentialdis_,lmpot,lmowner,lmstride);
         const double beta_pot = GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);

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
               break;
         }
         const DRT::UTILS::IntegrationPoints1D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           Epetra_SerialDenseVector  funct_pot(numnode_pot);
           Epetra_SerialDenseMatrix  deriv_pot(1,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

            // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;
           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
            //std::cout << "potderiv1 = " << potderiv1 << std::endl;
            //std::cout << "potderiv2 = " << potderiv2 << std::endl;

           const int numdof = 3;

           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 3; dim++)
               F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

           // computation of stiffness matrix (possibly with non-local values)
           for (int inode = 0;inode < numnode; ++inode)
             for(int dim = 0; dim < 3; dim++)
             {
               // k,ii
               for (int jnode = 0; jnode < numnode; ++jnode)
                 for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                   K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                     funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

               // k,ij
               for (int jnode = 0;jnode < numnode_pot; ++jnode)
                 for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                   K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

              }
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element
  return;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for surface elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::SurfacePotential::ComputeFactor(
    const DRT::Element*                     element,
    Epetra_SerialDenseVector&              funct,
    Epetra_SerialDenseMatrix&              deriv,
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

  Epetra_SerialDenseMatrix dXYZdrs(2,3);

  Epetra_SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);

  Epetra_SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_, element,x,3);


  dXYZdrs.Multiply('N','N',1.0,deriv,X,0.0);
  Epetra_SerialDenseMatrix  metrictensor(2,2);
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
 |  compute gauss points for tree approx search              u.may 08/09|
 *----------------------------------------------------------------------*/
std::vector< LINALG::Matrix<3,1> > POTENTIAL::SurfacePotential::ComputeGP(
    const DRT::Element*                     element,
    const DRT::UTILS::IntegrationPoints2D&  intpoints)
{
  std::vector< LINALG::Matrix<3,1> > gaussPoints;
  const int numnode = element->NumNode();

  Epetra_SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_,element,x,3);

  // compute Gauss points
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    const double e0 = intpoints.qxg[gp][0];
    const double e1 = intpoints.qxg[gp][1];

    // get shape functions and derivatives of the element
    Epetra_SerialDenseVector   funct(numnode);
    DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());

    LINALG::Matrix<3,1> x_gp(true);
    // compute gauss point in physical coordinates
    for (int inode = 0; inode < numnode; inode++)
      for(int dim = 0; dim < 3; dim++)
        x_gp(dim) += funct(inode)*x(inode,dim);

    gaussPoints.push_back(x_gp);
  }
  return gaussPoints;
}



/*-------------------------------------------------------------------*
|   VON MIR GESCHRIEBEN                                              |
*--------------------------------------------------------------------*/
double POTENTIAL::SurfacePotential::ComputeFactorApprox(
    const DRT::Element*                     element,
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
    Epetra_SerialDenseMatrix&               FInvers,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    LINALG::Matrix<3,1>&                    X_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  //§1 und §2 Komponenter des Gauspunktes wird gespeichert
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  //werte der shape funktions am gp werden in funct gespeichert
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  //werte der ersten Ableitung der shape functions in gp wird in Matrix deriv gespeichert
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  Epetra_SerialDenseMatrix dXYZdrs(2,3); //2 Reihen 3 Spalten
  //TODO check rows and columns

  Epetra_SerialDenseMatrix X(numnode,3);
  //Die räumlichen Knotenkoordinaten des aktiven Elements werden gespeichert
  ReferenceConfiguration(element,X,3);

  Epetra_SerialDenseMatrix x(numnode,3);
  //Die materiellen Knotenkoordinaten des aktiven Elements werden gespeichert
  SpatialConfiguration(currentpositions_, element,x,3);

  dXYZdrs.Multiply('N','N',1.0,deriv,X,0.0);
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  //Jakobi-Matrix Referenzkonfig in Parameterraum
  metrictensor.Multiply('N','T',1.0,dXYZdrs,dXYZdrs,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = sqrt(  metrictensor(0,0)*metrictensor(1,1)
      -metrictensor(0,1)*metrictensor(1,0));
  double factor = intpoints.qwgt[gp] * detA * curve_fac;

  //Berechnung Inverser Deformationsgradient
  Epetra_SerialDenseMatrix dudrs(3,2);
  Epetra_SerialDenseMatrix dxyzdrs(3,2);
  Epetra_SerialDenseMatrix dxyzdrsInvers(2,3);
  Epetra_SerialDenseMatrix u(numnode,3);
  for(int i=0; i<numnode; i++)
  {
    for (int j=0; j<3; j++)
      u(i,j)=x(i,j)-X(i,j);
  }
  //Ableitung von u nach Parameterraum, keine quadratische Matrix
  dudrs.Multiply('T','T',1.0,u,deriv,0.0);

  //Berechnung Pseudoinverse
  dxyzdrs.Multiply('T','T',1.0,x,deriv,0.0);

  Epetra_SerialDenseMatrix	W (2,1);
  Epetra_SerialDenseMatrix	V (2,2);
  // Epetra_SerialDenseMatrix dxyzdrsInvers(2,3);

  //sigular value decomposition
  GEO::svdcmpSerialDense(dxyzdrs,W,V);

  /*
    //Test ob dxyzdrs richtig zerlegt wurde
    Epetra_SerialDenseMatrix  Test (3,2);
    Epetra_SerialDenseMatrix  Zwischenergebnis (2,2);

      for (int i=0;i<2;i++)
      {
    	  for (int j=0;j<2;j++)
    		  Zwischenergebnis(i,j)=VSD(j,i)*WSD(i,1);
      }
      Test.Multiply('N','N',1.0,dxyzdrs,Zwischenergebnis,0.0);
      std::cout<<"Test"<<std::endl;
      Test.Print(std::cout);
   */

  //V*WInvers ausrechnen und in V Speichern
  for(int i=0;i<2;i++)
  {
    //Verhindern, dass durch Null geteilt wird
    if(W(i,0)!=0.0)
      W(i,0)=(1/W(i,0));
    for (int j=0;j<2;j++)
      V(j,i)*=W(i,0);

  }

  dxyzdrsInvers.Multiply('N','T',1.0,V,dxyzdrs,0.0);
  //Berechung Pseudoinverse Ende

  //Berechnen des Inversen Deformationsgradieneten FInvers

  FInvers.Multiply('N','N',-1.0,dudrs,dxyzdrsInvers,0.0);
  for(int i=0;i<3;i++)
    FInvers(i,i)+=1.0;

  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  X_gp = 0.0;
  // Gausspunkt in Referenzkonfiguration berechnen
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      X_gp(dim) += funct(inode)*X(inode,dim);
  return factor;
}


/*-------------------------------------------------------------------*
|   VON MIR GESCHRIEBEN                                              |
*--------------------------------------------------------------------*/
double POTENTIAL::SurfacePotential::ComputeFactorApprox_new(
    const DRT::Element*                     element,
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    LINALG::Matrix<3,1>&                    n_xsi,
    const double                            curve_fac)
{
  const int    numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  Epetra_SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  Epetra_SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_, element,x,3);

  Epetra_SerialDenseMatrix dXYZdrs(2,3); //2 Reihen 3 Spalten
  dXYZdrs.Multiply('N','N', 1.0, deriv, X, 0.0);
  Epetra_SerialDenseMatrix dxyzdrs(2,3); //2 Reihen 3 Spalten
  dxyzdrs.Multiply('N','N', 1.0, deriv, x, 0.0);

  // compute n_xsi = x,r x x,s
  n_xsi(0) = dxyzdrs(0,1)*dxyzdrs(1,2) - dxyzdrs(0,2)*dxyzdrs(1,1) ;
  n_xsi(1) = dxyzdrs(0,2)*dxyzdrs(1,0) - dxyzdrs(0,0)*dxyzdrs(1,2) ;
  n_xsi(2) = dxyzdrs(0,0)*dxyzdrs(1,1) - dxyzdrs(0,1)*dxyzdrs(1,0) ;


  // compute J_x_xsi = || n_xsi ||
  const double J_x_xsi = n_xsi.Norm2();
  n_xsi.Scale(1.0/J_x_xsi);

  // compute N_xsi = X,r x X,s
  LINALG::Matrix<3,1> N_xsi(true);
  N_xsi(0) = dXYZdrs(0,1)*dXYZdrs(1,2) - dXYZdrs(0,2)*dXYZdrs(1,1) ;
  N_xsi(1) = dXYZdrs(0,2)*dXYZdrs(1,0) - dXYZdrs(0,0)*dXYZdrs(1,2) ;
  N_xsi(2) = dXYZdrs(0,0)*dXYZdrs(1,1) - dXYZdrs(0,1)*dXYZdrs(1,0) ;

  // compute J_X_xsi = || N_xsi ||
  const double J_X_xsi = N_xsi.Norm2();
  // compute J_x_X accounts for beta/Jx_X = beta_0
  const double J_x_X = J_x_xsi/J_X_xsi;

  // compute factor
  // (1.0/J_x_X) scales beta
  double factor = intpoints.qwgt[gp] * (1.0/J_x_X) * curve_fac*J_x_xsi;

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
double POTENTIAL::SurfacePotential::ComputeFactor(
    const DRT::Element*                     element,
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
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

  Epetra_SerialDenseMatrix dXYZdr(1,3);
  Epetra_SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  Epetra_SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_, element,x,3);
  dXYZdr.Multiply('N','N',1.0,deriv,X,0.0);
  Epetra_SerialDenseMatrix  metrictensor(1,1);
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
 |  get Gauss rule and integraion points                    u.may 12/09 |
 *----------------------------------------------------------------------*/
std::map< int, std::map<int, GEO::NearestObject> > POTENTIAL::SurfacePotential::GetPotentialObjectsAtGP(
    const DRT::Element*                 element,
    const DRT::UTILS::GaussRule2D&      gaussrule,
    const int                           label,
    const double                        cutoff)
{
  // compute Gauss point for element
  const std::vector< LINALG::Matrix<3,1> > gaussPoints = ComputeGP(element, DRT::UTILS::IntegrationPoints2D(gaussrule));

  // compute XAABB for element
  LINALG::Matrix<3,2> eleXAABB =  GEO::computeFastXAABB(element->Shape(),
                                  GEO::getCurrentNodalPositions(element, currentpositions_),
                                  GEO::HIGHERORDER);

  // enlarge box by cut off radius
  for(int dim = 0; dim < 3; dim++)
  {
    eleXAABB(dim,0) = eleXAABB(dim,0) - cutoff;
    eleXAABB(dim,1) = eleXAABB(dim,1) + cutoff;
  }

  // search for each gauss point the nearest object per structure
  std::map< int, std::map<int, GEO::NearestObject> >   potentialObjectsAtGP;

  int num = 1;
  searchTree_->queryPotentialElements_Approx2(*potentialdis_, currentpositions_, eleXAABB, gaussPoints,
                                               potentialObjectsAtGP, cutoff,label,  num);

  return potentialObjectsAtGP;
}


/*-------------------------------------------------------------------*
|   compute unit distance vector                         u.may 10/09 |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::computeDistanceVector(
    const LINALG::Matrix<3,1>&  x,
    const LINALG::Matrix<3,1>&  y,
    LINALG::Matrix<3,1>&        dist_unit)
{
  LINALG::Matrix<3,1>         dist_vec;
  double                      distance;

  dist_vec.Update(1.0, x, -1.0, y);
  distance = dist_vec.Norm2();
  dist_unit = dist_vec;
  dist_unit.Scale(1.0/distance);
  return;
}



/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 04/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::GetGaussRule2D(
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
      break;
  }
  return;
}


/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 12/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::ComputeNormalOnPotentialElement(
    const DRT::Element*                     element,
    const LINALG::Matrix<3,1>&              xp,
    LINALG::Matrix<3,1>&                    np)
{
  //Berechnung der Normalen
  const Epetra_SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element, currentpositions_));
  LINALG::Matrix<2,1> elecoord_xp(true);
  GEO::CurrentToSurfaceElementCoordinates(element->Shape(), xyze, xp, elecoord_xp);
  //Berechnet np in physikalischen Koordinaten?
  GEO::computeNormalToSurfaceElement(element->Shape(), xyze, elecoord_xp, np);
  return;
}




/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 12/09|
 *----------------------------------------------------------------------*/
double POTENTIAL::SurfacePotential::ComputeNormalAndDetFinXp(
    const DRT::Element*                     element,
    const LINALG::Matrix<3,1>&              xp,
    LINALG::Matrix<3,1>&                    np)
{
  //Berechnung der Normalen
  const Epetra_SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element, currentpositions_));
  LINALG::Matrix<2,1> elecoord_xp(true);
  GEO::CurrentToSurfaceElementCoordinates(element->Shape(), xyze, xp, elecoord_xp);
  //Berechnet np in physikalischen Koordinaten?
  GEO::computeNormalToSurfaceElement(element->Shape(), xyze, elecoord_xp, np);

  //Berechnung der Determinante
  const double e0 = elecoord_xp(0);
  const double e1 = elecoord_xp(1);

  const int numnode = element->NumNode();
  Epetra_SerialDenseVector   funct(numnode);
  Epetra_SerialDenseMatrix   deriv(2,numnode);
  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  Epetra_SerialDenseMatrix dXYZdrs(3,2);
  Epetra_SerialDenseMatrix dudrs(3,2);
  Epetra_SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_, element,x,3);
  //Epetra_SerialDenseMatrix dxyzdrs(3,2);
  //dxyzdrs.Multiply('T','T',1.0, x, deriv, 0.0);
  Epetra_SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);

  //Berechnung der discreten Knotenverschiebungen
  Epetra_SerialDenseMatrix u(numnode,3);
  for(int i=0; i<numnode; i++)
    {
      for (int j=0; j<3; j++)
        u(i,j)=x(i,j)-X(i,j);
    }

  dXYZdrs.Multiply('T','T',1.0,X,deriv,0.0);
  dudrs.Multiply('T','T',1.0,u,deriv,0.0);

  //Berechnung der Pseudoinversen von dXYZdrs
  Epetra_SerialDenseMatrix W (2,1);
  Epetra_SerialDenseMatrix V (2,2);
  Epetra_SerialDenseMatrix dXYZdrsInvers(2,3);

  //sigular value decomposition
  GEO::svdcmpSerialDense(dXYZdrs,W,V);

  //V*WInvers ausrechnen und in V Speichern
  for(int i=0;i<2;i++)
    {
        //Verhindern, dass durch Null geteilt wird
        if(W(i,0)!=0.0)
          W(i,0)=(1/W(i,0));
        for (int j=0;j<2;j++)
          V(j,i)*=W(i,0);
      }

   // TODO does not compile from here on
   dXYZdrsInvers.Multiply('N','T',1.0,V,dXYZdrs,0.0);
   //Berechung Pseudoinverse Ende

   //Berchnung des Deformationsgradienten
   Epetra_SerialDenseMatrix F (3,3);
   F.Multiply('N','N',1.0,dudrs,dXYZdrsInvers,0.0);
   for(int i=0;i<3;i++)
       F(i,i)+=1.0;

   //F.Print(std::cout);

   const double detF =
     F(0,0)*(F(1,1)*F(2,2)-F(1,2)*F(2,1))
    -F(1,0)*(F(0,1)*F(2,2)-F(0,2)*F(2,1))
    +F(2,0)*(F(0,1)*F(1,2)-F(0,2)*F(1,1));

  return detF;
}




/*----------------------------------------------------------------------*
 |  compute surface normal in gauss point                    u.may 04/09|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> POTENTIAL::SurfacePotential::ComputeNormalInGP(
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

/*-------------------------------------------------------------------*
| 	VON MIR GESCHRIEBEN                                              |
*--------------------------------------------------------------------*/

// Methode zur Berechnung der Normalen im Gausspunkt in der Anfangskonfiguration

LINALG::Matrix<3,1> POTENTIAL::SurfacePotential::ComputeNormalInGP_Initialconf(
    const DRT::Element*                     element,
    const LINALG::Matrix<3,1>&              X_gp)
{
  const LINALG::SerialDenseMatrix XYZe(GEO::InitialPositionArray(element));
  LINALG::Matrix<2,1> elecoord_gp(true);
  GEO::CurrentToSurfaceElementCoordinates(element->Shape(), XYZe, X_gp, elecoord_gp);
  LINALG::Matrix<3,1> Normal_gp(true);
  GEO::computeNormalToSurfaceElement(element->Shape(), XYZe, elecoord_gp, Normal_gp);

  return Normal_gp;
}


/*----------------------------------------------------------------------*
 |  compute surface normal in gauss point                    u.may 04/09|
 *----------------------------------------------------------------------*/
bool POTENTIAL::SurfacePotential::DetermineValidContribution(
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

  if((radius.Norm2()-radius_pot.Norm2()) > GEO::TOL7)
    dserror("radius is not the same");


  // compare normals
  const double scalarproduct = radius(0)*normal_gp(0) + radius(1)*normal_gp(1) + radius(2)*normal_gp(2);

  const double scalarproduct_pot =  radius_pot(0)*normal_pot_gp(0) +
                                    radius_pot(1)*normal_pot_gp(1) +
                                    radius_pot(2)*normal_pot_gp(2);



  // if valid contribution
  if(scalarproduct > (-1)*GEO::TOL7 && scalarproduct_pot > (-1)*GEO::TOL7 )
    return true;

  //std::cout << "scalarproduct =" <<  scalarproduct <<  std::endl;
  //std::cout << "scalarproduct_pot =" <<  scalarproduct_pot <<  std::endl;
  //std::cout << "radius.Norm2() ="<< radius.Norm2() << std::endl;

  return false;
}




///////////////////////// test potential //////////////////////////////////////
/*-------------------------------------------------------------------*
| (public)                                                 umay 01/10|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::TestEvaluatePotential(
  Teuchos::ParameterList&                      p,
  RCP<Epetra_Vector>          disp,
  RCP<Epetra_Vector>          fint,
  RCP<LINALG::SparseMatrix>   stiff,
  const double                        time,
  const int                           step)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);

  // TODO compute elements by label for volume elements of the spheres
  EvaluateSurfacePotentialCondition(p,stiff,Teuchos::null,fint,Teuchos::null,Teuchos::null,"Potential");

  std::map<int, std::set<int> > elementsByLabel_Vol = computeEleByLabelVol(discret_.GetState("displacement") ,elementsByLabel_);
  // compute test results
  if( p.get<std::string>("solution type") == "Sphere" )
    computeTestVanDerWaalsSpheres(potentialdis_, elementsByLabel_Vol, elementsByLabel_,
                                  discret_.GetState("displacement")  , fint, time, step,
                                  p.get("vdw_radius", 0.0), p.get("n_offset", 0.0));
  else if( p.get<std::string>("solution type") == "Sphere" )
  computeTestVanDerWaalsMembranes(potentialdis_, elementsByLabel_Vol, elementsByLabel_,
                                  discret_.GetState("displacement")  , fint, time, step,
                                  p.get("vdw_radius", 0.0), p.get("n_offset", 0.0),
                                  p.get("thickness", 0.0));

  else
    dserror("specify proper solution type");
  return;
}


/*-------------------------------------------------------------------*
| (public)                                                 umay 01/10|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
std::map<int, std::set<int> > POTENTIAL::SurfacePotential::computeEleByLabelVol(
                                                RCP<const Epetra_Vector>    disp,
                                                std::map<int, std::set<int> > elementList)
{
 //TODO check for 2 dims
 std::map<int, std::set<int> > elementsByLabel_Vol;

 std::map<int, LINALG::Matrix<3,2> > label_AABBs;
 // get bounding boxes of two spheres
 for(std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
 {
   // initialize xaabb_label with box around first element
   const int eleId = *((labelIter->second).begin());
   const DRT::Element* element_in = potentialdis_->gElement(eleId);
   LINALG::SerialDenseMatrix xyze(3, element_in->NumNode());
   getPhysicalEleCoords(potentialdis_, disp,element_in,xyze);
   LINALG::Matrix<3,2> xaabb_label = GEO::computeFastXAABB(element_in->Shape(), xyze, GEO::EleGeoType(GEO::LINEAR));

   // run over set elements
   for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
   {
     const DRT::Element* element = potentialdis_->gElement(*eleIter);
     getPhysicalEleCoords(potentialdis_,disp,element,xyze);
     LINALG::Matrix<3,2> xaabbEle = GEO::computeFastXAABB(element->Shape(), xyze, GEO::EleGeoType(GEO::LINEAR));
     xaabb_label = GEO::mergeAABB(xaabb_label, xaabbEle);
   }
   label_AABBs[labelIter->first] = xaabb_label;
 }


 for(int i_rowele = 0; i_rowele < discretRCP_->NumMyRowElements(); i_rowele++)
 {
   const DRT::Element* element_row = discretRCP_->lRowElement(i_rowele);
   LINALG::SerialDenseMatrix xyze(3, element_row->NumNode());
   getPhysicalEleCoords(discretRCP_, disp,element_row,xyze);
   LINALG::Matrix<3,2> xaabb_ele = GEO::computeFastXAABB(element_row->Shape(), xyze, GEO::EleGeoType(GEO::LINEAR));

   for(std::map<int, LINALG::Matrix<3,2> >::const_iterator mapIter = label_AABBs.begin() ; mapIter != label_AABBs.end(); ++mapIter)
     if( GEO::intersectionOfXAABB<3>(mapIter->second, xaabb_ele) )
       elementsByLabel_Vol[mapIter->first].insert(element_row->Id());
 }
 return elementsByLabel_Vol;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 11/09|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
/*
void POTENTIAL::SurfacePotential::TestEvaluatePotential(Teuchos::ParameterList& p,
                                                      RCP<Epetra_Vector> disp,
                                                      RCP<Epetra_Vector> fint,
                                                      RCP<LINALG::SparseMatrix> stiff,
                                                      const double time,
                                                      const int                           step)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);

  EvaluateSurfacePotentialCondition(p,stiff,Teuchos::null,fint,Teuchos::null,Teuchos::null,"Potential");

  // center of gravity

  //Speichert Gesammtkraft auf einzelnen Körper
  LINALG::Matrix<3,1> fint_sum_Body1 (true);
  LINALG::Matrix<3,1> fint_sum_Body2 (true);
  //Speichert den Schwerpunkt der einzelnen Körper
  LINALG::Matrix<3,1> Schwerpunkt_Body1 (true);
  LINALG::Matrix<3,1> Schwerpunkt_Body2 (true);
  //Distanzvektor zwischen den Schwerpunkten
  LINALG::Matrix<3,1> DistanceVector (true);

  //Schleife über die Elemente des ersten Körpers, muss je nach Geometrie angepasst werde
  //FintSumAndCenterOfGravityVector(fint_sum_Body1, Schwerpunkt_Body1, fint, disp, 0, 1080 , 0, 3645);
  FintSumAndCenterOfGravityVector(fint_sum_Body1, Schwerpunkt_Body1, fint, disp, 0, 2048 , 0, 6819);
  //FintSumAndCenterOfGravityVector(fint_sum_Body2, Schwerpunkt_Body2, fint, disp, 1080, 2160 , 3645, 7290);
  FintSumAndCenterOfGravityVector(fint_sum_Body2, Schwerpunkt_Body2, fint, disp, 2048, 4096 , 6819, 13638);

  //Abstandsvektor zeigt von Körper 1 nach Körper 2
  for(int i=0; i<3; i++)
    DistanceVector(i) = Schwerpunkt_Body2(i) - Schwerpunkt_Body1(i);



  double dDistance = DistanceVector.Norm2();
  double dForce1 = fint_sum_Body1.Norm2();
  double dForce2 = fint_sum_Body2.Norm2();
  std::cout<<std::endl<<std::endl<<"Abstand:"<<dDistance<<std::endl;


  fint_sum_Body1.Print(std::cout);
  std::cout<<"Kraft1:"<<dForce1<<std::endl;
  std::cout<<"Kraft2:"<<dForce2<<std::endl;
  AusgabedateiSchreiben(dDistance, dForce1, dForce2, time);

  return;
}

//*/

/*-------------------------------------------------------------------*
| (public)                                                 umay 11/09|
|                                                                    |
| compute potential forces and center of gravity                     |
*--------------------------------------------------------------------*/
//Rechnet Schwerpunkt des Körpers und die im Schwerpunkt wirkende Kraft aus
void POTENTIAL::SurfacePotential::FintSumAndCenterOfGravityVector(  LINALG::Matrix<3,1>& fint_sum_Body,
                                  LINALG::Matrix<3,1>& Schwerpunkt_Body,
                                  RCP<Epetra_Vector> fint,
                                  RCP<Epetra_Vector> disp,
                                  int gidBegin,
                                  int gidEnd,
                                  int dofBegin,
                                  int dofEnd)
{
  double Volumen_Body=0;

  for(int gid=gidBegin; gid < gidEnd ; gid++)
    {
      const DRT::Element* element = discret_.gElement(gid);
      int iNode = element->NumNode();
      int dim = 3;
      LINALG::SerialDenseMatrix xyze(dim , iNode);

      UpdateDisplacements(disp, element, xyze);

      //xsi = [0,0,0] weil Elementschwerpunkt im Parameterraum im Ursprung
      double VolElement=0;
      LINALG::Matrix<3,1> xsi(true);
      LINALG::Matrix<3,1> x(true);

      GEO::elementToCurrentCoordinatesT<DRT::Element::hex8>(xyze, xsi, x);
      VolElement = GEO::ElementVolumeT<DRT::Element::hex8>(xyze);
      //Schwerpunkt multipliziert mit Elementvolumen
      x.Scale(VolElement);

      //Aufsummieren des Volumens und (Elementvolumen x Schwerpunktelement)
      Volumen_Body += VolElement;
      Schwerpunkt_Body += x;
    }

  std::cout<<"Volumen:"<<Volumen_Body<<std::endl;

  Schwerpunkt_Body.Scale(1/Volumen_Body);

  // total force in center of gravity

  for(int j=0; j<3; j++)
    for(int i=j+dofBegin; i< dofEnd; i+=3)
        {
            fint_sum_Body(j) += fint->operator[] (i);
          }

   return;
}


/*-------------------------------------------------------------------*
| (public)                                                 umay 11/09|
|                                                                    |
| TEST update displacements                                          |
*--------------------------------------------------------------------*/
//Rechnet die aktuellen Kontenpositionen eines Elements aus und speichet sie in xyze
void POTENTIAL::SurfacePotential::UpdateDisplacements(
    Teuchos::RCP<Epetra_Vector>     idisp_solid,
    const DRT::Element*             element,
    LINALG::SerialDenseMatrix&      xyze)
{
  const DRT::Node*const* node = element->Nodes();

  for (int i=0; i< element->NumNode(); i++)
  {
    std::vector<int> lm;
    lm.reserve(3);
    discretRCP_->Dof(node[i], lm);
    std::vector<double> mydisp(3);

    //Updaten der Knotenpositonen und speichern in xyze
    DRT::UTILS::ExtractMyValues(*idisp_solid,mydisp,lm);
    xyze(0,i) = node[i]->X()[0] + mydisp[0];
    xyze(1,i) = node[i]->X()[1] + mydisp[1];
    xyze(2,i) = node[i]->X()[2] + mydisp[2];
  }

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                 umay 11/09|
|                                                                    |
| test output                                                        |
*--------------------------------------------------------------------*/
void POTENTIAL::SurfacePotential::AusgabedateiSchreiben(double dDistance,
                            double dForce1,
                            double dForce2,
                            const double time)
{
  std::ofstream AusgabeDatei("NumerischeLoesung.txt",std::ios_base::app);

  if(AusgabeDatei.good())
  {
    AusgabeDatei<<time<<"\t\t\t"<<dForce1<<"\t\t\t"<<dDistance<<std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/* TEST evaluate _certain_ potential forces and stiffness
 * evaluation happens internal-force like */
/*
void STR::TimIntImpl::TestForceStiffPotential
(
  const double time,
  const Teuchos::RCP<Epetra_Vector> dis
)
{
  // potential force loads (but on internal force vector side)
  if (potman_ != Teuchos::null)
  {
    Teuchos::ParameterList p; // create the parameters for manager
    p.set("pot_man", potman_);
    p.set("total time", time);

    Teuchos::RCP<LINALG::SparseMatrix> stiff_test=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false, LINALG::SparseMatrix::FE_MATRIX));
    Teuchos::RCP<Epetra_Vector> fint_test=LINALG::CreateVector(*dofrowmap_, true);
    fint_test->PutScalar(0.0);
    stiff_test->Zero();

    potman_->TestEvaluatePotential(p, dis, fint_test, stiff_test, time);
  }
  // wooop
  return;
}
*/

