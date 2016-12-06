/*----------------------------------------------------------------------*/
/*!
\file xfem_neumann.cpp

\brief base xfem Neumann boundary conditions

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#include "../drt_xfem/xfem_neumann.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../linalg/linalg_utils.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"




/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                    schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumann(                    Teuchos::RCP<GEO::CutWizard>         wizard,
                                               Teuchos::ParameterList&              params,
                                               Teuchos::RCP<DRT::Discretization>    discret,
                                               Teuchos::RCP<Epetra_Vector>          systemvector,
                                               Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (systemmatrix==Teuchos::null)
    EvaluateNeumann(wizard,params,discret,*systemvector);
  else
    EvaluateNeumann(wizard,params,discret,*systemvector,systemmatrix.get());
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                    schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumann(  Teuchos::RCP<GEO::CutWizard>         wizard,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<DRT::Discretization>    discret,
                             Epetra_Vector&                       systemvector,
                             LINALG::SparseOperator*              systemmatrix)
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 5) EvaluateNeumann" );


  if (!discret->Filled()) dserror("FillComplete() was not called");
  if (!discret->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  bool assemblemat = (systemmatrix != NULL);

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  std::multimap<std::string,DRT::Condition* >::iterator fool;
  std::multimap<std::string,DRT::Condition* > condition;

  // vector for conditions of one special type
  std::vector<DRT::Condition *> condition_vec;

  //================================================
  // Load Neumann conditions from discretization
  // REMARK: ° standard volume Neumann conditions are not loaded -> evaluated in Evaluate
  //         ° for XFEM Neumann boundaries: we assumme only XFEM Surface(!) Neumann conditions
  //================================================

  // get standard Point Neumann conditions
  condition_vec.clear();
  discret->GetCondition("PointNeumann", condition_vec);
  // copy conditions to a condition multimap
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("PointNeumann"),condition_vec[i]));
  }

  // get standard Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("LineNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("LineNeumann"),condition_vec[i]));
  }

  // get standard Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("SurfaceNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("SurfaceNeumann"),condition_vec[i]));
  }

  // get XFEM Point Neumann conditions
  condition_vec.clear();
  discret->GetCondition("PointXFEMNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("PointXFEMNeumann"),condition_vec[i]));
  }

  // get XFEM Line Neumann conditions
  condition_vec.clear();
  discret->GetCondition("LineXFEMNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("LineXFEMNeumann"),condition_vec[i]));
  }

  // get XFEM Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("SurfaceXFEMNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("SurfaceXFEMNeumann"),condition_vec[i]));
  }

  // get XFEM Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("VolXFEMNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("VolXFEMNeumann"),condition_vec[i]));
  }

  // evaluate standard Neumann conditions
  EvaluateNeumannStandard(condition,
      usetime,
      time,
      assemblemat,
      params,
      discret,
      systemvector,
      systemmatrix);

  // evaluate XFEM Neumann conditions
  EvaluateNeumannXFEM(   wizard,
      condition,
      usetime,
      time,
      assemblemat,
      params,
      discret,
      systemvector,
      systemmatrix);


  return;

}



/*----------------------------------------------------------------------*
 |  evaluate Neumann for standard conditions (public)       schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumannStandard( std::multimap<std::string,DRT::Condition* > &   condition,
                                    bool                                  usetime,
                                    const double                          time,
                                    bool                                  assemblemat,
                                    Teuchos::ParameterList&               params,
                                    Teuchos::RCP<DRT::Discretization>     discret,
                                    Epetra_Vector&                        systemvector,
                                    LINALG::SparseOperator*               systemmatrix)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::EvaluateNeumannStandard" );

  std::multimap<std::string,DRT::Condition* >::iterator fool;

  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition.begin(); fool!=condition.end(); ++fool)
  {
    if (fool->first != (std::string)"PointNeumann") continue;
    if (assemblemat && !systemvector.Comm().MyPID())
      std::cout << "WARNING: No linearization of PointNeumann conditions" << std::endl;
    DRT::Condition& cond = *(fool->second);
    const std::vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
    const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");
    const std::vector<double>* val    = cond.Get<std::vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      // call explicitly the main dofset, i.e. the first column
      std::vector<int> dofs = discret->Dof(0,actnode);
      const unsigned numdf = dofs.size();
      for (unsigned j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= curvefac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }

  //--------------------------------------------------------
  // loop through line/surface Neumann BCs and evaluate them
  // ATTENTION: VolumeNeumann conditions (bodyforces) are evaluated in Evaluate
  //--------------------------------------------------------
  for (fool=condition.begin(); fool!=condition.end(); ++fool)
    if (fool->first == (std::string)"LineNeumann" ||
        fool->first == (std::string)"SurfaceNeumann"
    )
    {
      DRT::Condition& cond = *(fool->second);
      std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        curr->second->LocationVector(*discret,lm,lmowner,lmstride);
        elevector.Size((int)lm.size());
        if (!assemblemat)
        {
          curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
        }
        else
        {
          const int size = (int)lm.size();
          if (elematrix.M() != size) elematrix.Shape(size,size);
          else memset(elematrix.A(),0,size*size*sizeof(double));
          curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector,&elematrix);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
          systemmatrix->Assemble(curr->second->Id(),lmstride,elematrix,lm,lmowner);
        }
      }
    }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Neumann for XFEM conditions (public)           schott 09/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumannXFEM( Teuchos::RCP<GEO::CutWizard>         wizard,
                                std::multimap<std::string,DRT::Condition* > &  condition,
                                bool                                 usetime,
                                const double                         time,
                                bool                                 assemblemat,
                                Teuchos::ParameterList&              params,
                                Teuchos::RCP<DRT::Discretization>    discret,
                                Epetra_Vector&                       systemvector,
                                LINALG::SparseOperator*              systemmatrix)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::EvaluateNeumannXFEM" );

  std::multimap<std::string,DRT::Condition* >::iterator fool;

  //--------------------------------------------------------
  // loop through line/surface Neumann BCs and evaluate them
  // ATTENTION: VolumeNeumann conditions (bodyforces) are evaluated in Evaluate
  //--------------------------------------------------------
  for (fool=condition.begin(); fool!=condition.end(); ++fool)
  {
    if (fool->first == (std::string)"PointXFEMNeumann" ||
        fool->first == (std::string)"LineXFEMNeumann" ||
        fool->first == (std::string)"VolXFEMNeumann") dserror(" XFEM Point/Line/Vol Neumann conditions are not implemented yet!");

    if (fool->first == (std::string)"SurfaceXFEMNeumann")
    {
      DRT::Condition& cond = *(fool->second);
      std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;

      // surfaces defined in conditions
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {

          //========================================================
          // use the following evaluate strategy
          // 1. get the parent volume element for this surface
          // 2. a) if(keine volume->cells) => standard element, standard evaluate ohne surface intersection check
          //    b) if(volume->cells && IsCut) => check for surface intersection -> yes: xfem assembly, getrennt für jede vc
          //                                                                    -> no:  vc outside: standard assembly
          //                                                                    -> no:  vc inside:  no assembly
          //    c) if(volume->cells && IsNotCut) => alle volume cells outside   -> standard assembly
          //                                     => alle volume cells inside    -> kein assembly(structure)
          //========================================================

          // 1. get the unique parent volume element of this surface condition
          int ele_id = getParentElementId(discret, curr->second);

          // get the element
          DRT::Element* parent_ele = discret->gElement(ele_id);

          DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( parent_ele );
          if ( ele==NULL ) dserror( "expect fluid element" );

          // ask wizard for this parent element
          GEO::CUT::ElementHandle * e = wizard->GetElement( parent_ele );


          // 2. DECIDE WHICH ASSEMBLY
          // decide if assembly has to be done and decide between standard or xfem
          bool eval_Neumann = true;         // false, if the whole surface lies within the structure (inside)
          bool xfem_eval_Neumann = false;   // true, if the surface is cut by the interface


          //a) standard element, standard evaluate without check for cut surface
          if(e==NULL) // cut does not know this element -> surface is uncut and always "outside/Fluid"
          {
             eval_Neumann = true;
             xfem_eval_Neumann = false;
          }
          if(e!=NULL)
          {
            if(e->IsCut())
            {
              //b) if(volume->cells && IsCut) => check for surface intersection -> yes: xfem assembly, for each vc
              //                                                                -> no:  vc outside: standard assembly
              //                                                                -> no:  vc inside:  no assembly

              //check for intersected Neumann surface and decide for evaluate Neumann
              CutNeumannSurf(curr->second, parent_ele, e, eval_Neumann, xfem_eval_Neumann);

            }
            if(!e->IsCut())
            {
              //c) if(volume->cells && IsNotCut) => all volume cells outside   -> standard assembly
              //                                 => all volume cells inside    -> no assembly(structure)

              GEO::CUT::plain_volumecell_set cells;
              e->GetVolumeCells( cells );

              int count = 0;
              for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
              {
                if(count>0) dserror("more than one volumecell_set in an uncut element");
                GEO::CUT::VolumeCell * vc = *i;
                if ( vc->Position()==GEO::CUT::Point::outside )
                {
                  eval_Neumann = true;
                  xfem_eval_Neumann = false;
                }
                else
                {
                  eval_Neumann = false;
                  xfem_eval_Neumann = false;
                }
                count++;
              }
            } // end if(!e->IsCut())

          } // end case b) and c)


          // standard Neumann assembly
          if(eval_Neumann && !xfem_eval_Neumann)
          {
              // get element location vector, dirichlet flags and ownerships
              std::vector<int> lm;
              std::vector<int> lmowner;
              std::vector<int> lmstride;
              curr->second->LocationVector(*discret,lm,lmowner,lmstride);
              elevector.Size((int)lm.size());
              if (!assemblemat)
              {
                curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector);
                LINALG::Assemble(systemvector,elevector,lm,lmowner);
              }
              else
              {
                const int size = (int)lm.size();
                if (elematrix.M() != size) elematrix.Shape(size,size);
                else memset(elematrix.A(),0,size*size*sizeof(double));
                curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector,&elematrix);
                LINALG::Assemble(systemvector,elevector,lm,lmowner);
                systemmatrix->Assemble(curr->second->Id(),lmstride,elematrix,lm,lmowner);
              }
          }

          // xfem Neumann assembly
          if(eval_Neumann && xfem_eval_Neumann)
          {
            std::cout << "/!\\ WARNING: XFEM_evaluate not implemented yet-> integrate integration cells" << std::endl;
          }

          // Print overview
          std::cout << "Neumann-surface at element : " << ele_id << "\tevalNeumann: " << eval_Neumann << "\txfem_eval: " << xfem_eval_Neumann << std::endl;


      } // end loop condition surfaces
    } // end if "SurfaceXFEMNeumann"
  } // end for conditions

  return;
}



int XFEM::getParentElementId(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<DRT::Element> surf_ele)
{
  int parent_ele_id = -1;

  // get the nodeIds of this surface element
  const int numnode = surf_ele->NumNode();
  const int* nodeids= surf_ele->NodeIds();

  //std::map(eid,#adjacent to nodes)
  std::map<int,int> ele_map;
  std::map<int,int>::iterator ele_map_it;

  // prepare finding the adjacent element for the current Neumann condition
  for(int node_it=0; node_it<numnode; node_it++)
  {
    // get element Ids of current node
    {
      DRT::Node* node = discret->gNode(nodeids[node_it]);

      const int numele = node->NumElement();
      DRT::Element* * elements = node->Elements();


      // find the element which is shared by all nodes
      for(int ele_it=0; ele_it<numele; ele_it++)
      {
        DRT::Element* actele = elements[ele_it];
        int act_eid = actele->Id();
        ele_map_it = ele_map.find(act_eid);
        if(ele_map_it!=ele_map.end()) ele_map_it->second++;
        else ele_map.insert(std::pair<int,int>(act_eid,1));
      }

    }
  }

  // find the adjacent element for the current Neumann condition
  int tmp = -1;
  for(ele_map_it=ele_map.begin(); ele_map_it!=ele_map.end(); ele_map_it++)
  {
    if(ele_map_it->second > tmp)
    {
      tmp = ele_map_it->second;
      parent_ele_id= ele_map_it->first;
    }
  }

  if(parent_ele_id == -1) dserror("the ele-Id of the parent volume element could not be found!");

  return parent_ele_id;
}

/*------------------------------------------------------------------------------*
 |  check for cut Neumann surface and decide for evaluate strategy  schott 09/11|
 *-----------------------------------------------------------------------------*/
void XFEM::CutNeumannSurf(Teuchos::RCP<DRT::Element> neumann_surface, DRT::Element* parentele, GEO::CUT::ElementHandle * parentele_handle, bool & eval_Neumann, bool & xfem_eval_Neumann)
{

  // -> yes: xfem assembly, separated for each vc
  // -> no:  vc outside: standard assembly
  // -> no:  vc inside:  no assembly

  // get the parent element surface which is equal to the Neumann-surface
  // REMARK: this step is needed, because the Neumann condition surfaces do not know their parent element


  // get surfaces of parent element
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(parentele->Surfaces());


  // compare the Neumann surface with the surfaces of the element
  const int  numnode_Neum = neumann_surface->NumNode();
  const int* nodes_Neum   = neumann_surface->NodeIds();

  //  // find local surface id of parent element
  //  int psurf_id = -1;
  //  for(size_t id=0; id<surfaces.size(); id++)
  //  {
  //    const int  numnode_act_surf = surfaces[id]->NumNode();
  //    const int* nodes_act_surf   = surfaces[id]->NodeIds();
  //
  //    if(numnode_act_surf == numnode_Neum)
  //    {
  //      bool psurf_found =  false;
  //      for(int j_node=0; j_node< numnode_Neum; j_node++)
  //      {
  //        // compare the nodeids
  //        if(nodes_act_surf[j_node] == nodes_Neum[j_node]) psurf_found= true;
  //        else
  //        {
  //          // at least one node not identical
  //          psurf_found = false;
  //          continue;
  //        }
  //      }
  //      if(psurf_found) psurf_id = id;
  //    }
  //  }
  //
  //  if(psurf_id == -1) dserror("Neumann surface in element surfaces not found!");
  //
  //  // get the surface of the element identical to the Neumann surface
  //  Teuchos::RCP<DRT::Element> parent_surf = surfaces[psurf_id];


  // find the cut_side, there must be one cut side identical to the Neumann boundary side
  GEO::CUT::plain_element_set ele_set;
  parentele_handle->CollectElements(ele_set);

  int count_collectElements=0;

  // iterator on collected elements for current parent element
  for(GEO::CUT::plain_element_set::iterator it=ele_set.begin(); it!=ele_set.end(); it++)
  {
    if(count_collectElements > 0) dserror("more than one GEO::CUT::Element! -> Quadratic case?");

    // get all cut sides
    const std::vector<GEO::CUT::Side*> & sides = (*it)->Sides();

    int side_id = -1;
    for(size_t i=0; i< sides.size(); i++)
    {
      // get all nodes and there Ids
      std::vector<GEO::CUT::Node*> act_side_nodes = sides[i]->Nodes();
      std::vector<int> act_side_nodeids;

      for(std::vector<GEO::CUT::Node*>::iterator it_node=act_side_nodes.begin(); it_node!=act_side_nodes.end(); it_node++)
      {
        act_side_nodeids.push_back( (*it_node)->Id() );
      }


      // compare surface_nodeids with side_nodeids
      if(numnode_Neum == (int)act_side_nodeids.size())// check same number of nodes
      {
        for(int i=0; i<numnode_Neum; i++)
        {
          bool side_found =  false;

          for(int j_node=0; j_node< numnode_Neum; j_node++)
          {
            // compare the nodeids
            if(act_side_nodeids[j_node] == nodes_Neum[j_node]) side_found= true;
            else
            {
              // at least one node not identical
              side_found = false;
              continue;
            }
          }

          if(side_found) side_id = i;
        }
      }

    }// end sides

    if(side_id==-1) dserror("no side equal to the Neumann-boundary surface found!");



    // decide betwenn assembly strategies
    if(sides[side_id]->IsCut()) // b) -> yes: xfem assembly, getrennt für jede vc
    {
      eval_Neumann = true;
      xfem_eval_Neumann = true;
    }
    if(!(sides[side_id]->IsCut()))
    {
      // all facets outside: standard assembly
      // all facets inside:  no assembly

      const std::vector<GEO::CUT::Facet*> & facets = sides[side_id]->Facets();

      int outside_facets = 0;
      int inside_facets = 0;

      for(int f_index=0; f_index< (int)facets.size(); f_index++)
      {
        GEO::CUT::Facet* f = facets[f_index];

        if(f->Position()==GEO::CUT::Point::outside) outside_facets++;
        if(f->Position()!=GEO::CUT::Point::outside) inside_facets++;
      }

      if(outside_facets>0 && inside_facets<=0)
      {
        // standard assembly
        eval_Neumann = true;
        xfem_eval_Neumann = false;
      }
      else if(outside_facets<=0 && inside_facets>0)
      {
        // no Neumann evaluate
        eval_Neumann = false;
        xfem_eval_Neumann = false;
      }
      else dserror("this can not happen!");

    }

    count_collectElements++;
  }

//  {
//
//    GEO::CUT::plain_volumecell_set cells;
//    parentele_handle->VolumeCells(cells);
//
//    for(GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end();i++)
//    {
//      GEO::CUT::VolumeCell * vc = *i;
//      if(vc->Position() == GEO::CUT::Point::outside)
//      {
//        GEO::CUT::plain_integrationcell_set int_cells;
//        vc->GetIntegrationCells(int_cells);
//
//        for(GEO::CUT::plain_integrationcell_set::iterator j=int_cells.begin(); j!=int_cells.end(); j++)
//        {
//          GEO::CUT::IntegrationCell* intcell = *j;
//
//          if(intcell->Shape() == DRT::Element::tet4) std::cout << "tet4 integration cell" << std::endl;
//
//          if(intcell->Shape() == DRT::Element::hex8) std::cout << "hex8 integration cell" << std::endl;
//
//          std::vector<GEO::CUT::Point*> points = intcell->Points();
//          //  for(int i=0; i<points.size(); i++)
//          //  {
//          //    GEO::CUT::Point* point = points[i];
//          //
//          //  }
//          GEO::CUT::Facet* facet;
//
//          const std::vector<GEO::CUT::Facet*> & facets = sides[side_id]->Facets();
//
//          for(int f_index=0; f_index< (int)facets.size(); f_index++)
//          {
//            GEO::CUT::Facet* f = facets[f_index];
//
//            if(f->Position()==GEO::CUT::Point::outside) facet = f;
//          }
//
//          if(intcell->Shape() == DRT::Element::tet4)
//          {
//            //  vc->NewBoundaryCell(wizard_.CutWizard().Mesh().NormalMesh(), DRT::Element::tri3, facet, points );
//            //  vc->NewTri3Cell(wizard_.CutWizard().Mesh().NormalMesh(), facet, points );
//            std::cout << "new tri3 cell created" << std::endl;
//            vc->NewIntegrationCell(wizard_.CutWizard().Mesh().NormalMesh(), DRT::Element::tet4, points );
//
//          }
//          if(intcell->Shape() == DRT::Element::hex8)
//          {
//            //  vc->NewBoundaryCell(wizard_.CutWizard().Mesh().NormalMesh(), DRT::Element::quad4, facet, points);
//            std::cout << "new quad4 cell created" << std::endl;
//          }
//        }
//
//
//
//
//      }
//
//    }
//
//
//  }

  return;
}

