/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret_evaluate.cpp
 *
 * \brief overloaded evaluation methods for meshfree discretisations
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"
#include "drt_meshfree_utils.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

void DRT::MESHFREE::MeshfreeDiscretization::EvaluateDirichlet(
  Teuchos::ParameterList& params,
  Teuchos::RCP<Epetra_Vector> sysvec,
  Teuchos::RCP<Epetra_Vector> sysvecd,
  Teuchos::RCP<Epetra_Vector> sysvecdd,
  Teuchos::RCP<Epetra_Vector> toggle,
  Teuchos::RCP<LINALG::MapExtractor> inputmapextractor)
{
  //----------------------------------------------------------------------------
  // check that no time derivatives prescribed
  //----------------------------------------------------------------------------
  if (sysvecd!=Teuchos::null or sysvecdd!=Teuchos::null)
    dserror("Meshfree Dirichlet BC are not tested for prescribing time derivatives. Remove at own risk.");

  //----------------------------------------------------------------------------
  // get values at nodes for DBC
  //----------------------------------------------------------------------------

  // create map extractor to always (re)build dbcmapextractor which is needed later
  Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor = Teuchos::rcp(new LINALG::MapExtractor());
  // call base class EvaluateDirichlet
  DRT::Discretization::EvaluateDirichlet(params,sysvec,sysvecd,sysvecdd,toggle,dbcmapextractor);
  // pass dbcmapextractor if rebuild was required
  if (inputmapextractor!=Teuchos::null)
    *inputmapextractor = *dbcmapextractor;

  //----------------------------------------------------------------------------
  // compute nodal values for DBC
  //----------------------------------------------------------------------------

  // iterator over conditions
  std::multimap<std::string,Teuchos::RCP<Condition> >::iterator fool;
  // number of degrees of freedom of Dirichlet BCs
  int numdof = -1;
  // flag for constant DBC
  bool isconst = true;
  // constant value of DBC
  std::vector<const double*> constvalue;
  // create set to memorize dof-gids of lower-dimensional faces
  std::set<int> nodegids;
  // create matrix for dependencies
  const Teuchos::RCP<const Epetra_Map> dbcdofmap = dbcmapextractor->CondMap();
  Teuchos::RCP<LINALG::SparseMatrix>   dbcmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*dbcdofmap,81,true,true));

  //----------------------------------------------------------------------------
  // loop through Dirichlet conditions and and fill dbcmatrix
  //----------------------------------------------------------------------------
  // Due to the weak Kronecker-delta property, basis function may only be
  // evaluate in a relative interior or at vertices (where they are 1). Thus
  // we need to evaluate in reverse order, starting from vertices to volumes:
  //                PointDirichlet
  //                LineDirichlet
  //                SurfaceDirichlet
  //                VolumeDirichlet
  // Furthermore, points, lines, surfaces, and volumes have to coincide with
  // the faces of the convex hull of nodes. Not checked - own responsibility!!

  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    // select PointDirichlet conditions only
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    if (numdof==-1)
    {
      numdof = fool->second->GetInt("numdof");
      constvalue = std::vector<const double*>(numdof,NULL);
    }
    if (numdof!=fool->second->GetInt("numdof"))
      dserror("NUMDOF must be the same for all meshfree Dirichlet BCs! Check dat-file!");
    FillDBCMatrix(*(fool->second), 0, numdof, nodegids, isconst, constvalue, dbcmatrix, dbcdofmap);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    {
      numdof = fool->second->GetInt("numdof");
      constvalue = std::vector<const double*>(numdof,NULL);
    }
    if (numdof!=fool->second->GetInt("numdof"))
      dserror("NUMDOF must be the same for all meshfree Dirichlet BCs! Check dat-file!");
    FillDBCMatrix(*(fool->second), 1, numdof, nodegids, isconst, constvalue, dbcmatrix, dbcdofmap);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    {
      numdof = fool->second->GetInt("numdof");
      constvalue = std::vector<const double*>(numdof,NULL);
    }
    if (numdof!=fool->second->GetInt("numdof"))
      dserror("NUMDOF must be the same for all meshfree Dirichlet BCs! Check dat-file!");
    FillDBCMatrix(*(fool->second), 2, numdof, nodegids, isconst, constvalue, dbcmatrix, dbcdofmap);
  }
  // Do VolumeDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    {
      numdof = fool->second->GetInt("numdof");
      constvalue = std::vector<const double*>(numdof,NULL);
    }
    if (numdof!=fool->second->GetInt("numdof"))
      dserror("NUMDOF must be the same for all meshfree Dirichlet BCs! Check dat-file!");
    FillDBCMatrix(*(fool->second), 3, numdof, nodegids, isconst, constvalue, dbcmatrix, dbcdofmap);
  }

  // make dbcmatrix ready for solver
  dbcmatrix->Complete();

  // if necessary, solve for nodal values and write them to respective vector
  if (sysvec!=Teuchos::null and !isconst)
  {
    // extract dbc-values at nodes from sysvec and copy to nadalvalues as first guess
    Teuchos::RCP<Epetra_Vector> valuesatnodes = dbcmapextractor->ExtractCondVector(sysvec);
    Teuchos::RCP<Epetra_Vector> nodalvalues = Teuchos::rcp(new Epetra_Vector(*(valuesatnodes)));

    // solve for nodal values
    if (dbcsolver_==Teuchos::null) dserror("No DBC_SOLVER defined for non-constant meshfree DBC. Check MESHFREE section in dat-file!");
    dbcsolver_->Solve(dbcmatrix->EpetraMatrix(),nodalvalues,valuesatnodes,true,true);

    // insert nodal values to sysvec
    dbcmapextractor->InsertCondVector(nodalvalues,sysvec);
  }
  if (sysvecd!=Teuchos::null and !isconst)
  {
    // extract dbc-values at nodes from sysvecd and copy to nadalvalues as first guess
    Teuchos::RCP<Epetra_Vector> valuesatnodes = dbcmapextractor->ExtractCondVector(sysvec);
    Teuchos::RCP<Epetra_Vector> nodalvalues = Teuchos::rcp(new Epetra_Vector(*(valuesatnodes)));

    // solve for nodal values
    if (dbcsolver_==Teuchos::null) dserror("No DBC_SOLVER defined for non-constant meshfree DBC. Check MESHFREE section in dat-file!");
    dbcsolver_->Solve(dbcmatrix->EpetraMatrix(),nodalvalues,valuesatnodes,true,true);

    // insert nodal values to sysvecd
    dbcmapextractor->InsertCondVector(nodalvalues,sysvecd);
  }
  if (sysvecdd!=Teuchos::null and !isconst)
  {
    // extract dbc-values at nodes from sysvecdd and copy to nadalvalues as first guess
    Teuchos::RCP<Epetra_Vector> valuesatnodes = dbcmapextractor->ExtractCondVector(sysvec);
    Teuchos::RCP<Epetra_Vector> nodalvalues = Teuchos::rcp(new Epetra_Vector(*(valuesatnodes)));

    // solve for nodal values
    if (dbcsolver_==Teuchos::null) dserror("No DBC_SOLVER defined for non-constant meshfree DBC. Check MESHFREE section in dat-file!");
    dbcsolver_->Solve(dbcmatrix->EpetraMatrix(),nodalvalues,valuesatnodes,true,true);

    // insert nodal values to sysvecd
    dbcmapextractor->InsertCondVector(nodalvalues,sysvecdd);
  }

  return;
}

//------------------------------------------------------------------------------
// fill Dirichlet BC to compute nodal values from values at nodes if necessary
//------------------------------------------------------------------------------
void DRT::MESHFREE::MeshfreeDiscretization::FillDBCMatrix(
  const DRT::Condition& cond,
  const int             facedim,
  const int             numdof,
  std::set<int>&        nodegids,
  bool&                 isconst,
  std::vector<const double*>&          constvalue,
  Teuchos::RCP<LINALG::SparseMatrix>   dbcmatrix,
  const Teuchos::RCP<const Epetra_Map> dbcdofmap)
{
  //----------------------------------------------------------------------------
  // get condition and relevant information
  //----------------------------------------------------------------------------

  const std::vector<int>*    nodeids = cond.Nodes();
  const std::vector<int>*    onoff   = cond.Get<std::vector<int> >("onoff");
  const std::vector<int>*    funct   = cond.Get<std::vector<int> >("funct");
  const std::vector<double>* val     = cond.Get<std::vector<double> >("val");

  //----------------------------------------------------------------------------
  // declare variables needed in loop
  //----------------------------------------------------------------------------

  // prepare auxiliary variables
  int nneighbour = -1;
  bool haveneighbourhood = false;

  // node positions
  Teuchos::RCP<LINALG::SerialDenseMatrix> nxyz = Teuchos::null;
  // node positions the way a searchtree needs it
  std::map<int,LINALG::Matrix<3,1> > nxyz_searchtree;
  // reduced node positions on face
  LINALG::SerialDenseMatrix nxyz_face;

  // prepare vectors of node gids and basis function values for non-constant DBC
  std::vector<double> distnn;
  LINALG::SerialDenseMatrix distnn_sdm;
  std::vector<int> ngids;
  LINALG::SerialDenseVector basisfunct;
  Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::null;

  // prepare vectors for constant (or point) DBC
  std::vector<int> const_ngids(1,-1);
  LINALG::SerialDenseVector const_basisfunct(1);
  const_basisfunct(0) = 1;

  // pointers to node gids and basis function values
  std::vector<int>* temp_ngids;
  LINALG::SerialDenseVector* temp_basisfunct;

  // prepare vector for dof-gids
  std::vector<int> gids;

  //----------------------------------------------------------------------------
  // loop over all nodes of this condition
  //----------------------------------------------------------------------------
  const int nnode = (*nodeids).size();
  for (int inode=0; inode<nnode; ++inode)
  {
    // do only nodes in my row map
    const int ngid = (*nodeids)[inode];
    if (!(this->NodeRowMap()->MyGID(ngid))) continue;

    // do only nodes that are not done yet
    if (nodegids.find(ngid)!=nodegids.end()) continue;

    // get active node and list node gid
    // (the latter not done for surfaces/volumes in 2D/3D, respectively)
    DRT::Node* actnode = this->gNode(ngid);
    if (facedim!=DRT::Problem::Instance()->NDim()) nodegids.insert(ngid);

    // loop over all dofs of this node
    for (int idof=0; idof<numdof; ++idof)
    {
      if ((*onoff)[idof]!=0)
      {
        // determine node gids and basis functions necessary to fill dbcmatrix:
        // ====================================================================
        // this rather awkward implementation with temp_x, const_x, and
        // x-variables was done to only compute once the neighbourhood and
        // basis function values for each node and this only for non-constant
        // DBCs. For constant DBCs, the diagonal entry is set to 1 and nodal
        // values are equal to values at nodes.

        // check for constant DBC
        if (constvalue[idof]==NULL)
          constvalue[idof] = &(*val)[idof];
        if ((*(constvalue[idof]) != (*val)[idof] or (*funct)[idof]!=0) and facedim!=0)
        {
          //----------------------------------------------------------------
          // non-constant DBC:
          //----------------------------------------------------------------

          isconst = false;

          // get range of basis solution functions - declared outside for dserror
          double range = solutionapprox_->GetRange();

          // compute nodal neighbourhood and basis functions if not already done
          if (!haveneighbourhood)
          {
            distnn.clear();
            ngids.clear();
            haveneighbourhood = true;

            if (facedim==1)
            {
              //----------------------------------------------------------------
              // Line DBC - non-constant
              //----------------------------------------------------------------

              // create SerialDenseMatrix with all node positions if not already done
              if (nxyz == Teuchos::null)
              {
                nxyz = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,nnode));
                for (int i=0; i<nnode; ++i)
                  for (int j=0; j<3; ++j)
                    (*nxyz)(j,i) = this->gNode((*nodeids)[i])->X()[j];

                // get reduced node positions
                nxyz_face.LightShape(facedim,nnode);
                ReduceDimensionOfFaceNodes(*nxyz, nxyz_face);
              }

              // brute force neighbourhood search
              const double actnodexyz = nxyz_face(0,inode);
              for (int i=0; i<nnode; ++i)
              {
                const double dist = nxyz_face(0,i)-actnodexyz;
                if (std::abs(dist)<range)
                {
                  distnn.push_back(dist);
                  ngids.push_back((*nodeids)[i]);
                }
              }

              nneighbour = ngids.size();
              if ((int)(distnn.size())!=nneighbour) dserror("Something went seriously wrong!");

              distnn_sdm.LightShape(1,nneighbour);
              distnn_sdm = LINALG::SerialDenseMatrix(Copy,distnn.data(),1,1,nneighbour);
            }
            else
            {
              //----------------------------------------------------------------
              // Surface and Volume DBC - non-constant
              //----------------------------------------------------------------

              // initialize searchtree if not already done
              if (searchTree == Teuchos::null)
              {
                // get node coordinates searchtree-style
                nxyz_searchtree.clear();
                for (int i=0; i<nnode; ++i)
                  nxyz_searchtree.insert( std::pair<int,LINALG::Matrix<3,1> >((*nodeids)[i],LINALG::Matrix<3,1>(this->gNode((*nodeids)[i])->X())));

                // set up searchtree - always OCTREE although for surfaces a QUADTREE would be optimal
                searchTree = Teuchos::rcp(new GEO::SearchTree(8));
                searchTree->initializePointTree(GEO::getXAABBofPositions(nxyz_searchtree), nxyz_searchtree, GEO::TreeType(GEO::OCTTREE));
              }

              // find neighbours of active node
              const LINALG::Matrix<3,1> actnodexyz(this->gNode((*nodeids)[inode])->X());
              ngids = searchTree->searchPointsInRadius(nxyz_searchtree, actnodexyz, range);
              nneighbour = ngids.size();

              // setup SerialDenseMatrix with positions of neighbours
              nxyz = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,nneighbour));
              for (int i=0; i<nneighbour; ++i)
              {
                const double* neighbourxyz = this->gNode(ngids[i])->X();
                for (int j=0; j<3; ++j)
                  (*nxyz)(j,i) = neighbourxyz[j];
              }

              std::vector<int> dims(facedim);
              // get reduced node positions if necessary
              if (facedim!=3)
              {
                nxyz_face.LightShape(facedim,nneighbour);
                dims = ReduceDimensionOfFaceNodes(*nxyz, nxyz_face);
              }
              else
              {
                dserror("Meshfree Dirichlet BCs not yet tested for volumes. Should work. Remove dserror at own risk.");

                for (int i=0; i<facedim; ++i)
                  dims[i] = i;
                nxyz_face = LINALG::SerialDenseMatrix(View,nxyz->A(),1,facedim,nneighbour);
              }

              // get distnn-matrix
              distnn_sdm.LightShape(facedim, nneighbour);
              for (int i=0; i<nneighbour; ++i)
                for (int j=0; j<facedim; ++j)
                  distnn_sdm(j,i) = nxyz_face(j,i) - actnodexyz(dims[j]);
            }
          }

          // check whether enough neighbours were found
          nneighbour = ngids.size(); // needs to be repeated - trust me!
          if(nneighbour < (facedim+1))
            dserror("Only %i node(s) found in neighbourhood (range=%g) of node#%i/dof#%i on %i-dimensional face.",nneighbour,ngid,idof,facedim,range);

          // get basis function values of neighbours
          basisfunct.LightSize(nneighbour);
          LINALG::SerialDenseMatrix dummy(0,0);
          this->GetMeshfreeSolutionApprox()->GetMeshfreeBasisFunction(basisfunct,dummy,distnn_sdm,facedim);

          // use neighbourhood-variables for node gids and basis functions
          temp_ngids = &ngids;
          temp_basisfunct = &basisfunct;
        }
        else
        {
          //----------------------------------------------------------------
          // constant DBC or point DBC:
          //----------------------------------------------------------------

          // use const_x-variables for node gids and basis functions
          nneighbour = 1;
          const_ngids[0] = ngid;
          temp_ngids = &const_ngids;
          temp_basisfunct = &const_basisfunct;
        }

        // get dof-gid of current dof of active node
        const int gid = this->Dof(0,actnode)[idof];
        if (gid==-1) dserror("Something went wrong: dof-gid is not in dbcdofmap.");

        // get dof-gids of current dofs of neighbouring nodes
        gids.resize(nneighbour);
        for(int ineighbour=0; ineighbour<nneighbour; ++ineighbour)
          gids[ineighbour] = this->Dof(0,this->gNode((*temp_ngids)[ineighbour]))[idof];
        // add values in dbcmatrix
        dbcmatrix->EpetraMatrix()->InsertGlobalValues(gid,nneighbour,temp_basisfunct->Values(),gids.data());
      }
    }

    // reset nneighbour to compute neighbourhood
    haveneighbourhood = false;
  }

  return;
}

