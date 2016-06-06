/*----------------------------------------------------------------------*/
/*!
\file drt_nurbs_discret.cpp

\brief a class to manage one nurbs discretization

\maintainer Anh-Tu Vuong

\level 1

*/
/*----------------------------------------------------------------------*/

#include <Epetra_Vector.h>
#include <Epetra_Time.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "drt_nurbs_utils.H"
#include "drt_nurbs_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::NurbsDiscretization(
  const std::string             name,
  Teuchos::RCP<Epetra_Comm> comm)
  :
  DRT::Discretization::Discretization(name,comm    ),
  npatches_                          (            0),
  knots_                             (Teuchos::null)
{
//  dbcsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(1),
//                                               this->Comm(),
//                                               DRT::Problem::Instance()->ErrorFile()->Handle()));
//  this->ComputeNullSpaceIfNecessary(dbcsolver_->Params());

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::~NurbsDiscretization()
{
  return;
}


/*----------------------------------------------------------------------*
 |  add a knotvector to the discretization (public)          gammi 05/08|
 *----------------------------------------------------------------------*/
void
DRT::NURBS::NurbsDiscretization::SetKnotVector
(Teuchos::RCP<DRT::NURBS::Knotvector> knots)
{

  if(knots==Teuchos::null)
  {
    dserror("trying to set invalid knotvector (%s)\n",(this->Name()).c_str());
  }

  knots_=knots;
  return;
}

/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |                                                           gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::NURBS::Knotvector>
DRT::NURBS::NurbsDiscretization::GetKnotVector
()
{
  if(knots_==Teuchos::null)
  {
    dserror("knotvector invalid (%s)\n",(this->Name()).c_str());
  }
  return knots_;
}


/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |  (const version, read-only)                               gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const DRT::NURBS::Knotvector>
DRT::NURBS::NurbsDiscretization::GetKnotVector() const
{
  if(knots_==Teuchos::null)
  {
    dserror("knotvector invalid (%s)\n",(this->Name()).c_str());
  }
  return knots_;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
void DRT::NURBS::NurbsDiscretization::EvaluateDirichlet(Teuchos::ParameterList& params,
                                            Teuchos::RCP<Epetra_Vector> systemvector,
                                            Teuchos::RCP<Epetra_Vector> systemvectord,
                                            Teuchos::RCP<Epetra_Vector> systemvectordd,
                                            Teuchos::RCP<Epetra_Vector> toggle,
                                            Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{


  // call base class EvaluateDirichlet
  //Discretization::EvaluateDirichlet(params,systemvector,systemvectord,systemvectordd,toggle,dbcmapextractor);

  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::null;
  if (dbcmapextractor != Teuchos::null) dbcgids = Teuchos::rcp(new std::set<int>());

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  std::multimap<std::string,Teuchos::RCP<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in systemvector.
  // For this reason, Dirichlet BCs are evaluated hierarchical meaning
  // in this order:
  //                VolumeDirichlet
  //                SurfaceDirichlet
  //                LineDirichlet
  //                PointDirichlet
  // This way, lower entities override higher ones which is
  // equivalent to inheritance of dirichlet BCs as done in the old
  // ccarat discretization with design          (mgee 1/07)

  // Do VolumeDirichlet first
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    if (fool->first == "Dirichlet")
      Discretization::DoDirichletCondition(
                           *(fool->second),
                           usetime,
                           time,
                           systemvector,
                           systemvectord,
                           systemvectordd,
                           toggle,
                           dbcgids);

    else if (fool->first == "NurbsLSDirichlet")
    {
      FindDBCgidAndToggle(
          *(fool->second),
          systemvector,
          systemvectord,
          systemvectordd,
          toggle,
          dbcgids);

      DoNurbsLSDirichletCondition(
          *(fool->second),
          usetime,
          time,
          systemvector,
          systemvectord,
          systemvectordd);
    }
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    if (fool->first == "Dirichlet")
      Discretization::DoDirichletCondition(
                           *(fool->second),
                           usetime,
                           time,
                           systemvector,
                           systemvectord,
                           systemvectordd,
                           toggle,
                           dbcgids);

    else if (fool->first == "NurbsLSDirichlet")
    {
      FindDBCgidAndToggle(
          *(fool->second),
          systemvector,
          systemvectord,
          systemvectordd,
          toggle,
          dbcgids);

      DoNurbsLSDirichletCondition(
          *(fool->second),
          usetime,
          time,
          systemvector,
          systemvectord,
          systemvectordd);
    }
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    if (fool->first == "Dirichlet")
      Discretization::DoDirichletCondition(
                           *(fool->second),
                           usetime,
                           time,
                           systemvector,
                           systemvectord,
                           systemvectordd,
                           toggle,
                           dbcgids);

    else if (fool->first == "NurbsLSDirichlet")
    {
      FindDBCgidAndToggle(
          *(fool->second),
          systemvector,
          systemvectord,
          systemvectordd,
          toggle,
          dbcgids);

      DoNurbsLSDirichletCondition(
          *(fool->second),
          usetime,
          time,
          systemvector,
          systemvectord,
          systemvectordd);
    }
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    if (fool->first == "Dirichlet")
      Discretization::DoDirichletCondition(
                           *(fool->second),
                           usetime,
                           time,
                           systemvector,
                           systemvectord,
                           systemvectordd,
                           toggle,
                           dbcgids);

    else if (fool->first == "NurbsLSDirichlet")
    {
      FindDBCgidAndToggle(
          *(fool->second),
          systemvector,
          systemvectord,
          systemvectordd,
          toggle,
          dbcgids);

      DoNurbsLSDirichletCondition(
          *(fool->second),
          usetime,
          time,
          systemvector,
          systemvectord,
          systemvectordd);
    }
  }

  // create DBC and free map and build their common extractor
  if (dbcmapextractor != Teuchos::null)
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids->size() > 0)
    {
      dbcgidsv.reserve(dbcgids->size());
      dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap
      = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DofRowMap()->IndexBase(), DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *dbcmapextractor = LINALG::MapExtractor(*(DofRowMap()), dbcmap);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
void DRT::NURBS::NurbsDiscretization::DoNurbsLSDirichletCondition(
  DRT::Condition&              cond,
  const bool                   usetime,
  const double                 time,
  Teuchos::RCP<Epetra_Vector>  systemvector,
  Teuchos::RCP<Epetra_Vector>  systemvectord,
  Teuchos::RCP<Epetra_Vector>  systemvectordd)
{
  Epetra_Time timer(Comm());

  // get the processor ID from the communicator
  const int myrank  = Comm().MyPID();
  if(myrank==0)
    std::cout << "calculating least squares Dirichlet condition in ... ";

  Teuchos::RCP<std::set<int> > nurbslsdbcgids = Teuchos::rcp(new std::set<int>());
  //for integration over elements with DBC we need the column map
  Teuchos::RCP<std::set<int> > nurbslsdbccolgids = Teuchos::rcp(new std::set<int>());

  FindDBCgidAndToggle(cond,
      systemvector,systemvectord,systemvectordd,
      Teuchos::null,nurbslsdbcgids,nurbslsdbccolgids);

  // create map extractor to always (re)build dbcmapextractor which is needed later
  Teuchos::RCP<LINALG::MapExtractor> auxdbcmapextractor = Teuchos::rcp(new LINALG::MapExtractor());
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (nurbslsdbcgids->size() > 0)
    {
      dbcgidsv.reserve(nurbslsdbcgids->size());
      dbcgidsv.assign(nurbslsdbcgids->begin(),nurbslsdbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap
      = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DofRowMap()->IndexBase(), DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *auxdbcmapextractor = LINALG::MapExtractor(*(DofRowMap()), dbcmap);
  }

  //column map of all DOFs subjected to a least squares Dirichlet condition
  Teuchos::RCP<Epetra_Map> dbccolmap=Teuchos::null;
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (nurbslsdbccolgids->size() > 0)
    {
      dbcgidsv.reserve(nurbslsdbccolgids->size());
      dbcgidsv.assign(nurbslsdbccolgids->begin(),nurbslsdbccolgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    dbccolmap
      = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DofRowMap()->IndexBase(), DofRowMap()->Comm()));
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Teuchos::RCP<const Epetra_Map> dofrowmap = auxdbcmapextractor->CondMap();

  if(dofrowmap->NumGlobalElements()==0)
    return;//no dbc gids ->leave

  //read information from condition
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
  const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
  const std::vector<double>* val    = cond.Get<std::vector<double> >("val");


  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");


  // -------------------------------------------------------------------
  // create empty mass matrix
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix
    = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));

  // -------------------------------------------------------------------
  // create empty right hand side vector
  // -------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector>  dbcvector = LINALG::CreateVector(*dofrowmap,true);

  Teuchos::RCP<Epetra_Vector> rhsd=Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectord=Teuchos::null;
  if (systemvectord != Teuchos::null)
  {
    rhsd = LINALG::CreateVector(*dofrowmap,true);
    dbcvectord = LINALG::CreateVector(*dofrowmap,true);
  }

  Teuchos::RCP<Epetra_Vector> rhsdd=Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectordd=Teuchos::null;
  if (systemvectord != Teuchos::null)
  {
    rhsdd = LINALG::CreateVector(*dofrowmap,true);
    dbcvectordd = LINALG::CreateVector(*dofrowmap,true);
  }

  const bool assemblevecd = rhsd       !=Teuchos::null;
  const bool assemblevecdd = rhsdd       !=Teuchos::null;

  // -------------------------------------------------------------------
  // call elements to calculate massmatrix and righthandside
  // -------------------------------------------------------------------
  {
    // call elements and assemble
    if (!Filled())   dserror("FillComplete() was not called");
    if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    // see what we have for input
    bool assemblemat = massmatrix!=Teuchos::null;
    bool assemblevec = rhs       !=Teuchos::null;

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elemass;
    //Epetra_SerialDenseVector elerhs;
    std::vector<Epetra_SerialDenseVector> elerhs(deg+1);

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;

    std::vector<int> lm_full;
    std::vector<int> lmowner_full;
    std::vector<int> lmstride_full;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      Teuchos::RCP<DRT::Element> actele = curr->second;

      static const int probdim = DRT::Problem::Instance()->NDim();
      const DRT::Element::DiscretizationType distype = actele->Shape();
      const int dim = DRT::UTILS::getDimension(distype);
      const bool isboundary = (dim!=probdim);
      const int nen = DRT::UTILS::getNumberOfElementNodes(distype);

      // access elements knot span
      std::vector<Epetra_SerialDenseVector> eleknots(dim);
      Epetra_SerialDenseVector weights(nen);

      bool zero_size = false;
      if(isboundary)
      {
        Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(actele,true);
        double normalfac = 0.0;
        std::vector<Epetra_SerialDenseVector> pknots(probdim);
        zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
            actele.get(), faceele->FaceMasterNumber(), faceele->ParentElement()->Id(), *this, pknots, eleknots, weights, normalfac);
      }
      else
        zero_size =DRT::NURBS::GetMyNurbsKnotsAndWeights(*this,actele.get(),eleknots,weights);

      // nothing to be done for a zero sized element
      if(zero_size)
      {
        continue;
      }

      // get element full location vector, dirichlet flags and ownerships
      lm_full.clear();
      lmowner_full.clear();
      lmstride_full.clear();
      actele->LocationVector(*this,lm_full,lmowner_full,lmstride_full);

      //we are only interested in DOFs with dirichlet condition, hence we compare the location vector with the
      // drichlet condition map
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      for(unsigned j=0;j<lm_full.size();++j)
      {
        int gid =lm_full[j];
        if(dbccolmap->MyGID(gid))
        {
          lm.push_back(gid);
          lmowner.push_back(lmowner_full[j]);
          lmstride.push_back(lmstride_full[j]);
        }
      }

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      if (assemblemat)
      {
        if (elemass.M()!=eledim or elemass.N()!=eledim)
        elemass.Shape(eledim,eledim);
        else
        memset(elemass.A(),0,eledim*eledim*sizeof(double));
      }
      if (assemblevec)
      {
        for (unsigned i=0; i<deg+1; ++i)
        {
          if (elerhs[i].Length()!=eledim)
          elerhs[i].Size(eledim);
          else
          memset(elerhs[i].Values(),0,eledim*sizeof(double));
        }
      }

      if(isboundary)
        switch(distype)
        {
        case DRT::Element::nurbs2:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs2>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs3:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs3>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs4:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs4>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs9:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs9>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs8:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs8>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs27:
          FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs27>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        default:
          dserror("invalid element shape for least squares dirichlet evaluation: %s",DistypeToString(distype).c_str());
          break;
        }
      else
        switch(distype)
        {
        case DRT::Element::nurbs2:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs2>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs3:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs3>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs4:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs4>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs9:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs9>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs8:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs8>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        case DRT::Element::nurbs27:
          FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs27>(actele,&eleknots,lm,funct,curve,val,usetime,deg,time,elemass,elerhs);
          break;
        default:
          dserror("invalid element shape for least squares dirichlet evaluation: %s",DistypeToString(distype).c_str());
          break;
        }

      int eid = actele->Id();
      if (assemblemat) massmatrix->Assemble(eid,elemass,lm,lmowner);
      if (assemblevec) LINALG::Assemble(*rhs,elerhs[0],lm,lmowner);
      if (assemblevecd) LINALG::Assemble(*rhsd,elerhs[1],lm,lmowner);
      if (assemblevecdd) LINALG::Assemble(*rhsdd,elerhs[2],lm,lmowner);

    }
  }


  // -------------------------------------------------------------------
  // finalize the system matrix
  // -------------------------------------------------------------------
  massmatrix->Complete();

  // -------------------------------------------------------------------
  // solve system
  // -------------------------------------------------------------------

  // always refactor and reset the matrix before a single new solver call
  bool refactor=true;
  bool reset   =true;

  // Owing to experience a very accurate solution has to be enforced here!
  // Thus, we allocate an own solver with VERY strict tolerance!
  // One could think of verifiying an extra solver in the input file...
  Teuchos::ParameterList p = DRT::Problem::Instance()->UMFPACKSolverParams();
//  const double origtol = p.get<double>("AZTOL");
//  const double newtol  = 1.0e-11;
//  p.set("AZTOL",newtol);

//  if(myrank==0)
//    cout<<"\nSolver tolerance for least squares problem set to "<<newtol<<"\n";

  Teuchos::RCP<LINALG::Solver> solver =
      Teuchos::rcp(new LINALG::Solver(p,
          Comm(),
          DRT::Problem::Instance()->ErrorFile()->Handle()));
  ComputeNullSpaceIfNecessary(solver->Params());

  //solve for control point values
  solver->Solve(massmatrix->EpetraOperator(),
               dbcvector                 ,
               rhs                         ,
               refactor                    ,
               reset                       );

  //solve for first derivatives in time
  if (assemblevecd)
    solver->Solve(massmatrix->EpetraOperator(),
                 dbcvectord                 ,
                 rhsd                         ,
                 refactor                    ,
                 reset                       );

  //solve for second derivatives in time
  if (assemblevecdd)
    solver->Solve(massmatrix->EpetraOperator(),
                 dbcvectordd                 ,
                 rhsdd                         ,
                 refactor                    ,
                 reset                       );

  // perform resets for solver and matrix
  solver->Reset();
  massmatrix->Reset();

  // insert nodal values to sysvec
  auxdbcmapextractor->InsertCondVector(dbcvector,systemvector);
  if (assemblevecd) auxdbcmapextractor->InsertCondVector(dbcvectord,systemvectord);
  if (assemblevecdd) auxdbcmapextractor->InsertCondVector(dbcvectordd,systemvectordd);

  if(myrank==0)
    std::cout << timer.ElapsedTime() << " seconds \n\n";

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
void DRT::NURBS::NurbsDiscretization::FindDBCgidAndToggle(
    DRT::Condition&              cond,
    const Teuchos::RCP<Epetra_Vector>  systemvector,
    const Teuchos::RCP<Epetra_Vector>  systemvectord,
    const Teuchos::RCP<Epetra_Vector>  systemvectordd,
    Teuchos::RCP<Epetra_Vector>        toggle,
    Teuchos::RCP<std::set<int> >       dbcgids,
    Teuchos::RCP<std::set<int> >       dbccolgids)
{
  const bool findcolgids = (dbccolgids!=Teuchos::null);
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    bool iscol=false;

    DRT::Node* actnode = NULL;

    // do only nodes in my row map
    int nlid = this->NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0)
    {
      if(not findcolgids) continue;  //not in row map and column dofs not needed -> next node
      //check if node is in col map
      else
      {
        // do nodes in my col map
        nlid = this->NodeColMap()->LID((*nodeids)[i]);
        if (nlid < 0) continue;   //node not on this processor -> next node
        iscol =true;
        //get node from col node list
        actnode = this->lColNode( nlid );
      }
    }
    else
      //get node from row node list
      actnode = this->lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = this->Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // only continue if there are node dofs
    if (total_numdf == 0) continue;

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    if(not iscol)
      for (unsigned j=0; j<total_numdf; ++j)
      {
        int onesetj = j % numdf;
        if ((*onoff)[onesetj]==0)
        {
          const int lid = (*systemvectoraux).Map().LID(dofs[j]);
          if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
          if (toggle!=Teuchos::null)
            (*toggle)[lid] = 0.0;
          // get rid of entry in DBC map - if it exists
          if (dbcgids != Teuchos::null)
            (*dbcgids).erase(dofs[j]);;
          continue;
        }
        const int gid = dofs[j];

        // assign value
        const int lid = (*systemvectoraux).Map().LID(gid);
        // set toggle vector
        if (toggle != Teuchos::null)
          (*toggle)[lid] = 1.0;
        // amend vector of DOF-IDs which are Dirichlet BCs
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }  // loop over nodal DOFs

    if(findcolgids)
      for (unsigned j=0; j<total_numdf; ++j)
      {
        int onesetj = j % numdf;
        if ((*onoff)[onesetj]==0)
        {
          // get rid of entry in DBC map - if it exists
          (*dbccolgids).erase(dofs[j]);
          continue;
        }
        const int gid = dofs[j];

        // assign value
        (*dbccolgids).insert(gid);
      }  // loop over nodal DOFs
  }  // loop over nodes

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::NURBS::NurbsDiscretization::FillMatrixAndRHSForLSDirichletBoundary(
    Teuchos::RCP<DRT::Element>              actele,
    const std::vector<Epetra_SerialDenseVector>* knots,
    const  std::vector<int> &               lm,
    const std::vector<int>*                 funct,
    const std::vector<int>*                 curve,
    const std::vector<double>*              val,
    const bool                              usetime,
    const unsigned                          deg,
    const double                            time,
    Epetra_SerialDenseMatrix&               elemass,
    std::vector<Epetra_SerialDenseVector>&  elerhs)
{
  if(deg+1!=elerhs.size())
    dserror("given degree of time derivative does not match number or rhs vectors!");

  static const int dim= DRT::UTILS::DisTypeToDim<distype>::dim;

  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs/nen;

  // get node coordinates of element
  LINALG::Matrix<dim+1,nen> xyze;
  DRT::Node** nodes = actele->Nodes();

  for(int inode=0;inode<nen;inode++)
  {
    const double* x = nodes[inode]->X();
    for(int idim=0;idim<dim+1;++idim)
    {
      xyze(idim,inode)=x[idim];
    }
  }

  // aquire weights from nodes
  Epetra_SerialDenseVector weights(nen);

  for (int inode=0; inode<nen; ++inode)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

    weights(inode) = cp->W();
  }

  //shape functions
  LINALG::Matrix<nen,1> shpfunct;
  //coordinates of integration points in parameter space
  LINALG::Matrix<dim,1> xsi;
  //first derivative of shape functions
  LINALG::Matrix<dim,nen>  deriv;
  //coordinates of integration point in physical space
  Epetra_SerialDenseVector  position(3); // always three-dimensional coordinates for function evaluation!
  //auxiliary date container for dirichlet evaluation
  std::vector<Epetra_SerialDenseVector> value(deg+1,dofblock);
  //unit normal on boundary element
  LINALG::Matrix<dim+1,1> unitnormal;

  // gaussian points
  const DRT::UTILS::IntPointsAndWeights<dim> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {
    //integration factor
    double fac=0.0;
    double drs =0.0;

    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(
        shpfunct,
        deriv,
        fac,
        unitnormal,
        drs,
        xsi,
        xyze,
        intpoints,
        iquad,
        knots,
        &weights,
        true);

    // get real physical coordinates of integration point
    /*
    //              +-----
    //               \
    //    pos (x) =   +      N (x) * x
    //               /        j       j
    //              +-----
    //              node j
    */
    for (int rr=0;rr<dim+1;++rr)
    {
      position(rr)=shpfunct(0)*xyze(rr,0);
      for (int mm=1;mm<nen;++mm)
      {
        position(rr)+=shpfunct(mm)*xyze(rr,mm);
      }
    }
    // if dim < 3, ensure we define a valid z-coordinate!
    for (int rr=dim+1;rr<3;++rr)
      position(rr)=0.0;

    for(int rr=0;rr<dofblock;++rr)
    {
      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[rr];
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      double functfac = 1.0;
      const int funct_num = (*funct)[rr];
      if (funct_num>0)
        // important: position has to have always three components!!
        functfac =DRT::Problem::Instance()->Funct((*funct)[rr]-1).Evaluate(rr,position.Values(),0.0,NULL);

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i](rr) = (*val)[rr]*functfac * curvefac[i];
      }
    }

    for (int vi=0; vi<nen; ++vi) // loop rows  (test functions)
    {
      const int fvi=dofblock*vi;

      for (int ui=0; ui<nen; ++ui) // loop columns  (test functions)
      {
        const int fui=dofblock*ui;

        const double diag=fac*shpfunct(ui)*shpfunct(vi);

        for(int rr=0;rr<dofblock;++rr)
        {
          elemass(fvi+rr,fui+rr)+=diag;
        }
      }
      for(int rr=0;rr<dofblock;++rr)
      {
        for (unsigned i=0; i<deg+1; ++i)
          elerhs[i](fvi+rr)+=fac*shpfunct(vi)*value[i](rr);
      }
    }
  } // end gaussloop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::NURBS::NurbsDiscretization::FillMatrixAndRHSForLSDirichletDomain(
    Teuchos::RCP<DRT::Element>              actele,
    const std::vector<Epetra_SerialDenseVector>* knots,
    const  std::vector<int> &               lm,
    const std::vector<int>*                 funct,
    const std::vector<int>*                 curve,
    const std::vector<double>*              val,
    const bool                              usetime,
    const unsigned                          deg,
    const double                            time,
    Epetra_SerialDenseMatrix&               elemass,
    std::vector<Epetra_SerialDenseVector>&  elerhs)
{
  if(deg+1!=elerhs.size())
    dserror("given degree of time derivative does not match number or rhs vectors!");

  static const int dim= DRT::UTILS::DisTypeToDim<distype>::dim;
  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs/nen;

  // get node coordinates of element
  LINALG::Matrix<dim,nen> xyze;
  DRT::Node** nodes = actele->Nodes();

  for(int inode=0;inode<nen;inode++)
  {
    const double* x = nodes[inode]->X();
    for(int idim=0;idim<dim;++idim)
    {
      xyze(idim,inode)=x[idim];
    }
  }

  // aquire weights from nodes
  Epetra_SerialDenseVector weights(nen);

  for (int inode=0; inode<nen; ++inode)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

    weights(inode) = cp->W();
  }

  //shape functions
  LINALG::Matrix<nen,1> shpfunct;
  //coordinates of integration points in parameter space
  LINALG::Matrix<dim,1> xsi;
  // transposed jacobian "dx/ds"
  LINALG::Matrix<dim,dim> xjm;
  //first derivative of shape functions
  LINALG::Matrix<dim,nen>  deriv;
  //coordinates of integration point in physical space
  Epetra_SerialDenseVector  position(3); // always three-dimensional coordinates for function evaluation!
  //auxiliary date container for dirichlet evaluation
  std::vector<Epetra_SerialDenseVector> value(deg+1,dofblock);

  // gaussian points
  const DRT::UTILS::IntPointsAndWeights<dim> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {

    for (int idim=0;idim<dim;idim++)
    {
       xsi(idim) = intpoints.IP().Point(iquad)[idim];
    }

    //evaluate shape function and derivatevs at integration point
    DRT::NURBS::UTILS::nurbs_get_funct_deriv
    (shpfunct  ,
      deriv  ,
      xsi    ,
      *knots,
      weights,
      distype );

    xjm.MultiplyNT(deriv,xyze);
    double det = xjm.Determinant();

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

    // compute integration factor
    double fac = intpoints.IP().qwgt[iquad]*det;

    // get real physical coordinates of integration point
    /*
    //              +-----
    //               \
    //    pos (x) =   +      N (x) * x
    //               /        j       j
    //              +-----
    //              node j
    */
    for (int rr=0;rr<dim;++rr)
    {
      position(rr)=shpfunct(0)*xyze(rr,0);
      for (int mm=1;mm<nen;++mm)
      {
        position(rr)+=shpfunct(mm)*xyze(rr,mm);
      }
    }
    // if dim < 3, ensure we define a valid z-coordinate!
    for (int rr=dim;rr<3;++rr)
      position(rr)=0.0;

    for(int rr=0;rr<dofblock;++rr)
    {
      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[rr];
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      double functfac = 1.0;
      const int funct_num = (*funct)[rr];
      if (funct_num>0)
        // important: position has to have always three components!!
        functfac =DRT::Problem::Instance()->Funct((*funct)[rr]-1).Evaluate(rr,position.Values(),0.0,NULL);

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i](rr) = (*val)[rr]*functfac * curvefac[i];
      }
    }

    for (int vi=0; vi<nen; ++vi) // loop rows  (test functions)
    {
      const int fvi=dofblock*vi;

      for (int ui=0; ui<nen; ++ui) // loop columns  (test functions)
      {
        const int fui=dofblock*ui;

        const double diag=fac*shpfunct(ui)*shpfunct(vi);

        for(int rr=0;rr<dofblock;++rr)
        {
          elemass(fvi+rr,fui+rr)+=diag;
        }
      }
      for(int rr=0;rr<dofblock;++rr)
      {
        for (unsigned i=0; i<deg+1; ++i)
          elerhs[i](fvi+rr)+=fac*shpfunct(vi)*value[i](rr);
      }
    }
  } // end gaussloop

  return;
}
