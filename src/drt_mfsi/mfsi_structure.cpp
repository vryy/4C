#ifdef CCADISCRET

#include "mfsi_structure.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::Structure::~Structure()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::StructureAdapter::StructureAdapter(Teuchos::RCP<Teuchos::ParameterList> params,
                                         Teuchos::RCP<DRT::Discretization> dis,
                                         Teuchos::RCP<LINALG::Solver> solver,
                                         Teuchos::RCP<IO::DiscretizationWriter> output)
  : structure_(*params, *dis, *solver, *output),
    interface_(dis),
    dis_(dis),
    params_(params),
    solver_(solver),
    output_(output)
{
  interface_.SetupCondDofMap("FSICoupling");
  interface_.SetupOtherDofMap();
  sumdisi_ = Teuchos::rcp(new Epetra_Vector(Dispm()->Map()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> MFSI::StructureAdapter::InitialGuess()
{
  return Teuchos::rcp(&structure_.Getdu(),false);
  //return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::StructureAdapter::RHS()
{
  return structure_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::StructureAdapter::Dispm()
{
  return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> MFSI::StructureAdapter::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> MFSI::StructureAdapter::SysMat()
{
  return structure_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> MFSI::StructureAdapter::Discretization()
{
  return structure_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MFSI::StructureAdapter::DispIncrFactor()
{
  return structure_.DispIncrFactor();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::StructureAdapter::FluidCondRHS()
{
  // structure part of the rhs to enforce
  // u(n+1) dt = d(n+1) - d(n)

  // extrapolate d(n+1) at the interface and substract d(n)

  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp ());

  double alphaf = structure_.AlphaF();
  idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf)-1.);
  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::StructureAdapter::MeshCondRHS()
{
  // structure part of the rhs to enforce
  // d(G,n+1) = d(n+1)

  // extrapolate d(n+1) at the interface

  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp ());

  double alphaf = structure_.AlphaF();
  idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf));
  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> MFSI::StructureAdapter::InterfaceDisp()
// {
//   Teuchos::RefCountPtr<Epetra_Vector> idis  = rcp(new Epetra_Vector(*idispmap_));

//   int err = idis->Import(*structure_.Disp(),*extractor_,Insert);
//   if (err)
//     dserror("Import using importer returned err=%d",err);

//   return idis;
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::StructureAdapter::PrepareTimeStep()
{
  sumdisi_->PutScalar(0.);

  std::string pred = params_->get<string>("predictor","consistent");
  if (pred=="constant")
  {
    structure_.ConstantPredictor();
  }
  else if (pred=="consistent")
  {
    structure_.ConsistentPredictor();
  }
  else
    dserror("predictor %s unknown", pred.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::StructureAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (disp!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    disi->Update(-1.0,*sumdisi_,1.0);
    structure_.Evaluate(disi);
    sumdisi_->Update(1.0,*disp,0.0);
  }
  else
  {
    structure_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::StructureAdapter::Update()
{
  structure_.UpdateandOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::StructureAdapter::Output()
{
  // noop
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::StructureAdapter::SetInterfaceMap(Teuchos::RCP<Epetra_Map> im)
{
#if 0
  // build inner displacement map
  // dofs at the interface are excluded
  // we use only dofs without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map> dispmap = DofRowMap();
  Teuchos::RCP<Epetra_Vector> dirichtoggle = structure_.Dirichlet();

  int numids = dispmap->NumMyElements();
  std::vector<int> ids;
  ids.reserve(numids);
  for (int i=0; i<numids; ++i)
  {
    int gid = dispmap->GID(i);
    if (not interface_.CondDofMap()->MyGID(gid) and (*dirichtoggle)[i]==0.)
    {
      ids.push_back(gid);
    }
  }

  innerdispmap_ = Teuchos::rcp(new Epetra_Map(-1,ids.size(), &ids[0], 0, dispmap->Comm()));
#endif
}


#endif
