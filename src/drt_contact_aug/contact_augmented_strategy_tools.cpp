/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_strategy_tools.cpp

\brief Tools for the augmented contact solving strategy with
       standard Lagrangian multipliers.

\level 3

\maintainer Michael Hiermeier

\date Jun 3, 2014

*/
/*---------------------------------------------------------------------*/

#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_node.H"

#include "../linalg/linalg_utils.H"

#include "../drt_fsi/fsi_matrixtransform.H"

//#define CONTACTFD_DLMGAPLINMATRIX      /* flag for global FD-check of the weighted gap gradient
// w.r.t. displ. */ #define CONTACTFD_DGLMLINMATRIX        /* flag for global FD-check of the
// dGLmLinMatrix w.r.t. displ. */ #define CONTACTFD_DGGLINMATRIX         /* flag for global FD-check
// of the dGGLinMatrix w.r.t. displ. */

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AugFDCheckGP(CONTACT::ParamsInterface& cparams)
{
  if (not IsInContact()) return;

  for (auto& interface_ptr : interface_)
  {
    Teuchos::RCP<Interface> aug_interface_ptr =
        Teuchos::rcp_dynamic_cast<Interface>(interface_ptr, true);

    aug_interface_ptr->FD_Debug_GPCheck1stOrder(cparams);
    aug_interface_ptr->FD_Debug_GPCheck2ndOrder(cparams);

    aug_interface_ptr->FD_DebugVec_GPCheck1stOrder(cparams);
    aug_interface_ptr->FD_DebugVec_GPCheck2ndOrder(cparams);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AugFDCheckGlobal(CONTACT::ParamsInterface& cparams)
{
  // *** linearization w.r.t. displ. *********************************
#ifdef CONTACTFD_DLMGAPLINMATRIX
  static bool first_attempt = true;
  if (first_attempt)
  {
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << "| Finite Difference check of [g_{N}]_{i,j}              |\n";
    std::cout << "*-------------------------------------------------------*\n";
    first_attempt = false;
  }

  FD_Debug::Instance(this)->Evaluate(Data().DLmNWGapLinMatrixPtr(), Data().WGapPtr(), cparams);
#endif

  /*---------------------------------*
   | Std force terms                 |
   *---------------------------------*/
#ifdef CONTACTFD_DGLMLINMATRIX
  static Teuchos::RCP<Epetra_Vector> dGLm = Teuchos::null;
  dGLm = Teuchos::rcp(new Epetra_Vector(*Data().GSlMaDofRowMapPtr(), true));

  LINALG::AssembleMyVector(0.0, *dGLm, 1.0, Data().SlForceLm());
  LINALG::AssembleMyVector(1.0, *dGLm, 1.0, Data().MaForceLm());

  double nrm2 = 0.0;
  dGLm->Norm2(&nrm2);
  if (nrm2 == 0.0) return;

  static bool first_attempt = true;
  if (first_attempt)
  {
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << "| Finite Difference check of [g_{N}]_{i,jk} [l_{N}]^{i} |\n";
    std::cout << "*-------------------------------------------------------*\n";
    first_attempt = false;
  }

  FD_Debug::Instance(this)->Evaluate(Data().DGLmLinMatrixPtr(), dGLm, cparams);
#endif
  /*---------------------------------*
   | Regularization term             |
   *---------------------------------*/
#ifdef CONTACTFD_DGGLINMATRIX
  static Teuchos::RCP<Epetra_Vector> dGG = Teuchos::null;
  dGG = Teuchos::rcp(new Epetra_Vector(*Data().GSlMaDofRowMapPtr(), true));

  LINALG::AssembleMyVector(0.0, *dGG, 1.0, Data().SlForceG());
  LINALG::AssembleMyVector(1.0, *dGG, 1.0, Data().MaForceG());

  static bool first_attempt = true;
  if (first_attempt)
  {
    std::cout << "*------------------------------------------------------------*\n";
    std::cout << "| Finite Difference check of ([g_{N}]_{i} [g_{N}]^{i})_{,jk} |\n";
    std::cout << "*------------------------------------------------------------*\n";
    first_attempt = false;
  }
  FD_Debug::Instance(this)->Evaluate(Data().DGGLinMatrixPtr(), dGG, cparams);
#endif

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Strategy::FD_Debug* CONTACT::AUG::Strategy::FD_Debug::Instance(
    Strategy* strat, const double delta, const bool delete_me)
{
  static FD_Debug* instance = NULL;

  if (delete_me)
  {
    if (instance) delete instance;

    return NULL;
  }

  if (not instance)
  {
    instance = new FD_Debug();
    instance->Init(strat, delta);
  }

  return instance;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::FD_Debug::Done() { Instance(NULL, 0.0, true); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::FD_Debug::Init(Strategy* strat, const double delta)
{
  if (not strat) dserror("NULL pointer!");

  strat_ = strat;
  delta_ = delta;
  is_fd_check_ = false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::FD_Debug::Evaluate(
    const Teuchos::RCP<LINALG::SparseMatrix>& derivMatrixPtr,
    Teuchos::RCP<Epetra_Vector>& rhsVector, CONTACT::ParamsInterface& cparams)
{
  if (is_fd_check_)
    return;
  else
    is_fd_check_ = true;

  if (strat_->Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  LINALG::SparseMatrix derivMatrix(*derivMatrixPtr);

  const Epetra_Map rowMap = derivMatrix.RowMap();
  const Epetra_Map colMap = derivMatrix.ColMap();

  const Epetra_BlockMap rhsMap = rhsVector->Map();

  LINALG::SparseMatrix fdMatrixRef = LINALG::SparseMatrix(rowMap, 100);
  LINALG::SparseMatrix fdMatrixNew = LINALG::SparseMatrix(rowMap, 100);
  int dim = strat_->Dim();

  // create reference Matrix:
  // loop over all columns
  for (int c = 0; c < colMap.NumMyElements(); ++c)
  {
    int colId = colMap.GID(c);
    // loop over all rows of the right hand side vector
    for (int r = 0; r < rhsMap.NumMyElements(); ++r)
    {
      int rowId = rhsMap.GID(r);
      if (rowMap.LID(rowId) == -1)
        dserror(
            "ERROR: Row gids of the corresponding rhs-vector and the "
            "derivative matrix do not match!");

      double val = (*rhsVector)[r];

      fdMatrixRef.Assemble(val, rowId, colId);
    }
  }
  fdMatrixRef.Complete(colMap, rowMap);

  ref_x_.clear();
  // loop over all columns of the reference Matrix
  for (int fd = 0; fd < colMap.NumMyElements(); ++fd)
  {
    const int colId = colMap.GID(fd);
    const int gid = colId / dim;
    const int dof = colId % dim;
    const char xyz[3] = {'x', 'y', 'z'};

    std::cout << "Performing FD perturbation of column " << (fd + 1) << "/"
              << colMap.NumMyElements() << " corresponding to DOF " << colId << "(" << xyz[dof]
              << ") of node " << gid << "..." << std::flush;

    // do the finite difference step
    DoPerturbation(gid, dof);

    // Update matrix and rhs
    for (auto& iptr : strat_->Interfaces()) iptr->SetElementAreas();
    strat_->EvalForceStiff(cparams);

    // loop over all rows of the updated right hand side vector
    // and save the values in a new matrix
    for (int r = 0; r < rhsMap.NumMyElements(); ++r)
    {
      const int rowId = rhsMap.GID(r);
      if (rowId == -1) dserror("Couldn't find the GID for lid %d!", r);
      const double val = (*rhsVector)[r];
      fdMatrixNew.Assemble(val, rowId, colId);
    }

    // Undo finite difference step
    UndoPerturbation(gid, dof);

    // Update matrix and rhs
    for (auto& iptr : strat_->Interfaces()) iptr->SetElementAreas();
    strat_->EvalForceStiff(cparams);

    std::cout << "done!" << std::endl;
  }

  // calculate the finite difference
  const double delta_inv = 1.0 / delta_;
  fdMatrixNew.Add(fdMatrixRef, false, -delta_inv, delta_inv);
  fdMatrixNew.Complete(colMap, rowMap);

  // loop over all rows
  for (int r = 0; r < rowMap.NumMyElements(); ++r)
  {
    int rowId = rowMap.GID(r);
    // check if the row belongs to the slave or master side
    std::string rSlMa = "(S)";
    if (strat_->SlaveRowDofs()->LID(rowId) == -1) rSlMa = "(M)";

    int w = 0;

    // *** finite differences ***
    // get all non-zero values and the corresponding ids of the current row
    int rLengthFD = fdMatrixNew.EpetraMatrix()->NumGlobalEntries(rowId);
    int numEntriesFD = 0;
    std::vector<double> rValFD(rLengthFD);
    std::vector<int> cIdsFD(rLengthFD);
    fdMatrixNew.EpetraMatrix()->ExtractGlobalRowCopy(
        rowId, rLengthFD, numEntriesFD, &rValFD[0], &cIdsFD[0]);

    // *** analytical solution ***
    // get all non-zero values and the corresponding ids of the current row
    int rLengthAna = derivMatrix.EpetraMatrix()->NumGlobalEntries(rowId);
    int numEntriesAna = 0;
    std::vector<double> rValAna(rLengthAna);
    std::vector<int> cIdsAna(rLengthAna);
    derivMatrix.EpetraMatrix()->ExtractGlobalRowCopy(
        rowId, rLengthAna, numEntriesAna, &rValAna[0], &cIdsAna[0]);

    /*-------------------------------------------------------------*
     |   Compare analytical and finite difference solution         |
     *-------------------------------------------------------------*/
    std::cout << "\nFINITE DIFFERENCE CHECK FOR ROW-ID # " << rowId << rSlMa << std::endl;
    for (int c = 0; c < rLengthFD; ++c)
    {
      // Do the finite difference check only for values which are greater than some threshold
      if (std::abs(rValFD[c]) < 1e-12) continue;

      // check if the column belongs to the slave or master side
      std::string cSlMa = "(S)";
      if (strat_->SlaveRowDofs()->LID(cIdsFD[c]) == -1) cSlMa = "(M)";

      // search for entry in the analytical solution
      int anaId = -1;
      for (int cc = 0; cc < rLengthAna; ++cc)
        if (cIdsFD[c] == cIdsAna[cc])
        {
          anaId = cc;
          break;
        }

      if (std::abs(rValFD[c]) < delta_ and anaId == -1)
        continue;
      else if (anaId == -1)
        std::cout << "*** WARNING: Global column #" << cIdsFD[c] << " in global row #" << rowId
                  << " could not be found in the analytical derivative matrix! (fd= " << rValFD[c]
                  << ") ***" << std::endl;
      else
      {
        double dev = rValFD[c] - rValAna[anaId];

        std::cout << cIdsFD[c] << cSlMa
                  << ":"
                     "   fd="
                  << std::setw(14) << std::setprecision(5) << std::scientific << rValFD[c]
                  << "   ana=" << std::setw(14) << std::setprecision(5) << std::scientific
                  << rValAna[anaId] << "   DEVIATION=" << std::setw(14) << std::setprecision(5)
                  << std::scientific << dev << "   REL-ERROR [%]=" << std::setw(14)
                  << std::setprecision(5) << std::scientific << abs(dev / rValFD[c]) * 100;

        if (abs(dev) > 1.0e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1.0e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }
    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // deactivate global finite difference indicator
  is_fd_check_ = false;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::FD_Debug::DoPerturbation(const int gid, const int dof)
{
  DRT::Node* node = FindINode(gid);

  CoNode* cnode = dynamic_cast<CoNode*>(node);

  // store current position
  LINALG::Matrix<3, 1>& x = ref_x_[gid];
  std::copy(cnode->xspatial(), cnode->xspatial() + 3, x.A());

  // change forward step to backward step
  switch (dof)
  {
    case 0:
      cnode->xspatial()[0] += delta_;
      break;
    case 1:
      cnode->xspatial()[1] += delta_;
      break;
    case 2:
      cnode->xspatial()[2] += delta_;
      break;
    default:
      dserror("Unsupported case!");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::FD_Debug::UndoPerturbation(const int gid, const int dof) const
{
  DRT::Node* node = FindINode(gid);

  CoNode* cnode = dynamic_cast<CoNode*>(node);

  // get stored position
  const LINALG::Matrix<3, 1>& x = ref_x_.at(gid);
  std::copy(x.A(), x.A() + 3, cnode->xspatial());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Node* CONTACT::AUG::Strategy::FD_Debug::FindINode(const int gid) const
{
  DRT::Node* node = NULL;

  // do the finite difference step
  for (auto& interface : strat_->Interfaces())
  {
    node = interface->Discret().gNode(gid);
    if (node) return node;
  }
  dserror("Node %d not found!", gid);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::ExtractMatrix(const LINALG::SparseMatrix& source,
    const Epetra_Map& target_range_map, const Epetra_Map& target_domain_map)
{
  if (not source.Filled()) dserror("The source matrix must be filled!");

  const int maxnumentries = source.MaxNumEntries();
  Teuchos::RCP<LINALG::SparseMatrix> target_ptr =
      Teuchos::rcp(new LINALG::SparseMatrix(target_range_map, maxnumentries));
  LINALG::SparseMatrix& target = *target_ptr;

  if (target_range_map.NumGlobalElements())
  {
    FSI::UTILS::MatrixLogicalSplitAndTransform extractor;
    extractor(source, target_range_map, target_domain_map, 1.0, NULL, NULL, target);
  }

  target.Complete(target_domain_map, target_range_map);

  return target_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::MultiplyElementwise(const Epetra_Vector& source,
    const Epetra_Map& source2targetMap, Epetra_Vector& target, const bool inverse)
{
  // consistency check
  if (source2targetMap.NumMyElements() != target.Map().NumMyElements())
  {
    const Epetra_Comm& comm = source.Comm();
    dserror(
        "The number of local elements of the source2targetMap and the "
        "target.Map() have to be equal! \n"
        ".........................................................\n"
        "source2targetMap.NumMyElements() = %d on proc %i \n"
        "target.Map().NumMyElements() = %d on proc %i \n"
        ".........................................................\n",
        source2targetMap.NumMyElements(), comm.MyPID(), target.Map().NumMyElements(), comm.MyPID());
  }

  // nothing to do, if the target map size is equal zero
  if (source2targetMap.NumGlobalElements() == 0) return;

  Teuchos::RCP<Epetra_Vector> source_exp_ptr = LINALG::ExtractMyVector(source, source2targetMap);

  Epetra_Vector& source_exp = *source_exp_ptr;
  int error = source_exp.ReplaceMap(target.Map());

  if (error) dserror("The source map couldn't be replaced by the target map! (error=%d)", error);

  if (inverse)
    error = target.ReciprocalMultiply(1.0, source_exp, target, 0.0);
  else
    error = target.Multiply(1.0, source_exp, target, 0.0);

  if (error) dserror("The element-wise multiplication failed! (error=%d)", error);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::RedistributeRowMap(const Epetra_Map& ref_map, Epetra_Map& red_map)
{
  const Epetra_Comm& comm = ref_map.Comm();

  const int nummyelements = ref_map.NumMyElements();
  int* mygids = ref_map.MyGlobalElements();

  int count = 0;
  std::vector<int> myGids(nummyelements, -1);

  const Teuchos::RCP<Epetra_Map> allreducedMap = LINALG::AllreduceEMap(red_map);

  for (int i = 0; i < nummyelements; ++i)
  {
    const int gid = mygids[i];
    if (allreducedMap->LID(gid) >= 0)
    {
      myGids[count] = gid;
      ++count;
    }
  }

  myGids.resize(count);
  int gCount = 0;
  comm.SumAll(&count, &gCount, 1);
  red_map = Epetra_Map(gCount, count, &myGids[0], 0, comm);
}
