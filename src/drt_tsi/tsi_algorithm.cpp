/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.cpp

\brief  Basis of all TSI algorithms that perform a coupling between the linear
        momentum equation and the heat conduction equation

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_contact/contact_abstract_strategy.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(
  Epetra_Comm& comm
  )
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    tempincnp_(rcp(new Epetra_Vector(*(ThermoField().Tempnp()))))
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();
  // what kind of one-way coupling should be applied
  displacementcoupling_
    = tsidyn.get<std::string>("COUPVARIABLE") == "Displacement";
    
  if (displacementcoupling_) // (temperature change due to deformation)
  {
    // build a proxy of the structure discretization for the temperature field
    Teuchos::RCP<DRT::DofSet> structdofset
      = StructureField().Discretization()->GetDofSetProxy();
    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in thermo field");

    // build matrices again and now consider displacements in the thermo field
    ThermoField().TSIMatrix();
  }
  
  else // temperature as coupling variable (deformation due to temperature change)
  {
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSet> thermodofset
      = ThermoField().Discretization()->GetDofSetProxy();
    // check if StructField has 2 discretizations, so that coupling is possible
    if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
      dserror("unexpected dof sets in structure field");
    // now build the matrices again and consider the temperature dependence
    StructureField().TSIMatrix();
  }

  //  // now check if the two dofmaps are available and then bye bye
  //  cout << "structure dofmap" << endl;
  //  cout << *StructureField().DofRowMap(0) << endl;
  //  cout << "thermo dofmap" << endl;
  //  cout << *StructureField().DofRowMap(1) << endl;
  //  cout << "thermo dofmap" << endl;
  //  cout << *ThermoField().DofRowMap(0) << endl;
  //  cout << "structure dofmap" << endl;
  //  cout << *ThermoField().DofRowMap(1) << endl;
  //  exit(0);

}


/*----------------------------------------------------------------------*
 | destructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::ReadRestart(int step)
{
  ThermoField().ReadRestart(step);
  StructureField().ReadRestart(step);
  SetTimeStep(ThermoField().GetTime(),step);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TSI::Algorithm::TimeLoop()
{
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = tsidyn.get<int>("ITEMAX"); // default: =1
  ittol_ = tsidyn.get<double>("CONVTOL"); // default: =1e-6

  // get an idea of the temperatures (like in partitioned FSI)
  itempn_ = ThermoField().ExtractTempnp();
  idispn_ = StructureField().ExtractDispn();


  // ==================================================================

  // time loop
  while (NotFinished())
  {
    // get active nodes from structural contact simulation
    RCP<MORTAR::ManagerBase> cmtman = StructureField().ContactManager();

    // tsi with or without contact
    // only tsi
    if (cmtman == Teuchos::null)
    {
      // counter and print header
      IncrementTimeAndStep();
      PrintHeader();

      // what kind of coupling should be applied
      displacementcoupling_
        = tsidyn.get<std::string>("COUPVARIABLE") == "Displacement";

      if (displacementcoupling_) // (temperature change due to deformation)
      {
        // predict the thermal field without influence of structure
        StructureField().PrepareTimeStep();

        // Begin Nonlinear Solver / Outer Iteration ******************************

        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;

        // Outer Iteration loop starts
        if (Comm().MyPID()==0)
        {
          cout<<"\n";
          cout<<"**************************************************************\n";
          cout<<"      OUTER ITERATION LOOP \n";
          printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
            ThermoField().GetTimeNumStep());
          cout<<"**************************************************************\n";
        }

        // Start OUTER ITERATION
        while (stopnonliniter == false)
        {
          itnum ++;

          // store temperature from first solution for convergence check
          tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

          // get current displacement due to Solve Structure Step
          const Teuchos::RCP<Epetra_Vector> idisp = DoStructureStep();

          // solve thermo system
          // and therefore use the u-solution calculated in DoStructureStep
          // including: - ApplyDisplacements()
          //            - PrepareTimeStep()
          //            - Solve()
          DoThermoStep(idisp);

          // check convergence of temperature field for "partitioned scheme"
          stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        } // end OUTER ITERATION

        // End Nonlinear Solver **************************************

        // update all single field solvers
        Update();

        // extract final temperatures,
        // since we did update, this is very easy to extract
        idispn_ = StructureField().ExtractDispnp();

      } // displacement coupling

      else // temperature coupling (deformation due to heating)
      {
        // predict the thermal field without influence of structure
        ThermoField().PrepareTimeStep();

        // Begin Nonlinear Solver / Outer Iteration ******************************

        // iterate between the two fields
        int  itnum = 0;
        bool stopnonliniter = false;

        // Outer Iteration loop starts
        if (Comm().MyPID()==0)
        {
          cout<<"\n";
          cout<<"**************************************************************\n";
          cout<<"      OUTER ITERATION LOOP \n";
          printf("      Time Step %3d/%3d \n",ThermoField().GetTimeStep(),
            ThermoField().GetTimeNumStep());
          cout<<"**************************************************************\n";
        }

        // Start OUTER ITERATION
        while (stopnonliniter == false)
        {
          itnum ++;

          // store temperature from first solution for convergence check
          tempincnp_->Update(1.0,*ThermoField().Tempnp(),0.0);

          // get current temperatures due to Solve ThermoStep
          const Teuchos::RCP<Epetra_Vector> itemp = DoThermoStep();

          // solve structure system
          // and therefore use the temperature solution calculated in DoThermoStep
          // including: - ApplyTemperatures()
          //            - PrepareTimeStep()
          //            - Solve()
          DoStructureStep(itemp);

          // check convergence of temperature field for "partitioned scheme"
          stopnonliniter = ConvergenceCheck(itnum,itmax_,ittol_);

        } // end OUTER ITERATION

        // End Nonlinear Solver **************************************

        // update all single field solvers
        Update();

        // extract final temperatures,
        // since we did update, this is very easy to extract
        itempn_ = ThermoField().ExtractTempnp();

      }  // temperature coupling

      // write output to screen and files
      Output();
    } // tsi

    // thermo-structure interaction with contact in seperate routine
    else
    {
      //****************************************************************//
      // this algorithm consists of two parts within one time step      //
      // 1) solution of nonlinear structural system without influence   //
      //    of temperature.                                             //
      // 2) solution of linear temperature problem                      //
      //    with heat transfer over the contacting surface (as a result //
      //    from the structural problem). This linear problem is solved //
      //    within one Newton step.                                     //
      //******************************************************************

      // only for frictional contact so far
      // as information are needed from the friction node
      if((cmtman->GetStrategy().Friction())==false)
        dserror ("Thermo-Structure interaction only for frictional contact so far");

      // reset thermo field
      // FIXGIT: Is there another possibility to reset
      // I had to omit the "const"
      (ThermoField().Tempn())->PutScalar(0.0);

      // counter and print header
      IncrementTimeAndStep();
      PrintHeader();

      // predict ans solve structural system
      StructureField().PrepareTimeStep();
      StructureField().Solve();

      // predict and evaluate the thermal field without influence of structure
      // evaluate assembles the stiffness matrix and the rhs of
      // the thermo field
      ThermoField().PrepareTimeStep();
      ThermoField().Evaluate(Teuchos::null);

      // get the linear system of equations form the thermo field
      Teuchos::RCP<const LINALG::SparseMatrix> stiff = ThermoField().SystemMatrix();
      Teuchos::RCP<const Epetra_Vector> rhs = ThermoField().RHS();

      // copy the linear system of equations to new matrix/vector
      RCP<LINALG::SparseMatrix> stiffthermcont = rcp(new LINALG::SparseMatrix(*stiff));
      RCP<Epetra_Vector> rhsthermcont = rcp(new Epetra_Vector(*rhs));

      // make modifications towards thermo contact to the stiffness matrix
      // and the right hand side
      ApplyThermoContact(stiffthermcont,rhsthermcont,cmtman);

      // solution vector
      // FIXGIT: perhaps this should be written in existing vector
      Teuchos::RCP<Epetra_Vector> inc = LINALG::CreateVector(*(ThermoField().DofRowMap()),true);
      inc = LINALG::CreateVector(*(ThermoField().DofRowMap()),true);

      // Create the solver and solve modified system of equations
      Teuchos::RCP<LINALG::Solver> solver = ThermoField().LinearSolver();
      solver->Solve(stiffthermcont->EpetraMatrix(),inc,rhsthermcont,true,true);

      // write solution
      ThermoField().UpdateNewton(inc);

      // update all single field solvers
      Update();

      // write output to screen and files
      Output();
    }
  } // time loop

  // ==================================================================

}

/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 03/10 |
 | apply temperatures, coupling variable in structure field             |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoStructureStep(Teuchos::RCP<Epetra_Vector> itemp)
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }
  // call the current temperatures
  StructureField().ApplyTemperatures(itemp);

  // call the predictor here, because the temperature is considered too
  StructureField().PrepareTimeStep();

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  return;
}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 05/10 |
 | extract end temperatures for coupling to displacement field          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::DoThermoStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }
  // call the predictor here, because the displacements influences the thermo
  // solution (2nd discretization is available here)
  ThermoField().PrepareTimeStep();
  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();
  // now extract the current temperatures and pass it to the structure
  return ThermoField().ExtractTempnp();
}


/*----------------------------------------------------------------------*
 | Solve the structure system (protected)                    dano 05/10 |
 | extract current displacements needed for coupling to thermo field    |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>  TSI::Algorithm::DoStructureStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"************************\n";
    cout<<"    STRUCTURE SOLVER    \n";
    cout<<"************************\n";
  }

  // solve structure system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  StructureField().Solve();
  // now extract the current temperatures and pass it to the structure
  return StructureField().ExtractDispnp();
}


/*----------------------------------------------------------------------*
 | Solve the thermo system (protected)                       dano 05/10 |
 | coupling of displacements to thermo field                            |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::DoThermoStep(Teuchos::RCP<Epetra_Vector> idisp)
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n";
    cout<<"*********************\n";
    cout<<"    THERMO SOLVER    \n";
    cout<<"*********************\n";
  }
  // the displacement -> velocity conversion
  const Teuchos::RCP<Epetra_Vector> ivel = CalcVelocity(idisp);

  // call the current displacements and velocities
  ThermoField().ApplyStructVariables(idisp,ivel);

  // call the predictor here, because displacement field is considered too
  ThermoField().PrepareTimeStep();

  /// solve temperature system
  /// Do the nonlinear solve for the time step. All boundary conditions have
  /// been set.
  ThermoField().Solve();
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Update()
{
  StructureField().Update();
  ThermoField().Update();
  return;
}


/*----------------------------------------------------------------------*
 | (protected)                                               dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField().Output();

  // write the thermo output (temperatures at the moment) to the the structure output
  // get disc writer from structure field
  Teuchos::RCP<IO::DiscretizationWriter> output = StructureField().DiscWriter();

  // get the temperature and the noderowmap of thermo discretization
  Epetra_Vector temperature = *(ThermoField().Tempn());
  const Epetra_Map* temprowmap = ThermoField().Discretization()->NodeRowMap();

  // replace map and write it to output
  temperature.ReplaceMap(*temprowmap);
  RCP<Epetra_Vector> temp = rcp(new Epetra_Vector(temperature));
  output->WriteVector("temperature",temp);

  ThermoField().Output();
}


/*----------------------------------------------------------------------*
 | convergence check for the temperature field               dano 02/10 |
 | originally by vg 01/09                                               |
 *----------------------------------------------------------------------*/
bool TSI::Algorithm::ConvergenceCheck(
  int itnum,
  const int itmax,
  const double ittol
  )
{
  // convergence check based on the temperature increment
  bool stopnonliniter = false;

  //    | temperature increment |_2
  //  -------------------------------- < Tolerance
  //     | temperature_n+1 |_2

  // Variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double tempincnorm_L2(0.0); // pot
  double tempnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  tempincnp_->Update(1.0,*(ThermoField().Tempnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  ThermoField().Tempnp()->Norm2(&tempnorm_L2);

//  cout << "\ntempincnorm_L2\n" << tempincnorm_L2 << endl;
//  cout << "\ntempnorm_L2\n" << tempnorm_L2 << endl;

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  // if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    cout<<"\n";
    cout<<"****************************\n";
    cout<<"    OUTER ITERATION STEP    \n";
    cout<<"****************************\n";
    printf("+------------+-------------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|- temp-inc   -|\n");
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |",
         itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2);
    printf("\n");
    printf("+------------+-------------------+--------------+\n");
  }

  // Converged
  if (tempincnorm_L2/tempnorm_L2 <= ittol)
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d ! |\n", itnum,itmax);
      printf("+-----------------------------------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax) and (tempincnorm_L2/tempnorm_L2 > ittol))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf("|     >>>>>> not converged in itemax steps!     |\n");
      printf("+-----------------------------------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Algorithm::ConvergenceCheck


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 05/10 |
 | like InterfaceVelocity(idisp) in FSI::DirichletNeumann               |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Algorithm::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> idispnp
  ) const
{
  Teuchos::RCP<Epetra_Vector> ivel = Teuchos::null;
  // copy D_n onto V_n+1
  ivel = rcp(new Epetra_Vector(*idispn_));
  cout << "CalcVel kopier idispn_" << *ivel << endl;
  // calculate velocity with timestep Dt()
  //  V_n+1 = (D_n+1 - D_n) / Dt
  ivel->Update(1./Dt(), *idispnp, -1./Dt());

  return ivel;
}


/*----------------------------------------------------------------------*
 | apply modifications towards thermal contact               mgit 04/10 |
 *----------------------------------------------------------------------*/

void TSI::Algorithm::ApplyThermoContact(RCP<LINALG::SparseMatrix>& kteff,
                                        RCP<Epetra_Vector>& feff,
                                        RCP<MORTAR::ManagerBase> cmtman)
{

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();

  // convert maps (from structure discretization to thermo discretization)
  // slave-, active-, inactive-, master-, activemaster- and ndofs
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs;
  ConvertMaps (sdofs,adofs,mdofs,cmtman);

  // map of active and master dofs
  amdofs = LINALG::MergeMap(adofs,mdofs,false);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);

  // row map of thermal problem
  RCP<const Epetra_Map> problemrowmap = ThermoField().DofRowMap();

  // split problemrowmap in n+am
  ndofs = LINALG::SplitMap(*problemrowmap,*amdofs);

  // modifications only for active nodes
  if (adofs->NumGlobalElements()==0)
    return;

  // assemble Mortar Matrices D and M in thermo dofs for active nodes
  RCP<LINALG::SparseMatrix> dmatrix = rcp(new LINALG::SparseMatrix(*adofs,10));
  RCP<LINALG::SparseMatrix> mmatrix = rcp(new LINALG::SparseMatrix(*adofs,100));
  AssembleDM(*dmatrix,*mmatrix,cmtman);

  // FillComplete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*mdofs,*adofs);

  // matrices from linearized thermal contact condition
  RCP<LINALG::SparseMatrix>  thermcontLM = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<LINALG::SparseMatrix>  thermcontTEMP = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<Epetra_Vector>         thermcontRHS = LINALG::CreateVector(*adofs,true);

  // assemble thermal contact contition
  AssembleThermContCondition(*thermcontLM,*thermcontTEMP,*thermcontRHS,adofs,cmtman);

  // complete the matrices
  thermcontLM->Complete(*adofs,*adofs);
  thermcontTEMP->Complete(*amdofs,*adofs);

  /**********************************************************************/
  /* Modification of the stiff matrix and rhs towards thermo contact    */
  /**********************************************************************/

  /**********************************************************************/
  /* Create inv(D)                                                      */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*adofs,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  RCP<LINALG::SparseMatrix> mhatmatrix;
  mhatmatrix = LINALG::MLMultiply(*invd,false,*mmatrix,false,false,false);

  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  RCP<LINALG::SparseMatrix> kteffmatrix = rcp(new LINALG::SparseMatrix(*kteff));

  LINALG::SplitMatrix2x2(kteffmatrix,amdofs,ndofs,amdofs,ndofs,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,adofs,mdofs,adofs,mdofs,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,adofs,mdofs,ndofs,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,ndofs,tempmap,adofs,mdofs,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*problemrowmap,*feff,amdofs,fsm,ndofs,fn);
  LINALG::SplitVector(*amdofs,*fsm,adofs,fs,mdofs,fm);

  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,adofs,idofs,adofs,idofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,adofs,idofs,ndofs,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,adofs,idofs,mdofs,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,mdofs,tempmap,adofs,idofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*adofs));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*idofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*sdofs,*fs,adofs,fa,idofs,fi);

  // abbreviations for active and inactive set
  int aset = adofs->NumGlobalElements();
  int iset = idofs->NumGlobalElements();

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmnmod->Add(*kmn,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhatmatrix,true,*ksn,false,false,false,true);
  kmnmod->Add(*kmnadd,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmmmod->Add(*kmm,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhatmatrix,true,*ksm,false,false,false,true);
  kmmmod->Add(*kmmadd,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmimod->Add(*kmi,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhatmatrix,true,*kai,false,false,false,true);
    kmimod->Add(*kmiadd,false,1.0,1.0);
    kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  }

  // kmi: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmamod->Add(*kma,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhatmatrix,true,*kaa,false,false,false,true);
    kmamod->Add(*kmaadd,false,1.0,1.0);
    kmamod->Complete(kma->DomainMap(),kma->RowMap());
  }

  // kan: thermcontlm*invd*kan
  RCP<LINALG::SparseMatrix> kanmod;
  if (aset)
  {
    kanmod = LINALG::MLMultiply(*thermcontLM,false,*invd,false,false,false,true);
    kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
    kanmod->Complete(kan->DomainMap(),kan->RowMap());
  }

  // kam: thermcontlm*invd*kam
  RCP<LINALG::SparseMatrix> kammod;
  if (aset)
  {
    kammod = LINALG::MLMultiply(*thermcontLM,false,*invd,false,false,false,true);
    kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
    kammod->Complete(kam->DomainMap(),kam->RowMap());
  }

  // kai: thermcontlm*invd*kai
  RCP<LINALG::SparseMatrix> kaimod;
  if (aset && iset)
  {
    kaimod = LINALG::MLMultiply(*thermcontLM,false,*invd,false,false,false,true);
    kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
    kaimod->Complete(kai->DomainMap(),kai->RowMap());
  }

  // kaa: thermcontlm*invd*kaa
  RCP<LINALG::SparseMatrix> kaamod;
  if (aset)
  {
    kaamod = LINALG::MLMultiply(*thermcontLM,false,*invd,false,false,false,true);
    kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
    kaamod->Complete(kaa->DomainMap(),kaa->RowMap());
  }

  // Modifications towards rhs
  // FIXGIT: pay attention to genalpha
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*mdofs));
  mhatmatrix->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  // fa: mutliply with thermcontlm
  RCP<Epetra_Vector> famod;
  {
    famod = rcp(new Epetra_Vector(*adofs));
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*thermcontLM,false,*invd,false,false,false,true);
    temp->Multiply(false,*fa,*famod);
  }

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap,81,true,false,kteffmatrix->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap);

  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kmamod,false,1.0,1.0);
  if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);

  // add i submatrices to kteffnew
  if (iset) kteffnew->Add(*kin,false,1.0,1.0);
  if (iset) kteffnew->Add(*kim,false,1.0,1.0);
  if (iset) kteffnew->Add(*kii,false,1.0,1.0);
  if (iset) kteffnew->Add(*kia,false,1.0,1.0);

  // add a submatrices to kteffnew
  if (aset) kteffnew->Add(*kanmod,false,1.0,1.0);
  if (aset) kteffnew->Add(*kammod,false,1.0,1.0);
  if (aset && iset) kteffnew->Add(*kaimod,false,1.0,1.0);
  if (aset) kteffnew->Add(*kaamod,false,1.0,1.0);

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp;
  if (iset)
  {
    fiexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*fi,*fiexp);
    feffnew->Update(1.0,*fiexp,1.0);
  }

  // add a subvector to feffnew
  RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*famod,*famodexp);
    feffnew->Update(1.0,*famodexp,+1.0);
  }

  // add linearized thermo contact condition
  kteffnew->Add(*thermcontTEMP,false,-1.0,+1.0);

  // rhs zero till now
  //RCP<Epetra_Vector> linslipRHSexp = rcp(new Epetra_Vector(*problemrowmap_));
  //LINALG::Export(*linslipRHS_,*linslipRHSexp);
  //feffnew->Update(-1.0,*linslipRHSexp,1.0);

  // FillComplete kteffnew (square)
  kteffnew->Complete();

  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  kteff = kteffnew;
  feff = feffnew;

  return;
}

/*----------------------------------------------------------------------*
 | convert maps form structure dofs to thermo dofs            mgit 04/10 |
 *----------------------------------------------------------------------*/

void TSI::Algorithm::ConvertMaps(RCP<Epetra_Map>& slavedofs,
                                 RCP<Epetra_Map>& activedofs,
                                 RCP<Epetra_Map>& masterdofs,
                                 RCP<MORTAR::ManagerBase> cmtman)
{

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes/dofs
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // define local variables
    int slavecountnodes = 0;
    vector<int> myslavegids(slavenodes->NumMyElements());

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node = (StructureField().Discretization())->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleDM: Node ownership inconsistency!");

      myslavegids[slavecountnodes] = ((StructureField().Discretization())->Dof(1,node))[0];
      ++slavecountnodes;
    }

    // resize the temporary vectors
    myslavegids.resize(slavecountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gslavecountnodes;
    Comm().SumAll(&slavecountnodes,&gslavecountnodes,1);

    // create active node map and active dof map
    slavedofs = rcp(new Epetra_Map(gslavecountnodes,slavecountnodes,&myslavegids[0],0,Comm()));

    // active nodes/dofs
    const RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // define local variables
    int countnodes = 0;
    vector<int> mynodegids(activenodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node = (StructureField().Discretization())->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleDM: Node ownership inconsistency!");

      mynodegids[countnodes] = ((StructureField().Discretization())->Dof(1,node))[0];
      ++countnodes;
    }

    // resize the temporary vectors
    mynodegids.resize(countnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gcountnodes;
    Comm().SumAll(&countnodes,&gcountnodes,1);

    // create active node map and active dof map
    activedofs = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));

    // master nodes/dofs
    const RCP<Epetra_Map> masternodes = interface[m]->MasterRowNodes();

    // define local variables
    int mastercountnodes = 0;
    vector<int> mymastergids(masternodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = (StructureField().Discretization())->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleDM: Node ownership inconsistency!");

      mymastergids[mastercountnodes] = ((StructureField().Discretization())->Dof(1,node))[0];
      ++mastercountnodes;
    }

    // resize the temporary vectors
    mymastergids.resize(mastercountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gmastercountnodes;
    Comm().SumAll(&mastercountnodes,&gmastercountnodes,1);

    // create active node map and active dof map
    masterdofs = rcp(new Epetra_Map(gmastercountnodes,mastercountnodes,&mymastergids[0],0,Comm()));
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemnle mortar matrices in thermo dofs (active nodes)     mgit 04/10 |
 *----------------------------------------------------------------------*/

void TSI::Algorithm::AssembleDM(LINALG::SparseMatrix& dmatrix,
                                LINALG::SparseMatrix& mmatrix,
                                RCP<MORTAR::ManagerBase> cmtman)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // active
    const RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = (StructureField().Discretization())->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      /************************************************** D-matrix ******/
      if ((cnode->MoData().GetD()).size()>0)
      {
        vector<map<int,double> > dmap = cnode->MoData().GetD();
        int rowtemp = (StructureField().Discretization())->Dof(1,nodeges)[0];
        int rowdisp = cnode->Dofs()[0];
        double val = (dmap[0])[rowdisp];
        dmatrix.Assemble(val, rowtemp, rowtemp);
      }

      /************************************************** M-matrix ******/
      if ((cnode->MoData().GetM()).size()>0)
      {
        vector<map<int,double> > mmap = cnode->MoData().GetM();

        int rowtemp = (StructureField().Discretization())->Dof(1,nodeges)[0];

        set<int> mnodes = cnode->Data().GetMNodes();
        set<int>::iterator mcurr;

        // loop over all according master nodes
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = (interface[m]->Discret()).gNode(gid);
          DRT::Node* mnodeges = (StructureField().Discretization())->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          const int* mdofs = cmnode->Dofs();
          int coltemp = (StructureField().Discretization())->Dof(1,mnodeges)[0];
          double val = (mmap[0])[mdofs[0]];

          if(abs(val)>1e-12)
            mmatrix.Assemble(val, rowtemp, coltemp);
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble the thermal contact conditions                    mgit 04/10 |
 *----------------------------------------------------------------------*/

void TSI::Algorithm::AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                                LINALG::SparseMatrix& thermcontTEMP,
                                                Epetra_Vector& thermcontRHS,
                                                RCP<Epetra_Map> activedofs,
                                                RCP<MORTAR::ManagerBase> cmtman)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {

    // heat transfer coefficient, should be greater than zero
    double heattrancoeff = interface[m]->IParams().get<double>("HEATTRANSFERCOEFF");
    if (heattrancoeff <= 0)
      dserror("Error: Choose realistic heat transfer coefficient");

    // active
    const RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = (StructureField().Discretization())->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // with respect to Lagrange multipliers
      int row = (StructureField().Discretization())->Dof(1,nodeges)[0];
      int val = 1;
      thermcontLM.Assemble(val,row,row);

      // with resprect to temperatur dof
      /************************************************** D-matrix ******/
      if ((cnode->MoData().GetD()).size()>0)
      {
        vector<map<int,double> > dmap = cnode->MoData().GetD();
        int rowtemp = (StructureField().Discretization())->Dof(1,nodeges)[0];
        int rowdisp = cnode->Dofs()[0];
        double val = -heattrancoeff*(dmap[0])[rowdisp];
        if (abs(val)>1.0e-12) thermcontTEMP.Assemble(val,rowtemp,rowtemp);
      }

      /************************************************** M-matrix ******/
      if ((cnode->MoData().GetM()).size()>0)
      {
        vector<map<int,double> > mmap = cnode->MoData().GetM();

        int rowtemp = (StructureField().Discretization())->Dof(1,nodeges)[0];

        set<int> mnodes = cnode->Data().GetMNodes();
        set<int>::iterator mcurr;

        // loop over all according master nodes
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = (interface[m]->Discret()).gNode(gid);
          DRT::Node* mnodeges = (StructureField().Discretization())->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          const int* mdofs = cmnode->Dofs();
          int coltemp = (StructureField().Discretization())->Dof(1,mnodeges)[0];
          double val = +heattrancoeff*(mmap[0])[mdofs[0]];

          if(abs(val)>1e-12)
            thermcontTEMP.Assemble(val, rowtemp, coltemp);
        }
      }
    }
  }
  return;
}
/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
