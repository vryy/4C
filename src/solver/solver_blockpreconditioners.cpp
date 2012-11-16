/*
 * solver_blockpreconditioners.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#include "../drt_lib/drt_dserror.H"

#include "solver_blockpreconditioners.H"

// include header files for concrete implementation
#include "bgs2x2_operator.H"                   // Lena's BGS implementation
#include "saddlepointpreconditioner.H"         // Tobias' saddle point preconditioner
#include "solver_cheapsimplepreconditioner.H"  // Tobias' CheapSIMPLE

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::SimplePreconditioner::SimplePreconditioner( FILE * outfile,
                                                            Teuchos::ParameterList & params)
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::SimplePreconditioner::Setup( bool create,
                                                  Epetra_Operator * matrix,
                                                  Epetra_MultiVector * x,
                                                  Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    // SIMPLER does not need copy of preconditioning matrix to live
    // SIMPLER does not use the downwinding installed here, it does
    // its own downwinding inside if desired

    // free old matrix first
    P_ = Teuchos::null;

    // temporary hack: distinguish between "old" SIMPLER_Operator (for fluid
    // only) and "new" more general test implementation
    bool mt = params_.get<bool>("MESHTYING",false);
    bool co = params_.get<bool>("CONTACT",false);
    bool cstr = params_.get<bool>("CONSTRAINT",false);
    bool fl = params_.isSublist("SIMPLER") || params_.get<bool>("FLUID",false); //params_.get<bool>("FLUIDSIMPLE",false); // SIMPLE for fluids
    bool elch = params_.get<bool>("ELCH",false);
    if (mt || co || cstr)
    {
      // adapt ML null space for contact/meshtying/constraint problems
      RCP<BlockSparseMatrixBase> A = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp( matrix, false ));
      if (A==Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      // fix null space for "Inverse1"
      //      {
      //        const Epetra_Map& oldmap = A->FullRowMap();
      //        const Epetra_Map& newmap = A->Matrix(0,0).EpetraMatrix()->RowMap();
      //        LINALG::Solver::FixMLNullspace("Inverse1",oldmap, newmap, params_.sublist("Inverse1"));
      //      }

      // adapt null space for constraint equations
      Teuchos::ParameterList& inv2 = params_.sublist("Inverse2");
      if(inv2.isSublist("ML Parameters"))
      {
        // Schur complement system (1 degree per "node") -> standard nullspace
        inv2.sublist("ML Parameters").set("PDE equations",1);
        inv2.sublist("ML Parameters").set("null space: dimension",1);
        const int plength = (*A)(1,1).RowMap().NumMyElements();
        RCP<std::vector<double> > pnewns = Teuchos::rcp(new std::vector<double>(plength,1.0));
        //TODO: vector<double> has zero length for particular cases (e.g. no Lagrange multiplier on this processor)
        //      -> RCP for the null space is set to NULL in Fedora 12 -> dserror
        //      -> RCP points to a random memory field in Fedora 8 -> RCP for null space is not NULL
        // Temporary work around (ehrl, 21.12.11):
        // In the case of plength=0 the vector<double> is rescaled (size 0 -> size 1, initial value 0) in order to avoid problems with ML
        // (ML expects an RCP for the null space != NULL)
        if (plength==0)
          pnewns->resize(1,0.0);
        inv2.sublist("ML Parameters").set("null space: vectors",&((*pnewns)[0]));
        inv2.sublist("ML Parameters").remove("nullspace",false);
        inv2.sublist("Michael's secret vault").set<RCP<std::vector<double> > >("pressure nullspace",pnewns);
      }

      P_ = Teuchos::rcp(new LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A,params_.sublist("Inverse1"),params_.sublist("Inverse2"),outfile_));
    }
    else if(fl || elch) // CheapSIMPLE for pure fluid problems
    {
      // adapt nullspace for splitted pure fluid problem
      int nv = 0; // number of velocity dofs
      int np = 0; // number of pressure dofs
      int ndofpernode = 0; // dofs per node
      int nlnode;

      const Epetra_Map& fullmap = matrix->OperatorRangeMap();
      const int length = fullmap.NumMyElements();

      RCP<BlockSparseMatrixBase> A = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp( matrix, false ));
      if (A==Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

      // this is a fix for the old SIMPLER sublist
      if(!params_.isSublist("Inverse1") && params_.isSublist("SIMPLER"))
      {
        Teuchos::ParameterList& inv1 = params_.sublist("Inverse1");
        inv1 = params_;
        inv1.remove("SIMPLER");
        inv1.remove("Inverse1",false);
        Teuchos::ParameterList& inv2 = params_.sublist("Inverse2");
        inv2 = params_.sublist("SIMPLER");
        params_.remove("SIMPLER");
        params_.sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
        params_.set("FLUID",true);
      }

      // fix null spae for ML inverses
      Teuchos::ParameterList& inv1 = params_.sublist("Inverse1");
      if(inv1.isSublist("ML Parameters"))
      {
        ndofpernode = inv1.sublist("NodalBlockInformation").get<int>("numdf",0);
        nv = inv1.sublist("NodalBlockInformation").get<int>("nv",0);
        np = inv1.sublist("NodalBlockInformation").get<int>("np",0);
        if(ndofpernode==0) dserror("cannot read numdf from NodalBlockInformation");
        if (nv==0 || np==0) dserror("nv or np == 0?");
        nlnode = length/ndofpernode;

        inv1.sublist("ML Parameters").set("PDE equations",nv);
        inv1.sublist("ML Parameters").set("null space: dimension",nv);

        const int vlength = A->Matrix(0,0).RowMap().NumMyElements();
        Teuchos::RCP<std::vector<double> > vnewns = Teuchos::rcp(new std::vector<double>(nv*vlength,0.0));

        for (int i=0;i<nlnode; ++i)
        {
          (*vnewns)[i*nv] = 1.0;
          (*vnewns)[vlength+i*nv+1] = 1.0;
          if(nv > 2) (*vnewns)[2*vlength+i*nv+2] = 1.0;
        }
        inv1.sublist("ML Parameters").set("null space: vectors",&((*vnewns)[0]));
        inv1.sublist("ML Parameters").remove("nullspace",false); // necessary??
        inv1.sublist("Michael's secret vault").set<RCP<std::vector<double> > >("velocity nullspace",vnewns);
      }

      Teuchos::ParameterList& inv2 = params_.sublist("Inverse2");
      if(inv2.isSublist("ML Parameters"))
      {
        inv2.sublist("ML Parameters").set("PDE equations",1);
        inv2.sublist("ML Parameters").set("null space: dimension",1);
        const int plength = A->Matrix(1,1).RowMap().NumMyElements();
        Teuchos::RCP<std::vector<double> > pnewns = Teuchos::rcp(new std::vector<double>(plength,1.0));
        inv2.sublist("ML Parameters").set("null space: vectors",&((*pnewns)[0]));
        inv2.sublist("ML Parameters").remove("nullspace",false); // necessary?
        inv2.sublist("Michael's secret vault").set<RCP<std::vector<double> > >("pressure nullspace",pnewns);
      }

      P_ = Teuchos::rcp(new LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner(A,params_.sublist("Inverse1"),params_.sublist("Inverse2"),outfile_));
    }
    //else if(!params_.isSublist("Inverse1") || !params_.isSublist("Inverse2"))
    else
    {
      //cout << "************************************************" << endl;
      //cout << "WARNING: SIMPLE for Fluid? expect bugs..." << endl;
      //cout << "************************************************" << endl;
      // Michaels old CheapSIMPLE for Fluid
      // TODO replace me by CheapSIMPLE_BlockPreconditioner

      //P_ = Teuchos::rcp(new LINALG::SOLVER::SIMPLER_Operator(Teuchos::rcp( matrix, false ),params_,params_.sublist("SIMPLER"),outfile_));
      dserror("old SIMPLE not supported any more");
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::AMGBSPreconditioner::AMGBSPreconditioner( FILE * outfile,
    Teuchos::ParameterList & params )
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::AMGBSPreconditioner::Setup( bool create,
                                                 Epetra_Operator * matrix,
                                                 Epetra_MultiVector * x,
                                                 Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    // free old matrix first
    P_ = Teuchos::null;

    // Params().sublist("AMGBS") just contains the Fluid Pressure Solver block
    // from the dat file (not needed anymore)
    P_ = Teuchos::rcp(new LINALG::SaddlePointPreconditioner(Teuchos::rcp( matrix, false ),params_,outfile_));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BGSPreconditioner::BGSPreconditioner( FILE * outfile,
    Teuchos::ParameterList & params, Teuchos::ParameterList & bgslist )
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params ),
    bgslist_( bgslist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BGSPreconditioner::Setup( bool create,
                                               Epetra_Operator * matrix,
                                               Epetra_MultiVector * x,
                                               Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    P_ = Teuchos::null;

    int numblocks = bgslist_.get<int>("numblocks");

    if (numblocks == 2) // BGS2x2
    {
      // check whether sublists for individual block solvers are present
      bool haveprec1 = params_.isSublist("Inverse1");
      bool haveprec2 = params_.isSublist("Inverse2");
      if (!haveprec1 or !haveprec2)
        dserror("individual block solvers for BGS2x2 need to be specified");

      int global_iter = bgslist_.get<int>("global_iter");
      double global_omega = bgslist_.get<double>("global_omega");
      int block1_iter = bgslist_.get<int>("block1_iter");
      double block1_omega = bgslist_.get<double>("block1_omega");
      int block2_iter = bgslist_.get<int>("block2_iter");
      double block2_omega = bgslist_.get<double>("block2_omega");
      bool fliporder = bgslist_.get<bool>("fliporder");

      P_ = Teuchos::rcp(new LINALG::BGS2x2_Operator(Teuchos::rcp( matrix, false ),
                                           params_.sublist("Inverse1"),
                                           params_.sublist("Inverse2"),
                                           global_iter,
                                           global_omega,
                                           block1_iter,
                                           block1_omega,
                                           block2_iter,
                                           block2_omega,
                                           fliporder,
                                           outfile_));
    }
    else
      dserror("Block Gauss-Seidel BGS2x2 is currently only implemented for a 2x2 system. Use BGSnxn for a common block Gauss-Seidel implementation (based on Teko package in Trilinos).");
  }
}
