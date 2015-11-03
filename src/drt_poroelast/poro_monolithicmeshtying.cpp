/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicmeshtying.cpp

 \brief base for porous media monolithic meshtying method

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Christoph Ager
// Originates from poro_monolithic

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "poro_monolithicmeshtying.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_poro_mortar.H"

#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_blocksparsematrix.H"


/*----------------------------------------------------------------------*
 | constructor                                                  2015    |
 *----------------------------------------------------------------------*/
POROELAST::MonolithicMeshtying::MonolithicMeshtying(const Epetra_Comm& comm,
                                                              const Teuchos::ParameterList& timeparams)
  : Monolithic(comm, timeparams),
    normrhsfactiven_(0.0),
    tolfres_ncoup_(0.0)
{
  //Initialize mortar adapter for meshtying interface
  mortar_adapter_ = Teuchos::rcp(new ADAPTER::CouplingPoroMortar);

  const int ndim = DRT::Problem::Instance()->NDim();
  std::vector<int> coupleddof (ndim,1);//1,1,1 should be in coupleddof
  //coupleddof[ndim]=0; // not necessary because structural discretization is used
  mortar_adapter_->Setup(StructureField()->Discretization(),
                         StructureField()->Discretization(),
                         coupleddof,
                         "Mortar");

  fvelactiverowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);
}

/*----------------------------------------------------------------------*
 | setup system (called in poro_dyn.cpp)                         2015   |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::SetupSystem()
{
  Monolithic::SetupSystem();
} // SetupSystem()

/*----------------------------------------------------------------------*
 | Evaluate override for meshtying matrices                     2015    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::Evaluate(Teuchos::RCP<const Epetra_Vector> x, bool firstiter)
{
  //evaluate monolithic system for newton iterations
  Monolithic::Evaluate(x,firstiter);

  //get state vectors to store in contact data container
  Teuchos::RCP<Epetra_Vector> fvel = Teuchos::rcp(new Epetra_Vector(*FluidField()->ExtractVelocityPart(FluidField()->Velnp())));
  fvel = FluidStructureCoupling().SlaveToMaster(fvel);

  //in the case that the problem is no pure poroelastic but a structure - poroelastic one, the poro skeleton
  //velocities must be distinct from structural velocities. structural velocities and displacements are not needed
  //for structural-poro meshtying.
  //so this vector splitting procedure is here to not store structural quantities for poroelastic meshtying,
  //because it is not needed

  Teuchos::RCP<Epetra_Vector> tmpsvel = Teuchos::rcp(new Epetra_Vector(*StructureField()->Velnp()));
  Teuchos::RCP<Epetra_Vector> svel = Teuchos::rcp(new Epetra_Vector(*FluidStructureCoupling().MasterDofMap()));
  LINALG::Export(*tmpsvel, *svel);

  Teuchos::RCP<Epetra_Vector> tmpsdisp = Teuchos::rcp(new Epetra_Vector(*StructureField()->Dispnp()));
  Teuchos::RCP<Epetra_Vector> sdisp = Teuchos::rcp(new Epetra_Vector(*FluidStructureCoupling().MasterDofMap()));
  LINALG::Export(*tmpsdisp, *sdisp);

  Teuchos::RCP<const Epetra_Vector> fpres = FluidField()->ExtractPressurePart(FluidField()->Velnp());
  Teuchos::RCP<Epetra_Vector> modfpres = Teuchos::rcp(new Epetra_Vector(*FluidField()->VelocityRowMap(),true));

  const int ndim = DRT::Problem::Instance()->NDim();
  int* mygids = fpres->Map().MyGlobalElements();
  double* val = fpres->Values();
  for (int i = 0; i < fpres->MyLength() ; ++i)
  {
    int gid = mygids[i]-ndim;//should be ndim
    modfpres->ReplaceGlobalValues(1, &val[i], &gid);
  }
  modfpres = FluidStructureCoupling().SlaveToMaster(modfpres);

  Teuchos::RCP<LINALG::SparseMatrix> f = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(systemmatrix_->Matrix(1,1)));
  Teuchos::RCP<LINALG::SparseMatrix> k_fs = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(systemmatrix_->Matrix(1,0)));

  Teuchos::RCP<Epetra_Vector> frhs = Extractor()->ExtractVector(rhs_,1);

  mortar_adapter_->EvaluatePoroMt(fvel,svel,modfpres,sdisp,StructureField()->Discretization(),f,k_fs,frhs,FluidStructureCoupling(),FluidField()->DofRowMap());


  systemmatrix_->Assign(1,1,LINALG::Copy,*f);
  systemmatrix_->Assign(1,0,LINALG::Copy,*k_fs);

  Extractor()->InsertVector(*frhs, 1, *rhs_);

  //because the mesh tying interface stays the same, the map extractors for a separate convergence check
  //of the mesh tying fluid coupling condition is only build once
  if((iter_== 1) and (Step()==1))
    SetupExtractor();

  return;
} //Evaluate()

/*----------------------------------------------------------------------*
 | Update override for meshtying matrices and LMP               2015    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::Update()
{
  Monolithic::Update();
  mortar_adapter_->UpdatePoroMt();
return;
}//POROELAST::MonolithicMeshtying::Update()

/*----------------------------------------------------------------------*
| Recover the Lagrange multipliers after newton step            2015    |
*-----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<const Epetra_Vector> x)
{
  Monolithic::RecoverLagrangeMultiplierAfterNewtonStep(x);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  ExtractFieldVectors(x,sx,fx);

  //RecoverStructuralLM
  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*sx));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*fx));

  mortar_adapter_->RecoverFluidLMPoroMt(tmpsx,tmpfx);

  return;
}//POROELAST::MonolithicMeshtying::RecoverLagrangeMultiplierAfterNewtonStep

//mostly copied from Monolithic Method
/*----------------------------------------------------------------------*
 |   evaluate poroelastic meshtying specific constraint            2015 |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::BuildConvergenceNorms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = UTILS::CalculateVectorNorm(vectornormfres_,rhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_s;

  Teuchos::RCP<const Epetra_Vector> rhs_f;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel;
  //split velocity into part of normal coupling & tangential condition and velocities
  Teuchos::RCP<const Epetra_Vector> rhs_fvel_activen;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel_other;
  Teuchos::RCP<const Epetra_Vector> rhs_fpres;

  // process structure unknowns of the first field (structure)
  rhs_s = Extractor()->ExtractVector(rhs_, 0);

  // process fluid unknowns of the second field
  rhs_f = Extractor()->ExtractVector(rhs_, 1);
  rhs_fvel = FluidField()->ExtractVelocityPart(rhs_f);
  //now split it
  rhs_fvel_activen = FluidVelActiveDofExtractor()->ExtractVector(rhs_fvel,0);
  rhs_fvel_other = FluidVelActiveDofExtractor()->ExtractVector(rhs_fvel,1);
  //pressure is treated separately anyway
  rhs_fpres = FluidField()->ExtractPressurePart(rhs_f);

  if(porositydof_)
  {
    dserror("porosity dof not implemented for poro_monolithicmeshtying");
    //consult method of mother class for further hints how to do this
  }
  else
  {
    normrhsstruct_  = UTILS::CalculateVectorNorm(vectornormfres_,rhs_s);
  }

  normrhsfluid_     = UTILS::CalculateVectorNorm(vectornormfres_,rhs_f);
  normrhsfluidvel_  = UTILS::CalculateVectorNorm(vectornormfres_,rhs_fvel_other);
  //residual norm of normal coupling condition on poro-fluid
  normrhsfactiven_   = UTILS::CalculateVectorNorm(vectornormfres_,rhs_fvel_activen);

  normrhsfluidpres_ = UTILS::CalculateVectorNorm(vectornormfres_,rhs_fpres);


  //------------------------------------------------------------- build residual increment norms
  // can stay exactly the same because a monolithic scheme with the same increments as without meshtying is used
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincs;
  Teuchos::RCP<const Epetra_Vector> interincf;
  Teuchos::RCP<const Epetra_Vector> interincfvel;
  Teuchos::RCP<const Epetra_Vector> interincfpres;
  // process structure unknowns of the first field
  interincs = Extractor()->ExtractVector(iterinc_, 0);
  // process fluid unknowns of the second field
  interincf = Extractor()->ExtractVector(iterinc_, 1);
  interincfvel = FluidField()->ExtractVelocityPart(interincf);
  interincfpres = FluidField()->ExtractPressurePart(interincf);

  normincstruct_     = UTILS::CalculateVectorNorm(vectornorminc_,interincs);
  normincfluid_         = UTILS::CalculateVectorNorm(vectornorminc_,interincf);
  normincfluidvel_      = UTILS::CalculateVectorNorm(vectornorminc_,interincfvel);
  normincfluidpres_     = UTILS::CalculateVectorNorm(vectornorminc_,interincfpres);

  return;
}//POROELAST::MonolithicMeshtying::BuildConvergenceNorms()

/*----------------------------------------------------------------------*
 |   setup meshtying activedof extractors                          2015 |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::SetupExtractor()
{
  //some maps and vectors
  Teuchos::RCP<Epetra_Map> factivenmap;
  Teuchos::RCP<Epetra_Map> factivenmapcomplement;
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidveldofmapvec;

  //get activemap from poro lagrange strategy of the adapter
  factivenmap = mortar_adapter_->GetPoroStrategy()->FluidActiveNDofMap();

  //build the complement part of the map
  factivenmapcomplement = LINALG::SplitMap(*FluidField()->VelocityRowMap(),*factivenmap);

  //write things into the vector for ->Setup
  fluidveldofmapvec.push_back(factivenmap);
  fluidveldofmapvec.push_back(factivenmapcomplement);

  fvelactiverowdofmap_->Setup(*FluidField()->VelocityRowMap(),fluidveldofmapvec);

  return;
}//POROELAST::MonolithicMeshtying::SetupExtractor()

//mostly copied from Monolithic method
/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)                       |
 |         with separate check of the meshtying constraint       2015   |
 *----------------------------------------------------------------------*/
bool POROELAST::MonolithicMeshtying::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  //convinc can stay the same because the increments are the same as without meshtying
  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      convinc = norminc_ < tolinc_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      convinc = ( normincstruct_    < tolinc_struct_   and
                  normincfluidvel_  < tolinc_velocity_ and
                  normincfluidpres_ < tolinc_pressure_ and
                  normincporo_      < tolinc_porosity_
                );
      break;
    default:
      dserror("Cannot check for convergence of residual values!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      convfres = normrhs_ < tolfres_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      convfres = ( normrhsstruct_    < tolfres_struct_   and
                   normrhsfluidvel_  < tolfres_velocity_ and
                   normrhsfluidpres_ < tolfres_pressure_ and
                   normrhsporo_      < tolfres_porosity_ and
                   normrhsfactiven_  < tolfres_ncoup_
                );
      break;
    default:
    dserror("Cannot check for convergence of residual forces!");
    break;
  }

  // combine increments and forces
  bool conv = false;
  if (combincfres_==INPAR::POROELAST::bop_and)
    conv = convinc and convfres;
  else if (combincfres_==INPAR::POROELAST::bop_or)
    conv = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;
}  // POROELAST::MonolithicMeshtying::Converged()

/*----------------------------------------------------------------------*
 | setup solver for monolithic system                                   |
 |            with exisiting meshtying interface               2015     |
 *----------------------------------------------------------------------*/
bool POROELAST::MonolithicMeshtying::SetupSolver()
{
  Monolithic::SetupSolver();

  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroelastdyn = DRT::Problem::Instance()->PoroelastDynamicParams();

  tolfres_ncoup_ = poroelastdyn.get<double> ("TOLRES_NCOUP");

  return true;
}//POROELAST::Monolithic::SetupSolver()

/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file              |
 | originally by lw 12/07, tk 01/08                           2015      |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::PrintNewtonIterHeaderStream(std::ostringstream& oss)
{
  oss<<"------------------------------------------------------------"<<std::endl;
  oss<<"                   Newton-Raphson Scheme                    "<<std::endl;
  oss<<"                NormRES "<<VectorNormString(vectornormfres_);
  oss<<"     NormINC "<<VectorNormString(vectornorminc_)<<"                    "<<std::endl;
  oss<<"------------------------------------------------------------"<<std::endl;

  // enter converged state etc
  oss << "numiter";

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(15) << "abs-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_<<")";
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss <<std::setw(15)<< "abs-s-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_struct_<<")";
      if(porositydof_)
        oss <<std::setw(15)<< "abs-poro-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_porosity_<<")";
      oss <<std::setw(15)<< "abs-fvel-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_velocity_<<")";
      oss <<std::setw(15)<< "abs-fpres-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_pressure_<<")";
      oss <<std::setw(15)<< "abs-fncoup-res"<<"("<<std::setw(5)<<std::setprecision(2) <<tolfres_ncoup_<<")";
      break;
    default:
      dserror("Unknown or undefined convergence form for residual.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs_global :
      oss <<std::setw(15)<< "abs-inc"<<"("<<std::setw(5)<<std::setprecision(2) <<tolinc_<<")";
      break;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss <<std::setw(15)<< "abs-s-inc"<<"("<<std::setw(5)<<std::setprecision(2) <<tolinc_struct_<<")";
      if(porositydof_)
        oss <<std::setw(15)<< "abs-poro-inc"<<"("<<std::setw(5)<<std::setprecision(2) <<tolinc_porosity_<<")";
      oss <<std::setw(15)<< "abs-fvel-inc"<<"("<<std::setw(5)<<std::setprecision(2) <<tolinc_velocity_<<")";
      oss <<std::setw(15)<< "abs-fpres-inc"<<"("<<std::setw(5)<<std::setprecision(2) <<tolinc_pressure_<<")";
      break;
    default:
      dserror("Unknown or undefined convergence form for increment.");
    break;
  }

  return;
} //POROELAST::MonolithicMeshtying::PrintNewtonIterHeaderStream


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                             |
 | originally by lw 12/07, tk 01/08                          2015       |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicMeshtying::PrintNewtonIterTextStream(std::ostringstream& oss)
{
  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(22) << std::setprecision(5) << std::scientific
          << normrhs_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      break;
    default:
      dserror("Unknown or undefined convergence form for global residual.");
      break;
  }
  // increments
  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs_global :
      oss << std::setw(22) << std::setprecision(5) << std::scientific << norminc_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      break;
    default:
      dserror("Unknown or undefined convergence form for global increment.");
    break;
  }

  // --------------------------------------------------------single field test
  switch ( normtypefres_ )
  {
    case INPAR::POROELAST::convnorm_abs_singlefields :
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsstruct_;
      if(porositydof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidpres_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfactiven_;
    break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("Unknown or undefined convergence form for single field residual.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs_singlefields :
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincstruct_;
      if(porositydof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normincporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidpres_;
    break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("Unknown or undefined convergence form for single field increment.");
    break;
  }

}//POROELAST::MonolithicMeshtying::PrintNewtonIterTextStream()
