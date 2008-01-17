
#ifdef CCADISCRET

#include "../drt_lib/drt_dserror.H"
#include "mfsi_overlappreccondoperator.H"

#include <string>

#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include <Teuchos_dyn_cast.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::OverlappingPCOperator::OverlappingPCOperator(
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > structure_lows,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > fluid_lows,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > ale_lows,
  Teuchos::RCP<const Thyra::DefaultBlockedLinearOp<double> > blockFsiOp
  )
  : structure_lows_(structure_lows),
    fluid_lows_(fluid_lows),
    ale_lows_(ale_lows),
    blockFsiOp_(blockFsiOp)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > MFSI::OverlappingPCOperator::range() const
{
  return blockFsiOp_->domain();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > MFSI::OverlappingPCOperator::domain() const
{
  return blockFsiOp_->range();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlappingPCOperator::apply(Thyra::EConj conj,
                                        const Thyra::MultiVectorBase<double>& X,
                                        Thyra::MultiVectorBase<double>* Y,
                                        double alpha,
                                        double beta) const
{
  Teuchos::RCP<Thyra::MultiVectorBase<double> > T;
  if (beta==0.0)
  {
    T = Teuchos::rcp(Y,false);
  }
  else
  {
    T = Thyra::createMembers(Y->range(),Y->domain()->dim());
    Thyra::scale(beta,Y);
  }

  Teuchos::RCP<const Thyra::LinearOpBase<double> > structInnerOp = blockFsiOp_->getBlock(0,0);
  Teuchos::RCP<const Thyra::LinearOpBase<double> > fluidInnerOp  = blockFsiOp_->getBlock(1,1);
  Teuchos::RCP<const Thyra::LinearOpBase<double> > fluidBoundOp  = blockFsiOp_->getBlock(1,0);
  Teuchos::RCP<const Thyra::LinearOpBase<double> > aleInnerOp    = blockFsiOp_->getBlock(2,2);
  Teuchos::RCP<const Thyra::LinearOpBase<double> > aleBoundOp    = blockFsiOp_->getBlock(2,0);

  //
  // - vektor X in seine drei Bestandteile zerlegen
  // - Die Gleichungen nacheinander lösen, dabei die rechten Seite
  // entsprechend der Vorgaben herstellen
  // - die Lösungen wieder in den Vektor Y zusammensetzen
  //

  const Thyra::DefaultProductVector<double>& x =
    Teuchos::dyn_cast<const Thyra::DefaultProductVector<double> >(X);

  Thyra::DefaultProductVector<double>* y =
    dynamic_cast<Thyra::DefaultProductVector<double>*>(T.get());
  if (y==NULL)
  {
    dserror("illegal pointer type");
  }

  // Extract vector blocks

  Teuchos::RCP<const Thyra::VectorBase<double> > sx = x.getVectorBlock(0);
  Teuchos::RCP<const Thyra::VectorBase<double> > fx = x.getVectorBlock(1);
  Teuchos::RCP<const Thyra::VectorBase<double> > ax = x.getVectorBlock(2);

  Teuchos::RCP<Thyra::VectorBase<double> > sy = y->getNonconstVectorBlock(0);
  Teuchos::RCP<Thyra::VectorBase<double> > fy = y->getNonconstVectorBlock(1);
  Teuchos::RCP<Thyra::VectorBase<double> > ay = y->getNonconstVectorBlock(2);

  // Create the LOWSB objects that will be used to solve the linear systems

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invStruct
    = Thyra::linearOpWithSolve<double>(*structure_lows_,structInnerOp);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invFluid
    = Thyra::linearOpWithSolve<double>(*fluid_lows_,fluidInnerOp);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invAle
    = Thyra::linearOpWithSolve<double>(*ale_lows_,aleInnerOp);

  // Solve structure equations for sy with the rhs sx

  Thyra::SolveStatus<double> status = Thyra::solve<double>(*invStruct,Thyra::NOTRANS,*sx,&*sy);
  if (status.solveStatus!=Thyra::SOLVE_STATUS_CONVERGED)
  {
    dserror("structure solve failed");
  }

  // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy

  Teuchos::RCP<Thyra::VectorBase<double> > fb = Thyra::createMember(*fx->range());
  fluidBoundOp->apply(Thyra::NONCONJ_ELE,*sy,&*fb,-1.0);
  Thyra::update(1.0,*fx,&*fb);
  status = Thyra::solve<double>(*invFluid,Thyra::NOTRANS,*fb,&*fy);
  if (status.solveStatus!=Thyra::SOLVE_STATUS_CONVERGED)
  {
    dserror("fluid solve failed");
  }

  // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

  Teuchos::RCP<Thyra::VectorBase<double> > ab = Thyra::createMember(*ax->range());
  aleBoundOp->apply(Thyra::NONCONJ_ELE,*sy,&*ab,-1.0);
  Thyra::update(1.0,*ax,&*ab);
  status = Thyra::solve<double>(*invAle,Thyra::NOTRANS,*ab,&*ay);
  if (status.solveStatus!=Thyra::SOLVE_STATUS_CONVERGED)
  {
    dserror("ale solve failed");
  }

#if 0
  const SolveCriteria<double>* solveCriteria = fwdSolveCriteria_.get();

  Thyra::assign(T.get(), 0.0); // Have to initialize before solve!

  Thyra::SolveStatus<double> solveStatus = Thyra::solve<double>(*structure_lows_->getConstObj(),
                                                                Thyra::NOTRANS,
                                                                X,
                                                                &*T,
                                                                solveCriteria);
  if (solveCriteria && solveStatus.solveStatus!=SOLVE_STATUS_CONVERGED)
  {
    dserror("solve failed");
  }
#endif

  if (beta==0.0)
  {
    if (alpha!=1)
      Thyra::scale(alpha,Y);
  }
  else
  {
    Thyra::update(alpha, *T, Y);
  }
}

#endif
