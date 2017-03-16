/*----------------------------------------------------------------------*/
/*!
\file stru_aux.cpp

\brief structure-specific utils and auxiliary functions

\level 1

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238

*/
/*----------------------------------------------------------------------*/



#include "stru_aux.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/* Calculate vector norm */
double STR::AUX::CalculateVectorNorm
(
  const enum INPAR::STR::VectorNorm norm,
  const Teuchos::RCP<Epetra_Vector> vect,
  const int numneglect
)
{
  // L1 norm
  if (norm == INPAR::STR::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  else if (norm == INPAR::STR::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == INPAR::STR::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double) (vect->GlobalLength() - numneglect));
  }
  // infinity/maximum norm
  else if (norm == INPAR::STR::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::AUX::MapExtractor::Setup(const DRT::Discretization& dis, const Epetra_Map& fullmap, bool overlapping)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"BioGrCoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"AleWear",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FPSICoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"IMMERSEDCoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"ParticleWall",0,ndim)));

  mcs.SetupExtractor(dis,fullmap,*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > STR::AUX::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"StructAleCoupling");
  Teuchos::RCP<std::set<int> > condelements3 = DRT::UTILS::ConditionedElementMap(dis,"BioGrCoupling");
  Teuchos::RCP<std::set<int> > condelements4 = DRT::UTILS::ConditionedElementMap(dis,"AleWear");
  Teuchos::RCP<std::set<int> > condelements5 = DRT::UTILS::ConditionedElementMap(dis,"FPSICoupling");
  Teuchos::RCP<std::set<int> > condelements6 = DRT::UTILS::ConditionedElementMap(dis,"IMMERSEDCoupling");
  Teuchos::RCP<std::set<int> > condelements7 = DRT::UTILS::ConditionedElementMap(dis,"ParticleWall");


  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements3->begin(),condelements3->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements4->begin(),condelements4->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements5->begin(),condelements5->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements6->begin(),condelements6->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements7->begin(),condelements7->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

