/*----------------------------------------------------------------------*/
/*!
\file structure_utils.cpp

\brief structure-specific utils and auxiliary functions

<pre>
Maintainer: Thomas Klöppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
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
  const Teuchos::RCP<Epetra_Vector> vect
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
    return vectnorm/sqrt((double) vect->GlobalLength());
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
void STR::AUX::MapExtractor::Setup(const DRT::Discretization& dis, const Epetra_Map& fullmap)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,ndim)));
  mcs.SetupExtractor(dis,fullmap,*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > STR::AUX::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"StructAleCoupling");
  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

