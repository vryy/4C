
#include "drt_utils_polynomial.H"

template<int nsd_> DRT::UTILS::PolynomialSpaceCache<nsd_> * DRT::UTILS::PolynomialSpaceCache<nsd_> ::instance_;

template <int nsd_>
DRT::UTILS::PolynomialSpaceCache<nsd_> & DRT::UTILS::PolynomialSpaceCache<nsd_>::Instance()
{
  if ( instance_==NULL )
  {
    instance_ = new PolynomialSpaceCache<nsd_>;
  }
  return *instance_;
}

template <int nsd_>
void DRT::UTILS::PolynomialSpaceCache<nsd_>::Done()
{
  ps_cache_.clear();
  delete instance_;
  instance_ = NULL;
}

template <int nsd_>
Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_> > DRT::UTILS::PolynomialSpaceCache<nsd_>::Create(PolynomialSpaceParams params)
{
  typename std::map<PolynomialSpaceParams, Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_> > >::iterator
    i = ps_cache_.find(params);
  if ( i!=ps_cache_.end() )
  {
    return i->second;
  }
  std::cout<<"size "<<ps_cache_.size()<<std::endl;

  // this is expensive and should not be done too often
  Teuchos::RCP<PolynomialSpace<nsd_> > ps;
  ps = Teuchos::rcp( new PolynomialSpace<nsd_>(params) );

  ps_cache_[params] = ps;

  return ps;
}


// explicit instantations
template class DRT::UTILS::PolynomialSpaceCache<1>;
template class DRT::UTILS::PolynomialSpaceCache<2>;
template class DRT::UTILS::PolynomialSpaceCache<3>;

//template<> DRT::UTILS::PolynomialSpaceCache<1> * DRT::UTILS::PolynomialSpaceCache<1>::instance_ = NULL;
//template<> DRT::UTILS::PolynomialSpaceCache<2> * DRT::UTILS::PolynomialSpaceCache<2>::instance_ = NULL;
