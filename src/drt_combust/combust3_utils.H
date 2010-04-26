/*!
\file combust3_utils.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
#ifdef CCADISCRET
#ifndef COMBUST3_UTILS_H
#define COMBUST3_UTILS_H


namespace COMBUST
{

  template <size_t numparam>
  struct EnrViscs2
  {
    LINALG::Matrix<numparam,1> xx;
    LINALG::Matrix<numparam,1> xy;
    LINALG::Matrix<numparam,1> xz;
    LINALG::Matrix<numparam,1> yx;
    LINALG::Matrix<numparam,1> yy;
    LINALG::Matrix<numparam,1> yz;
    LINALG::Matrix<numparam,1> zx;
    LINALG::Matrix<numparam,1> zy;
    LINALG::Matrix<numparam,1> zz;
  };

template<class M1, class M2>
void fill_tau(
    const size_t         numparamtau,
    const M1&            shp_tau,
    const M2&            etau,
    LINALG::Matrix<3,3>& tau
    )
{
    tau.Clear();
    for (size_t i = 0; i < numparamtau; ++i)
    {
        const double fac = shp_tau(i);
        tau(0,0) += fac*etau(0,i);
        tau(1,1) += fac*etau(1,i);
        tau(2,2) += fac*etau(2,i);
        tau(0,1) += fac*etau(3,i);
        tau(0,2) += fac*etau(4,i);
        tau(1,2) += fac*etau(5,i);
    }
    tau(1,0) = tau(0,1);
    tau(2,0) = tau(0,2);
    tau(2,1) = tau(1,2);
}

  //! interpolate from nodal vector array to integration point vector using the shape function
  template <class M, class VS>
  static LINALG::Matrix<3,1> interpolateVectorFieldToIntPoint(
      const M&  eleVectorField,       ///< array with nodal vector values
      const VS& shp,                  ///< array with nodal shape function
      const int numparam              ///< number of parameters
      )
  {
    LINALG::Matrix<3,1> v(true);
    const int nsd = 3;
    for (int iparam = 0; iparam < numparam; ++iparam)
      for (int isd = 0; isd < nsd; ++isd)
        v(isd) += eleVectorField(isd,iparam)*shp(iparam);
    return v;
  }
} // namespace COMBUST

#endif
#endif  // #ifdef CCADISCRET
