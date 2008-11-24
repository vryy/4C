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
