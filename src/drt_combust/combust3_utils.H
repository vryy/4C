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
  /// hold arrays for all - potentially enriched - approximation functions
  /// and their derivatives with respect to the physical coordinates
  ///
  /// Let P be the approximation function and N the shape function.
  /// Then, for standard FEM, P = N
  ///
  /// Two examples:
  ///
  /// u^h = SUM_I P^I u^I
  ///  u^I are the free parameters (nodal unknowns in the standard FEM) and
  ///  normal and enriched parameters combined in the XFEM
  ///
  /// | P1 |   | N1     |
  /// | P2 |   | N2     |
  /// | P3 | = | N1*psi |
  /// | P4 |   | N2*psi |
  ///   ...
  /// similar
  /// du^h/dx = SUM_I dP^I/dx u^I
  ///   P is the approximation function
  ///
  /// | P1/dx |   | N1/dx                   |
  /// | P2/dx |   | N2/dx                   |
  /// | P3/dx | = | dN1/dx*psi + N1*dpsi/dx |
  /// | P4/dx |   | dN2/dx*psi + N2*dpsi/dx |
  ///   ...
  ///   Note the chain rule above!
  template <size_t numparam>
  struct ApproxFunc
  {
    LINALG::Matrix<numparam,1>  d0;
    LINALG::Matrix<numparam,1>  dx;
    LINALG::Matrix<numparam,1>  dy;
    LINALG::Matrix<numparam,1>  dz;
    LINALG::Matrix<numparam,1>  dxdx;
    LINALG::Matrix<numparam,1>  dxdy;
    LINALG::Matrix<numparam,1>  dxdz;
    LINALG::Matrix<numparam,1>  dydx;
    LINALG::Matrix<numparam,1>  dydy;
    LINALG::Matrix<numparam,1>  dydz;
    LINALG::Matrix<numparam,1>  dzdx;
    LINALG::Matrix<numparam,1>  dzdy;
    LINALG::Matrix<numparam,1>  dzdz;
  };

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
