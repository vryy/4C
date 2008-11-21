/*----------------------------------------------------------------------*/
/*!
\file inpar_structure.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_STRUCTURE_H
#define INPAR_STRUCTURE_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for structural problems
  namespace STR
  {

    /// Type of structural stress output
    /// (this enum represents the input file parameter STRUCT_STRESS)
    enum StressType
    {
      stress_none,    ///< no stress output
      stress_cauchy,  ///< output of Cauchy stresses
      stress_2pk      ///< output of 2nd Piola-Kirchhoff stresses
    };

    /// Type of structural strain output
    /// (this enum represents the input file parameter STRUCT_STRAIN)
    enum StrainType
    {
      strain_none,    ///< no strain output
      strain_ea,      ///< output of Euler-Almansi strains
      strain_gl       ///< output of Green-Lagrange strains
    };


  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STRUCTURE_H
#endif  // #ifdef CCADISCRET
