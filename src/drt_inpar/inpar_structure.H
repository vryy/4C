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

    /// @name Constraints
    //@{

    /// possible constraint solvers
    enum ConSolveAlgo
    {
      consolve_direct,  ///< ???
      consolve_iterative  ///< ???
    };

    //@}

    /// @name Output
    //@{

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

    //@}

    /// @name Time adapivity
    //@{

    /// type of adaptivity in time
    enum TimAdaKind
    {
      timada_kind_none,           ///< no time adaptivity
      timada_kind_zienxie,        ///< Zienkiewicz-Xie indicator
      timada_kind_ab2             ///< Adams-Bahsforth2 indicator
    };

    //@}

    /// type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_vague = 0,  ///< undetermined norm
      norm_l1,         ///< L1/linear norm
      norm_l2,         ///< L2/Euclidean norm
      norm_rms,        ///< root mean square (RMS) norm
      norm_inf         ///< Maximum/infinity norm
    };

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STRUCTURE_H
#endif  // #ifdef CCADISCRET
