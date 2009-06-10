/*----------------------------------------------------------------------*/
/*!
\file elast_summand.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELAST_SUMMAND_H
#define ELAST_SUMMAND_H

#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"

namespace MAT
{

  namespace ELASTIC
  {

    namespace PAR { }

    class Summand
    {
    public:

      /// standard constructor
      Summand() { ; }

      /// destructor
      virtual ~Summand() { ; }

      virtual INPAR::MAT::MaterialType MaterialType() const = 0;

      /// Create summand object by input parameter ID
      static Teuchos::RCP<Summand> Factory(
        int matnum  ///< material ID
        );

      virtual void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<8,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
        ) = 0;

      virtual void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< see above
        LINALG::Matrix<5,1>& delta,  ///< see above
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
        ) = 0;
      
      virtual void AddCoefficientsPrincipalAniso
      (
          bool& havecoefficients,  ///< non-zero coefficients of this type are provided
          LINALG::Matrix<3,1>& gamma,  ///< see above
          LINALG::Matrix<15,1>& delta,  ///< see above
          const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor         
      ) = 0;
      
      /// To be done and commented by Moritz :-)      
      virtual void AddCoefficientsModifiedAniso
      (
          bool& havecoefficients,  ///< non-zero coefficients of this type are provided
          LINALG::Matrix<3,1>& gamma,  ///< see above
          LINALG::Matrix<15,1>& delta,  ///< see above
          const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor         
      ) = 0;
      
    private:
      

    };  // class Summand

  }  // namespace ELAST

}

#endif

#endif
