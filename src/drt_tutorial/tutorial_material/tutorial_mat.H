/*!----------------------------------------------------------------------
\file tutorial_mat.H

\brief student's c++/baci tutorial materials

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef TUTORIAL_MAT_H
#define TUTORIAL_MAT_H

#include "../../linalg/linalg_fixedsizematrix.H"

namespace TUTORIAL {

/// FORWARD DECLARATION
/// Here we declare the class 'TutorialElement'.
/// This is needed because an pointer of this type is
/// provided to Evaluate(...). We could also include
/// the tutorial_ele.H-file but it should be avoided
/// to include headers in headers. This increases
/// dependencies in the code and slows down compilation
/// and linking.
namespace ELEMENTS {
class TutorialElement;
}

namespace MATERIAL {

class TutorialMat
{

public:

  /// constructor (ctor)
  TutorialMat(double young);

  /// empty destructor (dtor)
  virtual ~TutorialMat(){};

  /// pure virtual evaluation routine for material law and its linearization
  /// this pure virtual function ( with = 0 ) HAS TO be implemented by the
  /// subclass inheriting from TutorialMat (see TutorialMatStVenantKirchhoff).
  virtual void EvaluateMaterial(TUTORIAL::ELEMENTS::TutorialElement* ele, double* stress, LINALG::Matrix<2,2>* linearization, const double strain, const double defgrd, const double J, const LINALG::Matrix<2,1> B ) = 0;

  /// return young's modulus
  double YoungModulus(){return Young_;};

private:

  // young's modulus
  double Young_;

}; // class TutorialMat

/*----------------------------------------------------------------------*/

/// the St.-Venant-Kirchhoff Material inherits from the base class TutorialMat
class TutorialMatStVenantKirchhoff : public TutorialMat
{

public:

  /// constructor (ctor)
  TutorialMatStVenantKirchhoff(double young);

  /// empty destructor (dtor)
  virtual ~TutorialMatStVenantKirchhoff(){};

  /// here we implement the pure virtual
  virtual void EvaluateMaterial(TUTORIAL::ELEMENTS::TutorialElement* ele, double* stress, LINALG::Matrix<2,2>* linearization, const double strain, const double defgrd, const double J, const LINALG::Matrix<2,1> B);


}; // class TutorialMatStVenantKirchhoff


} // namespace MATERIAL
} // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef TUTORIAL_MAT_H*/
