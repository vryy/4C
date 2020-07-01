/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial materials


\level 2

*/
/*---------------------------------------------------------------------*/

// All header files should have #define guards to prevent multiple inclusion.
// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef TUTORIAL_MAT_H
#define TUTORIAL_MAT_H

#include "../../linalg/linalg_fixedsizematrix.H"

namespace TUTORIAL
{
  // FORWARD DECLARATION
  // Here we declare the class 'TutorialElement'.
  // This is needed because an pointer of this type is
  // provided to Evaluate(...). We could also include
  // the tutorial_ele.H-file but it should be avoided
  // to include headers in headers. This increases
  // dependencies in the code and slows down compilation
  // and linking.
  namespace ELEMENTS
  {
    class TutorialElement;
  }

  namespace MATERIAL
  {
    class TutorialMat
    {
     public:
      /*!
      \brief constructor (ctor)

      \param young (in) : Young's Modulus of this material

      \note Use the C++ keyword explicit for constructors with one argument.
            See LNM/Baci code style guide for further information.
      */
      explicit TutorialMat(const double young);

      /*!
      \brief destructor (dtor)

      \note This Destructor is declared virtual, since class NonlinearTruss contains virtual
      methods. This is an indicator, that this class is supposed to be part of some inheritance
      line. Check the internet on why it is then important to declare the destructor virtual.
      */
      virtual ~TutorialMat(){};

      /*!
      \brief evaluate stress and material linearization (called for each integration point of
      elements)

      \note Pure virtual evaluation routine for material law and its linearization.
            This pure virtual function ( with = 0 ) HAS TO be implemented by the
            subclass inheriting from TutorialMat (see TutorialMatStVenantKirchhoff).

       \param ele (in)            : element that is evaluated
       \param stress (out)        : stress at integration point (filled in this method)
       \param linearization (out) : element matrix contribution of integration point (filled in this
      method) \param strain (out)        : Green-Lagrange strain at integration point \param defgrd
      (in)         : deformation gradient at integration point \param J (in)              : Jacobian
      determinant at integration point \param B (in)              : B Operator at integration point

      */
      virtual void EvaluateMaterial(TUTORIAL::ELEMENTS::TutorialElement* ele, double* stress,
          LINALG::Matrix<2, 2>* linearization, const double strain, const double defgrd,
          const double J, const LINALG::Matrix<2, 1> B) = 0;

      /// return young's modulus
      double YoungModulus() { return Young_; };

     private:
      /// young's modulus
      double Young_;

    };  // class TutorialMat

    /*----------------------------------------------------------------------*/

    /// the St.-Venant-Kirchhoff Material inherits from the base class TutorialMat
    class TutorialMatStVenantKirchhoff : public TutorialMat
    {
     public:
      /*!
      \brief constructor (ctor)

      \param young (in) : Young's Modulus of this material

      \note Use the C++ keyword explicit for constructors with one argument.
            See LNM/Baci code style guide for further information.
      */
      TutorialMatStVenantKirchhoff(const double young);

      /*!
      \brief destructor (dtor)

      \note This Destructor is declared virtual, since class NonlinearTruss contains virtual
      methods. This is an indicator, that this class is supposed to be part of some inheritance
      line. Check the internet on why it is then important to declare the destructor virtual.
      */
      virtual ~TutorialMatStVenantKirchhoff(){};

      /*!
      \brief evaluate stress and material linearization (called for each integration point of
      elements)

      \note This is the implementation of EvaluateMaterial for the St. Venant-Kirchhoff material.
            Any other class that would inherit from 'class TutorialMat', would have ti implement
            this method, too. SEE 'class TutorialMat'.

       \param ele (in)            : element that is evaluated
       \param stress (out)        : stress at integration point (filled in this method)
       \param linearization (out) : element matrix contribution of integration point (filled in this
      method) \param strain (out)        : Green-Lagrange strain at integration point \param defgrd
      (in)         : deformation gradient at integration point \param J (in)              : Jacobian
      determinant at integration point \param B (in)              : B Operator at integration point

      */
      virtual void EvaluateMaterial(TUTORIAL::ELEMENTS::TutorialElement* ele, double* stress,
          LINALG::Matrix<2, 2>* linearization, const double strain, const double defgrd,
          const double J, const LINALG::Matrix<2, 1> B);


    };  // class TutorialMatStVenantKirchhoff


  }  // namespace MATERIAL
}  // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef TUTORIAL_MAT_H*/
