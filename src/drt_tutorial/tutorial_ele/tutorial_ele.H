/*!----------------------------------------------------------------------
\file tutorial_ele.H

\brief student's c++/baci tutorial element

\level 3

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/

// All header files should have #define guards to prevent multiple inclusion.
// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef TUTORIAL_ELE_H
#define TUTORIAL_ELE_H

#include <Teuchos_RCP.hpp>

#include "../../linalg/linalg_fixedsizematrix.H"


namespace TUTORIAL
{
  /// FORWARD DECLARATION
  /// Here we declare the class 'TutorialMat'.
  /// This is needed because a member of this type is
  /// defined in this header. We could also include
  /// the tutorial_mat.H-file but it should be avoided
  /// to include headers in headers. This increases
  /// dependencies in the code and slows down compilation
  /// and linking.
  namespace MATERIAL
  {
    class TutorialMat;
  }

  namespace ELEMENTS
  {
    class TutorialElement
    {
     public:
      /// constructor (ctor)
      ///
      /// construct ele with gid='id'
      ///         E
      ///  o-------------o
      ///  X0            X1
      /// ---> X
      ///
      explicit TutorialElement(int id, double E, double X0, double X1);

      /// this method just returns the element id
      int Id() const { return id_; };

      /// this method just returns the number of nodes of this element
      int NumNode() const { return numnode_; };

      /// this method just returns the number of dofs of this element
      int NumDof() const { return numdof_; };


      /// this method just returns the initial position of node='node' ('node' = 0 or 1)
      double X(int node) const { return X_[node]; };

      /// this method just returns the one dof (degree of freedom) per node='node' ('node' = 0 or 1)
      double Dof(int node) const { return dofs_[node]; };

      /// return the material of this element
      const Teuchos::RCP<TUTORIAL::MATERIAL::TutorialMat> Material() const { return material_; };

      /// enter the evaluation of this element
      void Evaluate(LINALG::Matrix<7, 1>* state, LINALG::Matrix<2, 2>* elematrix,
          LINALG::Matrix<2, 1>* elevector, bool eval_mat, bool eval_vec);


     private:
      /// the element definitely needs to have a unique id
      int id_;

      /// each element has a specific number of nodes
      int numnode_;

      /// each element has a specific number of dofs
      int numdof_;


      /// node coordinates in initial configuration
      std::vector<double> X_;

      /// dofs of this element
      std::vector<double> dofs_;

      /// pointer to the material of this element
      Teuchos::RCP<TUTORIAL::MATERIAL::TutorialMat> material_;

    };  // class TutorialEle


  }  // namespace ELEMENTS
}  // namespace TUTORIAL

/*----------------------------------------------------------------------*/
#endif /* #ifndef TUTORIAL_ELE_H*/
