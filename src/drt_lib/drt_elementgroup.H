#ifdef CCADISCRET

#ifndef DRT_ELEMENTGROUP_H
#define DRT_ELEMENTGROUP_H

#include <blitz/array.h>

#include <vector>
#include <list>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>


namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{
  class Discretization;
  class Element;

  /// Element groups for faster element evaluation
  namespace EGROUP
  {

    /// abstract element group base
    class Group
    {
    public:
      /// polymorph destruction
      virtual ~Group() {}

      void Add(DRT::Element* ele) { elements_.push_back(ele); }

      unsigned Size() const { return elements_.size(); }

      void Clear() { elements_.clear(); }

      DRT::Element* operator[](int i) { return elements_[i]; }

    private:
      std::vector<DRT::Element*> elements_;
    };


    /// group of all non-special elements
    class OtherElements : public Group
    {
    };


    /// hex8 elements aligned with axes
    class AlignedHex8 : public Group
    {
    public:

    private:

    };


    /// Evaluate element group in a specific way
    class Evaluator
    {
    public:
      virtual ~Evaluator() {}

      /// Evaluate generic elements
      virtual void Evaluate(OtherElements& elements) = 0;

      /// Evaluate axis aligned hex8 elements
      virtual void Evaluate(AlignedHex8& elements) = 0;

    };


    /// base class of "normal" assembling evaluators
    class SystemEvaluatorBase : public Evaluator
    {
    public:
      SystemEvaluatorBase(Teuchos::RCP<DRT::Discretization> dis,
                          const Teuchos::ParameterList& params,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                          Teuchos::RCP<Epetra_Vector> systemvector);

    protected:

      Teuchos::RCP<DRT::Discretization> dis_;
      const Teuchos::ParameterList& params_;
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
      Teuchos::RCP<Epetra_Vector> systemvector_;
    };


    /// manager of element groups
    /*!

      Elements are put into groups according to their geometry. The point is
      to find special elements (elements with particularly simple geometry)
      and evaluate these elements much faster by using their special
      characteristics.

      \author u.kue
      \date 03/08
     */
    class ElementGroupManager
    {
    public:

      /// Construct element group manager on given Discretization
      ElementGroupManager(const Discretization& dis);

      /// Evaluate all elements using the given Evaluator
      void Evaluate(Evaluator& evaluator);

    private:

      /// put elements of my Discretization in groups
      void GroupElements();

      /// find axis aligned hex8 elements and rotate those
      /*!
        There are six possible rotations of an axis aligned hex8
        element. However, we want to have the element local r,s,t axis be
        aligned with the global x,y,z axis. So we do an undercover rotation
        here.

        \warning This rotation can cause trouble in conjunction with
        non-standard boundary conditions or the like. The current stage needs
        rethinking.
       */
      bool FindHex8Align(DRT::Element* ele);

      /// whether given vector is parallel to x axis
      bool Xdir(blitz::Array<double,1>& e) const
      { return fabs(e(0)) > 1e-10 and fabs(e(1)) < 1e-10 and fabs(e(2)) < 1e-10; }

      /// whether given vector is parallel to y axis
      bool Ydir(blitz::Array<double,1>& e) const
      { return fabs(e(0)) < 1e-10 and fabs(e(1)) > 1e-10 and fabs(e(2)) < 1e-10; }

      /// whether given vector is parallel to z axis
      bool Zdir(blitz::Array<double,1>& e) const
      { return fabs(e(0)) < 1e-10 and fabs(e(1)) < 1e-10 and fabs(e(2)) > 1e-10; }

      /// rotate (cyclic swap) of three values
      template <class T>
      void Rotate(T& t1, T& t2, T& t3)
      {
        T t = t1;
        t1 = t2;
        t2 = t3;
        t3 = t;
      }

      /// rotate (cyclic swap) of four values
      template <class T>
      void Rotate(T& t1, T& t2, T& t3, T& t4)
      {
        T t = t1;
        t1 = t2;
        t2 = t3;
        t3 = t4;
        t4 = t;
      }

      /// my Discretization
      const Discretization& dis_;

      /// unspecified set of elements
      OtherElements others_;

      /// axis aligned hex8 elements
      AlignedHex8 aligned_;
    };
  }
}

#endif
#endif
