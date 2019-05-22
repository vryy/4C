/*!----------------------------------------------------------------------

\brief immersed specific background discretization

\level 1

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#ifndef DISCRET_IMMERSED_H
#define DISCRET_IMMERSED_H

#include "../drt_lib/drt_discret.H"


namespace DRT
{
  /*!
  \brief immersed specific background discretization

  Specialization of class \ref DRT::Discretization for background discr.
  in immersed simulations which require the so called 'selective assembly'.

  <h3>Purpose</h3>
  This class is introduced to be used as discretization for the solid
  phase in immersed cell migration. Since the artificial solid phase,
  i.e., the part of the solid phase discr. which is covered by the cell,
  is physically occupied by the cell only, we cancel the mechanical
  effect of the artificial solid phase onto the surrounding physical
  structure by means of a 'selective assembly'.
  Therefore, this class overrides evaluate routines to allow us not to
  assemble the force contribution to specific nodes.


  \author rauch
  */
  class DiscretizationImmersed : public Discretization
  {
   public:
    /*!
    \brief Standard Constructor

    \param name (in): name of this discretization
    \param comm (in): An epetra comm object associated with this discretization
    */
    DiscretizationImmersed(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

    /*!
    \brief Destructor

    */
    virtual ~DiscretizationImmersed() { return; };

    /// Call elements to evaluate
    virtual void Evaluate(Teuchos::ParameterList& params, DRT::AssembleStrategy& strategy);

    /// Evaluate a specific condition using assemble strategy
    virtual void EvaluateCondition(Teuchos::ParameterList& params, DRT::AssembleStrategy& strategy,
        const std::string& condstring, const int condid = -1);

   private:
    /*!
    \brief Return the location vector of given element for selective assembly

    The method computes degrees of freedom the provided element addresses.
    This method works similar to \ref DRT::Element::LocationVector but with
    the difference, that certain degrees of freedom are omitted because we do
    not want to assemble the element distribution to these dofs. This is what
    selective assembly means.

    \note the LocationVector method usually belongs to the element. For the
    special purposes addressed in this class, this functionality is shifted
    to this specialized discretization class to combine the functionality in
    one concise implementation.

    \param ele (in)      : ptr to element for which we want to get the dofs locations
    \param dis (in)      : the discretization this element belongs to
    \param la (out)      : location data for all dofsets of the discretization

    */
    void LocationVectorSelectiveAssembly(
        const Element* ele, const Discretization& dis, DRT::Element::LocationArray& la) const;

  };  // class DiscretizationImmersed
}  // namespace DRT

#endif  // #ifndef DISCRET_IMMERSED_H
