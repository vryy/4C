/*!----------------------------------------------------------------------
\file deal_tria_wrapper.H
\brief Computes a deal.II Triangulation from DRT::Discretization

<pre>

\level 2
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*----------------------------------------------------------------------*/

#ifndef DEAL_TRIA_WRAPPER_H
#define DEAL_TRIA_WRAPPER_H

#ifdef HAVE_DEAL_II

#include <deal.II/grid/tria.h>
#include <deal.II/base/index_set.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>

#include "../drt_lib/drt_discret_hdg.H"

namespace dealii
{
  /**
   * This class acts like the dealii::Triangulation class, but it
   * distributes the mesh across a number of different processors when using
   * MPI using the parallel distribution of DRT::Discretization. The functionality
   * considerably reduced as compared to dealii::Triangulation, in particular,
   * we do not allow the grid to be refined. This is because the parallel distribution
   * of the mesh would change and there are no data communication facilities in deal.II
   * for this cases. However, this is not expected, either, because we have the mesh
   * available in DRT::Discretization and have no means to get reasonable data back from
   * a different deal.II mesh.
   */
  template <int dim>
  class DistributedTriangulation : public Triangulation<dim, dim>
  {
   public:
    /**
     * A typedef that is used to to identify cell iterators. See the deal.II
     * documentation for details.
     *
     * The current typedef identifies cells in a triangulation. You can find
     * the exact type it refers to in the base class's own typedef, but it
     * should be TriaIterator<CellAccessor<dim,spacedim> >. The TriaIterator
     * class works like a pointer that when you dereference it yields an
     * object of type CellAccessor. CellAccessor is a class that identifies
     * properties that are specific to cells in a triangulation, but it is
     * derived (and consequently inherits) from TriaAccessor that describes
     * what you can ask of more general objects (lines, faces, as well as
     * cells) in a triangulation.
     *
     * @ingroup Iterators
     */
    typedef typename dealii::Triangulation<dim, dim>::cell_iterator cell_iterator;

    /**
     * A typedef that is used to to identify "active cell iterators". See the
     * deal.II documentation for details.
     *
     * The current typedef identifies active cells in a triangulation. You
     * can find the exact type it refers to in the base class's own typedef,
     * but it should be TriaActiveIterator<CellAccessor<dim,spacedim> >. The
     * TriaActiveIterator class works like a pointer to active objects that
     * when you dereference it yields an object of type CellAccessor.
     * CellAccessor is a class that identifies properties that are specific
     * to cells in a triangulation, but it is derived (and consequently
     * inherits) from TriaAccessor that describes what you can ask of more
     * general objects (lines, faces, as well as cells) in a triangulation.
     *
     * @ingroup Iterators
     */
    typedef typename dealii::Triangulation<dim, dim>::active_cell_iterator active_cell_iterator;

    /**
     * Constructor.
     */
    DistributedTriangulation(Teuchos::RCP<DRT::DiscretizationHDG> discret);

    /**
     * Destructor.
     */
    ~DistributedTriangulation() {}

    /**
     * Overloads the respective method in dealii::Triangulation, disabling refinement.
     */
    virtual void execute_coarsening_and_refinement()
    {
      dserror("Refinement of deal.II triangulation is not supported.");
    }

    /**
     * Prints the global number of degrees of freedom
     */
    virtual types::global_dof_index n_global_active_cells() const
    {
      return discret->NumGlobalElements();
    }

    /**
     * Keeps a pointer to the discretization in order to access its data.
     */
    Teuchos::RCP<DRT::DiscretizationHDG> discret;
  };



  template <int dim>
  void assign_dg_dofs(const FiniteElement<dim> &fe, DoFHandler<dim> &dof_handler);
}  // namespace dealii

#endif

#endif
