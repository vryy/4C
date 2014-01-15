/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret.H
 *
 * \brief discretisation with additional reference point vector
 *        (meshfree analysis)
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#ifndef MESHFREE_DISCRET_H
#define MESHFREE_DISCRET_H

#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_meshfree.H"

namespace LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
}

namespace DRT
{

/*==========================================================================*/
/*!
 * \brief MESHFREE: namespace of all meshfree analysis tools
 *
 * blabladiblub
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2012
 *///                                                     (public) nis Jan12
/*==========================================================================*/
namespace MESHFREE
{
/*==========================================================================*/
// forward declaration
/*==========================================================================*/
  class MeshfreeNode;
  class MeshfreeApprox;

/*--------------------------------------------------------------------------*/
/*!
 * \brief A class to manage a meshfree discretization
 *
 * This is an extended discretisation for meshfree schemes like
 * max-ent. Compared to the standard discretization, it is has additional
 * knots describing the geometry of an intregration cell (=former
 * element). Furthermore, the function MakeNear() determines which nodes have
 * an influence on a certain integration cell.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeDiscretization : public DRT::Discretization
{
public:

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *
   * \param name (in): name of this meshfree discretization
   * \param comm (in): An epetra comm object associated with this
   *                   discretization
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/

  MeshfreeDiscretization(
    const std::string name,
    Teuchos::RCP<Epetra_Comm> comm,
    const Teuchos::ParameterList & params
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Desctructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/

  virtual ~MeshfreeDiscretization();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get map associated with the distribution of the ownership of
   *        knots (Filled()==true prerequisite).
   *
   * This map includes all knots stored on this proc and also owned by this
   * proc.
   *
   * This map is non-ambiguous, meaning that it is a non-overlapping map.
   *
   * \return NULL if Filled() is false. A call to FillComplete() is a
   *         prerequisite.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const Epetra_Map* KnotRowMap() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get map associated with the distribution of knots including
   *        ghosted knots (Filled()==true prerequisite).
   *
   * This map includes all knots stored on this proc including any ghosted
   * knots.
   *
   * This map is ambiguous, meaning that it is an overlapping map.
   *
   * \return NULL if Filled() is false. A call to FillComplete() is a
   *         prerequisite.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const Epetra_Map* KnotColMap() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get global number of knots (true number of total knots without
   *        ghosting) (Filled()==true prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumGlobalKnots() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get processor local number of knots owned by this processor
   *        (Filled()==true prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumMyRowKnots() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get processor local number of knots including ghost knots
   *        (Filled()==true NOT prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumMyColKnots() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Query whether a Knot with global id gid is stored on this proc
   *        (Filled()==true NOT prerequisite)
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool HaveGlobalKnot(int gid) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the knot with global row id gid (Filled()==true NOT prerequisite)
   *
   * Returns the knot with global row id gid if knot is on this proc.
   * Will return row or column knot, ghosted or  not.
   * This is an individual call
   *
   * \return Adress of knot if knot is stored on calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual MeshfreeNode* gKnot(int gid) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the knot with local row id lid (Filled()==true prerequisite)
   *
   * Returns the knot with local row index lid.
   * Will not return any ghosted knot.
   * This is an individual call and Filled()=true is a prerequisite
   *
   * \return Adress of knot if knot is owned and stored by calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline MeshfreeNode* lRowKnot(int lid) const
  {
  #ifdef DEBUG
    if (!Filled())
      dserror("DRT::MESHFREE::MeshfreeDiscretization::lRowKnot: Filled() != true");
  #endif
    return knotrowptr_[lid];
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the knot with local column id lid (Filled()==true prerequisite)
   *
   * Returns the knot with local column index lid.
   * Will return any knot stored on this proc.
   * This is an individual call and Filled()=true is a prerequisite
   *
   * \return Adress of knot if knot is stored by calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline MeshfreeNode* lColKnot(int lid) const
  {
    #ifdef DEBUG
      if (!Filled())
        dserror("DRT::MESHFREE::MeshfreeDiscretization::lColKnot: Filled() != true");
    #endif
    return knotcolptr_[lid];
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a knot to the discretization  (Filled()==true NOT prerequisite)
   *
   * The discretization takes ownership of the added knot.
   * Note that if a knot with the same Id() exists, it will be
   * deleted and replaced by the new one;
   * note furthermore that in this method Reset() is called only on the processor where
   * a knot has been added actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \note Sets Filled()=false
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AddKnot(Teuchos::RCP<MeshfreeNode> knot);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete an knot from the discretization (Filled()==true NOT prerequisite)
   *
   * Delete an knot from the discretization. Knot can either be ghosted or not.
   * If knot on calling processor is not found nothing is done and false is returned.
   * Note that this is not a fatal error and no message will be posted;
   * note furthermore that in this method Reset() is called only on the processor where
   * an knot has been deleted actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \return Returns true upon successful deletion or knot, returns false if knot
   *         was not found on calling proc.
   *
   * \note Sets Filled()=false and calls Reset() upon discretization.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool DeleteKnot(Teuchos::RCP<MeshfreeNode> knot);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete an knot with global id gid from the discretization
   *        (Filled()==true NOT prerequisite)
   *
   * Delete an knot from the discretization. Knot can either be ghosted or not.
   * If knot on calling processor is not found nothing is done and false is returned.
   * Note that this is not a fatal error and no message will be posted;
   * note furthermore that in this method Reset() is called only on the processor where
   * an knot has been deleted actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \return Returns true upon successful deletion or knot, returns false if knot
   *         was not found on calling proc.
   *
   * \note Sets Filled()=false and calls Reset() upon discretization.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool DeleteKnot(const int gid);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack local nodes (row map) into buffer.
   *
   * Call Pack on all local (row map) nodes and put the results into
   * a common vector. This is used to output the discretization.
   *
   * \note Filled()=true is a prerequisite
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  Teuchos::RCP<std::vector<char> > PackMyKnots() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack nodal buffer and create local nodes.
   *
   * Interprets the argument as packed nodes and unpacks them on the
   * local processor. Takes the ownership of the unpacked nodes.
   *
   * \param e (in): buffer of packed nodes
   *
   * \note Sets Filled()=false
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  void UnPackMyKnots(Teuchos::RCP<std::vector<char> > k);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate meshfree Dirichlet boundary conditions for
   *        non-interpolatory basis functions with weak Kronecker-delta and
   *        partion-of-unity property
   *
   * \param params (in): List of parameters
   * \param systemvector (out): Vector holding prescribed Dirichlet values
   * \param systemvectord (out): Vector holding 1st time derivative of prescribed Dirichlet values
   * \param systemvectordd (out): Vector holding 2nd time derivative prescribed Dirichlet values
   * \param toggle (out): Vector containing 1.0 for each Dirichlet dof and 0 for everything else
   * \param dbcmapextractor (out): Map extractor containing maps for the DOFs subjected to
   *                               Dirichlet boundary conditions and the remaining/free DOFs
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual void EvaluateDirichlet(Teuchos::ParameterList& params,
                                 Teuchos::RCP<Epetra_Vector> systemvector,
                                 Teuchos::RCP<Epetra_Vector> systemvectord,
                                 Teuchos::RCP<Epetra_Vector> systemvectordd,
                                 Teuchos::RCP<Epetra_Vector> toggle,
                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief assignment node<>{knot,cell} by neighbourhood search
   *
   *
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AssignNodesToKnotsAndCells();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Complete construction of a meshfree discretization
   *        (Filled()==true NOT prerequisite)
   *
   * After adding or deleting elements, nodes, or knots or redistributing them
   * in parallel, or adding/deleting boundary conditions, this method has to
   * be called to (re)construct pointer topologies.<br>
   * It builds in this order:<br>
   * - calls base class FillCompete(false,false,false)
   * - row map of knots
   * - column map of knots
   * - pointers from elements to knots
   * - pointers from knots to elements
   * - assigns nodes to knots and cells
   * - calls AssignDegreesOfFreedom() if desired
   * - calls InitializeElements() if desired
   * - calls BoundaryConditionsGeometry() if desried
   *
   * \param assigndegreesoffreedom (in) : if true, resets existing dofsets and
   *                                      performs assigning of degrees of
   *                                      freedoms to knots and elements.
   * \param initelements (in) : if true, build element register classes and
   *                            call Initialize() on each type of finite
   *                            element present
   * \param doboundaryconditions (in) : if true, build geometry of boundary
   *                                    conditions present.
   *
   * \note In order to receive a fully functional meshfree discretization,
   *       this method must be called with all parameters set to true (the
   *       default). The parameters though can be used to turn off specific
   *       tasks to allow for more flexibility in the construction of a
   *       meshfree discretization, were it is known that this method will be
   *       called more than once.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int FillComplete(bool assigndegreesoffreedom=true,
                           bool initelements = true,
                           bool doboundaryconditions = true);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get flag indicating whether discretization has been filled and assigned
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline bool Filled() const {return (assigned_ and filled_);}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Returns pointer to the meshfree approximation of the solution function
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual Teuchos::RCP<MeshfreeApprox> GetMeshfreeSolutionApprox(){return solutionapprox_;}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Returns pointer to the meshfree approximation of the weighting function
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual Teuchos::RCP<MeshfreeApprox> GetMeshfreeWeightingApprox(){return weightingapprox_;}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this discretization to os (Filled()==true NOT prerequisite)
   *       (std::ostream << also supported)
   *
   * \note This is a collective call
  *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Print(std::ostream& os) const;

/*--------------------------------------------------------------------------*/

private:

  /*------------------------------------------------------------------------*/
  //! Empty! No ctor(name,comm) wanted                    (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeDiscretization(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

  /*------------------------------------------------------------------------*/
  //! Empty! No = operator wanted.                        (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeDiscretization operator = (const MeshfreeDiscretization& old);

  /*------------------------------------------------------------------------*/
  //! Empty! No copy constructor wanted.                  (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeDiscretization(const MeshfreeDiscretization& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief commentation to be written...
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool InitAssignSingleNode(const double* const & xn,
                                    double const & range,
                                    const int & numeles,
                                    DRT::Element** ele,
                                    std::set<int>& elegid,
                                    std::set<int>& knotgid,
                                    int& init_knot
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief recursive assignment node<>{knot,cell}
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AssignSingleNode(const double* const & xn,         //!< (const): position of node
                                const int & nodeid,               //!< (const): global id of node
                                DRT::Node* & node,                //!< (in): pointer to node (element ptr is updated)
                                double const & range,
                                const int & numele,               //!< (const): number of elements for this recursion
                                DRT::Element** eles,              //!< (in): pointer to element pointers for this recursion
                                std::set<int>& elegid,                 //!< (in/out): global id of elements already busy
                                                                  //!        dealing or even dealt with (continuously
                                                                  //!        updated and common to all recursion states)
                                std::set<int>& knotgid,                //!< (in/out): global id of knots already busy
                                                                  //!        dealing or even dealt with (continuously
                                                                  //!        updated and common to all recursion states)
                                const int & myrank,               //!< (const): rank of proc
                                std::map<int, std::map<int,int> > & procmap //!< (out): map over ranks assigned to map over nodes
                                                                  //!         assigned to inital knots
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief unassigns all node information in all cells
   *
   * The cell/element information in the nodes is used for initilizing the
   * reassignment. Once the reassignment is initilized, this information is
   * deleted form the nodes, as well.
   *///                                                   (public) nis Apr12
  /*------------------------------------------------------------------------*/
  virtual void Unassign();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Reset all maps and set Filled()=false (Filled()==true NOT
   *        prerequisite)
   *
   * Resets all maps and sets flags filled_ and havedof_ to false.
   *
   * \param killdofs (in): if true reset existing dofsets in discretization
   *
   * \note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Reset(bool killdofs = true);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build knotrowmap_, knotcolmap_, knotrowptr, knotcolptr,
   *        and knotghoptr (Filled()==true NOT prerequisite)
   *
   * Build the parallel layout of nodes in this
   * discretization and store it as an Epetra_Map in noderowmap_
   * noderowmap_ is unique.
   * It considers nodes owned by a proc only.
   *
   * \note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildKnotMaps();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build nodecolmap_ (Filled()==true NOT prerequisite)
   *
   * Build the parallel layout of nodes in this
   * discretization and store it as an Epetra_Map in nodecolmap_
   * nodecolmap_ is potentially but not necessarily overlapping.
   * It considers nodes owned by a proc and its ghosted nodes
   *
   * \note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  //virtual void BuildKnotColMap();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointers elements -> Nodes (Filled()==true NOT prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildElementToKnotPointers();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointers Node -> Element (Filled()==true NOT prerequisite)
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildKnotToElementPointers();

protected:

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Fill Dirichlet BC to compute nodal values from values at nodes if
   *        necessary. Needs to be done since meshfree methods are usually not
   *        interpolatory. Prerequisite: weak Kronecker-delta property
   *///                                                (protected) nis Jan14
  /*------------------------------------------------------------------------*/
  void FillDBCMatrix(
    const DRT::Condition& cond,
    const int             facedim,
    const int             numdof,
    std::set<int>&        nodegids,
    bool&                 isconst,
    std::vector<const double*>&          constvalue,
    Teuchos::RCP<LINALG::SparseMatrix>   dbcmatrix,
    const Teuchos::RCP<const Epetra_Map> dbcdofmap
    );

  // internal data
  bool assigned_;                                  //!< status for assigned of nodes to cells/knots
  INPAR::MESHFREE::NodeAssignType nodeassigntype_; //!< way the neighbourhood search for nodes<>elements is managed

  Teuchos::RCP<Epetra_Map> knotrowmap_;            //!< unique distribution of nodal ownerships
  Teuchos::RCP<Epetra_Map> knotcolmap_;            //!< distribution of knots including ghost knots
  std::vector<MeshfreeNode*> knotrowptr_;          //!< vector of pointers to row knots for faster access
  std::vector<MeshfreeNode*> knotcolptr_;          //!< vector of pointers to column knots for faster access
  std::vector<MeshfreeNode*> knotghoptr_;          //!< vector of pointers to ghosted knots for faster access
  std::map<int,Teuchos::RCP<MeshfreeNode> > knot_; //!< map from knot Gid to knot pointers

  Teuchos::RCP<MeshfreeApprox> solutionapprox_;    //!< defines the meshfree approximation for the (basis) solution function
  Teuchos::RCP<MeshfreeApprox> weightingapprox_;   //!< defines the meshfree approximation for the (basis) weighting function

  Teuchos::RCP<LINALG::Solver> dbcsolver_;         //!< linear solver to strongly enforce non-constant DirichletBC

}; // class MeshfreeDiscretization

/*--------------------------------------------------------------------------*/

} // namespace MESHFREE
} // namespace DRT

/// << operator
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeDiscretization& dis);

#endif  // #ifndef MESHFREE_DISCRET_H
