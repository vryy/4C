/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret.H
 *
 * \brief discretisation for meshfree analysis
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
#include "../linalg/linalg_fixedsizematrix.H"

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
  class Cell;
  class MeshfreeNode;
  class MeshfreeApprox;

/*--------------------------------------------------------------------------*/
/*!
 * \brief A class to manage a meshfree discretization
 *
 * This is an extended discretisation for meshfree schemes like
 * max-ent. Compared to the standard discretization, it is has additional
 * points describing the geometry of an intregration cell (=former
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
   * \brief Subclass for meshfree discretization with the folowwing features:
   *        - holding nodal face information for 0-D to 3-D faces
   *        - unique information shared with respective boundary conditions
   *        - managing node to node assignements needed for coverting nodal values
   *          to values at nodes and vice versa. Used for
   *          - non-trivial Dirichlet boundary conditions
   *          - post processing
   *        - managing node to cell assignments for boundary elements
   *///                                                (protected) nis Jan14
  /*------------------------------------------------------------------------*/
  class Face
  {
    /*========================================================================*/
    //! @name Enums and friends
    /*========================================================================*/

    /*------------------------------------------------------------------------*/
    // \brief only MeshfreeDiscretization is allowed to deal with Face nis Jan14
    /*------------------------------------------------------------------------*/
    friend class MeshfreeDiscretization;

  public:
    /*========================================================================*/
    //! @name Face methods
    /*========================================================================*/

    /*------------------------------------------------------------------------*/
    // \brief standard constructor to be public to resize std::vector<Face>
    /*------------------------------------------------------------------------*/
    Face();

  private:
    /*------------------------------------------------------------------------*/
    // \brief constructor with global node ids and dimension         nis Jan14
    /*------------------------------------------------------------------------*/
    Face(std::vector<int>& nodeids, int dim, DRT::MESHFREE::MeshfreeDiscretization* discret);

    /*------------------------------------------------------------------------*/
    // \brief fills node pointers and assignes nodes on this face    nis Jan14
    /*------------------------------------------------------------------------*/
    void FillComplete(const double range, const int myrank, const bool setnodefacedim=false)
    {
      if (!filled_)   BuildFaceNodePointer();
      if (!assigned_) AssignNodesToNodes(setnodefacedim);
      return;
    }

    /*------------------------------------------------------------------------*/
    // \brief checks whether FillComplete() was called               nis Jan14
    /*------------------------------------------------------------------------*/
    bool Filled() {return (assigned_ and filled_);}

    /*------------------------------------------------------------------------*/
    // \brief returns pointer to global face node ids                nis Jan14
    /*------------------------------------------------------------------------*/
    Teuchos::RCP<std::vector<int> > Nodes() const {return nodeids_;}

    /*------------------------------------------------------------------------*/
    // \brief returns Teuchos::RCP to face geometry                   nis Jan14
    /*------------------------------------------------------------------------*/
    Teuchos::RCP<std::map<int,Teuchos::RCP<DRT::Element> > > Geometry(){return geometry_;}

    /*------------------------------------------------------------------------*/
    // \brief adds geometry to face and assigns nodes to its cells   nis Jan14
    /*------------------------------------------------------------------------*/
    void AddGeometry(Teuchos::RCP<std::map<int,Teuchos::RCP<DRT::Element> > > g){ geometry_ = g; }

    /*------------------------------------------------------------------------*/
    // \brief builds face node row and column pointers               nis Jan14
    /*------------------------------------------------------------------------*/
    void BuildFaceNodePointer();

    /*------------------------------------------------------------------------*/
    // \brief assigns column face nodes to row face nodes if within range  nis Jan14
    /*------------------------------------------------------------------------*/
    void AssignNodesToNodes(const bool setnodefacedim=false);

    /*------------------------------------------------------------------------*/
    // \brief assign nodes to cells on face (for Neumann conditions) nis Jan14
    /*------------------------------------------------------------------------*/
    void AssignNodesToCells();

    /*========================================================================*/
    //! @name Face members
    /*========================================================================*/

    Teuchos::RCP<std::vector<int> > nodeids_;  //!< GIDs of all (!!!) global face nodes
    std::vector<MeshfreeNode*> noderowptr_;    //!< vector of pointers to row nodes for faster access
    std::map<int,MeshfreeNode*> nodecolmap_;   //!< map of pointers to col nodes

    Teuchos::RCP<std::map<int,Teuchos::RCP<DRT::Element> > > geometry_; //!< elements associated with this face

    DRT::MESHFREE::MeshfreeDiscretization* discret_; //!< pointer to discretization

    int    dim_;      //!< dimension of face
    bool   filled_;   //!< flag whether node pointers are filled
    bool   assigned_; //!< flag whether nodes are assigned to nodes
  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get map associated with the distribution of the ownership of
   *        points (Filled()==true prerequisite).
   *
   * This map includes all points stored on this proc and also owned by this
   * proc.
   *
   * This map is non-ambiguous, meaning that it is a non-overlapping map.
   *
   * \return NULL if Filled() is false. A call to FillComplete() is a
   *         prerequisite.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const Epetra_Map* PointRowMap() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get map associated with the distribution of points including
   *        ghosted points (Filled()==true prerequisite).
   *
   * This map includes all points stored on this proc including any ghosted
   * points.
   *
   * This map is ambiguous, meaning that it is an overlapping map.
   *
   * \return NULL if Filled() is false. A call to FillComplete() is a
   *         prerequisite.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const Epetra_Map* PointColMap() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get global number of points (true number of total points without
   *        ghosting) (Filled()==true prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumGlobalPoints() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get processor local number of points owned by this processor
   *        (Filled()==true prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumMyRowPoints() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get processor local number of points including ghost points
   *        (Filled()==true NOT prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumMyColPoints() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Query whether a Point with global id gid is stored on this proc
   *        (Filled()==true NOT prerequisite)
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool HaveGlobalPoint(int gid) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the point with global row id gid (Filled()==true NOT prerequisite)
   *
   * Returns the point with global row id gid if point is on this proc.
   * Will return row or column point, ghosted or  not.
   * This is an individual call
   *
   * \return Adress of point if point is stored on calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual DRT::Node* gPoint(int gid) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the point with local row id lid (Filled()==true prerequisite)
   *
   * Returns the point with local row index lid.
   * Will not return any ghosted point.
   * This is an individual call and Filled()=true is a prerequisite
   *
   * \return Adress of point if point is owned and stored by calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline DRT::Node* lRowPoint(int lid) const
  {
  #ifdef DEBUG
    if (!Filled())
      dserror("DRT::MESHFREE::MeshfreeDiscretization::lRowPoint: Filled() != true");
  #endif
    return pointrowptr_[lid];
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get the point with local column id lid (Filled()==true prerequisite)
   *
   * Returns the point with local column index lid.
   * Will return any point stored on this proc.
   * This is an individual call and Filled()=true is a prerequisite
   *
   * \return Adress of point if point is stored by calling proc
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline DRT::Node* lColPoint(int lid) const
  {
    #ifdef DEBUG
      if (!Filled())
        dserror("DRT::MESHFREE::MeshfreeDiscretization::lColPoint: Filled() != true");
    #endif
    return pointcolptr_[lid];
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a point to the discretization  (Filled()==true NOT prerequisite)
   *
   * The discretization takes ownership of the added point.
   * Note that if a point with the same Id() exists, it will be
   * deleted and replaced by the new one;
   * note furthermore that in this method Reset() is called only on the processor where
   * a point has been added actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \note Sets Filled()=false
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AddPoint(Teuchos::RCP<DRT::Node> point);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete an point from the discretization (Filled()==true NOT prerequisite)
   *
   * Delete an point from the discretization. Point can either be ghosted or not.
   * If point on calling processor is not found nothing is done and false is returned.
   * Note that this is not a fatal error and no message will be posted;
   * note furthermore that in this method Reset() is called only on the processor where
   * an point has been deleted actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \return Returns true upon successful deletion or point, returns false if point
   *         was not found on calling proc.
   *
   * \note Sets Filled()=false and calls Reset() upon discretization.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool DeletePoint(Teuchos::RCP<DRT::Node> point);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete an point with global id gid from the discretization
   *        (Filled()==true NOT prerequisite)
   *
   * Delete an point from the discretization. Point can either be ghosted or not.
   * If point on calling processor is not found nothing is done and false is returned.
   * Note that this is not a fatal error and no message will be posted;
   * note furthermore that in this method Reset() is called only on the processor where
   * an point has been deleted actually; however, as such a modification may affect the
   * discretization as a whole it may be required to call Reset() on each processor
   * subsequently; to do so, please, call CheckFilledGlobally() if required
   *
   * \return Returns true upon successful deletion or point, returns false if point
   *         was not found on calling proc.
   *
   * \note Sets Filled()=false and calls Reset() upon discretization.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool DeletePoint(const int gid);

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
  Teuchos::RCP<std::vector<char> > PackMyPoints() const;

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
  void UnPackMyPoints(Teuchos::RCP<std::vector<char> > k);

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
  virtual void EvaluateDirichlet(
    Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector,
    Teuchos::RCP<Epetra_Vector> systemvectord,
    Teuchos::RCP<Epetra_Vector> systemvectordd,
    Teuchos::RCP<Epetra_Vector> toggle,
    Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Adds node set topology to faces of meshfree discretisation. Basis
   *        function are only evaluated within a face/node set.
   *
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual void AddNodeSetTopology(
    std::vector<std::vector<std::vector<int> >* > nodeset
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Overloading base class function to prevent storing condition node
   *        sets twice
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual void SetCondition(const std::string& name,Teuchos::RCP<Condition> cond);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief assignment node<>{point,cell} by neighbourhood search
   *
   *
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AssignNodesToPointsAndCells();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Complete construction of a meshfree discretization
   *        (Filled()==true NOT prerequisite)
   *
   * After adding or deleting elements, nodes, or points or redistributing them
   * in parallel, or adding/deleting boundary conditions, this method has to
   * be called to (re)construct pointer topologies.<br>
   * It builds in this order:<br>
   * - calls base class FillCompete(false,false,false)
   * - row map of points
   * - column map of points
   * - pointers from elements to points
   * - pointers from points to elements
   * - assigns nodes to points and cells
   * - calls AssignDegreesOfFreedom() if desired
   * - calls InitializeElements() if desired
   * - calls BoundaryConditionsGeometry() if desried
   *
   * \param assigndegreesoffreedom (in) : if true, resets existing dofsets and
   *                                      performs assigning of degrees of
   *                                      freedoms to points and elements.
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
  virtual int FillComplete(
    bool assigndegreesoffreedom=true,
    bool initelements = true,
    bool doboundaryconditions = true
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Computes values at nodes from nodal values for non-interpolatory
   *        basis function. Needs correct node-to-node assignment which is
   *        done when constructing faces.
   *        (Filled()==true prerequisite)
   *
   * \note This should be called and implemented as 'ComputeValuesAtPoints'
   *       since for post filtering, element topology is needed which is only
   *       provided by cells and points. This would require a second dof-set
   *       based on points. However, up to now nodes lie on points by
   *       construction, which is why we circumvent this problem by computing
   *       values at nodes instead of points, using the nodal dof-set.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual void ComputeValuesAtNodes(
    const Teuchos::RCP<const Epetra_Vector>& nodalvalues,
    const Teuchos::RCP<Epetra_Vector>      & valuesatnodes ) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Returns pointer to the meshfree approximation of the solution function
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual Teuchos::RCP<MeshfreeApprox> GetSolutionApprox(){return solutionapprox_;}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Returns pointer to the meshfree approximation of the weighting function
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual Teuchos::RCP<MeshfreeApprox> GetWeightingApprox(){return weightingapprox_;}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this discretization to os (Filled()==true NOT prerequisite)
   *       (std::ostream << also supported)
   *
   * \note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Print(std::ostream& os) const;

//  virtual void Redistribute(const Epetra_Map& noderowmap,
//                            const Epetra_Map& nodecolmap,
//                            bool assigndegreesoffreedom=true,
//                            bool initelements          =true,
//                            bool doboundaryconditions  =true,
//                            bool killdofs              =true,
//                            bool killcond              =true);

/*--------------------------------------------------------------------------*/

protected:

  /*------------------------------------------------------------------------*/
  //! Empty! No ctor(name,comm) wanted                    (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeDiscretization(const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DRT::Discretization(name,comm){};

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
   * \brief Build the geometry of lines for a certain line condition.
   *        Additonal node-to-cell-assignment compared to base class call.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual bool BuildLinesinCondition(
    const std::string name,
    Teuchos::RCP<DRT::Condition> cond);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build the geometry of surfaces for a certain surface condition
   *        Additonal node-to-cell-assignment compared to base class call.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual bool BuildSurfacesinCondition(
    const std::string name,
    Teuchos::RCP<DRT::Condition> cond);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief TODO commentation to be written...
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool InitAssignSingleNode(
    const double* const xn,
    const double range,
    const int numeles,
    DRT::Element** ele,
    std::set<int>& elegid,
    std::set<int>& pointgid,
    int& init_point
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief recursive assignment node<>{point,cell}
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AssignSingleNode(
    const double* const xn,                     //!< (const): position of node
    const int nodeid,                           //!< (const): global id of node
    DRT::Node* node,                            //!< (in): pointer to node (element ptr is updated)
    const double range,
    const int numele,                           //!< (const): number of elements for this recursion
    DRT::Element** eles,                        //!< (in): pointer to element pointers for this recursion
    std::set<int>& elegid,                      //!< (in/out): global id of elements already busy
                                                //!            dealing or even dealt with (continuously
                                                //!            updated and common to all recursion states)
    std::set<int>& pointgid,                    //!< (in/out): global id of points already busy
                                                //!        dealing or even dealt with (continuously
                                                //!        updated and common to all recursion states)
    const int   myrank,                         //!< (const): rank of proc
    std::map<int, std::map<int,int> > & procmap //!< (out): map over ranks assigned to map over nodes
                                                //!         assigned to inital points
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
  virtual void Reset(bool killdofs, bool killcond);
  virtual void Reset() { this->Reset(true, true); }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointrowmap_, pointcolmap_, pointrowptr and pointcolptr,
   *        (Filled()==true NOT prerequisite)
   *
   * Build the parallel layout of nodes in this
   * discretization and store it as an Epetra_Map in noderowmap_
   * noderowmap_ is unique.
   * It considers nodes owned by a proc only.
   *
   * \note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildPointMaps();

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
  //virtual void BuildPointColMap();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointers elements -> Nodes (Filled()==true NOT prerequisite)
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildElementToPointPointers();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointers Node -> Element (Filled()==true NOT prerequisite)
   *
   *\note This is a collective call
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void BuildPointToElementPointers();

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

protected:

  /*==========================================================================*/
  //! @name internal data additional to the one in parent Discretization
  /*==========================================================================*/
  INPAR::MESHFREE::NodeAssignType nodeassigntype_; //!< way the neighbourhood search for nodes<>elements is managed
  INPAR::MESHFREE::meshfreetype meshfreetype_;     //!< type of meshfree discretisation

  Teuchos::RCP<Epetra_Map> pointrowmap_;            //!< unique distribution of nodal ownerships
  Teuchos::RCP<Epetra_Map> pointcolmap_;            //!< distribution of points including ghost points
  std::vector<DRT::Node*> pointrowptr_;             //!< vector of pointers to row points for faster access
  std::vector<DRT::Node*> pointcolptr_;             //!< vector of pointers to column points for faster access
  std::map<int,Teuchos::RCP<DRT::Node> > point_;    //!< map from point Gid to point pointers

  Teuchos::RCP<MeshfreeApprox> solutionapprox_;     //!< defines the meshfree approximation for the (basis) solution function

  Teuchos::RCP<MeshfreeApprox> weightingapprox_;    //!< defines the meshfree approximation for the (basis) weighting function

  Teuchos::RCP<LINALG::Solver> dbcsolver_;          //!< linear solver to strongly enforce non-constant DirichletBC

  std::vector<std::vector<Face> > domaintopology_;  //!< the highest dimensional face, i.e. the domain

}; // class MeshfreeDiscretization

/*--------------------------------------------------------------------------*/

} // namespace MESHFREE
} // namespace DRT

/// << operator
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeDiscretization& dis);

#endif  // #ifndef MESHFREE_DISCRET_H
