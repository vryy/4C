/*---------------------------------------------------------------------*/
/*! \file

\brief Utils methods concerning the discretization


\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_UTILS_DISCRET_H_
#define SRC_DRT_LIB_DRT_UTILS_DISCRET_H_

#include <Teuchos_RCP.hpp>
// include due to enumerator...
#include "drt_condition.H"

class Epetra_Vector;
namespace LINALG
{
  class MapExtractor;
}  // namespace LINALG
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace DRT
{
  class AssembleStrategy;
  class DiscretizationInterface;
  class DiscretizationFaces;
  namespace UTILS
  {
    class Dbc;

    /** \brief Evaluate the elements of the given discretization and fill the
     *         system matrix and vector
     *
     *  This evaluate routine supports the evaluation of a subset of all column
     *  elements inside the given discretization. If the \c col_ele_map pointer
     *  is not set or set to \c NULL, this routine generates almost no overhead
     *  and is equivalent to the more familiar implementation in DRT::Discretization.
     *
     *  \param discret      (in)  : discretization containing the considered elements
     *  \param eparams      (in)  : element parameter list
     *  \param systemmatrix (out) : system-matrix which is supposed to be filled
     *  \param systemvector (out) : system-vector which is supposed to be filled
     *  \param col_ele_map  (in)  : column element map, which can be a subset of the
     *                              discretization column map ( optional )
     *
     *  \author hiermeier \date 12/16 */
    void Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& systemvector, const Epetra_Map* col_ele_map = NULL);

    /** \brief Evaluate the elements of the given discretization and fill the
     *         system matrices and vectors
     *
     *  This evaluate routine supports the evaluation of a subset of all column
     *  elements inside the given discretization. If the \c col_ele_map pointer
     *  is not set or set to \c NULL, this routine generates almost no overhead
     *  and is equivalent to the more familiar implementation in DRT::Discretization.
     *
     *  \param discret      (in)  : discretization containing the considered elements
     *  \param eparams      (in)  : element parameter list
     *  \param systemmatrix (out) : system-matrix vector which is supposed to be filled
     *  \param systemvector (out) : system-vector vector which is supposed to be filled
     *  \param col_ele_map  (in)  : column element map, which can be a subset of the
     *                              discretization column map ( optional )
     *
     *  \author hiermeier \date 12/16 */
    void Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
        std::vector<Teuchos::RCP<LINALG::SparseOperator>>& systemmatrix,
        std::vector<Teuchos::RCP<Epetra_Vector>>& systemvector,
        const Epetra_Map* col_ele_map = NULL);

    /** \brief Evaluate the elements of the given discretization and fill the
     *         system matrices and vectors
     *
     *  This evaluate routine supports the evaluation of a subset of all column
     *  elements inside the given discretization. If the \c col_ele_map pointer
     *  is not set or set to \c NULL, this routine generates almost no overhead
     *  and is equivalent to the more familiar implementation in DRT::Discretization.
     *
     *  \param discret      (in)  : discretization containing the considered elements
     *  \param eparams      (in)  : element parameter list
     *  \param strategy     (out) : assemble strategy containing all the system vectors
     *                              and matrices
     *  \param col_ele_map  (in)  : column element map, which can be a subset of the
     *                              discretization column map ( optional )
     *
     *  \author hiermeier \date 12/16 */
    void Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
        DRT::AssembleStrategy& strategy, const Epetra_Map* col_ele_map = NULL);

    /** \brief Evaluate Dirichlet boundary conditions
     *
     *  non-member functions to call the dbc public routines
     *
     *  \author hiermeier \date 10/16 */
    void EvaluateDirichlet(const DRT::DiscretizationInterface& discret,
        const Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_Vector>& systemvector,
        const Teuchos::RCP<Epetra_Vector>& systemvectord,
        const Teuchos::RCP<Epetra_Vector>& systemvectordd,
        const Teuchos::RCP<Epetra_Vector>& toggle,
        const Teuchos::RCP<LINALG::MapExtractor>& dbcmapextractor);

    /** \brief Evaluate Dirichlet boundary conditions
     *
     *  Call this variant, if you need no new dbc map extractor.
     *  See the corresponding called function for more detailed information.
     *
     *  \author hiermeier \date 10/16 */
    inline void EvaluateDirichlet(const DRT::DiscretizationInterface& discret,
        const Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_Vector>& systemvector,
        const Teuchos::RCP<Epetra_Vector>& systemvectord,
        const Teuchos::RCP<Epetra_Vector>& systemvectordd,
        const Teuchos::RCP<Epetra_Vector>& toggle)
    {
      EvaluateDirichlet(
          discret, params, systemvector, systemvectord, systemvectordd, toggle, Teuchos::null);
    }

    /*!
    \brief Evaluate a specified initial field (scalar or vector field)

    Loop all intial field conditions attached to the discretization and
    evaluate them if their names match the user-provided string fieldstring.
    Information on which local DOFs ids are addressed by the condition
    MUST be pre-defined and is represented by the locids vector.
    As an example, if we provide an initial velocity for a 3D structural
    dynamics simulation, locids must contain the local DOF ids {0,1,2}.
    Another example would be prescribing an initial pressure in a 3D
    fluid dynamics simulation, where locids would have to contain only
    the local pressure DOF id, namely {3}.

    */
    void EvaluateInitialField(const DRT::DiscretizationInterface& discret,
        const std::string& fieldstring, Teuchos::RCP<Epetra_Vector> fieldvector,
        const std::vector<int> locids);


    /*!
    \brief Evaluate a specified initial field (scalar or vector field)

    This is the actual evaluation method.

    */
    void DoInitialField(const DRT::DiscretizationInterface& discret, DRT::Condition& cond,
        Teuchos::RCP<Epetra_Vector> fieldvector, const std::vector<int> locids);

    /** \brief Build a Dbc object
     *
     *  The Dbc object is build in dependency of the given discretization.
     *
     *  \author hiermeier \date 10/16 */
    Teuchos::RCP<const Dbc> BuildDbc(const DRT::DiscretizationInterface* discret_ptr);

    /** \brief Default Dirchilet boundary condition evaluation class
     *
     *  \author hiermeier \date 10/16 */
    class Dbc
    {
     protected:
      enum DbcSet
      {
        set_row = 0,  ///< access the dbc row GID set
        set_col = 1   ///< access the dbc column GID set
      };

     public:
      /** \brief constructor
       *
       *  Intentionally left blank! */
      Dbc(){};

      /// destructor
      virtual ~Dbc(){};

      /** \brief Extract parameters and setup some temporal variables, before the actual
       *  evaluation process can start
       *
       *  \author hiermeier \date 10/16 */
      void operator()(const DRT::DiscretizationInterface& discret,
          const Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_Vector>& systemvector,
          const Teuchos::RCP<Epetra_Vector>& systemvectord,
          const Teuchos::RCP<Epetra_Vector>& systemvectordd,
          const Teuchos::RCP<Epetra_Vector>& toggle,
          const Teuchos::RCP<LINALG::MapExtractor>& dbcmapextractor) const;

     protected:
      /// create the toggle vector based on the given systemvector maps
      Teuchos::RCP<Epetra_Vector> CreateToggleVector(const Teuchos::RCP<Epetra_Vector> toggle_input,
          const Teuchos::RCP<Epetra_Vector>* systemvectors) const;

      /** \brief Evaluate Dirichlet boundary conditions
       *
       *  Loop all Dirichlet conditions attached to the discretization and evaluate
       *  them. This method considers all conditions in condition_ with the names
       *  "PointDirichlet", "LineDirichlet", "SurfaceDirichlet" and "VolumeDirichlet".
       *  It takes a current time from the parameter list params named "total time"
       *  and evaluates the appropiate time curves at that time for each
       *  Dirichlet condition separately. If "total time" is not included
       *  in the parameters, no time curves are used.
       *
       *  \note Opposed to the other 'Evaluate' method does this one NOT assembly
       *        but OVERWRITE values in the output vector systemvector. For this
       *        reason, dirichlet boundary conditions are evaluated in the
       *        following order: First "VolumeDirichlet", then "SurfaceDirichlet",
       *        then "LineDirichlet" and finally "PointDirichlet". This way, the
       *        lower entity dirichlet BCs override the higher ones and a point
       *        Dirichlet BCs has priority over other dirichlet BCs in the input
       *        file.
       *
       *  Parameters recognized by this method:
       *  \code
       *  params.set("total time",acttime); // current total time
       *  \endcode
       *
       *  \param discret          (in): discretization corresponding to the input
       *                                system vectors
       *  \param params           (in): List of parameters
       *  \param systemvector    (out): Vector holding prescribed Dirichlet values
       *  \param systemvectord   (out): Vector holding 1st time derivative of
       *                                prescribed Dirichlet values
       *  \param systemvectordd  (out): Vector holding 2nd time derivative prescribed
       *                                Dirichlet values
       *  \param toggle          (out): Vector containing 1.0 for each Dirichlet dof
       *                                and 0 for everything else
       *  \param dbcmapextractor (out): Map extractor containing maps for the DOFs
       *                                subjected to Dirichlet boundary conditions
       *                                and the remaining/free DOFs
       *
       *  \author mwgee (original) \date 02/08
       *  \author rauch \date 06/16 */
      virtual void Evaluate(const DRT::DiscretizationInterface& discret, const double& time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, Epetra_Vector& toggle,
          Teuchos::RCP<std::set<int>>* dbcgids) const;

      /** \brief loop through Dirichlet conditions and evaluate them
       *
       *  Note that this method does not sum up but 'sets' values in systemvector.
       *  For this reason, Dirichlet BCs are evaluated hierarchical meaning
       *  in this order:
       *                 VolumeDirichlet
       *                 SurfaceDirichlet
       *                 LineDirichlet
       *                 PointDirichlet
       *  This way, lower entities override higher ones which is
       *  equivalent to inheritance of dirichlet BCs as done in the old
       *  ccarat discretization with design          (mgee 1/07)
       *
       *  Lower entities MUST NOT set dof values in systemvector before
       *  we know if higher entities also prescribe/release Dirichlet BCs
       *  for the same dofs!
       *
       *  Therefore, we first have to assess the full hierarchy and
       *  set the toggle vector for a dof to 1 if an entity prescribes a
       *  dirichlet BC and we have to set it to 0 again if a higher entity
       *  does NOT prescribe a dirichlet BC for the same dof. This is done
       *  in ReadDirichletCondition(...). We do this for each type of entity,
       *  starting with volume DBCs.
       *
       *  Only then we call DoDirichletCondition(...) for each type of entity,
       *  starting with volume DBCs.
       *
       *  This way, it is guaranteed, that the highest entity defined in
       *  the input file determines if the systemvector for the corresponding
       *  dofs is actually touched, or not, irrespective of the dirichlet BC
       *  definition of a lower entity.
       *
       *  \author rauch \date 06/16 */
      void ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
          const std::vector<Teuchos::RCP<DRT::Condition>>& conds, Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;

      /// loop over the conditions and read the given type
      void ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
          const std::vector<Teuchos::RCP<DRT::Condition>>& conds, Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids,
          const enum DRT::Condition::ConditionType& type) const;

      /** \brief Determine dofs subject to Dirichlet condition from input file
       *
       *  \param discret  (in)  :  discretization corresponding to the input
       *                           system vectors
       *  \param cond     (in)  :  The condition object
       *  \param toggle   (out) :  Its i-th compononent is set 1 if it has a
       *                           DBC, otherwise this component remains untouched
       *  \param dbcgids  (out) :  Map containing DOFs subjected to Dirichlet
       *                           boundary conditions (row and optional column)
       *
       *  \remark If you want to be sure which Dirichlet values are set, look at
       *          the highest entity for a certain node in your input file.
       *
       *  The corresponding condition, e.g.:
       *  ---------DESIGN LINE DIRICH CONDITIONS
       *  DSURF  1
       *  // example_line
       *  E 1 - NUMDOF 6 ONOFF 0 1 1 VAL 0.0 1.0 1.0 FUNCT 1 0 0
       *
       *  tells you which dofs are actually conditioned. In the example given
       *  above, the highest entity which contains a certain node is a LINE.
       *  The first dof has an ONOFF-toggle of ZERO. This means, that the nodes
       *  in LINE 'E 1' definitely DO NOT HAVE a dirichlet BC on their first dof.
       *  No matter what is defined in a condition line of lower priority. The
       *  corresponding entries in the systemvectors remain untouched.
       *
       *  \author rauch \date 06/16 */
      virtual void ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
          const DRT::Condition& cond, Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;

      /** \brief Assignment of the values to the system vectors.
       *
       *  (1) Assign VolumeDirichlet DBC GIDs
       *  (2) Assign SurfaceDirichlet DBC GIDs
       *  (3) Assign LineDirichlet DBC GIDs
       *  (4) Assign PointDirichlet DBC GIDs
       *
       *  \author hiermeier \date 10/16 */
      void DoDirichletCondition(const DRT::DiscretizationInterface& discret,
          const std::vector<Teuchos::RCP<DRT::Condition>>& conds, const double& time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;

      /// loop over the conditions and assign the given type
      void DoDirichletCondition(const DRT::DiscretizationInterface& discret,
          const std::vector<Teuchos::RCP<DRT::Condition>>& conds, const double& time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids,
          const enum DRT::Condition::ConditionType& type) const;

      /** \brief Apply the Dirichlet values to the system vectors
       *
       *  \param discret         (in): discretization corresponding to the input
       *                               system vectors
       *  \param cond            (in): The condition object
       *  \param time            (in): Evaluation time
       *  \param systemvector   (out): Vector to apply DBCs to (eg displ. in
       *                               structure, vel. in fluids)
       *  \param systemvectord  (out): First time derivative of DBCs
       *  \param systemvectordd (out): Second time derivative of DBCs
       *  \param toggle          (in): Its i-th compononent is set 1 if it has
       *                               a DBC, otherwise this component remains
       *                               untouched
       *
       *  \remark If you want to be sure which Dirichlet values are set, look at
       *          the highest entity for a certain node in your input file.
       *
       *  The corresponding condition, e.g.:
       *  ---------DESIGN LINE DIRICH CONDITIONS
       *  DSURF  1
       *  // example_line
       *  E 1 - NUMDOF 6 ONOFF 0 1 1 VAL 0.0 1.0 1.0 FUNCT 1 0 0
       *
       *  tells you which dofs are actually conditioned. In the example given
       *  above, the highest entity which contains a certain node is a LINE.
       *  The first dof has an ONOFF-toggle of ZERO. This means, that the
       *  nodes in LINE 'E 1' definitely DO NOT HAVE a dirichlet BC on their
       *  first dof. No matter what is defined in a condition line of lower
       *  priority. The corresponding entries in the systemvectors remain
       *  untouched.
       *
       *  \version rauch 06/2016
       *  Shifted and rearranged parts of the former implementation to
       *  ReadDirichletCondition(...). This fixes inconsistency in hierarchy
       *  handling. Now, the Dirichlet conditions are first read by
       *  ReadDirichletCondition(...). Then, they are applied by
       *  DoDirichletCondition(...).
       *
       *  \author mwgee (original) \date 02/08
       *  \author rauch (modified) \date 06/2016 */
      virtual void DoDirichletCondition(const DRT::DiscretizationInterface& discret,
          const DRT::Condition& cond, const double& time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const;

      /** \brief Create a Dbc map extractor, if desired
       *
       *  \author hiermeier \date 10/16 */
      void BuildDbcMapExtractor(const DRT::DiscretizationInterface& discret,
          const Teuchos::RCP<const std::set<int>>& dbcrowgids,
          const Teuchos::RCP<LINALG::MapExtractor>& dbcmapextractor) const;

    };  // class Dbc
  }     // namespace UTILS
}  // namespace DRT


#endif /* SRC_DRT_LIB_DRT_UTILS_DISCRET_H_ */
