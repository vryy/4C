/*----------------------------------------------------------------------*/
/*! \file

\brief base class for all immersed algorithms

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_IMMERSED_PROBLEM_IMMERSED_BASE_HPP
#define FOUR_C_IMMERSED_PROBLEM_IMMERSED_BASE_HPP


#include "4C_config.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_cut_boundingbox.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_immersed_node.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fem_geometry_searchtree_service.hpp"
#include "4C_fluid_ele_immersed_base.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_mortar_calc_utils.hpp"
#include "4C_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}

namespace Adapter
{
  class FluidWrapper;
  class StructureWrapper;
}  // namespace Adapter


namespace Immersed
{
  class ImmersedBase
  {
   public:
    /*!
    \brief constructor

    */
    ImmersedBase();


    /*!
    \brief destructor

    */
    virtual ~ImmersedBase() = default;


    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(const Teuchos::ParameterList& params) = 0;

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup(){};


    /*!
    \brief create a volume condition on the fly

    \author rauch
    \date 05/14

    \param dis           (in) : discretization on which condition is to be constructed
    \param dvol_fenode   (in) : nodes to be conditioned
    \param condtype      (in) : type of condition
    \param condname      (in) : condition name
    \param buildgeometry (in) : build geometry for condition or not
    */
    void create_volume_condition(const Teuchos::RCP<Discret::Discretization>& dis,
        const std::vector<int> dvol_fenode, const Core::Conditions::ConditionType condtype,
        const std::string condname, bool buildgeometry);

    /*!
    \brief construct the dof row map for a given condition.

    \author rauch
    \date 02/17
    */
    void build_condition_dof_map(const Teuchos::RCP<const Discret::Discretization>& dis,
        const std::string condname, const Teuchos::RCP<const Epetra_Map>& cond_dofmap_orig,
        const int numdof, Teuchos::RCP<Epetra_Map>& cond_dofmap);


    /*!
    \brief Apply Dirichlet values to systemvector of conditioned field.

           Do not touch original dofs.

    \author rauch
    \date 02/17
    */
    void DoDirichletCond(const Teuchos::RCP<Epetra_Vector>& statevector,
        const Teuchos::RCP<const Epetra_Vector>& dirichvals,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_new,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_orig);


    /*!
    \brief Apply Dirichlet values to systemvector of conditioned field.

           Touch all dofs.

    \author rauch
    \date 02/17
    */
    void DoDirichletCond(const Teuchos::RCP<Epetra_Vector>& statevector,
        const Teuchos::RCP<const Epetra_Vector>& dirichvals,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_new);


    /*!
    \brief Apply Dirichlet boundary condition to a structure field.

    \author rauch
    \date 02/17
    */
    virtual void ApplyDirichlet(const Teuchos::RCP<Adapter::StructureWrapper>& field_wrapper,
        const Teuchos::RCP<Discret::Discretization>& dis, const std::string condname,
        Teuchos::RCP<Epetra_Map>& cond_dofrowmap, const int numdof,
        const Teuchos::RCP<const Epetra_Vector>& dirichvals);


    /*!
    \brief Apply Dirichlet boundary condition to a fluid field.

    \author rauch
    \date 02/17
    */
    void apply_dirichlet_to_fluid(const Teuchos::RCP<Adapter::FluidWrapper>& field_wrapper,
        const Teuchos::RCP<Discret::Discretization>& dis, const std::string condname,
        Teuchos::RCP<Epetra_Map>& cond_dofrowmap, const int numdof,
        const Teuchos::RCP<const Epetra_Vector>& dirichvals);


    /*!
    \brief Remove Dirichlet boundary condition from a structure field.

    \author rauch
    \date 02/17
    */
    void RemoveDirichlet(const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
        const Teuchos::RCP<Adapter::StructureWrapper>& field_wrapper);


    /*!
    \brief Remove Dirichlet boundary condition from a fluid field.

    \author rauch
    \date 02/17
    */
    void remove_dirichlet_from_fluid(const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
        const Teuchos::RCP<Adapter::FluidWrapper>& field_wrapper);


    /*!
    \brief specialized evaluate for immersed algorithm

    \author rauch
    \date 05/14
    */
    void EvaluateImmersed(Teuchos::ParameterList& params, Teuchos::RCP<Discret::Discretization> dis,
        Core::FE::AssembleStrategy* strategy, std::map<int, std::set<int>>* elementstoeval,
        Teuchos::RCP<Core::Geo::SearchTree> structsearchtree,
        std::map<int, Core::LinAlg::Matrix<3, 1>>* currpositions_struct, int action,
        bool evaluateonlyboundary = false);


    /*!
    \brief specialized evaluate for immersed algorithm without assembly

    \author rauch
    \date 05/14

    */
    void evaluate_immersed_no_assembly(Teuchos::ParameterList& params,
        Teuchos::RCP<Discret::Discretization> dis, std::map<int, std::set<int>>* elementstoeval,
        Teuchos::RCP<Core::Geo::SearchTree> structsearchtree,
        std::map<int, Core::LinAlg::Matrix<3, 1>>* currpositions_struct, int action);


    /*!
    \brief specialized evaluate for immersed algorithm allowing communication on element level

    \author rauch
    \date 05/14

    */
    void evaluate_sca_tra_with_internal_communication(Teuchos::RCP<Discret::Discretization> dis,
        const Teuchos::RCP<const Discret::Discretization> idis,
        Core::FE::AssembleStrategy* strategy, std::map<int, std::set<int>>* elementstoeval,
        Teuchos::RCP<Core::Geo::SearchTree> structsearchtree,
        std::map<int, Core::LinAlg::Matrix<3, 1>>* currpositions_struct,
        Teuchos::ParameterList& params, bool evaluateonlyboundary = false);


    /*!
    \brief specialized evaluate for immersed algorithm considering the need of communication inside
    the loop over all conditioned elements

    \author rauch
    \date 05/14

    */
    void evaluate_interpolation_condition(Teuchos::RCP<Discret::Discretization> evaldis,
        Teuchos::ParameterList& params, Core::FE::AssembleStrategy& strategy,
        const std::string& condstring, const int condid);


    /*!
    \brief search in search tree for background elements in given radius

    \author rauch
    \date 09/15

    */
    void search_potentially_covered_backgrd_elements(
        std::map<int, std::set<int>>* current_subset_tofill,
        Teuchos::RCP<Core::Geo::SearchTree> backgrd_SearchTree, const Discret::Discretization& dis,
        const std::map<int, Core::LinAlg::Matrix<3, 1>>& currentpositions,
        const Core::LinAlg::Matrix<3, 1>& point, const double radius, const int label);


    /*!
    \brief evaluate a given subset of elements

    \author rauch
    \date 09/15

    */
    void evaluate_subset_elements(Teuchos::ParameterList& params,
        Teuchos::RCP<Discret::Discretization> dis, std::map<int, std::set<int>>& elementstoeval,
        int action);


    /*!
    \brief write special output for immersed algorithms

           Writes vectors of values provided as argument to a text file with standard
           prefix and file ending which is also provided as argument to this method.
           The values are written in one row for the provided current time.

    \author rauch
    \date 03/17

    \param comm           (in) : communicator
    \param time           (in) : current time will be the first value in the row
    \param filenameending (in) : ending of file to write to
    \param valuetowrite   (in) : first set of values to write to row
    \param valuetowrite2  (in) : second set of values to wrtie to row
    \param valuetowrite3  (in) : third set of values to write to row
    */
    void WriteExtraOutput(const Epetra_Comm& comm, const double time,
        const std::string filenameending, const std::vector<double> valuetowrite,
        const std::vector<double> valuetowrite2, const std::vector<double> valuetowrite3);


    /*!
    \brief Calculate the global resulting force from an epetra vector.

    \author rauch
    \date 03/17
    */
    std::vector<double> calc_global_resultantfrom_epetra_vector(const Epetra_Comm& comm,
        const Teuchos::RCP<const Discret::Discretization>& dis,
        const Teuchos::RCP<const Epetra_Vector>& vec_epetra);

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool is_setup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool is_init() { return isinit_; };

    //! check if \ref Setup() was called
    void check_is_setup()
    {
      if (not is_setup()) FOUR_C_THROW("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void check_is_init()
    {
      if (not is_init()) FOUR_C_THROW("Init(...) was not called.");
    };

   protected:
    //! set flag true after setup or false if setup became invalid
    void set_is_setup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void set_is_init(bool trueorfalse) { isinit_ = trueorfalse; };


  };  // class ImmersedBase


  /*!
  \brief interpolate quantity from the background field to a given immersed point

  \author rauch
  \date 09/15

  */
  template <Core::FE::CellType sourcedistype, Core::FE::CellType targetdistype>
  int InterpolateToImmersedIntPoint(const Teuchos::RCP<Discret::Discretization> sourcedis,
      const Teuchos::RCP<Discret::Discretization> targetdis, Core::Elements::Element& targetele,
      const std::vector<double>& targetxi, const std::vector<double>& targetedisp, int action,
      std::vector<double>& targetdata)
  {
    //////////////////////////////////////////////////////////////////////////////////
    ///////
    ///////  SOURCE  -->  TARGET
    ///////
    ///////  maps source values to target discretization.
    ///////
    ///////  comments on parallelism:
    ///////
    ///////  every proc needs to enter this routine to make communication work.
    ///////  in the extreme scenario where a proc owns no target elements at all
    ///////  no evaluate will be called and this proc doesn't enter this routine
    ///////  which means in turn, that a Barrier() would never be reached. in
    ///////  this case the attempt to communicate would definitely fail.
    ///////  so make sure, that each proc owns at least one target element and/or
    ///////  it is also to be assured by external measures, that each proc
    ///////  reenters this routine until the last element on every proc is
    ///////  evaluated (e.g. dummy call of Evaluate(...)).
    ///////
    /////////////////////////////////////////////////////////////////////////////////

    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();

    // get the global problem
    Global::Problem* problem = Global::Problem::Instance();

    // true in case we have to deal wit ha deformable background mesh
    static int isALE =
        (problem->immersed_method_params().get<std::string>("DEFORM_BACKGROUND_MESH") == "yes");

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    Core::Elements::Element* sourceele;

    // sourceele displacements
    std::vector<double> myvalues;

    // declare vector of global current coordinates
    std::vector<double> x(globdim);

    // dimension of source discretization
    static constexpr int source_dim = Core::FE::dim<sourcedistype>;
    static constexpr int target_dim = Core::FE::dim<targetdistype>;

    // declarations and definitions for round robin loop
    const int numproc = comm.NumProc();
    const int myrank = comm.MyPID();                        // myrank
    const int torank = (myrank + 1) % numproc;              // sends to
    const int fromrank = (myrank + numproc - 1) % numproc;  // recieves from

    Core::Communication::Exporter exporter(comm);

    // get current global coordinates of the given point xi of the target dis
    Mortar::UTILS::LocalToCurrentGlobal<targetdistype>(
        targetele, globdim, targetxi.data(), targetedisp, x.data());
    std::vector<double> xvec(globdim);
    xvec[0] = x[0];
    xvec[1] = x[1];
    xvec[2] = x[2];

    // vector to fill by sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> vectofill =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(targetdata.size()));
    (*vectofill)(0) = -1234.0;

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> xi_dense =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(source_dim));
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> targetxi_dense =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(target_dim));

    Teuchos::RCP<Core::LinAlg::SerialDenseVector> normal_at_targetpoint =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(3));
    std::vector<double> normal_vec(3);

    Core::LinAlg::SerialDenseMatrix dummy1;
    Core::LinAlg::SerialDenseMatrix dummy2;
    Core::LinAlg::SerialDenseVector dummy3;

    Teuchos::ParameterList params;

    params.set<std::string>("action", "calc_cur_normal_at_point");
    for (int i = 0; i < target_dim; ++i) (*targetxi_dense)(i) = targetxi[i];

    Core::Elements::Element::LocationArray targetla(1);
    targetele.LocationVector(*targetdis, targetla, false);

    targetele.Evaluate(params, *targetdis, targetla, dummy1, dummy2, *normal_at_targetpoint,
        *targetxi_dense, dummy3);
    normal_at_targetpoint->scale(1.0 / (Core::LinAlg::Norm2(*normal_at_targetpoint)));
    normal_vec[0] = (*normal_at_targetpoint)(0);
    normal_vec[1] = (*normal_at_targetpoint)(1);
    normal_vec[2] = (*normal_at_targetpoint)(2);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->length();

    // get possible elements being intersected by immersed structure
    Core::Conditions::Condition* searchbox = sourcedis->GetCondition("ImmersedSearchbox");
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& searchboxgeom = searchbox->Geometry();
    int mysearchboxgeomsize = searchboxgeom.size();

    // round robin loop
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all background elements
      //
      // determine background element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator curr;
      Core::LinAlg::Matrix<source_dim, 1> xi(true);
      if (!matched and mysearchboxgeomsize > 0)
      {
        // every proc that has searchboxgeom elements
        for (curr = searchboxgeom.begin(); curr != searchboxgeom.end(); ++curr)
        {
          // only proc who owns background searchbox ele is supposed to interpolate and is set
          // BoundaryIsImmersed(true)
          if (curr->second->Owner() == myrank)
          {
            bool converged = false;
            double residual = -1234.0;

            if (isALE)
            {
              Teuchos::RCP<const Epetra_Vector> state = sourcedis->GetState("dispnp");

              Core::Elements::Element::LocationArray la(1);
              curr->second->LocationVector(*sourcedis, la, false);
              // extract local values of the global vectors
              myvalues.resize(la[0].lm_.size());
              Core::FE::ExtractMyValues(*state, myvalues, la[0].lm_);
              double sourceeledisp[24];
              for (int node = 0; node < 8; ++node)
                for (int dof = 0; dof < 3; ++dof)
                  sourceeledisp[node * 3 + dof] = myvalues[node * 4 + dof];

              // node 1  and node 7 coords of current source element (diagonal points)
              const auto& X1 = curr->second->Nodes()[1]->X();
              double x1[3];
              x1[0] = X1[0] + sourceeledisp[1 * 3 + 0];
              x1[1] = X1[1] + sourceeledisp[1 * 3 + 1];
              x1[2] = X1[2] + sourceeledisp[1 * 3 + 2];
              const auto& X7 = curr->second->Nodes()[7]->X();
              double diagonal =
                  sqrt(pow(X1[0] - X7[0], 2) + pow(X1[1] - X7[1], 2) + pow(X1[2] - X7[2], 2));

              // calc distance of current target point to arbitrary node (e.g. node 1) of curr
              // source element
              double distance =
                  sqrt(pow(x1[0] - xvec[0], 2) + pow(x1[1] - xvec[1], 2) + pow(x1[2] - xvec[2], 2));

              // get parameter space coords xi in source element of global point xvec of target
              // element NOTE: if the target point xvec is very far away from the source element
              // curr
              //       it is unnecessary to jump into this functon and invoke a newton iteration.
              // Therfore: only call GlobalToCurrentLocal if distance is smaller than
              // fac*characteristic element length
              if (distance < 1.5 * diagonal)
              {
                Mortar::UTILS::GlobalToCurrentLocal<sourcedistype>(
                    *(curr->second), sourceeledisp, xvec.data(), &xi(0), converged, residual);
                if (converged == false)
                {
                  //                  std::cout<<"Warning!! GlobalToCurrentLocal() in
                  //                  InterpolateToImmersedIntPoint() did not converge! res
                  //                  ="<<std::setprecision(12)<<residual<<"\n " "Source Element ID
                  //                  = "<<curr->second->Id()<< "  d="<<distance<<" < "<<"
                  //                  diagonal="<<diagonal<<"\n" "pos target point = ["<<xvec[0]<<"
                  //                  "<<xvec[1]<<" "<<xvec[2]<<"] "<< "xi in source
                  //                  ele="<<xi<<std::endl;
                  xi(0) = 2.0;
                  xi(1) = 2.0;
                  xi(2) = 2.0;
                }
              }
              else
              {
                xi(0) = 2.0;
                xi(1) = 2.0;
                xi(2) = 2.0;
              }
            }
            else
            {
              // node 1  and node 7 coords of current source element (diagonal points)
              const auto& X1 = curr->second->Nodes()[1]->X();
              const auto& X7 = curr->second->Nodes()[7]->X();
              double diagonal =
                  sqrt(pow(X1[0] - X7[0], 2) + pow(X1[1] - X7[1], 2) + pow(X1[2] - X7[2], 2));

              // calc distance of current target point to arbitrary node (e.g. node 1) of curr
              // source element
              double distance =
                  sqrt(pow(X1[0] - xvec[0], 2) + pow(X1[1] - xvec[1], 2) + pow(X1[2] - xvec[2], 2));

              if (distance < 1.5 * diagonal)
              {
                // get parameter space coords xi in source element of global point xvec of target
                // element
                Mortar::UTILS::GlobalToLocal<sourcedistype>(
                    *(curr->second), xvec.data(), &xi(0), converged);
                if (converged == false)
                {
                  xi(0) = 2.0;
                  xi(1) = 2.0;
                  xi(2) = 2.0;
                }
              }
              else
              {
                xi(0) = 2.0;
                xi(1) = 2.0;
                xi(2) = 2.0;
              }
            }

            /////////////////////////////////////////////////
            ////
            ////     MATCH
            ////
            /////////////////////////////////////////////////
            double tol = 1e-13;
            sourceele = curr->second.getRawPtr();
            bool validsource = false;
            // given point lies in element curr
            if (abs(xi(0)) < (1.0 - tol) and abs(xi(1)) < (1.0 - tol) and abs(xi(2)) < (1.0 - tol))
            {
              validsource = true;

              (*xi_dense)(0) = xi(0);
              (*xi_dense)(1) = xi(1);
              (*xi_dense)(2) = xi(2);
            }
            // does targetpoint lie on any element edge?
            // targetpoint lies within tolerance away from any element edge
            // only the sourceele located outward from the boundary should interpolate
            else if ((abs((abs(xi(0)) - 1.0)) <= tol)     // at xi_0 surface
                     or (abs((abs(xi(1)) - 1.0)) <= tol)  // at xi_1 surface
                     or (abs((abs(xi(2)) - 1.0)) <= tol)  // at xi_2 surface
            )
            {
              double scalarproduct = 0.0;
              Teuchos::RCP<Core::LinAlg::SerialDenseVector> vector;

              // get nodal coords of source ele
              vector = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(3));
              std::vector<double> distances(sourceele->num_node());

              // loop over source element nodes and calc distance from targetpoint to those nodes
              for (int sourcenode = 0; sourcenode < sourceele->num_node(); sourcenode++)
              {
                if (isALE)
                {
                  const auto& X = sourceele->Nodes()[sourcenode]->X();

                  distances[sourcenode] =
                      sqrt(pow(xvec[0] - (X[0] + myvalues[sourcenode * 4 + 0]), 2) +
                           pow(xvec[1] - (X[1] + myvalues[sourcenode * 4 + 1]), 2) +
                           pow(xvec[2] - (X[2] + myvalues[sourcenode * 4 + 2]), 2));
                }
                else
                {
                  const auto& X = sourceele->Nodes()[sourcenode]->X();

                  distances[sourcenode] = sqrt(
                      pow(xvec[0] - X[0], 2) + pow(xvec[1] - X[1], 2) + pow(xvec[2] - X[2], 2));
                }
              }

              // get max. distance and corresponding index of src node
              int maxdistanceindex = -1;
              double maxvalue = 0.0;
              for (int i = 0; i < (int)distances.size(); ++i)
              {
                if (distances[i] > maxvalue)
                {
                  maxvalue = distances[i];
                  maxdistanceindex = i;
                }
              }

              double x[3];
              if (isALE)
              {
                const auto& X = sourceele->Nodes()[maxdistanceindex]->X();
                x[0] = X[0] + myvalues[maxdistanceindex * 4 + 0];
                x[1] = X[1] + myvalues[maxdistanceindex * 4 + 1];
                x[2] = X[2] + myvalues[maxdistanceindex * 4 + 2];
              }
              else
              {
                x[0] = sourceele->Nodes()[maxdistanceindex]->X()[0];
                x[1] = sourceele->Nodes()[maxdistanceindex]->X()[1];
                x[2] = sourceele->Nodes()[maxdistanceindex]->X()[2];
              }

              // build vector directed from targetele point (xvec) to sourceele node (X)
              // xvec should not lie at same ele surface as sourcenode
              (*vector)(0) = x[0] - xvec[0];
              (*vector)(1) = x[1] - xvec[1];
              (*vector)(2) = x[2] - xvec[2];
              vector->scale(Core::LinAlg::Norm2(*vector));

              // build scalar product between normal and vector
              scalarproduct = (*vector)(0) * normal_vec[0] + (*vector)(1) * normal_vec[1] +
                              (*vector)(2) * normal_vec[2];

              if (scalarproduct > tol)
                validsource = true;
              else
                validsource = false;

              std::cout << "WARNING !! Immersed point is lying on element edge! validsource="
                        << validsource << "for eleid=" << sourceele->Id() << std::endl;

#ifdef FOUR_C_ENABLE_ASSERTIONS
              if (abs(scalarproduct) < 1e-13 and validsource)
              {
                std::cout << "scalarproduct = " << scalarproduct << std::endl;
                FOUR_C_THROW(
                    "normal vec and vec from point to node perpendicular.\n"
                    "no robust determination of outward lying src ele possible");
              }
              if (abs(scalarproduct) < -1e-12 and validsource)
              {
                std::cout << "scalarproduct = " << scalarproduct << std::endl;
                FOUR_C_THROW(
                    "normal vec and vec from point to node have no sharp angle.\n"
                    "no robust determination of outward laying src ele possible");
              }
#endif

              // count++;
            }  // if on element edge


            // evaluate and perform action handed in to this function
            if (action == 81311)
            {  // experimental pre cmbe
              params.set<std::string>("action", "interpolate_porosity_to_given_point");
            }
            else
              params.set<int>("action", action);
            params.set<int>("Physical Type", Inpar::FLUID::poro_p1);
            if (validsource)
            {
              //              std::cout<<"\n ---------- \n PROC "<<comm.MyPID()<<" ; sourcele :
              //              "<<sourceele->Id()<<std::endl; std::cout<<"POINT: ["<<xvec[0]<<"
              //              "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl; std::cout<<"normal:
              //              ["<<normal_vec[0]<<" "<<normal_vec[1]<<"
              //              "<<normal_vec[2]<<"]"<<std::endl; std::cout<<"vector:
              //              "<<*vector<<std::endl; std::cout<<"scalarproduct:
              //              "<<scalarproduct<<std::endl;

              // fill locationarray
              Core::Elements::Element::LocationArray sourcela(1);
              sourceele->LocationVector(*sourcedis, sourcela, false);
              sourceele->Evaluate(
                  params, *sourcedis, sourcela, dummy1, dummy2, *vectofill, *xi_dense, dummy3);
              matched = 1;

              // if element interpolates to structural intpoint set as "boundary is immersed"
              Teuchos::RCP<Discret::ELEMENTS::FluidImmersedBase> immersedele =
                  Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::FluidImmersedBase>(curr->second);
              // only possible and reasonable if current element is fluid element
              if (immersedele != Teuchos::null)
              {
                // std::cout<<"ID="<<immersedele->Id()<<"
                // scalarprod="<<std::setprecision(10)<<scalarproduct<<std::endl;
                immersedele->set_boundary_is_immersed(2);
                immersedele->ConstructElementRCP(
                    problem->immersed_method_params().get<int>("NUM_GP_FLUID_BOUND"));
              }

              // DEBUG output
              // std::cout<<"Tgt. Ele: "<<targetele.Id()<<" lies in src. ele:
              // "<<curr->second->Id()<<" at  pos. ["<<std::setprecision(5)<<xvec[0]<<"
              // "<<xvec[1]<<"
              // "<<xvec[2]<<"]"<<std::endl;

              break;  // break loop over all source eles
            }
            else
              matched = 0;

          }  // if myrank == owner of searchbox ele and proc has searchbox ele
        }    // loop over all seachbox elements
      }      // enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if (numproc > 1)
      {
        // ---- pack data for sending -----
        {
          Core::Communication::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->length());

            // point coordinate
            for (int dim = 0; dim < globdim; ++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if (matched == 1)
            {
              for (int dim = 0; dim < (int)vectofill->length(); ++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }

          // normal vector
          for (int dim = 0; dim < globdim; ++dim)
          {
            data.AddtoPack(normal_vec[dim]);
          }

          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.i_send(myrank, torank, sdata.data(), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          FOUR_C_THROW("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched = Core::Communication::ParObject::ExtractInt(position, rdata);
        owner = Core::Communication::ParObject::ExtractInt(position, rdata);
        datalength = Core::Communication::ParObject::ExtractInt(position, rdata);

        for (int i = 0; i < globdim; ++i)
          xvec[i] = Core::Communication::ParObject::ExtractDouble(position, rdata);

        if (matched == 1)
        {
          for (int dim = 0; dim < datalength; ++dim)
            (*vectofill)(dim) = Core::Communication::ParObject::ExtractDouble(position, rdata);
        }

        for (int i = 0; i < globdim; ++i)
          normal_vec[i] = Core::Communication::ParObject::ExtractDouble(position, rdata);

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      }  // only if num procs > 1

      if (mysearchboxgeomsize > 0)
      {
        if (curr == (searchboxgeom.end()--) and irobin == (numproc - 1) and matched == 0)
        {
          std::cout << "target position = [" << xvec[0] << " " << xvec[1] << " " << xvec[2] << "]"
                    << std::endl;
          FOUR_C_THROW("could not match given point on any proc. Element(0)=%f", (*vectofill)(0));
        }
      }

    }  // end for irobin

#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (owner != comm.MyPID())
      FOUR_C_THROW("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff

    for (int dim = 0; dim < datalength; ++dim) targetdata[dim] = (*vectofill)(dim);


    return err;
  };


  /*!
  \brief interpolate quantity from the immersed field to a given background point

  \author rauch
  \date 09/15

  */
  template <Core::FE::CellType sourcedistype, Core::FE::CellType targetdistype>
  int InterpolateToBackgrdPoint(std::map<int, std::set<int>>& curr_subset_of_structdis,
      const Teuchos::RCP<Discret::Discretization> sourcedis,
      const Teuchos::RCP<Discret::Discretization> targetdis, Core::Elements::Element& targetele,
      const std::vector<double>& targetxi, const std::vector<double>& targetedisp,
      const std::string action, std::vector<double>& targetdata, bool& match, bool vel_calculation,
      bool doCommunication = true)
  {
    //////////////////////////////////////////////////////////////////////////////////
    ///////
    /////// FAST version of InterpolateToBackgrdIntPoint
    ///////
    /////// employing search tree for immersed discretization.
    ///////
    ///////
    /////////////////////////////////////////////////////////////////////////////////

    // error flag
    int err = 0;
    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();
    // get the global problem
    Global::Problem* problem = Global::Problem::Instance();
    // dimension of global problem
    const int globdim = problem->NDim();
    // pointer to background element
    Core::Elements::Element* sourceele;
    // declare vector of global current coordinates
    std::vector<double> x(globdim);
    // dimension of source discretization
    static constexpr int source_dim = Core::FE::dim<sourcedistype>;
    // declarations and definitions for round robin loop
    int numproc = comm.NumProc();
    const int myrank = comm.MyPID();                        // myrank
    const int torank = (myrank + 1) % numproc;              // sends to
    const int fromrank = (myrank + numproc - 1) % numproc;  // recieves from

    if (numproc == 1) doCommunication = false;
    if (doCommunication == false) numproc = 1;

    Teuchos::RCP<const Epetra_Vector> dispnp = sourcedis->GetState("displacement");
    Core::Elements::Element::LocationArray la(1);

    // get current global coordinates of the given point xi of the target dis
    Mortar::UTILS::LocalToCurrentGlobal<targetdistype>(
        targetele, globdim, targetxi.data(), targetedisp, x.data());
    std::vector<double> xvec(globdim);
    xvec[0] = x[0];
    xvec[1] = x[1];
    xvec[2] = x[2];

    // vector to be filled in sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> vectofill =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(4));
    (*vectofill)(0) = -1234.0;
    if (match)  // was already matched, add no contribution
    {
      (*vectofill)(0) = 0.0;
      (*vectofill)(1) = 0.0;
      (*vectofill)(2) = 0.0;
      (*vectofill)(3) = 0.0;
    }

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> xi_dense =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(source_dim));

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    if (match == true) matched = 1;

    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->length();

    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all source elements
      //
      // determine source element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      Core::LinAlg::Matrix<source_dim, 1> xi(true);
      if (matched == 0)
      {
        for (std::map<int, std::set<int>>::const_iterator closele =
                 curr_subset_of_structdis.begin();
             closele != curr_subset_of_structdis.end(); closele++)
        {
          for (std::set<int>::const_iterator eleIter = (closele->second).begin();
               eleIter != (closele->second).end(); eleIter++)
          {
            bool converged = false;
            double residual = -1234.0;
            sourceele = sourcedis->gElement(*eleIter);

            // get parameter space coords xi in source element
            sourceele->LocationVector(*sourcedis, la, false);
            std::vector<double> mysourcedispnp(la[0].lm_.size());
            Core::FE::ExtractMyValues(*dispnp, mysourcedispnp, la[0].lm_);

            // construct bounding box around current source element
            Teuchos::RCP<Core::Geo::Cut::BoundingBox> bbside =
                Teuchos::rcp(Core::Geo::Cut::BoundingBox::Create());
            Core::LinAlg::Matrix<3, 1> nodalpos;
            for (int i = 0; i < (int)(sourceele->num_node()); ++i)
            {
              nodalpos(0, 0) = (sourceele->Nodes())[i]->X()[0] + mysourcedispnp[i * 3 + 0];
              nodalpos(1, 0) = (sourceele->Nodes())[i]->X()[1] + mysourcedispnp[i * 3 + 1];
              nodalpos(2, 0) = (sourceele->Nodes())[i]->X()[2] + mysourcedispnp[i * 3 + 2];
              bbside->AddPoint(nodalpos);
            }

            // check if given point lies within bounding box
            bool within = bbside->Within(1.0, xvec.data());

            // only try to match given target point and sourceele if within bounding box
            if (within)
            {
              Mortar::UTILS::GlobalToCurrentLocal<sourcedistype>(
                  *sourceele, mysourcedispnp.data(), xvec.data(), &xi(0), converged, residual);

              if (converged == false)
              {
                std::cout << "!! WARNING !! GlobalToCurrentLocal did not not converge (res="
                          << residual << ") though background point " << xvec[0] << " " << xvec[1]
                          << " " << xvec[2]
                          << "\n"
                             " is lying within bounding box of immersed element "
                          << sourceele->Id() << "." << std::endl;
                bbside->Print();
                // we assume not in bounding box. point is lying outside sourcele.
                xi(0) = 2.0;
                xi(1) = 2.0;
                xi(2) = 2.0;
              }
            }
            else
            {
              // not in bounding box. point is lying outside sourcele.
              xi(0) = 2.0;
              xi(1) = 2.0;
              xi(2) = 2.0;
            }

            /////////////////////////////////////////////////
            ////
            ////     MATCH
            ////
            /////////////////////////////////////////////////
            double tol = 1e-12;
            if (abs(xi(0)) <= (1.0 + tol) and abs(xi(1)) <= (1.0 + tol) and
                abs(xi(2)) <= (1.0 + tol))
            {  // given point lies in element curr
              matched = 1;
              match = true;

              // parameter list
              Teuchos::ParameterList params;
              params.set<std::string>("action", action);
              // fill locationarray
              Core::Elements::Element::LocationArray la(1);
              sourceele->LocationVector(*sourcedis, la, false);

              (*xi_dense)(0) = xi(0);
              (*xi_dense)(1) = xi(1);
              (*xi_dense)(2) = xi(2);

              Core::LinAlg::SerialDenseMatrix dummy1;
              Core::LinAlg::SerialDenseMatrix dummy2;
              Core::LinAlg::SerialDenseVector dummy3;

              // check whether velocity or divergence should be interpolated
              int vel_or_div = 1;
              if (not vel_calculation) vel_or_div = 0;
              params.set<int>("calculate_velocity", vel_or_div);

              if (action != "none")
                // evaluate and perform action handed in to this function
                sourceele->Evaluate(
                    params, *sourcedis, la, dummy1, dummy2, *vectofill, *xi_dense, dummy3);
              else
              {
                for (int dim = 0; dim < datalength; ++dim) (*vectofill)(dim) = 0.0;
              }

              break;  // break loop over all source elements
            }         // if match

            // given targetpoint point lies not underneath source discretization at all
            if (closele == (curr_subset_of_structdis.end()--) and
                eleIter == ((closele->second).end()--) and matched == 0 and irobin == (numproc - 1))
            {
              for (int dim = 0; dim < datalength; ++dim) (*vectofill)(dim) = -12345.0;
            }
          }  // loop over all seachbox elements
        }    // loop over closeele set
      }      // enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if (doCommunication)
      {
        Core::Communication::Exporter exporter(comm);

        // ---- pack data for sending -----
        {
          Core::Communication::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->length());

            // point coordinate
            for (int dim = 0; dim < globdim; ++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if (matched == 1)
            {
              for (int dim = 0; dim < (int)vectofill->length(); ++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }
          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.i_send(myrank, torank, sdata.data(), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          FOUR_C_THROW("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched = Core::Communication::ParObject::ExtractInt(position, rdata);
        owner = Core::Communication::ParObject::ExtractInt(position, rdata);
        datalength = Core::Communication::ParObject::ExtractInt(position, rdata);

        for (int i = 0; i < globdim; ++i)
          xvec[i] = Core::Communication::ParObject::ExtractDouble(position, rdata);

        if (matched == 1)
        {
          for (int dim = 0; dim < datalength; ++dim)
            (*vectofill)(dim) = Core::Communication::ParObject::ExtractDouble(position, rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      }  // only if numproc > 1

    }  // end for irobin

#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (owner != comm.MyPID())
      FOUR_C_THROW("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff
    for (int dim = 0; dim < datalength; ++dim) targetdata[dim] = (*vectofill)(dim);


    return err;
  };


  /*!
  \brief Find Closest Structure Surface Point

   Find the closest point on the structural surface for every fluid node of
   a cut fluid element that is not covered  by the structure. This structure
   point has to lie inside the element of the searching fluid node to
   enable later weighting via shape functions. Store velocity and local
   position of closest structure surface point.

  \author rauch
  \date 07/2015

   */
  template <Core::FE::CellType structdistype, Core::FE::CellType fluiddistype>
  int FindClosestStructureSurfacePoint(
      std::map<int, std::set<int>>& curr_subset_of_structdis,  //!< Input
      const Teuchos::RCP<Discret::Discretization> structdis,   //!< Input
      const Teuchos::RCP<Discret::Discretization> fluiddis,    //!< Input
      Core::Elements::Element& fluidele,                       //!< Input
      const std::vector<double>& fluidxi,                      //!< Input
      const std::vector<double>& fluideledisp,                 //!< Input
      const std::string action,                                //!< Input
      std::vector<double>& target_vel_data,                    //!< Output
      std::vector<double>& target_pos_data,                    //!< Output
      bool& match,                                             //!< Output
      bool doCommunication = true                              //!< Input
  )
  {
    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = structdis->Comm();

    // get the global problem
    Global::Problem* problem = Global::Problem::Instance();

    // true in case we have to deal wit ha deformable background mesh
    static int isALE =
        (problem->immersed_method_params().get<std::string>("DEFORM_BACKGROUND_MESH") == "yes");

    // dimension of global problem
    const int globdim = problem->NDim();

    // declare vector of global current coordinates of fluid node
    std::vector<double> x_fluid_node(globdim);

    // dimension of structure and fluid discretization
    static constexpr int struct_dim = Core::FE::dim<structdistype>;
    static constexpr int fluid_dim = Core::FE::dim<fluiddistype>;

    // get integration rule of structural surface element
    const Core::FE::IntPointsAndWeights<2> intpoints(Core::FE::GaussRule2D::quad_36point);

    // initialize structural shape functions matrix
    Core::LinAlg::Matrix<4, 1> shapefcts(true);

    // declarations and definitions for round robin loop
    int numproc = comm.NumProc();
    const int myrank = comm.MyPID();                        // myrank
    const int torank = (myrank + 1) % numproc;              // sends to
    const int fromrank = (myrank + numproc - 1) % numproc;  // recieves from

    if (numproc == 1) doCommunication = false;

    if (doCommunication == false) numproc = 1;

    Core::Nodes::Node** fluidnode = fluidele.Nodes();
    double char_fld_ele_length = sqrt(pow(fluidnode[1]->X()[0] - fluidnode[7]->X()[0], 2) +
                                      pow(fluidnode[1]->X()[1] - fluidnode[7]->X()[1], 2) +
                                      pow(fluidnode[1]->X()[2] - fluidnode[7]->X()[2], 2));

    // get current displacements and velocities of structure discretization
    Teuchos::RCP<const Epetra_Vector> dispnp = structdis->GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> velnp = structdis->GetState("velocity");
    Core::Elements::Element::LocationArray la(structdis->NumDofSets());

    // get current global coordinates of the given fluid node fluid_xi
    Mortar::UTILS::LocalToCurrentGlobal<fluiddistype>(
        fluidele, globdim, fluidxi.data(), fluideledisp, x_fluid_node.data());
    // get as vector
    std::vector<double> fluid_node_glob_coord(globdim);
    for (int idim = 0; idim < 3; ++idim) fluid_node_glob_coord[idim] = x_fluid_node[idim];

    // initialize vectors to fill
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> velnp_at_struct_point =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(3));
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> xi_pos_of_struct_point =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(3));
    (*velnp_at_struct_point)(0) = -1234.0;
    (*xi_pos_of_struct_point)(0) = -1234.0;

    // given point already matched?
    int matched = 0;
    if (match == true)
    {
      // do nothing if node is already matched
      matched = 1;
      match = false;
    }

    // store current smallest distance to check each point
    double stored_dist = 100000.0;
    // store vel and pos for current closest structure point
    std::vector<double> pos_to_store(globdim);
    std::vector<double> vel_to_store(globdim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////

    // owner of given point
    int owner = myrank;
    // length of velnp_at_struct_point
    int datalength = (int)velnp_at_struct_point->length();

    // loop over all procs
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      // local coordinates of structrual surface point in parameter space of the fluid element
      Core::LinAlg::Matrix<fluid_dim, 1> fluid_xi(true);

      // loop over points that are not yet matched
      if (matched == 0)
      {
        // loop over closeele set
        for (std::map<int, std::set<int>>::const_iterator closele =
                 curr_subset_of_structdis.begin();
             closele != curr_subset_of_structdis.end(); closele++)
        {
          // loop over all searchbox elements
          for (std::set<int>::const_iterator eleIter = (closele->second).begin();
               eleIter != (closele->second).end(); eleIter++)
          {
            // only surface elements (with immersed coupling condition) are relevant to find closest
            // structure point
            std::vector<Teuchos::RCP<Core::Elements::Element>> surface_eles =
                structdis->gElement(*eleIter)->Surfaces();

            // loop over surface element, find the elements with IMMERSEDCoupling condition
            for (std::vector<Teuchos::RCP<Core::Elements::Element>>::iterator surfIter =
                     surface_eles.begin();
                 surfIter != surface_eles.end(); ++surfIter)
            {
              // pointer to current surface element
              Core::Elements::Element* structele = surfIter->getRawPtr();

              // pointer to nodes of current surface element
              Core::Nodes::Node** NodesPtr = surfIter->getRawPtr()->Nodes();

              int numfsinodes = 0;

              // loop over all nodes of current surface element
              for (int surfnode = 0; surfnode < structele->num_node(); ++surfnode)
              {
                Core::Nodes::Node* checkNode = NodesPtr[surfnode];
                // check whether a IMMERSEDCoupling condition is active on this node
                if (checkNode->GetCondition("IMMERSEDCoupling") != nullptr) ++numfsinodes;
              }

              // only surface elements with IMMERSEDCoupling condition are relevant (all nodes of
              // this element have fsi condition)
              if (numfsinodes == structele->num_node())
              {
                // get displacements and velocities of structure dis
                structele->LocationVector(*structdis, la, false);
                std::vector<double> mydispnp(la[0].lm_.size());
                Core::FE::ExtractMyValues(*dispnp, mydispnp, la[0].lm_);
                std::vector<double> myvelnp(la[0].lm_.size());
                Core::FE::ExtractMyValues(*velnp, myvelnp, la[0].lm_);

                // 1.) check if closest point is a node
                // loop over all nodes of structural surface element
                for (int inode = 0; inode < structele->num_node(); ++inode)
                {
                  // get current global position of structure node
                  std::vector<double> struct_node_postion(3);
                  for (int idim = 0; idim < 3; ++idim)
                    struct_node_postion[idim] =
                        (structele->Nodes())[inode]->X()[idim] + mydispnp[inode * 3 + idim];

                  // distance between fluid node and given structure node
                  double dist = sqrt(pow(fluid_node_glob_coord[0] - struct_node_postion[0], 2) +
                                     pow(fluid_node_glob_coord[1] - struct_node_postion[1], 2) +
                                     pow(fluid_node_glob_coord[2] - struct_node_postion[2], 2));

                  // only store point information, if distance is smallest so far
                  if (dist < stored_dist)
                  {
                    // update distance, velocity and position
                    stored_dist = dist;
                    for (int idim = 0; idim < 3; ++idim)
                    {
                      vel_to_store[idim] = myvelnp[inode * 3 + idim];
                      pos_to_store[idim] = struct_node_postion[idim];
                    }
                  }
                }  // end loop over nodes

                // 2.) check if closest point is an integration point of the structural surface
                // element loop over all int points of structural surface element
                for (int igp = 0; igp < intpoints.IP().nquad; ++igp)
                {
                  // get local coordinates of gp
                  std::vector<double> struct_xsi(struct_dim);
                  struct_xsi[0] = intpoints.IP().qxg[igp][0];
                  struct_xsi[1] = intpoints.IP().qxg[igp][1];

                  // get local coordinates of gp as matrix
                  Core::LinAlg::Matrix<struct_dim, 1> struct_xi(true);
                  for (int i = 0; i < struct_dim; ++i) struct_xi(i) = struct_xsi[i];

                  // declare vector of global current coordinates
                  std::vector<double> structsurf_gp_glob_coord(globdim);

                  // get current global position of the given int point xi on the structural surface
                  Mortar::UTILS::LocalToCurrentGlobal<structdistype>(*(structele), globdim,
                      struct_xsi.data(), mydispnp, structsurf_gp_glob_coord.data());

                  // distance between fluid node and given structure gp
                  double dist =
                      sqrt(pow(fluid_node_glob_coord[0] - structsurf_gp_glob_coord[0], 2) +
                           pow(fluid_node_glob_coord[1] - structsurf_gp_glob_coord[1], 2) +
                           pow(fluid_node_glob_coord[2] - structsurf_gp_glob_coord[2], 2));

                  // only store point information, if distance is smallest so far
                  if (dist < stored_dist)
                  {
                    // update distance
                    stored_dist = dist;

                    // evaluate shape functions of structure dis at gp
                    Core::FE::shape_function<Core::FE::CellType::quad4>(struct_xi, shapefcts);
                    Core::LinAlg::Matrix<4, 3> myvelocitynp(true);
                    for (int node = 0; node < 4; ++node)
                      for (int dim = 0; dim < 3; ++dim)
                        myvelocitynp(node, dim) = myvelnp[node * 3 + dim];

                    // calculate velocity in gp
                    Core::LinAlg::Matrix<3, 1> vel_at_gp(true);
                    vel_at_gp.MultiplyTN(myvelocitynp, shapefcts);

                    // update velocity and position
                    for (int idim = 0; idim < 3; ++idim)
                    {
                      vel_to_store[idim] = vel_at_gp(idim, 0);
                      pos_to_store[idim] = structsurf_gp_glob_coord[idim];
                    }
                  }
                }  // end loop over int points
              }    // end if fsi coupling
            }      // end loop over all surface elements
          }        // loop over all seachbox elements
        }          // loop over closeele set

        // global position of closest point (that is now found) as vector
        std::vector<double> xvec(globdim);
        for (int idim = 0; idim < 3; ++idim) xvec[idim] = pos_to_store[idim];

        bool converged = false;
        double residual = -1234.0;

        if (stored_dist < 1.5 * char_fld_ele_length)
        {
          // get local coordinates of closest structural surface point in fluid element parameter
          // space
          if (isALE)
            Mortar::UTILS::GlobalToCurrentLocal<fluiddistype>(
                fluidele, fluideledisp.data(), xvec.data(), &fluid_xi(0), converged, residual);
          else
            Mortar::UTILS::GlobalToLocal<fluiddistype>(
                fluidele, xvec.data(), &fluid_xi(0), converged);

          if (converged == false)
          {
            std::cout << "WARNING!! GlobalToCurrentLocal() did not converge in "
                         "FindClosestStructureSurfacePoint( \n"
                         "xi="
                      << fluid_xi << std::endl;
            fluid_xi(0) = 2.0;
            fluid_xi(1) = 2.0;
            fluid_xi(2) = 2.0;
          }
        }
        else
        {
          fluid_xi(0) = 2.0;
          fluid_xi(1) = 2.0;
          fluid_xi(2) = 2.0;
        }

        // check if closest point lies in current fluid element
        double tol = 1e-13;
        if (abs(fluid_xi(0)) < (1.0 + tol) and abs(fluid_xi(1)) < (1.0 + tol) and
            abs(fluid_xi(2)) < (1.0 + tol))
        {  // given point lies current fluid element
          matched = 1;
          match = true;
          // give back velocity and local coordinates of closest point
          for (int dim = 0; dim < 3; ++dim)
          {
            (*velnp_at_struct_point)(dim) = vel_to_store[dim];
            (*xi_pos_of_struct_point)(dim) = fluid_xi(dim);
          }

        }  // if match
      }    // enter loop only if node not yet matched

      if (doCommunication)
      {
        Core::Communication::Exporter exporter(comm);
        // ---- pack data for sending -----
        {
          Core::Communication::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)velnp_at_struct_point->length());

            // point coordinate
            for (int dim = 0; dim < globdim; ++dim)
            {
              data.AddtoPack(fluid_node_glob_coord[dim]);
            }

            if (matched == 1)
            {
              for (int dim = 0; dim < (int)velnp_at_struct_point->length(); ++dim)
              {
                data.AddtoPack((*velnp_at_struct_point)(dim));
              }
            }
          }
          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.i_send(myrank, torank, sdata.data(), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          FOUR_C_THROW("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched = Core::Communication::ParObject::ExtractInt(position, rdata);
        owner = Core::Communication::ParObject::ExtractInt(position, rdata);
        datalength = Core::Communication::ParObject::ExtractInt(position, rdata);

        for (int i = 0; i < globdim; ++i)
          fluid_node_glob_coord[i] = Core::Communication::ParObject::ExtractDouble(position, rdata);

        if (matched == 1)
        {
          for (int dim = 0; dim < datalength; ++dim)
            (*velnp_at_struct_point)(dim) =
                Core::Communication::ParObject::ExtractDouble(position, rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      }  // only if numproc > 1

    }  // loop over all procs

#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (owner != comm.MyPID())
      FOUR_C_THROW("after round robin every proc should have received his original data");
#endif

    // now every proc should store the requested quantities in velnp_at_struct_point for his fluidxi
    // time to return that stuff
    for (int dim = 0; dim < datalength; ++dim)
    {
      target_vel_data[dim] = (*velnp_at_struct_point)(dim);
      target_pos_data[dim] = (*xi_pos_of_struct_point)(dim);
    }

    return err;
  };

}  // namespace Immersed

FOUR_C_NAMESPACE_CLOSE

#endif
