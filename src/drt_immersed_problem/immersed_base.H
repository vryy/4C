/*!----------------------------------------------------------------------
\file immersed_base.H

\brief base class for all immersed algorithms

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef IMMERSED_BASE_H
#define IMMERSED_BASE_H


#include "../drt_fluid_ele/fluid_ele_immersed_base.H"
#include "../drt_immersed_problem/immersed_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_calc_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_cut/cut_boundingbox.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

#include "../drt_inpar/inpar_fluid.H"

// experimental pre cmbe
#include "../drt_scatra_ele/scatra_ele_action.H"


namespace DRT{
  class Discretization;
}

namespace ADAPTER{
  class FluidWrapper;
  class StructureWrapper;
}


namespace IMMERSED{

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
  virtual ~ImmersedBase(){};


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

  \param dis         (in) : discretization on which condition is to be constructed
  \param dvol_fenode (in) : nodes to be conditioned
  \param condtype    (in) : type of condition
  \param condname    (in) : condition name

  */
  void CreateVolumeCondition(const Teuchos::RCP<DRT::Discretization>& dis,
                             const std::vector<int> dvol_fenode,
                             const DRT::Condition::ConditionType condtype,
                             const std::string condname);


  /*!
  \brief construct the dof row map for a given condition.

  \author rauch
  \date 02/17
  */
  void BuildConditionDofMap(
      const Teuchos::RCP<const DRT::Discretization>&  dis,
      const std::string                     condname,
      const Teuchos::RCP<const Epetra_Map>& cond_dofmap_orig,
      Teuchos::RCP<Epetra_Map> &            cond_dofmap);


  /*!
  \brief Apply Dirichlet values to systemvector of conditioned field.

  \author rauch
  \date 02/17
  */
  void DoDirichletCond(
      const Teuchos::RCP<Epetra_Vector>&       statevector,
      const Teuchos::RCP<const Epetra_Vector>& dirichvals,
      const Teuchos::RCP<const Epetra_Map>&    dbcmap);


  /*!
  \brief Apply Dirichlet boundary condition to a structure field.

  \author rauch
  \date 02/17
  */
  virtual void ApplyDirichlet(
      const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper,
      const Teuchos::RCP<DRT::Discretization>&  dis,
      const std::string                         condname,
      Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
      const Teuchos::RCP<const Epetra_Vector>&  dirichvals);


  /*!
  \brief Apply Dirichlet boundary condition to a fluid field.

  \author rauch
  \date 02/17
  */
  void ApplyDirichletToFluid(
      const Teuchos::RCP< ::ADAPTER::FluidWrapper>& field_wrapper,
      const Teuchos::RCP<DRT::Discretization>&  dis,
      const std::string                         condname,
      Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
      const Teuchos::RCP<const Epetra_Vector>&  dirichvals);


  /*!
  \brief Remove Dirichlet boundary condition from a structure field.

  \author rauch
  \date 02/17
  */
  void RemoveDirichlet(
      const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
      const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper);


  /*!
  \brief Remove Dirichlet boundary condition from a fluid field.

  \author rauch
  \date 02/17
  */
  void RemoveDirichletFromFluid(
      const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
      const Teuchos::RCP< ::ADAPTER::FluidWrapper>& field_wrapper);


  /*!
  \brief specialized evaluate for immersed algorithm

  \author rauch
  \date 05/14
  */
  void EvaluateImmersed(Teuchos::ParameterList& params,
                        Teuchos::RCP<DRT::Discretization> dis,
                        DRT::AssembleStrategy* strategy,
                        std::map<int,std::set<int> >* elementstoeval,
                        Teuchos::RCP<GEO::SearchTree> structsearchtree,
                        std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                        int action,
                        bool evaluateonlyboundary=false);


  /*!
  \brief specialized evaluate for immersed algorithm without assembly

  \author rauch
  \date 05/14

  */
  void EvaluateImmersedNoAssembly(Teuchos::ParameterList& params,
                                  Teuchos::RCP<DRT::Discretization> dis,
                                  std::map<int,std::set<int> >* elementstoeval,
                                  Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                  int action
                                   );


  /*!
  \brief specialized evaluate for immersed algorithm allowing communication on element level

  \author rauch
  \date 05/14

  */
  void EvaluateScaTraWithInternalCommunication(Teuchos::RCP<DRT::Discretization> dis,
                                               const Teuchos::RCP<const DRT::Discretization> idis,
                                               DRT::AssembleStrategy* strategy,
                                               std::map<int,std::set<int> >* elementstoeval,
                                               Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                               std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                               Teuchos::ParameterList& params,
                                               bool evaluateonlyboundary=false);


  /*!
  \brief specialized evaluate for immersed algorithm considering the need of communication inside the loop over all conditioned elements

  \author rauch
  \date 05/14

  */
  void EvaluateInterpolationCondition
  (
      Teuchos::RCP<DRT::Discretization> evaldis,
      Teuchos::ParameterList& params,
      DRT::AssembleStrategy & strategy,
      const std::string& condstring,
      const int condid
  );


  /*!
  \brief search in search tree for background elements in given radius

  \author rauch
  \date 09/15

  */
  void SearchPotentiallyCoveredBackgrdElements(
          std::map<int,std::set<int> >* current_subset_tofill,
          Teuchos::RCP<GEO::SearchTree> backgrd_SearchTree,
          const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
          const LINALG::Matrix<3, 1>& point,
          const double radius,
          const int label);


  /*!
  \brief evaluate a given subset of elements

  \author rauch
  \date 09/15

  */
  void EvaluateSubsetElements(Teuchos::ParameterList& params,
                              Teuchos::RCP<DRT::Discretization> dis,
                              std::map<int,std::set<int> >& elementstoeval,
                              int action);

private:

  //! flag indicating if class is setup
  bool issetup_;

  //! flag indicating if class is initialized
  bool isinit_;

protected:

  //! returns true if Setup() was called and is still valid
  bool IsSetup()
  {return issetup_;};

  //! returns true if Init(..) was called and is still valid
  bool IsInit()
  {return isinit_;};

  //! check if \ref Setup() was called
  void CheckIsSetup()
  {if(not IsSetup()) dserror("Setup() was not called.");};

  //! check if \ref Init() was called
  void CheckIsInit()
  {if(not IsInit()) dserror("Init(...) was not called.");};

protected:

  //! set flag true after setup or false if setup became invalid
  void SetIsSetup(bool trueorfalse){issetup_ = trueorfalse;};

  //! set flag true after init or false if init became invalid
  void SetIsInit(bool trueorfalse){isinit_ = trueorfalse;};


};// class ImmersedBase


/*!
\brief interpolate quantity from the background field to a given immersed point

\author rauch
\date 09/15

*/
  template < DRT::Element::DiscretizationType sourcedistype,
             DRT::Element::DiscretizationType targetdistype >
  int InterpolateToImmersedIntPoint(const Teuchos::RCP<DRT::Discretization> sourcedis,
                                    const Teuchos::RCP<DRT::Discretization> targetdis,
                                    DRT::Element& targetele,
                                    const std::vector<double>& targetxi,
                                    const std::vector<double>& targetedisp,
                                    int action,
                                    std::vector<double>& targetdata
                                                                    )
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
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    DRT::Element* sourceele;

    // sourceele displacements
    std::vector<double> myvalues;

    // declare vector of global current coordinates
    double x [globdim];

    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;
    static const int target_dim = DRT::UTILS::DisTypeToDim<targetdistype>::dim;

    // declarations and definitions for round robin loop
    const int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    DRT::Exporter exporter(comm);

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // vector to fill by sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Epetra_SerialDenseVector> vectofill = Teuchos::rcp(new Epetra_SerialDenseVector(targetdata.size()));
    (*vectofill)(0)=-1234.0;

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));
    Teuchos::RCP<Epetra_SerialDenseVector> targetxi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(target_dim));

    Teuchos::RCP<Epetra_SerialDenseVector> normal_at_targetpoint = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    std::vector<double> normal_vec(3);// normal_vec[0]=0.0;normal_vec[1]=0.0; normal_vec[2]=0.0;

    Epetra_SerialDenseMatrix dummy1;
    Epetra_SerialDenseMatrix dummy2;
    Epetra_SerialDenseVector dummy3;

    Teuchos::ParameterList params;

    params.set<std::string>("action", "calc_cur_normal_at_point");
    for(int i=0; i<target_dim;++i)
      (*targetxi_dense)(i) = targetxi[i];

    DRT::Element::LocationArray targetla(1);
    targetele.LocationVector(*targetdis,targetla,false);

    targetele.Evaluate(params,*targetdis,targetla,dummy1,dummy2,*normal_at_targetpoint,*targetxi_dense,dummy3);
    normal_at_targetpoint->Scale(1.0/(normal_at_targetpoint->Norm2()));
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
    int datalength = (int)vectofill->Length();

    // get possible elements being intersected by immersed structure
    DRT::Condition* searchbox = sourcedis->GetCondition("ImmersedSearchbox");
    std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();
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
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      LINALG::Matrix<source_dim,1> xi(true);
      if(!matched and mysearchboxgeomsize>0)
      {
        // every proc that has searchboxgeom elements
        for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
        {
          // only proc who owns background searchbox ele is supposed to interpolate and is set BoundaryIsImmersed(true)
          if(curr->second->Owner()==myrank)
          {
            bool converged = false;
            double residual = -1234.0;

            if(problem->ProblemType()==prb_immersed_ale_fsi)
            {
              Teuchos::RCP<const Epetra_Vector> state = sourcedis->GetState("dispnp");

              DRT::Element::LocationArray la(1);
              curr->second->LocationVector(*sourcedis,la,false);
              // extract local values of the global vectors
              myvalues.resize(la[0].lm_.size());
              DRT::UTILS::ExtractMyValues(*state,myvalues,la[0].lm_);
              double sourceeledisp[24];
              for(int node=0;node<8;++node)
                for(int dof=0; dof<3;++dof)
                  sourceeledisp[node*3+dof]=myvalues[node*4+dof];

              // node 1  and node 7 coords of current source element (diagonal points)
              const double* X1 = curr->second->Nodes()[1]->X();
              double x1[3];
              x1[0]=X1[0]+sourceeledisp[1*3+0];
              x1[1]=X1[1]+sourceeledisp[1*3+1];
              x1[2]=X1[2]+sourceeledisp[1*3+2];
              const double* X7 = curr->second->Nodes()[7]->X();
              double diagonal = sqrt(pow(X1[0]-X7[0],2)+pow(X1[1]-X7[1],2)+pow(X1[2]-X7[2],2));

              // calc distance of current target point to arbitrary node (e.g. node 1) of curr source element
              double distance = sqrt(pow(x1[0]-xvec[0],2)+pow(x1[1]-xvec[1],2)+pow(x1[2]-xvec[2],2));

              // get parameter space coords xi in source element of global point xvec of target element
              // NOTE: if the target point xvec is very far away from the source element curr
              //       it is unnecessary to jump into this functon and invoke a newton iteration.
              // Therfore: only call GlobalToCurrentLocal if distance is smaller than fac*characteristic element length
              if(distance < 1.5*diagonal)
              {
                MORTAR::UTILS::GlobalToCurrentLocal<sourcedistype>(*(curr->second),&sourceeledisp[0],&xvec[0],&xi(0),converged,residual);
                if(converged == false)
                {
//                  std::cout<<"Warning!! GlobalToCurrentLocal() in InterpolateToImmersedIntPoint() did not converge! res ="<<std::setprecision(12)<<residual<<"\n "
//                  "Source Element ID = "<<curr->second->Id()<<
//                  "  d="<<distance<<" < "<<"  diagonal="<<diagonal<<"\n"
//                  "pos target point = ["<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"] "<<
//                  "xi in source ele="<<xi<<std::endl;
                  xi(0)=2.0; xi(1)=2.0; xi(2)=2.0;
                }
              }
              else
              {
                xi(0)=2.0; xi(1)=2.0; xi(2)=2.0;
              }
            }
            else
            {
              // node 1  and node 7 coords of current source element (diagonal points)
              const double* X1 = curr->second->Nodes()[1]->X();
              const double* X7 = curr->second->Nodes()[7]->X();
              double diagonal = sqrt(pow(X1[0]-X7[0],2)+pow(X1[1]-X7[1],2)+pow(X1[2]-X7[2],2));

              // calc distance of current target point to arbitrary node (e.g. node 1) of curr source element
              double distance = sqrt(pow(X1[0]-xvec[0],2)+pow(X1[1]-xvec[1],2)+pow(X1[2]-xvec[2],2));

              if(distance<1.5*diagonal)
              {
                // get parameter space coords xi in source element of global point xvec of target element
                MORTAR::UTILS::GlobalToLocal<sourcedistype>(*(curr->second),&xvec[0],&xi(0),converged);
                if(converged==false)
                {
                  xi(0)=2.0; xi(1)=2.0; xi(2)=2.0;
                }
              }
              else
              {
                xi(0)=2.0; xi(1)=2.0; xi(2)=2.0;
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
            if (abs(xi(0))<(1.0-tol) and abs(xi(1))<(1.0-tol) and abs(xi(2))<(1.0-tol))
            {
              validsource = true;

              (*xi_dense)(0)=xi(0);
              (*xi_dense)(1)=xi(1);
              (*xi_dense)(2)=xi(2);
            }
            // does targetpoint lie on any element edge?
            // targetpoint lies within tolerance away from any element edge
            // only the sourceele located outward from the boundary should interpolate
            else if(( abs((abs(xi(0))-1.0))<=tol ) // at xi_0 surface
                or  ( abs((abs(xi(1))-1.0))<=tol ) // at xi_1 surface
                or  ( abs((abs(xi(2))-1.0))<=tol ) // at xi_2 surface
                   )
            {
              double scalarproduct= 0.0;
              Teuchos::RCP<Epetra_SerialDenseVector> vector;

              // get nodal coords of source ele
              vector = Teuchos::rcp(new Epetra_SerialDenseVector(3));
              std::vector<double> distances(sourceele->NumNode());

              // loop over source element nodes and calc distance from targetpoint to those nodes
              for(int sourcenode=0; sourcenode<sourceele->NumNode(); sourcenode++)
              {
                if(problem->ProblemType()==prb_immersed_ale_fsi)
                {
                  const double* X;
                  X = sourceele->Nodes()[sourcenode]->X();

                  distances[sourcenode] = sqrt(pow(xvec[0]-(X[0]+myvalues[sourcenode*4+0]),2)+pow(xvec[1]-(X[1]+myvalues[sourcenode*4+1]),2)+pow(xvec[2]-(X[2]+myvalues[sourcenode*4+2]),2));

                }
                else
                {
                  const double* X;
                  X = sourceele->Nodes()[sourcenode]->X();

                  distances[sourcenode] = sqrt(pow(xvec[0]-X[0],2)+pow(xvec[1]-X[1],2)+pow(xvec[2]-X[2],2));
                }
              }

              // get max. distance and corresponding index of src node
              int maxdistanceindex = -1;
              double maxvalue = 0.0;
              for(int i=0; i<(int)distances.size();++i)
              {
                if(distances[i]>maxvalue)
                {
                  maxvalue = distances[i];
                  maxdistanceindex = i;
                }
              }

              const double* X;
              double x[3];
              if(problem->ProblemType()==prb_immersed_ale_fsi)
              {
                X = sourceele->Nodes()[maxdistanceindex]->X();
                x[0] = X[0]+myvalues[maxdistanceindex*4+0];
                x[1] = X[1]+myvalues[maxdistanceindex*4+1];
                x[2] = X[2]+myvalues[maxdistanceindex*4+2];
              }
              else
              {
                x[0] = sourceele->Nodes()[maxdistanceindex]->X()[0];
                x[1] = sourceele->Nodes()[maxdistanceindex]->X()[1];
                x[2] = sourceele->Nodes()[maxdistanceindex]->X()[2];
              }

              // build vector directed from targetele point (xvec) to sourceele node (X)
              // xvec should not lie at same ele surface as sourcenode
              (*vector)(0) = x[0]-xvec[0];
              (*vector)(1) = x[1]-xvec[1];
              (*vector)(2) = x[2]-xvec[2];
              vector->Scale(vector->Norm2());

              // build scalar product between normal and vector
              scalarproduct = (*vector)(0) * normal_vec[0] + (*vector)(1) * normal_vec[1] + (*vector)(2) * normal_vec[2];

              if(scalarproduct > tol)
                validsource = true;
              else
                validsource = false;

              std::cout<<"WARNING !! Immersed point is lying on element edge! validsource="<<validsource<<"for eleid="<<sourceele->Id()<<std::endl;

#ifdef DEBUG
              if(abs(scalarproduct) < 1e-13 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node perpendicular.\n"
                    "no robust determination of outward lying src ele possible");
              }
              if(abs(scalarproduct) < -1e-12 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node have no sharp angle.\n"
                    "no robust determination of outward laying src ele possible");
              }
#endif

              //count++;
            } // if on element edge


            // evaluate and perform action handed in to this function
              if(action == 81311)
              { // experimental pre cmbe
                params.set<std::string>("action","interpolate_porosity_to_given_point");

              }
              else if(action == 813112)
              {
                params.set<int>("action",(int)SCATRA::calc_immersed_phi_at_given_point);
              }
              else
                params.set<int>("action", action);
            params.set<int>("Physical Type",INPAR::FLUID::poro_p1);
            if(validsource)
            {
//              std::cout<<"\n ---------- \n PROC "<<comm.MyPID()<<" ; sourcele : "<<sourceele->Id()<<std::endl;
//              std::cout<<"POINT: ["<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl;
//              std::cout<<"normal: ["<<normal_vec[0]<<" "<<normal_vec[1]<<" "<<normal_vec[2]<<"]"<<std::endl;
//              std::cout<<"vector: "<<*vector<<std::endl;
//              std::cout<<"scalarproduct: "<<scalarproduct<<std::endl;

              // fill locationarray
              DRT::Element::LocationArray sourcela(1);
              sourceele->LocationVector(*sourcedis,sourcela,false);
              sourceele->Evaluate(params,*sourcedis,sourcela,dummy1,dummy2,*vectofill ,*xi_dense,dummy3);
              matched = 1;

              // if element interpolates to structural intpoint set as "boundary is immersed"
              Teuchos::RCP<DRT::ELEMENTS::FluidImmersedBase> immersedele = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(curr->second);
              // only possible and reasonable if current element is fluid element
              if(immersedele != Teuchos::null)
              {
                //std::cout<<"ID="<<immersedele->Id()<<"  scalarprod="<<std::setprecision(10)<<scalarproduct<<std::endl;
                immersedele->SetBoundaryIsImmersed(2);
                immersedele->ConstructElementRCP(problem->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND"));
              }

              // DEBUG output
              //std::cout<<"Tgt. Ele: "<<targetele.Id()<<" lies in src. ele: "<<curr->second->Id()<<" at  pos. ["<<std::setprecision(5)<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl;

              break; // break loop over all source eles
            }
            else
              matched = 0;

          } // if myrank == owner of searchbox ele and proc has searchbox ele
        }// loop over all seachbox elements
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if(numproc>1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)vectofill->Length();++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }

          // normal vector
          for(int dim=0;dim<globdim;++dim)
          {
            data.AddtoPack(normal_vec[dim]);
          }

          std::swap(sdata, data());

        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*vectofill)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        for(int i=0;i<globdim;++i)
          normal_vec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if num procs > 1

      if(mysearchboxgeomsize>0)
      {
        if(curr==(searchboxgeom.end()--) and irobin == (numproc-1) and matched == 0)
        {
          std::cout<<"target position = ["<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl;
          dserror("could not match given point on any proc. Element(0)=%f",(*vectofill)(0));
        }
      }

    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff

    for(int dim=0;dim<datalength;++dim)
      targetdata[dim]=(*vectofill)(dim);


    return err;
  };


  /*!
  \brief interpolate quantity from the immersed field to a given background point

  \author rauch
  \date 09/15

  */
  template < DRT::Element::DiscretizationType sourcedistype,
  DRT::Element::DiscretizationType targetdistype >
  int InterpolateToBackgrdPoint(
      std::map<int,std::set<int> >& curr_subset_of_structdis,
      const Teuchos::RCP<DRT::Discretization> sourcedis,
      const Teuchos::RCP<DRT::Discretization> targetdis,
      DRT::Element& targetele,
      const std::vector<double>& targetxi,
      const std::vector<double>& targetedisp,
      const std::string action,
      std::vector<double>& targetdata,
      bool& match,
      bool vel_calculation,
      bool doCommunication = true
  )
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
    DRT::Problem* problem = DRT::Problem::Instance();
    // dimension of global problem
    const int globdim = problem->NDim();
    // pointer to background element
    DRT::Element* sourceele;
    // declare vector of global current coordinates
    double x [globdim];
    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;
    // declarations and definitions for round robin loop
          int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    if(numproc == 1)
      doCommunication=false;
    if(doCommunication==false)
      numproc=1;

    Teuchos::RCP<const Epetra_Vector> dispnp = sourcedis->GetState("displacement");
    DRT::Element::LocationArray la(1);

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // vector to be filled in sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Epetra_SerialDenseVector> vectofill = Teuchos::rcp(new Epetra_SerialDenseVector(4));
    (*vectofill)(0)=-1234.0;
    if(match) // was already matched, add no contribution
    {(*vectofill)(0)=0.0;
     (*vectofill)(1)=0.0;
     (*vectofill)(2)=0.0;
     (*vectofill)(3)=0.0; }

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));

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
    if(match == true)
      matched = 1;

    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->Length();

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
      LINALG::Matrix<source_dim,1> xi(true);
      if(matched == 0)
      {
        for(std::map<int, std::set<int> >::const_iterator closele = curr_subset_of_structdis.begin(); closele != curr_subset_of_structdis.end(); closele++)
        {
          for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
          {
            bool converged = false;
            double residual = -1234.0;
            sourceele = sourcedis->gElement(*eleIter);

            // get parameter space coords xi in source element
            sourceele->LocationVector(*sourcedis,la,false);
            std::vector<double> mysourcedispnp(la[0].lm_.size());
            DRT::UTILS::ExtractMyValues(*dispnp,mysourcedispnp,la[0].lm_);

            // construct bounding box around current source element
            Teuchos::RCP<GEO::CUT::BoundingBox> bbside = Teuchos::rcp( GEO::CUT::BoundingBox::Create() );
            LINALG::Matrix<3, 1> nodalpos;
            for ( int i=0; i<(int)(sourceele->NumNode()); ++i )
            {
              nodalpos(0,0) = (sourceele->Nodes())[i]->X()[0] + mysourcedispnp[i*3+0];
              nodalpos(1,0) = (sourceele->Nodes())[i]->X()[1] + mysourcedispnp[i*3+1];
              nodalpos(2,0) = (sourceele->Nodes())[i]->X()[2] + mysourcedispnp[i*3+2];
              bbside->AddPoint( nodalpos );
            }

            // check if given point lies within bounding box
            bool within = bbside->Within(1.0,&xvec[0]);

            // only try to match given target point and sourceele if within bounding box
            if(within)
            {
              MORTAR::UTILS::GlobalToCurrentLocal<sourcedistype>(*sourceele,&mysourcedispnp[0],&xvec[0],&xi(0),converged,residual);

              if(converged==false)
              {
                std::cout<<"!! WARNING !! GlobalToCurrentLocal did not not converge (res="<<residual<<") though background point " <<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"\n"
                           " is lying within bounding box of immersed element "<<sourceele->Id()<<"."<<std::endl;
                bbside->Print();
                // we assume not in bounding box. point is lying outside sourcele.
                xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
              }
            }
            else
            {
              // not in bounding box. point is lying outside sourcele.
              xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
            }

            /////////////////////////////////////////////////
            ////
            ////     MATCH
            ////
            /////////////////////////////////////////////////
            double tol = 1e-12;
            if ( abs(xi(0))<=(1.0+tol) and abs(xi(1))<=(1.0+tol) and abs(xi(2))<=(1.0+tol))
            { // given point lies in element curr
              matched = 1;
              match = true;

              // parameter list
              Teuchos::ParameterList params;
              params.set<std::string>("action", action);
              // fill locationarray
              DRT::Element::LocationArray la(1);
              sourceele->LocationVector(*sourcedis,la,false);

              (*xi_dense)(0)=xi(0);
              (*xi_dense)(1)=xi(1);
              (*xi_dense)(2)=xi(2);

              Epetra_SerialDenseMatrix dummy1;
              Epetra_SerialDenseMatrix dummy2;
              Epetra_SerialDenseVector dummy3;

              // check whether velocity or divergence should be interpolated
              int vel_or_div = 1;
              if (not vel_calculation)
                vel_or_div = 0;
              params.set<int>("calculate_velocity",vel_or_div);

              if(action != "none")
                // evaluate and perform action handed in to this function
                sourceele->Evaluate(params,*sourcedis,la,dummy1,dummy2,*vectofill ,*xi_dense,dummy3);
              else
              {
                 for(int dim=0;dim<datalength;++dim)
                   (*vectofill)(dim)=0.0;
              }

              break; // break loop over all source elements
            }// if match

            // given targetpoint point lies not underneath source discretization at all
            if(closele == (curr_subset_of_structdis.end()--) and eleIter == ((closele->second).end()--) and matched == 0 and irobin == (numproc-1))
            {
              for(int dim=0;dim<datalength;++dim)
                (*vectofill)(dim)=-12345.0;
            }
          }// loop over all seachbox elements
        }// loop over closeele set
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if(doCommunication)
      {
        DRT::Exporter exporter(comm);

        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)vectofill->Length();++dim)
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
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*vectofill)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if numproc > 1

    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff
    for(int dim=0;dim<datalength;++dim)
      targetdata[dim]=(*vectofill)(dim);


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
  template < DRT::Element::DiscretizationType structdistype,
  DRT::Element::DiscretizationType fluiddistype >
  int FindClosestStructureSurfacePoint(
      std::map<int,std::set<int> >& curr_subset_of_structdis,   //!< Input
      const Teuchos::RCP<DRT::Discretization> structdis,        //!< Input
      const Teuchos::RCP<DRT::Discretization> fluiddis,         //!< Input
      DRT::Element& fluidele,                                   //!< Input
      const std::vector<double>& fluidxi,                       //!< Input
      const std::vector<double>& fluideledisp,                  //!< Input
      const std::string action,                                 //!< Input
      std::vector<double>& target_vel_data,                     //!< Output
      std::vector<double>& target_pos_data,                     //!< Output
      bool& match,                                              //!< Output
      bool doCommunication = true                               //!< Input
  )
  {

    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = structdis->Comm();

    // get the global problem
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // declare vector of global current coordinates of fluid node
    double x_fluid_node[globdim];

    // dimension of structure and fluid discretization
    static const int struct_dim = DRT::UTILS::DisTypeToDim<structdistype>::dim;
    static const int fluid_dim = DRT::UTILS::DisTypeToDim<fluiddistype>::dim;

    // get integration rule of structural surface element
    const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_quad_36point);

    // initialize structural shape functions matrix
    LINALG::Matrix<4,1> shapefcts(true);

    // declarations and definitions for round robin loop
    int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    if(numproc == 1)
      doCommunication=false;

    if(doCommunication==false)
      numproc=1;

    DRT::Node** fluidnode = fluidele.Nodes();
    double char_fld_ele_length = sqrt(pow(fluidnode[1]->X()[0]-fluidnode[7]->X()[0],2)+pow(fluidnode[1]->X()[1]-fluidnode[7]->X()[1],2)+pow(fluidnode[1]->X()[2]-fluidnode[7]->X()[2],2));

    // get current displacements and velocities of structure discretization
    Teuchos::RCP<const Epetra_Vector> dispnp = structdis->GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> velnp = structdis->GetState("velocity");
    DRT::Element::LocationArray la(structdis->NumDofSets());

    // get current global coordinates of the given fluid node fluid_xi
    MORTAR::UTILS::LocalToCurrentGlobal<fluiddistype>(fluidele,globdim,&fluidxi[0],fluideledisp,&x_fluid_node[0]);
    // get as vector
    std::vector<double> fluid_node_glob_coord(globdim);
    for(int idim=0; idim<3;++idim)
      fluid_node_glob_coord[idim] = x_fluid_node[idim];

    // initialize vectors to fill
    Teuchos::RCP<Epetra_SerialDenseVector> velnp_at_struct_point = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    Teuchos::RCP<Epetra_SerialDenseVector> xi_pos_of_struct_point = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    (*velnp_at_struct_point)(0) = -1234.0;
    (*xi_pos_of_struct_point)(0) = -1234.0;

    // given point already matched?
    int matched = 0;
    if(match == true)
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
    int datalength = (int)velnp_at_struct_point->Length();

    // loop over all procs
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      // local coordinates of structrual surface point in parameter space of the fluid element
      LINALG::Matrix<fluid_dim,1> fluid_xi(true);

      // loop over points that are not yet matched
      if(matched == 0)
      {
        // loop over closeele set
        for(std::map<int, std::set<int> >::const_iterator closele = curr_subset_of_structdis.begin(); closele != curr_subset_of_structdis.end(); closele++)
        {
          // loop over all searchbox elements
          for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
          {
            // only surface elements (with immersed coupling condition) are relevant to find closest structure point
            std::vector<Teuchos::RCP<DRT::Element> > surface_eles = structdis->gElement(*eleIter)->Surfaces();

            // loop over surface element, find the elements with IMMERSEDCoupling condition
            for (std::vector<Teuchos::RCP<DRT::Element> >::iterator surfIter = surface_eles.begin(); surfIter != surface_eles.end(); ++surfIter)
            {
              // pointer to current surface element
              DRT::Element* structele  = surfIter->getRawPtr();

              // pointer to nodes of current surface element
              DRT::Node** NodesPtr = surfIter->getRawPtr()->Nodes();

              int numfsinodes = 0;

              // loop over all nodes of current surface element
              for (int surfnode = 0; surfnode < structele->NumNode(); ++surfnode)
              {
                DRT::Node* checkNode = NodesPtr[surfnode];
                // check whether a IMMERSEDCoupling condition is active on this node
                if (checkNode->GetCondition("IMMERSEDCoupling") != NULL)
                  ++numfsinodes;
              }

              // only surface elements with IMMERSEDCoupling condition are relevant (all nodes of this element have fsi condition)
              if (numfsinodes==structele->NumNode())
              {
                // get displacements and velocities of structure dis
                structele->LocationVector(*structdis,la,false);
                std::vector<double> mydispnp(la[0].lm_.size());
                DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,la[0].lm_);
                std::vector<double> myvelnp(la[0].lm_.size());
                DRT::UTILS::ExtractMyValues(*velnp,myvelnp,la[0].lm_);

                // 1.) check if closest point is a node
                // loop over all nodes of structural surface element
                for (int inode = 0; inode < structele->NumNode(); ++inode)
                {
                  // get current global position of structure node
                  std::vector<double> struct_node_postion(3);
                  for(int idim=0; idim<3;++idim)
                    struct_node_postion[idim]=(structele->Nodes())[inode]->X()[idim] + mydispnp[inode*3+idim];

                  // distance between fluid node and given structure node
                  double dist = sqrt(pow(fluid_node_glob_coord[0]-struct_node_postion[0],2)+pow(fluid_node_glob_coord[1]-struct_node_postion[1],2)+pow(fluid_node_glob_coord[2]-struct_node_postion[2],2));

                  // only store point information, if distance is smallest so far
                  if (dist<stored_dist)
                  {
                    // update distance, velocity and position
                    stored_dist = dist;
                    for (int idim=0; idim<3;++idim)
                    {
                      vel_to_store[idim] = myvelnp[inode*3+idim];
                      pos_to_store[idim] = struct_node_postion[idim];
                    }
                  }
                } // end loop over nodes

                // 2.) check if closest point is an integration point of the structural surface element
                // loop over all int points of structural surface element
                for (int igp = 0; igp <intpoints.IP().nquad; ++igp)
                {
                  // get local coordinates of gp
                  std::vector<double> struct_xsi(struct_dim);
                  struct_xsi[0] = intpoints.IP().qxg[igp][0];
                  struct_xsi[1] = intpoints.IP().qxg[igp][1];

                  // get local coordinates of gp as matrix
                  LINALG::Matrix<struct_dim,1> struct_xi(true);
                  for(int i=0;i<struct_dim;++i)
                    struct_xi(i) = struct_xsi[i];

                  // declare vector of global current coordinates
                  double structsurf_gp_glob_coord[globdim];

                  // get current global position of the given int point xi on the structural surface
                  MORTAR::UTILS::LocalToCurrentGlobal<structdistype>(*(structele),globdim,&struct_xsi[0],mydispnp,&structsurf_gp_glob_coord[0]);

                  // distance between fluid node and given structure gp
                  double dist = sqrt(pow(fluid_node_glob_coord[0]-structsurf_gp_glob_coord[0],2)+pow(fluid_node_glob_coord[1]-structsurf_gp_glob_coord[1],2)+pow(fluid_node_glob_coord[2]-structsurf_gp_glob_coord[2],2));

                  // only store point information, if distance is smallest so far
                  if (dist<stored_dist)
                  {
                    // update distance
                    stored_dist = dist;

                    // evaluate shape functions of structure dis at gp
                    DRT::UTILS::shape_function<DRT::Element::quad4>(struct_xi,shapefcts);
                    LINALG::Matrix<4,3> myvelocitynp(true);
                    for(int node=0;node<4;++node)
                      for(int dim=0;dim<3;++dim)
                        myvelocitynp(node,dim) = myvelnp[node*3+dim];

                    // calculate velocity in gp
                    LINALG::Matrix<3,1> vel_at_gp(true);
                    vel_at_gp.MultiplyTN(myvelocitynp,shapefcts);

                    // update velocity and position
                    for(int idim=0; idim<3;++idim)
                    {
                      vel_to_store[idim] = vel_at_gp(idim,0);
                      pos_to_store[idim] = structsurf_gp_glob_coord[idim];
                    }

                  }
                } // end loop over int points
              } // end if fsi coupling
            } // end loop over all surface elements
          }// loop over all seachbox elements
        }// loop over closeele set

        // global position of closest point (that is now found) as vector
        std::vector<double> xvec(globdim);
        for(int idim=0; idim<3;++idim)
          xvec[idim] = pos_to_store[idim];

        bool converged = false;
        double residual = -1234.0;

        if(stored_dist < 1.5*char_fld_ele_length)
        {
          // get local coordinates of closest structural surface point in fluid element parameter space
          if(problem->ProblemType()==prb_immersed_ale_fsi)
            MORTAR::UTILS::GlobalToCurrentLocal<fluiddistype>(fluidele,&fluideledisp[0],&xvec[0],&fluid_xi(0),converged,residual);
          else
            MORTAR::UTILS::GlobalToLocal<fluiddistype>(fluidele,&xvec[0],&fluid_xi(0),converged);

          if(converged==false)
          {
            std::cout<<"WARNING!! GlobalToCurrentLocal() did not converge in FindClosestStructureSurfacePoint( \n"
                "xi="<<fluid_xi<<std::endl;
            fluid_xi(0) = 2.0; fluid_xi(1) = 2.0; fluid_xi(2) = 2.0;
          }
        }
        else
        {
          fluid_xi(0) = 2.0; fluid_xi(1) = 2.0; fluid_xi(2) = 2.0;
        }

        // check if closest point lies in current fluid element
        double tol = 1e-13;
        if ( abs(fluid_xi(0))<(1.0+tol) and abs(fluid_xi(1))<(1.0+tol) and abs(fluid_xi(2))<(1.0+tol) )
        { // given point lies current fluid element
          matched = 1;
          match = true;
          // give back velocity and local coordinates of closest point
          for (int dim=0;dim<3;++dim)
          {
            (*velnp_at_struct_point)(dim) = vel_to_store[dim];
            (*xi_pos_of_struct_point)(dim) = fluid_xi(dim);
          }

        }// if match
      }// enter loop only if node not yet matched

      if(doCommunication)
      {
        DRT::Exporter exporter(comm);
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)velnp_at_struct_point->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(fluid_node_glob_coord[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)velnp_at_struct_point->Length();++dim)
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
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          fluid_node_glob_coord[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*velnp_at_struct_point)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if numproc > 1

    } // loop over all procs

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have received his original data");
#endif

    // now every proc should store the requested quantities in velnp_at_struct_point for his fluidxi
    // time to return that stuff
    for(int dim=0;dim<datalength;++dim)
    {
      target_vel_data[dim] = (*velnp_at_struct_point)(dim);
      target_pos_data[dim] = (*xi_pos_of_struct_point)(dim);
    }

    return err;
  };


  /*!
  \brief Find the intersection of a line normal to a given surface point with the same surface but on the opposite side

  \author rauch
  \date 2015

   */
  template < DRT::Element::DiscretizationType sourcedistype,
  DRT::Element::DiscretizationType targetdistype >
  double IntersectLineWithBoundary
  (
      Teuchos::RCP<DRT::Discretization> evaldis,
      const Teuchos::RCP<const Epetra_Vector> dispnp,
      Teuchos::ParameterList& params,
      LINALG::Matrix<3,1>& searchdirection,
      LINALG::Matrix<3,1>& spacepoint,
      const std::string& condstring,
      const int condid
  )
  {
    // is to be set 1 if intersection is found
    int matched = 0;

    // result variable
    double tempresult = -1234.0;

# ifdef DEBUG
    if (!(evaldis->Filled()) ) dserror("FillComplete() was not called");
    if (!(evaldis->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");
# endif

    // get communicator
    const Epetra_Comm& comm = evaldis->Comm();

    // declarations and definitions for round robin loop
    int numproc  = comm.NumProc();
    static const int myrank   = comm.MyPID();                     // myrank
    static const int torank   = (myrank + 1) % numproc;           // sends to
    static const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    int owner = myrank;

    DRT::Element::LocationArray la(evaldis->NumDofSets());

    std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator fool;

    //----------------------------------------------------------------------
    // loop through conditions and evaluate them if they match the criterion
    //----------------------------------------------------------------------
    for (fool=evaldis->GetAllConditions().begin(); fool!=evaldis->GetAllConditions().end(); ++fool)
    {
      if (fool->first == condstring)
      {
        DRT::Condition& cond = *(fool->second);
        if (condid == -1 || condid ==cond.GetInt("ConditionID"))
        {
          // get geometry
          std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
          // DEBUG check
          if (geom.empty()) dserror("evaluation of condition with empty geometry on proc %d",evaldis->Comm().MyPID());

          // define iterator over surface geometry elements
          std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;

          // loop over all geometry elements and check for intersection with
          // boundary along search direction
          LINALG::Matrix<3,1> p;
          LINALG::Matrix<3,1> m1;
          LINALG::Matrix<3,1> m2;
          LINALG::Matrix<3,3> A;
          LINALG::Matrix<3,1> b;
          LINALG::Matrix<3,1> parameters;
          double determinant = -1234.0;

          // round robin loop
          for (int irobin = 0; irobin < numproc; ++irobin)
          {
            std::vector<char> sdata;
            std::vector<char> rdata;

            // if the spacepoint on the calling proc has not yet been matched
            // to a point on the opposite surface (in inward normal direction) ...
            // if matched == 1 then tempresult is already correct
            if (matched == 0)
            {
              // ... loop over all geometry elements on this proc
              for (curr=geom.begin(); curr!=geom.end(); ++curr)
              {
                // get nodes
                DRT::Node** bdrynodes = curr->second->Nodes();
                const double* X0 = bdrynodes[0]->X();
                const double* X1 = bdrynodes[1]->X();
                const double* X3 = bdrynodes[3]->X();

                // get location vector
                curr->second->LocationVector(*evaldis,la,false,condstring,params);
                // extract  bdry ele displacements
                std::vector<double> bdryeledisp(la[0].lm_.size());
                DRT::UTILS::ExtractMyValues(*dispnp,bdryeledisp,la[0].lm_);

                // target ele space point p is curr. pos. of node 0
                // from this point a aux plane is built
                // warping of the element is neglected
                p(0,0) = X0[0] + bdryeledisp[0*3 + 0];
                p(1,0) = X0[1] + bdryeledisp[0*3 + 1];
                p(2,0) = X0[2] + bdryeledisp[0*3 + 2];

                // direction vector m1 is vector from node 0 to node 1
                // first vector to span the aux plane
                m1(0,0) = X1[0]+bdryeledisp[1*3 + 0] - p(0,0);
                m1(1,0) = X1[1]+bdryeledisp[1*3 + 1] - p(1,0);
                m1(2,0) = X1[2]+bdryeledisp[1*3 + 2] - p(2,0);

                // direction vector m1 is vector from node 0 to node 3
                // second vector to span the aux plane
                m2(0,0) = X3[0]+bdryeledisp[3*3 + 0] - p(0,0);
                m2(1,0) = X3[1]+bdryeledisp[3*3 + 1] - p(1,0);
                m2(2,0) = X3[2]+bdryeledisp[3*3 + 2] - p(2,0);

                // fill matrix A = [m1,m2,-searchdirection]
                A(0,0) = m1(0);
                A(1,0) = m1(1);
                A(2,0) = m1(2);

                A(0,1) = m2(0);
                A(1,1) = m2(1);
                A(2,1) = m2(2);

                A(0,2) = -searchdirection(0);
                A(1,2) = -searchdirection(1);
                A(2,2) = -searchdirection(2);

                determinant=A.Invert();

                // vector b = spacepoint - p
                b.Update(1.0,spacepoint,-1.0,p,0.0);

                // solve for parameters
                parameters.Multiply(A,b);
                tempresult = parameters(2,0);

                double eps = 1.0e-06;

                // parameters(2,0) is the length of the line intersecting the surface vis-a-vis the spacepoint
                if( tempresult>1e-04 and (parameters(0,0)<(1.0+eps) and parameters(0,0)>(0.0-eps)) and (parameters(1,0)<(1.0+eps) and parameters(1,0)>(0.0-eps)) )
                {
                  matched = 1;
                  if(abs(determinant)<1e-08)
                    dserror("determinant close to zero or negative! det=%f",determinant);
                  break;
                }

              } // loop over all geometry elements
            } // only if not yet matched

            // send data around to the other procs if there are any
            if(comm.NumProc() > 1)
            {
              DRT::Exporter exporter(comm);

              // ---- pack data for sending -----
              {
                DRT::PackBuffer data;
                data.StartPacking();
                {
                  data.AddtoPack(matched); //< owner already found matching point vis-a-vis spacepoint
                  data.AddtoPack(owner);   //< myrank = owner of data
                  DRT::ParObject::AddtoPack<3,1>(data,searchdirection); //< searchdirection of owner
                  DRT::ParObject::AddtoPack<3,1>(data, spacepoint);     //< spacepoint of owner

                  data.AddtoPack(tempresult);
                }

                std::swap(sdata, data());
              }
              //---------------------------------

              // ---- send ----
              MPI_Request request;
              exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

              // ---- receive ----
              int length = rdata.size();
              int tag = -1;
              int from = -1;
              exporter.ReceiveAny(from,tag,rdata,length);
              if (tag != 1234 or from != fromrank)
                dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

              // ---- unpack data -----
              std::vector<char>::size_type position = 0;
              matched    = DRT::ParObject::ExtractInt(position,rdata);
              owner      = DRT::ParObject::ExtractInt(position,rdata);
              DRT::ParObject::ExtractfromPack<3,1>(position,rdata,searchdirection);
              DRT::ParObject::ExtractfromPack<3,1>(position,rdata,spacepoint);

              tempresult = DRT::ParObject::ExtractDouble(position,rdata);

              // wait for all communication to finish
              exporter.Wait(request);
              comm.Barrier();
            } // only if numproc > 1

//            if(irobin == numproc-1 and matched == 0 )
//            {
//              // DEBUG output
//              std::cout<<"parameters:"<<parameters<<std::endl;
//              std::cout<<"m1"<<m1<<std::endl;
//              std::cout<<"m2"<<m2<<std::endl;
//              std::cout<<"searchdirection"<<searchdirection<<std::endl;
//              std::cout<<"spacepoint"<<spacepoint<<std::endl;
//              std::cout<<"p"<<p<<std::endl;
//              b(0,0)=p(0,0)+parameters(0,0)*m1(0,0)+parameters(1,0)*m2(0,0);
//              b(1,0)=p(1,0)+parameters(0,0)*m1(1,0)+parameters(1,0)*m2(1,0);
//              b(2,0)=p(2,0)+parameters(0,0)*m1(2,0)+parameters(1,0)*m2(2,0);
//              std::cout<<"intersectionpoint: "<<b<<std::endl;
//              b.Update(-1.0,spacepoint,1.0);
//              std::cout<<"connectingvector "<<b<<"\n with length: "<<b.Norm2()<<std::endl;
//                std::cout<<"Warning! Found no intersection for evaluation of cell diameter!"<<std::endl;
//            }
          } // robin loop

          // if not matched after round robin loop, something went wrong
          if(matched == 0)
            std::cout<<"No intersection for evaluation of cell diameter on rank "<<myrank<<std::endl;

        } // if check of condid successful
      } // if condstring found
    } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)

    return tempresult;
  };

  /*!
  \brief get submatrix of element least squares matrix and vector for immersed structural nodes and gp

  \author fuchs
  \date 11/16

  */
  template < DRT::Element::DiscretizationType structdistype,
  DRT::Element::DiscretizationType fluiddistype >
  int GetInFluidEleImmersedPoints(
      std::map<int,std::set<int> >& curr_subset_of_structdis,   // relevant structural elements
      const Teuchos::RCP<DRT::Discretization> structdis,        // structure discretization
      DRT::Element& fluid_ele,                                  // current fluid element
      Epetra_SerialDenseMatrix& elemat1,                        // element least squares matrix
      Epetra_SerialDenseVector& elevec1,                        // element least squares vector
      bool fillelemat                                           // true: fill least squares matrix and vector, false: check if element is immersed
  )
  {
    // error flag
    int err = 0;

    Teuchos::RCP<const Epetra_Vector> struct_dispnp = structdis->GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> struct_velnp = structdis->GetState("velocity");

    DRT::Element::LocationArray struct_la(structdis->NumDofSets());

    // pointer to structural element
    DRT::Element* struct_ele;

    // space dimension
    static const int space_dim = 3;

    // nodal degree of freedom
    static const int node_dof = 4;

    // position and velocity of node or gp of structural element
    std::vector<double> struct_pos(space_dim);
    std::vector<double> struct_vel(space_dim);

    // and stored in linalg matrix
    LINALG::Matrix<space_dim,1> struct_pos_m(true);
    LINALG::Matrix<space_dim,1> struct_vel_m(true);

    // number of nodes of fluid element and structure element
    static const int fluid_numnode = DRT::UTILS::DisTypeToNumNodePerEle<fluiddistype>::numNodePerElement;
    static const int struct_numnode = DRT::UTILS::DisTypeToNumNodePerEle<structdistype>::numNodePerElement;

    // current configuration and velocity of structural element
    LINALG::Matrix<struct_numnode,space_dim> struct_ele_curr_config(true);
    LINALG::Matrix<struct_numnode,space_dim> struct_ele_curr_vel(true);

    // reference configuration of fluid element
    LINALG::Matrix<fluid_numnode,space_dim> fluid_ele_refe_config(true);

    for (int i=0; i<fluid_numnode; ++i)
    {
      fluid_ele_refe_config(i,0) = (fluid_ele.Nodes())[i]->X()[0];
      fluid_ele_refe_config(i,1) = (fluid_ele.Nodes())[i]->X()[1];
      fluid_ele_refe_config(i,2) = (fluid_ele.Nodes())[i]->X()[2];
    }

    // structural position in local coordinates of current fluid element
    LINALG::Matrix<space_dim,1> curr_xi(true);

    // shape function for fluid and structure
    LINALG::Matrix<fluid_numnode,1> fluid_shapefcts(true);
    LINALG::Matrix<struct_numnode,1> struct_shapefcts(true);

    // shape function derivatives for structure
    LINALG::Matrix<2,struct_numnode> struct_derivs(true);

    // get nodal parameter space coordinates of structural element
    Epetra_SerialDenseMatrix struct_nodpoints = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(structdistype);

    // get gauss point parameter space coordinates of structural element
    DRT::UTILS::IntegrationPoints2D struct_intpoints = DRT::UTILS::IntegrationPoints2D( DRT::UTILS::intrule_quad_9point );

    // number of search points in structural element
    int num_search_pts = struct_numnode + struct_intpoints.nquad;

    // initialize vector of search point coordinates in parameter space
    std::vector<double> search_pts_xi(num_search_pts);
    std::vector<double> search_pts_eta(num_search_pts);

    // fill nodal parameter space coordinates
    for (int nd=0; nd<struct_numnode; ++nd)
    {
      search_pts_xi[nd] = struct_nodpoints[nd][0];
      search_pts_eta[nd] = struct_nodpoints[nd][1];
    }

    // fill gauss point parameter space coordinates
    for (int gp=0; gp<struct_intpoints.nquad; ++gp)
    {
      search_pts_xi[struct_numnode+gp] = struct_intpoints.qxg[gp][0];
      search_pts_eta[struct_numnode+gp] = struct_intpoints.qxg[gp][1];
    }

    // convergence flag of inverse mapping
    bool converged = false;

    // loop over all (potentially) immersed structural elements
    for (std::map<int, std::set<int> >::const_iterator closele = curr_subset_of_structdis.begin(); closele != curr_subset_of_structdis.end(); closele++)
    {
      for (std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
      {
        // get current structural element
        struct_ele = structdis->gElement(*eleIter);

        // get current structural displacements and velocities
        struct_ele->LocationVector(*structdis,struct_la,false);

        std::vector<double> mydispnp(struct_la[0].lm_.size());
        DRT::UTILS::ExtractMyValues(*struct_dispnp,mydispnp,struct_la[0].lm_);

        std::vector<double> myvelnp(struct_la[0].lm_.size());
        DRT::UTILS::ExtractMyValues(*struct_velnp,myvelnp,struct_la[0].lm_);

        // get current configuration and velocity of current structural element
        for (int i=0; i<struct_numnode; ++i)
        {
          struct_ele_curr_config(i,0) = (struct_ele->Nodes())[i]->X()[0]+mydispnp[i*space_dim+0];
          struct_ele_curr_config(i,1) = (struct_ele->Nodes())[i]->X()[1]+mydispnp[i*space_dim+1];
          struct_ele_curr_config(i,2) = (struct_ele->Nodes())[i]->X()[2]+mydispnp[i*space_dim+2];

          struct_ele_curr_vel(i,0) = myvelnp[i*space_dim+0];
          struct_ele_curr_vel(i,1) = myvelnp[i*space_dim+1];
          struct_ele_curr_vel(i,2) = myvelnp[i*space_dim+2];
        }

        // loop over all search points of current structural element
        for (int i=0; i<num_search_pts; ++i)
        {
          converged = false;

          // evaluate shape functions of the structural element at current search point
          DRT::UTILS::shape_function_2D(struct_shapefcts,search_pts_xi[i],search_pts_eta[i],struct_ele->Shape());

          // position of current search point
          struct_pos_m.MultiplyTN(1.0,struct_ele_curr_config,struct_shapefcts,0.0);

          struct_pos[0] = struct_pos_m(0);
          struct_pos[1] = struct_pos_m(1);
          struct_pos[2] = struct_pos_m(2);

          // velocity of current search point
          struct_vel_m.MultiplyTN(1.0,struct_ele_curr_vel,struct_shapefcts,0.0);

          struct_vel[0] = struct_vel_m(0);
          struct_vel[1] = struct_vel_m(1);
          struct_vel[2] = struct_vel_m(2);

          // get structural position in local coordinates of current fluid element
          MORTAR::UTILS::GlobalToLocal<fluiddistype>(fluid_ele,&struct_pos[0],&curr_xi(0),converged);

          // check next search point
          if (converged==false)
            continue;

          // current structural point should lie purely inside current fluid element
          // attention: a fluid element with just one structural point on the element surface/edge would result in a singular least squares matrix
          if (abs(curr_xi(0))<1.0 and abs(curr_xi(1))<1.0 and abs(curr_xi(2))<1.0)
          {
            // just check if element is immersed
            if (not fillelemat)
            {
              // set the flags for the immersed boundary fluid element
              DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(&fluid_ele);

              immersedele->SetIsImmersed(1);
              immersedele->SetHasProjectedDirichlet(1);

              // element is marked so we can leave this method
              return err;
            }

            // evaluate shape functions of current fluid element at current structural point
            DRT::UTILS::shape_function<fluiddistype>(curr_xi,fluid_shapefcts);

            // double loop over fluid element nodes
            for (int k=0; k<fluid_numnode; ++k)
            {
              // add contribution to least squares vector
              elevec1((node_dof*k)+0) += 2.0*fluid_shapefcts(k)*struct_vel[0];
              elevec1((node_dof*k)+1) += 2.0*fluid_shapefcts(k)*struct_vel[1];
              elevec1((node_dof*k)+2) += 2.0*fluid_shapefcts(k)*struct_vel[2];

              for (int l=0; l<fluid_numnode; ++l)
              {
                // add contribution to least squares matrix
                elemat1( (node_dof*k)+0 , (node_dof*l)+0 ) += 2.0*fluid_shapefcts(k)*fluid_shapefcts(l);
                elemat1( (node_dof*k)+1 , (node_dof*l)+1 ) += 2.0*fluid_shapefcts(k)*fluid_shapefcts(l);
                elemat1( (node_dof*k)+2 , (node_dof*l)+2 ) += 2.0*fluid_shapefcts(k)*fluid_shapefcts(l);
              }
            }
          }
        } // end loop over all search points of current structural element
      }
    } // end loop over all (potentially) immersed structural elements

    return err;
  }

}// namespace IMMERSED

#endif

