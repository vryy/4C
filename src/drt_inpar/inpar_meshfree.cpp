/*-----------------------------------------------------------*/
/*!
\file inpar_meshfree.cpp

\brief inpar meshfree

\maintainer Keijo Nissen

\level 2

*/
/*-----------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_meshfree.H"



void INPAR::MESHFREE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& meshfree = list->sublist("MESHFREE",false,"");
  setStringToIntegralParameter<int>("TYPE","MaxEnt","Type of meshfree discretisation.",
                                    tuple<std::string>("MaxEnt","Particle","GeoDecoupled"),
                                    tuple<int>(maxent,
                                               particle,
                                               geo_decoupled),
                                    &meshfree);
  setStringToIntegralParameter<int>("NODEPOINTASSIGNMENT","procwise","Type of assignment of nodes to cells/knots.",
                                    tuple<std::string>("procwise","blockwise"),
                                    tuple<int>(procwise,
                                               blockwise),
                                    &meshfree);
  DoubleParameter("NEWTON_TOL",1e-6,"Tolerance at which Newton is considered to be converged.",&meshfree);
  DoubleParameter("NEWTON_MAXITER",10,"Maximum number of Newton steps.",&meshfree);
  IntParameter("DBC_SOLVER",-1,"Solver number for solving non-constant Dirichlet BC if necessary.",&meshfree);
  BoolParameter("PARTITION_OF_UNITY","Yes","Enforcement of partition of unity constraint",&meshfree);
  DoubleParameter("NEGATIVITY",0.0,"Decides if and to which degree negativity is allowed",&meshfree);
  DoubleParameter("VARIANCE",1,"Variance of the basis solution function prior.",&meshfree);
  DoubleParameter("RANGE_TOL",1e-6,"Threshhold at which basis solution function prior is considered nmuerically zero.",&meshfree);
  DoubleParameter("CUTOFF_RADIUS",-1.0,"Cutoff radius for influence of meshfree points on each other.",&meshfree);
  setNumericStringParameter("BIN_PER_DIR","-1 -1 -1",
                            "Number of bins per direction (x, y, z) in particle simulations.",
                            &meshfree);
  setNumericStringParameter("BOUNDINGBOX","-1e12 -1e12 -1e12 1e12 1e12 1e12",
                            "Bounding box for binning strategy in particle simulations.",
                            &meshfree);
  setStringToIntegralParameter<int>("WRITEBINS","none","Visualize none, row or column bins",
                                    tuple<std::string>("none","rows","cols"),
                                    tuple<int>(
                                      none,
                                      rows,
                                      cols),
                                    &meshfree);
  setStringToIntegralParameter<int>("PRIOR","Gauss","Defines the prior type of the basis solution function.",
                                    tuple<std::string>("Gauss"),
                                    tuple<int>(p_gauss),
                                    &meshfree);
  setStringToIntegralParameter<int>("S_COMPLIANCE","linear","Defines the compliance type enforced for max-ent scheme of the basis solution functions.",
                                    tuple<std::string>("linear","stream","freespace"),
                                    tuple<int>(
                                      c_linear,
                                      c_stream,
                                      c_freesp),
                                    &meshfree);

  setStringToIntegralParameter<int>("W_COMPLIANCE","linear","Defines the compliance type enforced for max-ent scheme of the basis solution functions.",
                                    tuple<std::string>("linear","stream","freespace"),
                                    tuple<int>(
                                      c_linear,
                                      c_stream,
                                      c_freesp),
                                    &meshfree);
}
