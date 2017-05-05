/*-----------------------------------------------------------*/
/*!
\file str_dbc_biopolynet.H

\brief dbc management for biopolymer network simulations

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_DBC_BEAMINTERACTION_H_
#define STR_DBC_BEAMINTERACTION_H_

#include <Teuchos_RCP.hpp>
#include "str_dbc.H"


namespace GEO{
namespace MESHFREE{
  class BoundingBox;
}
}
namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
  class SerialDenseMatrix;
  class MapExtractor;
}

namespace STR {
namespace TIMINT {
  class MeshFreeData;
}

  class DbcBioPolyNet : public Dbc
  {
  public:
    //! Constructor
    DbcBioPolyNet();

    //! Destructor
    virtual ~DbcBioPolyNet() {};


    /*! \brief Evaluate and apply the DBC
     *
     * Note: Stay in the global coordinate system (Rotation: global-->local-->global).*/
    virtual void ApplyDirichletBC(
        const double& time,
        Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> vel,
        Teuchos::RCP<Epetra_Vector> acc,
        bool recreatemap);

    //! Initialize class variables
    virtual void Init(const Teuchos::RCP<DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_Vector>& freact,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

  protected:

    /*! \brief adjusted DBC evaluation for periodic case
     *  \param params          (in)     a parameter list
     *  \param dis             (in)     displacement vector
     *  \param vel             (in)     velocity vector
     *  \param dbcmapextractor (in/out) map extractor for Dirichlet DOFs (new method)*/
    void EvaluateDirichletPeriodic(Teuchos::ParameterList&            params,
                                   Teuchos::RCP<Epetra_Vector>        dis,
                                   Teuchos::RCP<Epetra_Vector>        vel,
                                   Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

    /*! \brief Determine Dirichlet condition at given time and apply its values to a system vector
     *  \param nodeids      (in)      vector containing DBC Node IDs
     *  \param onoff        (in)      vector containing onoff values for DOFs
     *  \param dis          (in)      displacement vector
     *  \param deltadbc     (in)      vector holding the increment in Dirichlet values
     *  \param dbcgids      (in/out)  set holding dbc DOF GIDs (new method)*/
    void DoDirichletConditionPeriodic(std::vector<int>*            nodeids,
                                      std::vector<int>*            onoff,
                                      Teuchos::RCP<Epetra_Vector>  dis,
                                      Teuchos::RCP<Epetra_Vector>  deltadbc,
                                      Teuchos::RCP<std::set<int> > dbcgids);

    /*! \brief Dirichlet conditions needed for viscoelastic measurements*/
    void DBCOscillatoryMotion(Teuchos::ParameterList&     params,
                              Teuchos::RCP<Epetra_Vector> dis,
                              Teuchos::RCP<Epetra_Vector> vel,
                              Teuchos::RCP<Epetra_Vector> deltadbc);

    /*! \brief Apply affine shear deformation to network
     *  \param params   (in)  Parameter List
     *  \param dis      (in)  displacement vector
     *  \param vel      (in)  velocity vector
     *  \param dbcgids  (in/out)  set of DOF GIDs on which Dirichlet values are imposed*/
    void DBCAffineShear(Teuchos::ParameterList&     params,
                        Teuchos::RCP<Epetra_Vector> dis,
                        Teuchos::RCP<Epetra_Vector> vel,
                        Teuchos::RCP<Epetra_Vector> deltadbc);

    /*! \brief pin down and release nodes at their location by Dirichlet Values*/
    void DBCPinNodes();

    /*! \brief Dirichlet condition accompanying constant Neumann shear force
     *  \param params (in)  parameter list
     *  \param dis    (in)  row displacement vector*/
    void DBCMovableSupport1D(Teuchos::ParameterList&     params,
                             Teuchos::RCP<Epetra_Vector> dis);

    /*! \brief Get DBCs defined in Input file and ad the DBCs DOFs to dbcgids
     *  \param params   (in)  Parameter List
     *  \param dis      (in)  displacement vector
     *  \param vel      (in)  velocity vector
     *  \param dbcgids  (in/out)  set of DOF GIDs on which Dirichlet values are imposed*/
    void DBCGetPredefinedConditions(Teuchos::ParameterList&       params,
                                    Teuchos::RCP<Epetra_Vector>   dis,
                                    Teuchos::RCP<Epetra_Vector>   vel,
                                    Teuchos::RCP<std::set<int> >  dbcgids);

    /*! \brief set Dirichlet values according to given Dirichlet map)
     *  \param disn            (in)     displacement vector
     *  \param deltadbc        (in)     Dirichlet increments for DOFs
     *  \param  dbcgids        (in/out)  set holding dbc DOF GIDs (new method)*/
    void DBCSetValues(Teuchos::RCP<Epetra_Vector>        dis,
                      Teuchos::RCP<Epetra_Vector>        deltadbc,
                      Teuchos::RCP<std::set<int> >       dbcgids);

    /*! \brief pin down a node at its location by Dirichlet Values
     *  \param dbcgids               (in/out)  set holding dbc DOF GIDs (new method)
     *  \param dbcmapextractor       (in)      vector containing onoff values for free DOFs */
    void DBCCreateMap(Teuchos::RCP<std::set<int> >       dbcgids,
                      Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor);

    /*!\brief return the current coordinates and DOF LIDs of all nodes of given element
       \param element (in)      an element
       \param discol  (in)      current displacement vector (col map)
       \param coord   (in/out)  current node coordinates of element
       \param lids    (in/out)  nodal DOF LIDs of the element*/
    void GetElementNodeCoords(DRT::Element*               element,
                              Teuchos::RCP<Epetra_Vector> discol,
                              LINALG::SerialDenseMatrix&  coord,
                              std::vector<int>*           lids=NULL);

    /*!\brief check whether element is "broken" due to periodic boundary conditions
      \param coord  (in)  set of nodal coordinates
      \param cut    (in/out)  in case of a broken element, cut holds component and information about the location of the nodes of the element*/
    bool CheckForBrokenElement(LINALG::SerialDenseMatrix& coord,
                               LINALG::SerialDenseMatrix& cut);

  protected:

    //! periodic bounding box
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> periodic_boundingbox_;

    //! \brief vector of node sets for DBC application
    std::vector<std::vector<int> > dbcnodesets_;

    //! \brief flag signaling the use of the initially found DBC node set throughout the simulation
    bool useinitdbcset_;

  };
}


#endif /* STR_DBC_BEAMINTERACTION_H_ */
