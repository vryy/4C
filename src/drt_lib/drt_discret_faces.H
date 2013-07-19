/*!----------------------------------------------------------------------
\file drt_discret_faces.H

\brief a class to manage an enhanced discretization including the faces between elements


</pre>

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/
#ifndef DISCRET_FACES_H
#define DISCRET_FACES_H

#include "drt_discret.H"
#include "drt_dserror.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include <string>
#include "Teuchos_ParameterList.hpp"

namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
}


/*!
\brief DRT: namespace of the baci discretization module

*/
namespace DRT
{

class DiscretizationFaces : public Discretization
{
public:

  /*!
   * \brief internal class that holds the information used to create internal faces elements
   *
   * \author Benedikt Schott (schott@lnm.mw.tum.de)
   */
  class InternalFacesData
  {
  public:

    /*!
    \brief Standard Constructor

    \param master_peid (in): element id of master parent element
    \param slave_peid (in): element id of slave parent element
    \param lsurface_master (in): local index of surface w.r.t master parent element
    \param nodes (in): vector of nodes building the surface element
    */
    InternalFacesData(int master_peid, std::vector<DRT::Node*> nodes, int lsurface_master)
    {
      master_peid_ = master_peid;
      slave_peid_ = -1;
      lsurface_master_ = lsurface_master;
      lsurface_slave_ = -1;
      nodes_ = nodes;
    }

    /*!
    \brief Destructor

    */
    virtual ~InternalFacesData()
    {
      return;
    };


    /*--- set ------------------------------------------*/

    //! set the parent element id for slave parent element
    void SetSlavePeid(int eid){slave_peid_ = eid;}

    //! set the local surface number w.r.t slave parent element
    void SetLSurfaceSlave(int lsurface_slave){lsurface_slave_ = lsurface_slave;}

    /*!
    \brief set the map for the face's nodes between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system
           and the slave element's face's coordinate system
    */
    void SetLocalNumberingMap(std::vector<int> localtrafomap){localtrafomap_ = localtrafomap;}


    /*--- get ------------------------------------------*/

    //! get the master parent element id
    int GetMasterPeid(){return master_peid_;}

    //! get the slave parent element id
    int GetSlavePeid(){return slave_peid_;}

    //! get the local surface number w.r.t master parent element
    int GetLSurfaceMaster(){return lsurface_master_;}

    //! get the local surface number w.r.t slave parent element
    int GetLSurfaceSlave(){return lsurface_slave_;}

    //! get the transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system
    //! and the slave element's face's coordinate system
    std::vector<int> GetLocalNumberingMap(){return localtrafomap_;}

    //! get surface's nodes (unsorted, original)
    std::vector<DRT::Node*> GetNodes(){return nodes_;}

  private:

    int master_peid_;          //!< master parent element id
    int slave_peid_;           //!< slave parent element id

    int lsurface_master_;      //!< local surface number w.r.t master parent element
    int lsurface_slave_;       //!< local surface number w.r.t slave parent element

    std::vector<DRT::Node*> nodes_; //!< vector of surface nodes, order w.r.t master parent element

    /*!
     \brief map for the face's nodes between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system
            and the slave element's face's coordinate system
     */
    std::vector<int> localtrafomap_;

  };




  /*!
  \brief Standard Constructor

  \param name (in): name of this discretization
  \param comm (in): An epetra comm object associated with this discretization
  */
  DiscretizationFaces(const string name, RCP<Epetra_Comm> comm)
  :Discretization(name, comm) { extension_filled_ = false; }; // use base class constructor

  /*!
  \brief Destructor

  */
  virtual ~DiscretizationFaces()
  {
    return;
  };



  /*!
  \brief Compute the nullspace of the discretization

  This method looks in the solver parameters whether algebraic multigrid (AMG)
  is used as preconditioner. AMG desires the nullspace of the
  system of equations which is then computed here if it does not already exist
  in the parameter list.

  \note This method is supposed to go away and live somewhere else soon....

  \param solveparams (in): List of parameters
  \param recompute (in)  : force method to recompute the nullspace
  */
  virtual void ComputeNullSpaceIfNecessary(Teuchos::ParameterList& solveparams,
                                           bool recompute = false)
  {
    // remark: the null space is not computed correctly for XFEM discretizations, since the number of
    //         degrees of freedom per node is not fixed.
    //         - it is not clear what happens with respect to the Krylov projection
    //           (having XFEM dofs seems to render the system non-singular, but it should be singular
    //           so the null space has a non-zero dimension)
    //         - the ML preconditioner also relies on a fixed number of dofs per node
    DRT::Discretization::ComputeNullSpaceIfNecessary(solveparams, recompute);
  }

  /*!
  \brief Complete construction of a discretization  (Filled()==true NOT prerequisite)

  After adding or deleting nodes or elements or redistributing them in parallel,
  or adding/deleting boundary conditions, this method has to be called to (re)construct
  pointer topologies.<br>
  It builds in this order:<br>
  Standard FillComplete of base class
  - row map of nodes
  - column map of nodes
  - row map of elements
  - column map of elements
  - pointers from elements to nodes
  - pointers from nodes to elements
  - assigns degrees of freedoms
  - map of element register classes
  - calls all element register initialize methods
  - build geometries of all Dirichlet and Neumann boundary conditions

  Additional features
  - build internal faces elements
  - build maps and pointers for internal faces

  \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                       assigning of degrees of freedoms to nodes and
                                       elements.
  \param initelements (in) : if true, build element register classes and call Initialize()
                             on each type of finite element present
  \param doboundaryconditions (in) : if true, build geometry of boundary conditions
                                     present.
  \param createinternalfacdes (in) : if true, build geometry of internal faces.

  \note In order to receive a fully functional discretization, this method must be called
        with all parameters set to true (the default). The parameters though can be
        used to turn off specific tasks to allow for more flexibility in the
        construction of a discretization, were it is known that this method will
        be called more than once.

  \note Sets Filled()=true
  \author schott 03/12
  */
  int FillCompleteFaces(Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave,
                        bool assigndegreesoffreedom = true,
                        bool initelements           = true,
                        bool doboundaryconditions   = true,
                        bool createinternalfaces    = false);

  /*!
  \brief Get flag indicating whether CreateInternalFacesExtension() has been called
  */
  virtual inline bool FilledExtension() const {return extension_filled_;}

  /*!
  \brief Get map associated with the distribution of the ownership of internal faces
         (Filled()==true prerequisite)

  This map includes all internal faces stored on this proc and also owned by this proc.
  This map is non-ambiguous, meaning that it is a non-overlapping map.

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* IntFacesRowMap() const;
  /*!
  \brief Get map associated with the distribution of elements including ghosted internal faces
         (Filled()==true prerequisite)

  This map includes all internal faces stored on this proc including any ghosted internal faces
  This map is ambiguous, meaning that it is an overlapping map

  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite.
  */
  virtual const Epetra_Map* IntFacesColMap() const;

  /*!
  \brief Get global number of internal faces (true number of total elements)
         (Filled()==true prerequisite)

  This is a collective call
  */
  virtual int NumGlobalIntFaces() const;

  /*!
  \brief Get processor local number of internal faces owned by this processor
         (Filled()==true prerequisite)
  */
  virtual int NumMyRowIntFaces() const;

  /*!
  \brief Get processor local number of internal faces including ghost elements
         (Filled()==true NOT prerequisite)
  */
  virtual int NumMyColIntFaces() const;

  /*!
  \brief Get the internal face element with local row id lid (Filled()==true prerequisite)

  Returns the internal face element with local row index lid.
  Will not return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite

  \return Adress of internal face element if element is owned by calling proc
  */
  virtual inline DRT::Element* lRowIntFace(int lid) const
  {
    #ifdef DEBUG
    if (!Filled()) dserror("DRT::DiscretizationFaces::lRowIntFace: Filled() != true");
    #endif
    return intfacerowptr_[lid];
  }

  /*!
  \brief Get the element with local column id lid (Filled()==true prerequisite)

  Returns the internal face element with local column index lid.
  Will also return any ghosted element.
  This is an individual call and Filled()=true is a prerequisite

  \return Adress of internal face element if element is stored by calling proc
  */
  virtual inline DRT::Element* lColIntFace(int lid) const
  {
  #ifdef DEBUG
    if (!Filled()) dserror("DRT::DiscretizationFaces::lColIntFace: Filled() != true");
  #endif
    return intfacecolptr_[lid];
  }

  /*!
  \brief Build internal faces extension
  */
  void CreateInternalFacesExtension(Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave,
                                    const bool verbose = false);

  /*!
  \brief Evaluate edge-based integrals
  */
  void EvaluateEdgeBased(
       Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
       Teuchos::RCP<Epetra_Vector>          systemvector1
  );

  /*!
  \brief Complete construction of a internal faces elements
  */
  void BuildInternalFaces(Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave,
                          const bool verbose = false);

  /*!
  \brief Build intfacerowmap_ (Filled()==true NOT prerequisite)

  Build the parallel layout of internal faces in this
  discretization and store it as an Epetra_Map in intfacerowmap_
  intfacerowmap_ is unique.
  It considers internal faces owned by a proc only

  \note This is a collective call

  */
  virtual void BuildIntFaceRowMap();

  /*!
  \brief Build intfacecolmap_ (Filled()==true NOT prerequisite)

  Build the potentially overlapping parallel layout of internal faces in this
  discretization and store it as an Epetra_Map in intfacecolmap_
  intfacecolmap_ includes ghosted internal faces and is potentially overlapping.

  \note This is a collective call

  */
  virtual void BuildIntFaceColMap();

  /*!
  \brief Print Print internal faces discretization to os (Filled()==true NOT prerequisite)
         (ostream << also supported)

  \note This is a collective call
  */
  void PrintIntFaces(ostream& os) const;


protected:

  bool                                  extension_filled_;      ///< flag indicating whether intfaces extension has been filled

  Teuchos::RCP<Epetra_Map>              intfacerowmap_;         ///< unique distribution of element ownerships
  Teuchos::RCP<Epetra_Map>              intfacecolmap_;         ///< distribution of elements including ghost elements
  std::vector<DRT::Element*>            intfacerowptr_;         ///< vector of pointers to row elements for faster access
  std::vector<DRT::Element*>            intfacecolptr_;         ///< vector of pointers to column elements for faster access
  std::map<int,RCP<DRT::Element> >      faces_;                 ///< map of internal faces elements


}; // class DiscretizationXFEM
} // namespace DRT

/// << operator
ostream& operator << (ostream& os, const DRT::DiscretizationFaces& dis);


#endif  // #ifndef DISCRET_FACES_H
