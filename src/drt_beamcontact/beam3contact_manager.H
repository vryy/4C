/*!----------------------------------------------------------------------
\file beam3contact_manager.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifndef BEAM3C_MANAGER_H
#define BEAM3C_MANAGER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RCP.hpp"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

#include "../linalg/linalg_fixedsizematrix.H"

namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
  class Element;
}

class Beam3ContactOctTree;

namespace CONTACT
{



/*!
 \brief contact manager for contact in a beam3 discretization

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author popp (popp@lnm.mw.tum.de), cyron (cyron@lnm.mw.tum.de)
 */
class Beam3contact;

class Beam3cmanager
{
public:

  //! @name Friends

  //no fried classes defined

  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param discret (in): A discretization containing beam elements

  */
  Beam3cmanager(DRT::Discretization& discret, double alphaf);

  /*!
  \brief Destructor

  */
  virtual ~Beam3cmanager() {};

  //@}

  //! @name Access methods

  /*!
  \brief Print this beam3 contact manager

  */
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Get problem discretization

  */
  inline const DRT::Discretization& ProblemDiscret() const { return pdiscret_; }

  /*!
  \brief Get beam contact discretization

  */
  inline DRT::Discretization& ContactDiscret() { return *cdiscret_; }

  /*!
  \brief Get offset of dofs between cdiscret and pdiscret

  */
  inline int& DofOffset() { return dofoffset_; }

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return pdiscomm_; }

  /*!
  \brief Get different node or element maps

  */
  inline const Teuchos::RCP<Epetra_Map> RowNodes() const { return noderowmap_; }
  inline const Teuchos::RCP<Epetra_Map> ColNodes() const { return nodecolmap_; }
  inline const Teuchos::RCP<Epetra_Map> FullNodes() const { return nodefullmap_; }
  inline const Teuchos::RCP<Epetra_Map> RowElements() const { return elerowmap_; }
  inline const Teuchos::RCP<Epetra_Map> ColElements() const { return elecolmap_; }
  inline const Teuchos::RCP<Epetra_Map> FullElements() const { return elefullmap_; }
  inline const std::vector<Teuchos::RCP<Beam3contact> >& Pairs() const { return outputpairs_; }
  inline const Teuchos::RCP<Beam3ContactOctTree> OcTree() const {return tree_; }

  /*!
  \brief Get list of beam contact input parameters

  */
  inline const Teuchos::ParameterList& BeamContactParameters() { return sbeamcontact_; }

  /*!
  \brief Get list of general contact input parameters

  */
  inline const Teuchos::ParameterList& GeneralContactParameters() { return scontact_; }

  /*!
  \brief Get current constraint norm

  */
  const double GetConstrNorm() { return constrnorm_; }

  /*!
  \brief Get current Uzawa iteration index

  */
  const int GetUzawaIter() { return uzawaiter_; }

  // \brief Get current penalty parameter
  const double GetCurrentpp() { return currentpp_; }

  //@}

  //! @name Public evaluation methods

  /*!
  \brief Evaluate beam contact

  First, we search for potential beam element pairs coming into contact.
  For each pair, a temporary Beam3contact object is generated, which handles
  penalty force and stiffness computation. Then, this method calls each beam
  contact pair to compute its contact forces and stiffness. Finally, all entries
  are assembles into global force resdiual and global stiffness matrix.

  */
  void Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fres,
                const Epetra_Vector& disrow, bool newsti = false);

  /*!
  \brief Update beam contact

  Stores fc_ into fcold_ and clears fc_ as needed for generalized alpha time
  integration scheme at the end of each time step. ASCII output files for
  visualization in GMSH will be written. Also some output to screen is done.

  */
  void Update(const Epetra_Vector& disrow, const int& timestep, const int& newtonstep);

  /*!
  \brief Initialize second, third, ... Uzawa step

  This method is called at the beginning of the second, third, ... Uzawa
  iterarion in order to create an of an out-of-balance force again. First,
  the contact force and stiffness terms are removed from stiff and fres.
  Then the new contact forces and contact stiffness terms for the updated LM
  (AugmentedLagrange) are created by calling Evaluate().

  */
  void InitializeUzawa(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fres,
                       const Epetra_Vector& disrow, bool newsti = false);

  /*!
  \brief Rest Lagrange Multipliers

  Reset all Lagrange multipliers of each contact pair to zero at the beginning
  of each new time step. This is only needed for an augmented Lagrange / Uzawa scheme.

  */
  void ResetAlllmuzawa();

  /*!
  \brief Update constraint norm

  In the penalty case, we only compute the constraint norm and print to screen,
  whereas in the augmented Lagrange case, an empiric update formula for the
  penalty parameter is also applied to speed up Uzawa convergence.

  \param cnorm (in/out) read-only (no norm update) if cnorm is not the null pointer
  \param convstatus (in) status of convergence, default value NULL for converged

  */
  void UpdateConstrNorm(double* cnorm = NULL);

  /*!
  \brief Shift current normal "normal_" vector to old normal vector "normal_old_"

  The converged normal vector of the last time step is stored as "normal_old_" to serve
  as a reference for the modified gap function definition

  */
  void ShiftAllNormal();

  /*!
  \brief Update Lagrange multipliers

  Update all Lagrange multipliers of each contact pair at the beginning of each
  new Uzawa time step. This is only needed for an augmented Lagrange / Uzawa scheme.

   */
  void UpdateAlllmuzawa();

  /*!
  \brief Reset penalty parameter

  In an augmented Lagrange strategy the penalty parameter might be increased to
  accelerate convergence of the Uzawa iteration. At the beginning of each time step,
  it has to be reset to the initial value given in the input file.

  */
  void ResetCurrentpp();

  /*!
  \brief Reset Uzawa iteration index

  */
  void ResetUzawaIter();

  /*!
  \brief Reset contact Pairs (e.g. if time step is redone)

  */
  void ResetPairs();

  /*!
  \brief Update Uzawa iteration index

  */
  void UpdateUzawaIter();

  /*!
  \brief Create output files for GMSH visualization

  Create ASCII-files to visualize beams with GMSH. The circular cross section will
  be approximated by prisms, which are rotated around the element's axis. This output
  method only works safely for the serial case, the parallel case is not yet implemented!

  */
  void GmshOutput(const Epetra_Vector& disrow, const int& timestep,
                  const int& newtonstep, bool endoftimestep = false);

  /*!
  \brief Print active set

  Print some output data to screen at the end of each time step.
  Interesting values are:
    a) IDs of current pairs and their elements
    b) the residual gap of this pair
    c) the current (augmented part) Lagrange multiplier of this pair
    d) the current element coordinates of the contact point

  NOTE: This method can also be called after each newton-step (e.g. if you want to check
  convergence problems). In this case, you have to uncomment the GMSHNEWTONSTEP preprocessor
  flag in 'beam3contact_defines.h'.

  */
  void ConsoleOutput();

  /*!
  \brief Write Bearing Reactions into a csv-file

  */
  void Reactions(const Epetra_Vector& fint,
                 const Epetra_Vector& dirichtoggle,
                 const int& timestep);

  //@}

private:

  // don't want = operator and cctor
  Beam3cmanager operator = (const Beam3cmanager& old);
  Beam3cmanager(const Beam3cmanager& old);

  //! @name member variables

  //! problem discretizaton
  DRT::Discretization& pdiscret_;

  //! contact discretization (basically a copy)
  Teuchos::RCP<DRT::Discretization> cdiscret_;

  //! the Comm interface of the problem discretization
  const Epetra_Comm& pdiscomm_;

  //! dof offset between pdicsret and cdiscret
  int dofoffset_;

  //! node and element maps
  Teuchos::RCP<Epetra_Map> noderowmap_;
  Teuchos::RCP<Epetra_Map> nodecolmap_;
  Teuchos::RCP<Epetra_Map> nodefullmap_;
  Teuchos::RCP<Epetra_Map> elerowmap_;
  Teuchos::RCP<Epetra_Map> elecolmap_;
  Teuchos::RCP<Epetra_Map> elefullmap_;
  
  //! occtree for contact search
  Teuchos::RCP<Beam3ContactOctTree> tree_;

  //! vector of contact pairs (pairs of elements, which might get in contact)
  std::vector<Teuchos::RCP<Beam3contact> > pairs_;
  //! vector of contact pairs stored for output (output is called after Update() which erases all pairs!)
  std::vector<Teuchos::RCP<Beam3contact> > outputpairs_;

  //! 2D-map with pointers on the contact pairs_. This map is necessary, to call a contact pair directly by the two element-iD's of the pair.
  // It is not needed at the moment due to the direct neigbour determination in the constructor but may be useful for future operations
  std::map<std::pair<int,int>, Teuchos::RCP<Beam3contact> > contactpairmap_;
  //std::map<int, Teuchos::RCP<Beam3contact> > contactpairmap_;

  //! parameter list for beam contact options
  Teuchos::ParameterList sbeamcontact_;

  //! parameter list for general contact options
  Teuchos::ParameterList scontact_;

  //! search radius
  double searchradius_;

  //! contact forces of current time step
  Teuchos::RCP<Epetra_Vector> fc_;

  //! contact forces of previous time step (for generalized alpha)
  Teuchos::RCP<Epetra_Vector> fcold_;

  //! contact stiffness matrix of current time step
  Teuchos::RCP<LINALG::SparseMatrix> stiffc_;

  //! time integration parameter (0.0 for statics)
  double alphaf_;

  //! current constraint norm (violation of non-penetration condition)
  double constrnorm_;

  //! current penalty parameter (might be modified within augmented Lagrange strategy)
  double currentpp_;

  //! current Uzawa iteration index (for augmented Lagrange strategy)
  int uzawaiter_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Search contact pairs

  We search pairs of elements that might get in contact and then create the
  corresponding beam3contact objects. Pairs of element that are direct neighbours,
  i.e. share one node, will be rejected.

  */
  void SearchPossibleContactPairs(std::map<int,LINALG::Matrix<3,1> >& currentpositions);

  /*!
  \brief Compute the search radius

  The search radius will be computed by analyzing the chracteristic length of each
  element. To guarantee, that each possible contact pair will be detected some
  empiric criterion will define the search radius, taking into account:

    a) the maximum element radius
     b) the maximum element length

  These two characteric lengths will be compared, the larger one is the characteristic
  length for this processor. Then via a communication among all procs the largest
  characteristic length in the whole discretization is found. Using this global
  characteristic length, we can compute a searchradius by multiplying with a constant factor.
  This method is called only once at the beginning of the simulation. If axial deformation
  of beam elements was high, it would have to be called more often!

  */
  void ComputeSearchRadius();

  /*!
  \brief Get maximum element radius

  Finds maximum element radius in the whole discretization for circular cross
  sections. Stores the maximum radius to 'max_ele_radius'.

  */
  void GetMaxEleRadius(double& max_ele_radius);

  /*!
  \brief Get maximum element length

  Finds maximum element radius in the whole discretization for circular cross
  sections. Stores the maximum radius to 'max_ele_length'. For higher-order-elements
  an approximation of the true element length is introduced, as only the direct distance
  of the two end nodes is computed. Yet, this is assumed to be accurate enough.

  */
  void GetMaxEleLength(double& max_ele_length);

  /*!
  \brief Compute rotation matrix R from given angle theta in 3D

  This function computes from a three dimensional rotation angle theta
  (which gives rotation axis and absolute value of rotation angle) the related
  rotation matrix R. Note that theta is given in radiant.

  */
  void TransformAngleToTriad(Epetra_SerialDenseVector& theta,
                             Epetra_SerialDenseMatrix& R);

  /*!
  \brief Compute spin

  Compute spin matrix according to Crisfield Vol. 2, equation (16.8)

  */
  void ComputeSpin(Epetra_SerialDenseMatrix& spin, Epetra_SerialDenseVector& rotationangle);

  /*!
  \brief Update penalty parameter

  In an augmented Lagrangian scheme, the penalty paramter might be increased to accelerate
  convergence of the Uzawa iteration. For doing this, the constraint norms of the previous
  and the current Uzawa step are compared. If the improvemenet is too small, the penalty
  paramter will be increased. Returns bool indicating update. Note that the rule when and
  by how much the penalty paramter is increased is NOT analytically  derived. It's just
  experience. It moght also be problem-dependent.

  */
  bool IncreaseCurrentpp(const double& GlobalConstrNorm);

  /*!
  \brief Set displacment state

  In this context this means that the map of matrices currentpositions is filled with
  the current nodal positions taken from the global displacement vector.

  */
  void SetState(std::map<int,LINALG::Matrix<3,1> >& currentpositions, const Epetra_Vector& disrow);

  /*!
  \brief Compute coordinates for GMSH-Output for two-noded-elements

  */
  void GMSH_2_noded(const int& n, const Epetra_SerialDenseMatrix& coord,
                    const DRT::Element* thisele, std::stringstream& gmshfilecontent);

  /*!
  \brief Compute coordinates for GMSH-Output for three-noded-elements

  */
  void GMSH_3_noded(const int& n, const Epetra_SerialDenseMatrix& allcoord,
                    const DRT::Element* thisele, std::stringstream& gmshfilecontent);

  //@}

}; // class Beam3cmanager
} // namespace CONTACT

#endif  // #ifndef BEAM3C_MANAGER_H
