/*!----------------------------------------------------------------------
\file beam3contact_manager.H

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BEAM3C_MANAGER_H
#define BEAM3C_MANAGER_H

#include "beam3contact.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
}

namespace CONTACT
{

using namespace std;
using namespace Teuchos;

/*!
 \brief contact manager for contact in a beam3 discretization

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author popp (popp@lnm.mw.tum.de), cyron (cyron@lnm.mw.tum.de)
 */
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
  Beam3cmanager(DRT::Discretization& discret);

  /*!
  \brief Destructor

  */
  virtual ~Beam3cmanager() {};

  //@}

  //! @name Access methods

  /*!
  \brief Print this beam3 contact manager

  */
  virtual void Print(ostream& os) const;

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
  virtual const Epetra_Comm& Comm() const { return pdiscret_.Comm(); }

  /*!
  \brief Get different node or element maps

  */
  inline const RCP<Epetra_Map> RowNodes() const { return noderowmap_; }
  inline const RCP<Epetra_Map> ColNodes() const { return nodecolmap_; }
  inline const RCP<Epetra_Map> FullNodes() const { return nodefullmap_; }
  inline const RCP<Epetra_Map> RowElements() const { return elerowmap_; }
  inline const RCP<Epetra_Map> ColElements() const { return elecolmap_; }
  inline const RCP<Epetra_Map> FullElements() const { return elefullmap_; }

  /*!
	\brief Get list of contact input parameters

	*/
	inline const Teuchos::ParameterList& InputParameters() { return scontact_; }

	/*!
  \brief Get current constraint norm

  */
  const double GetConstrNorm() { return constrnorm_; }

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
                const Epetra_Vector& disrow, double alphaf);

  /*!
  \brief Update beam contact

  Stores fc_ into fcold_ and clears fc_ as needed for generalized alpha time
  integration scheme at the end of each time step. ASCII output files for
  visualization in GMSH will be written. Also some output to screen is done.

  */
  void Update(const Epetra_Vector& disrow, const int& timestep,
              const int& uzawastep, const int& newtonstep);

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

  */
  void UpdateConstrNorm(const int uzawaiter = 0);

  /*!
  \brief Update Lagrange multipliers

  Update all Lagrange multipliers of each contact pair at the beginning of each
  new Uzawa time step. This is only needed for an augmented Lagrange / Uzawa scheme.

  */
  void UpdateAlllmuzawa(const int& uzawaiter);

	/*!
  \brief Reset penalty parameter

  In an augmented Lagrange strategy the penalty parameter might be increased to
  accelerate convergence of the Uzawa iteration. At the beginning of each time step,
  it has to be reset to the initial value given in the input file.

  */
	void ResetCurrentpp();

	/*!
  \brief Create output files for GMSH visualization

  Create ASCII-files to visualize beams with GMSH. The circular cross section will
  be approximated by prisms, which are rotated around the element's axis. This output
  method only works safely for the serial case, the parallel case is not yet implemented!

  */
	void GmshOutput(const Epetra_Vector& disrow, const int& timestep,
	                const int& uzawastep, const int& newtonstep);

	/*!
	\brief Print active set

	Print some output data to screen at the end of each time step.
	Interesting values are:
		a) IDs of current pairs and their elements
		b) the residual gap of this pair
		c) the current (augmented part) Lagrange multiplier of this pair
		d) the current element coordinates of the contact point

	NOTE: This method can also be called after each newton-step (e.g. if you want to check
	convergence problems). In this case, you have to uncomment the #define GMSHNEWTONSTEP
	flag in 'beam3contactstrugenalpha.h'.

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
  RCP<DRT::Discretization> cdiscret_;

  //! dof offset between pdicsret and cdiscret
  int dofoffset_;

  //! node and element maps
  RCP<Epetra_Map> noderowmap_;
  RCP<Epetra_Map> nodecolmap_;
  RCP<Epetra_Map> nodefullmap_;
  RCP<Epetra_Map> elerowmap_;
  RCP<Epetra_Map> elecolmap_;
  RCP<Epetra_Map> elefullmap_;

  //! vector of contact pairs (pairs of elements, which might get in contact)
  vector<RCP<Beam3contact> > pairs_;

  //! parameter list for contact
  Teuchos::ParameterList scontact_;

  //! search radius
  double searchradius_;

  //! contact forces of current time step
  RCP<Epetra_Vector> fc_;

  //! contact forces of previous time step (for generalized alpha)
  RCP<Epetra_Vector> fcold_;

  //! current constraint norm (violation of non-penetration condition)
  double constrnorm_;

  //! current penalty parameter (might be modified within augmented Lagrange strategy)
  double currentpp_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Search contact pairs

  We search pairs of elements that might get in contact and then create the
  corresponding beam3contact objects. Pairs of element that are direct neighbours,
  i.e. share one node, will be rejected.

  */
  void SearchPossibleContactPairs(map<int,LINALG::Matrix<3,1> >& currentpositions);

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
  by how much the penalty paramter is increased is NOT analytically	derived. It's just
  experience. It moght also be problem-dependent.

  */
	bool UpdateCurrentpp(const double& GlobalConstrNorm, const int& uzawaiter);

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
#endif  // #ifdef CCADISCRET
