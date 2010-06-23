/*!----------------------------------------------------------------------
\file beam3contactstrugenalpha.H

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BEAM3CONTACTSTRUGENALPHA_H
#define BEAM3CONTACTSTRUGENALPHA_H

#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_ana.H"
#include "../drt_io/io.H"
#include "../drt_structure/strugenalpha.H"

#ifdef D_BEAM3
#include "beam3contact_manager.H"
#endif // #ifdef D_BEAM3

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

// forward declarations

/*!
\brief Generalized Alpha time integration for structural problems with 3D beam contact

\author popp (popp@lnm.mw.tum.de), cyron (cyron@lnm.mw.tum.de)
*/
class Beam3ContactStruGenAlpha : public virtual StruGenAlpha
{
public:

  /*!
  \brief Standard Constructor

  */
  explicit Beam3ContactStruGenAlpha(ParameterList& params,
                                    DRT::Discretization& dis,
                                    LINALG::Solver& solver,
                                    IO::DiscretizationWriter& output);

  /*!
  \brief Destructor

  */
  virtual ~Beam3ContactStruGenAlpha() {};

  /*!
  \brief do constant predictor step

  */
  void ConstantPredictor();

  /*!
  \brief do consistent constant predictor step

  */
  void ConsistentPredictor();

  /*!
  \brief do full Newton iteration

  */
  void FullNewton();

  /*!
  \brief initialize full Newton iteration for Uzawa

  This method is called at the beginning of the Newton scheme in the second, third
  and so on Uzawa iterarion in order to create an of an out-of-balance force again.

  */
  void InitializeNewtonUzawa();

  /*!
  \brief integration control routine

  */
  void Integrate();

  /*!
  \brief do update at the end of time step

  */
  void Update();

protected:

  // don't want = operator and cctor
  Beam3ContactStruGenAlpha operator = (const Beam3ContactStruGenAlpha& old);
  Beam3ContactStruGenAlpha(const Beam3ContactStruGenAlpha& old);

  // uzawa iteration counter
  int uzawaiter_;

  // beam contact manager
#ifdef D_BEAM3
  RCP<CONTACT::Beam3cmanager> beamcmanager_;
#endif // #ifdef D_BEAM3

}; // class Beam3ContactStruGenAlpha
} // namespace CONTACT

#endif  // #ifndef BEAM3CONTACTSTRUGENALPHA_H
#endif  // #ifdef CCADISCRET
