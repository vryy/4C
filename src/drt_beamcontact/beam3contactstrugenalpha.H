/*!----------------------------------------------------------------------
\file beam3contactstrugenalpha.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BEAM3CONTACTSTRUGENALPHA_H
#define BEAM3CONTACTSTRUGENALPHA_H

#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "../drt_structure/strugenalpha.H"
#include "beam3contact_manager.H"


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

  // beam contact manager
  RCP<CONTACT::Beam3cmanager> beamcmanager_;

}; // class Beam3ContactStruGenAlpha
} // namespace CONTACT

#endif  // #ifndef BEAM3CONTACTSTRUGENALPHA_H
#endif  // #ifdef CCADISCRET
