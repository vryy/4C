/*!----------------------------------------------------------------------
\file statmech.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_H
#define STATMECH_H

#include "../drt_structure/strugenalpha.H"


// forward declarations

/*!
\brief time integration for structural problems with statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/
class StatMech : public virtual StruGenAlpha
{
public:

  //! \brief Standard Constructor
  explicit StatMech(ParameterList& params,
                   DRT::Discretization& dis,
                   LINALG::Solver& solver,
                   IO::DiscretizationWriter& output);

  //! \brief Destructor
  virtual ~StatMech() {};

  //! \brief Do time integration
  void Integrate();

  //! \brief do full Newton iteration
  void FullNewton();
    
  //! \brief update configuration
  virtual void Update();

  //! \brief output to file/screen
  virtual void Output();



protected:

  // don't want = operator and cctor
  StatMech operator = (const StatMech& old);
  StatMech(const StatMech& old);


}; // class StatMech

#endif  // #ifndef STATMECH_H
#endif  // #ifdef CCADISCRET
