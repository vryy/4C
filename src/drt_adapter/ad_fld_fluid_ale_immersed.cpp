/*----------------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ale_immersed.cpp

\brief

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            089 - 289 10240
</pre>
*/
/*----------------------------------------------------------------------------*/
#include "ad_fld_fluid_ale_immersed.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::FluidAleImmersed::FluidAleImmersed(const Teuchos::ParameterList& prbdyn,
                                            std::string condname)
:FluidAle(prbdyn,condname)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAleImmersed::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  FluidField()->AddDirichCond(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAleImmersed::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  FluidField()->RemoveDirichCond(maptoremove);
}

