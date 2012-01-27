/*!------------------------------------------------------------------------------------------------*
\file topopt_utils.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_UTILS_CPP_
#define TOPOPT_UTILS_CPP_

#include <iostream>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_scatra/scatra_element.H"
#include "topopt_utils.H"



void TOPOPT::printTopOptLogo()
{
  std::cout << "      _________            " << std::endl;
  std::cout << "     /         \\          " << std::endl;
  std::cout << "    /   _____   \\         " << " Das ist das       " << std::endl;
  std::cout << "   |   /     \\@ @\\        " << " Gebiets-          " << std::endl;
  std::cout << "   |   |      \\__/        " << " Optimierungsmodul " << std::endl;
  std::cout << "   |   \\         \\        " << " in BACI           " << std::endl;
  std::cout << "    \\   \\         \\        " << "                   " << std::endl;
  std::cout << "     \\   \\                " << " Die Schlange      " << std::endl;
  std::cout << "      \\   \\_________      " << " wird sich bald    " << std::endl;
  std::cout << "       \\             \\    " << " teilen und wieder " << std::endl;
  std::cout << "        \\__________   \\   " << " zusammenwachsen   " << std::endl;
  std::cout << "                   \\  |   " << " können!           " << std::endl;
  std::cout << "        _          |  |    " << std::endl;
  std::cout << "       | \\        /  /    " << std::endl;
  std::cout << "        \\ \\______/  /     " << std::endl;
  std::cout << "         \\_________/      " << std::endl;
  std::cout << "                           " << std::endl;
}

#endif /* TOPOPT_UTILS_CPP_ */
#endif  // #ifdef CCADISCRET
