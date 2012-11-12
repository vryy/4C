/*!---------------------------------------------------------------------*
\file combust_utils.cpp

\brief collection of functions in namespace COMBUST

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "combust_utils.H"

#include "../drt_io/io_pstream.H"

#include <iostream>


/*------------------------------------------------------------------------------------------------*
 | print COMBUST module logo on screen                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::printCombustLogo()
{
    IO::cout << "     ___            ___    \n"
             << "    /   \\          /   \\ \n"
             << "    \\_   \\        /  __/ \n"
             << "     _\\   \\      /  /__  " << " Das ist               \n"
             << "     \\___  \\____/   __/  " << " das Verbrennungsmodul \n"
             << "         \\_       _/      " << " in BACI               \n"
             << "           | @ @  \\_      " << "                       \n"
             << "           |               " << " Der Elch wird bald    \n"
             << "         _/     /\\        " << " ein feuerspeiender    \n"
             << "        /o)  (o/\\ \\_     " << " Drache sein!          \n"
             << "        \\_____/ /         \n"
             << "          \\____/          \n"
             << "                           " << IO::endl;
}

