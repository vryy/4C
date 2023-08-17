.. _free_surface_flow:

old_Free Surface Flow
=====================

Free surface flow is a special problem type:

::

   ----------------- PROBLEM TYP
   PROBLEMTYP  Fluid_FreeSurface

.. note::

   The rest of this is probably outdated in the meantime.

Genkinger – Veröffentlichung mit 2d und 3d Beispiel Fluid Fall mit 2
Feldern Verfahren partioned implicit local Lagrange vertical
heightfunction Oberflächenspannung nur 2d implementiert, wahrscheinlich
kaputt durch mass-rhs (Im Rahmen der Umstellung auf die Verwendung der
Massen-RHS wird etforce abgeschafft. Die korrekte Integration der
Oberflaechenspannungseinfluesse muss ueberprueft werden! chfoe 02/05)

Tests
~~~~~

There are a number of working tests:


- f2_freeosz20x20_drt.dat
- fs_kugel_tet.dat
- fs_quadrat.dat
- fs_wuerfel_tet.dat