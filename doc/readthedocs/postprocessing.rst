
.. _postprocessing:

Postprocessing
----------------

BACI supports traditional text file output and more recent binary
output. The later needs to be postprocessed by some conversion filter before it can
be viewed.


Conversion to readable formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to make the simulation results accessible, they need to be converted using the ``post_processor`` script located in your BACI build directory. To have the executable ``post_processor`` ready, one must compile BACI using ``make post`` in addition to ``make`` or ``make full`` command. As per default, running this script on the result (.control) file of your simulation, converts the results to the Ensight format. Alternatively, one can use additional option ``--filter=ensight``.

The ``post_processor`` script has a lot of options to specify the post processing details. A complete list of these options is printed to the screen when executing ``./post_processor --help`` . Most of these options are self-explanatory. Other ouput format filters are available via the ``--filter`` option. Here the vtu filter is of particular interest when you want to convert the results to VTK file format.

Note: When using the BACI post processing script ``post_processor``, be aware that derived quantities like e.g. stresses and strains are not automatically extracted from the simulation results! Presumed you set the proper flags in the \*.dat file (i.e. the quantities were actually calculated), you still have to set the ``--stress`` etc. options of the post processing script to extract these quantities from the simulation results.

Available post-processing options (output from ``post_processor --help``) are given in the following table:

+----------------+----------+---------+-------------------------------------------------+
| Parameter      |value type|default  |Explanation                                      |
+================+==========+=========+=================================================+
|--filter        |string    |"ensight"|"ensight"|"gid"|"vtu"                            |
+----------------+----------+---------+-------------------------------------------------+
|--start         |int       |0        |                                                 |
+----------------+----------+---------+-------------------------------------------------+
|--end           |int       |-1       |                                                 |
+----------------+----------+---------+-------------------------------------------------+
|--step          |int       |1        |                                                 |
+----------------+----------+---------+-------------------------------------------------+
|--output        |string    |"xxx"    |output file name                                 |
+----------------+----------+---------+-------------------------------------------------+
|--stress        |string    |""       |cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123\_nd123|
+----------------+----------+---------+-------------------------------------------------+
|--strain        |string    |""       |see stress                                       |
+----------------+----------+---------+-------------------------------------------------+
|--mortar        |string    |"no"     |"no" \|"yes"                                     |
+----------------+----------+---------+-------------------------------------------------+
|--optquantity   |string    |""       |cxyz, ndxyz, cxyz\_ndxyz                         |
+----------------+----------+---------+-------------------------------------------------+
|--heatflux      |string    |""       |see stress                                       |
+----------------+----------+---------+-------------------------------------------------+
|--tempgrad      |string    |""       |see stress                                       |
+----------------+----------+---------+-------------------------------------------------+
|--structvelacc  |string    |"no"     |"no" \| "yes"                                    |
+----------------+----------+---------+-------------------------------------------------+
|--rotation      |string    |"no"     |"no" \| "yes"                                    |
+----------------+----------+---------+-------------------------------------------------+
|--structumatdisp|string    |"no"     |"no" \| "yes"                                    |
+----------------+----------+---------+-------------------------------------------------+
|--outputtype    |string    |         |bin  \| ascii (only for vtu)                     |
+----------------+----------+---------+-------------------------------------------------+
 
For the respective filters, there are wrappers available:

+----------------+-------------------------------+
|wrapper         |is equal to                    |
+================+===============================+
|post_drt_ensight|post_processor --filter=ensight|
+----------------+-------------------------------+
|post_drt_gid    |post_processor --filter=gid    |
+----------------+-------------------------------+
|post_drt_vti    |post_processor --filter=vti    |
+----------------+-------------------------------+
|post_drt_vtu    |post_processor --filter=vtu    |
+----------------+-------------------------------+


Process output steps from ``start`` to ``end`` every ``step``. Works on
real time steps, steps not written by BACI are counted, too. Both
``start`` and ``end`` can be empty, in which case the filter will
process from the first and to the last step, respectively.

post_monitor
~~~~~~~~~~~~

generate gnuplot files for values at selected dofs:

::

   ./post_monitor [options] control-file monitor-file

with a monitor-file in control-file format consisting of at least one
block:

::

   monitor:

      field = 'fluid‘

      field_pos = 0

      discretization = 0

      node = 440

      group = ‘velocity‘

      dof = 0

      dof = 1

Here the three variables ``field``, ``field_pos`` and ``discretization``
specify the discretization that is to be used. ``node`` gives the global
node id. ``group`` gives the result name that should be monitored. And
any number of ``dof`` variables (in the example two dofs are set, no
mistake!) gives the dofs at the node those values are to be plotted.

ParaView
~~~~~~~~~~

Paraview <https://www.paraview.org> can read various post processing data formats. 
The most interesting ones are the *ensight* format and the *vtk/vtu* format. 
Since BACI is writing its data in a proprietary hdf5 format, 
we have to convert the data first, for which scripts are readily available.

By applying the filter ``post_drt_ensight`` to your binary BACI result
data, a *.case* file is created. This result format can then be
loaded into *paraview* for visualization. Start paraview with the
command ``paraview`` A number of ParaView tutorials are available in the internet:

::

   http://www.paraview.org/Wiki/The_ParaView_Tutorial

.. ifconfig:: institution in ("lnm", )

   Ensight
   --------

   By applying the filter ``post_drt_ensight`` (see next section) to your
   binary BACI result data, a *\*.case* file is created. This result
   format can then be loaded into *Ensight* for visualization. 

   Start Ensight with the command ``/lnm/programs/CEI/bin/ensight8``



Animations
~~~~~~~~~~

The ultimate goal of scientific research is a beautiful movie!

There are several way to create animations using BACI output files.
Movies should be playable accross platforms (at least Linux and
Windows?) and embeddable inside MS Powerpoint presentations without
the need of having different movie versions in different formats. My
newest finding: it seems that ffmpeg (available for Debian through apt)
can simplify that process without the need to install any additional
codecs on Windows (nice for conferences when other peoples laptops/PCs?
have to be used). The previous guide using XviD will become obsolete but
remains here, until enough experience with the new encoding process
could be gathered.

**Animations from GiD**

Using GiD Postprocessing, one is able to create MPEG2 Movies but these
are very large, the quality is not good and don’t play in PowerPoint
(They play on Windows in the MS MediaPlayer, though). Another simple
way is to use avi/mjpeg. The quality is generally good, but they are of
huge size. recompress them as described below. The best solution when
using GiD is to create individual pictures and encode them afterwards.
Output the images in the animation dialog using the uncompressed TIFF
format.


**Encoding animations using the ffmpeg encoder**

*Encoding an MPEG2 movie from a different format using ffmpeg*

ffmpeg can read a lot of video sources, so most likely, it will read
your in.avi or in.mpg just fine.

::

   ffmpeg -i in.avi -sameq -b 6400 out.mpg
   ffmpeg -i in.mpg -sameq -b 6400 out.mpg

If the initial .avi file has a framerate lower than 25 (see error
message), use the -r option to force 25 frames/s in the output MPEG2
movie (MPEG2 standard is 25 frames/s) with

::

   ffmpeg -i in.avi -r 25 -b 6400 out.mpg

Note that the movie speed won’t change.

The bitrate option -b is described below.

*Animations from several image files using ffmpeg*

Postprocessing from GiD, post_visual2 can provide a series of image
files, hopefully numbered in a consistent order (For weird and stupid
GiD numbering: Axel has a python script to start from)

Providing the image files

Produce a series of images consistently numbered as test0001.jpg,
test0002.jpg, ..., test0152.jpg. If you don’t have the leading zeros,
the order of the images in the movie will be wrong (1, 11, 12, 13, 14,
15, 16, 17, 18, 19, 2, 20, 21... you get the idea?).

ffmpeg can encode directly from PNG images, consequently, they are
prefered because of their lossless image compression. To convert other
formats into the PNG format, use a shell script such as

::

           for i in *.tiff ; do
           echo $i
           convert $i -depth 24 `basename $i .tiff`.png
           done

Encoding the MPEG2 movie from PNG files:

::

   ffmpeg -i output_%05d.png -sameq -b 6400 out.mpg

using bitrate 6400 results in high quality movies (note the quality
indicator q= output during encoding. q=2.0 seems to be the highest
possible value here). In practice, ffmpeg reduces bitrate when q=2.0 is
reached and a lower bitrate is used depending on the images content.

The -sameq flag here is important! It tells ffmpeg to use the same
quality as the input, which means highest possible quality if the input
is looseless png. This is the way to create movies!

*Encoding the MPEG2 movie from the PNG files at a lower speed*

This is achieved by using less frames per second (e.g. 12.5 frames/s)
for the input. Note that the low framerate is given before! the input
files which means that the input has 12.5 frames/s. An MPEG2 movies
always has 25 frames/s, which now has to be given explicitly

::

   ffmpeg -r 12.5 -i output_%05d.png -r 25 -b 6400 out.mpg

The quality of the resulting movie strongly depends on allowed bitrate,
the quality of the initial image files/movie file and the content of the
images. Pictures with lots of features, e.g. showing the FE-mesh, a more
likely to become blury. Read about the ffmpeg parameters to improve the
quality as needed.

Behaviour of MPEG2 movies created by ffmpeg on Windows and in
PowerPoint

The created MPEG2 movie files will play on Windows and in PowerPoint
without any additional codec installation.

*Encoding animations using mencoder and XviD*

Encoding an XviD movie from a different format

::

   mencoder old.avi -ovc xvid -oac mp3lame -o new.avi

*Animations from several image files using Mencoder*

Postprocessing from GiD or post_visual2 can also provide a series of
image files, hopefully numbered in a consistent order. (For weird and
stupid GiD numbering: Axel has a python script)

Providing the image files

The procedure for numbering is the same as above. However, mencoder only
takes JPG files which can be produced with a shell script as

::

           for i in *.tif ; do
           echo $i
           convert $i -quality 100 -depth 24 `basename $i .tif`.jpg
           done

Encoding the XviD movie from the JPEG files:

::

   mencoder "mf://*.jpg" -o new.avi -ovc xvid -xvidencopts fixed_quant=4

*Encoding the XviD movie from the JPEG files at a lower speed*

This is achieved by using less frames per second

::

   mencoder "mf://*.jpg" -mf fps=12.5 -o new.avi -ovc xvid -xvidencopts fixed_quant=4

Possible options to improve quality (see "man mencoder" or search the
web for more details)

::

   -xvidencopts fixed_quant=4
   -xvidencopts me_quality=0
   -xvidencopts quant_type=mpeg
   -xvidencopts hq_ac
   -xvidencopts vhq=4
   -xvidencopts notrellis

The quality of the resulting movie strongly depends on the above
parameters, the quality of the initial JPEG files and, of course the
content of the images. Read about the XviD parameters to improve the
quality as needed.

XviD movies on Windows and in PowerPoint

If the steps give above are followed, the XviD encoded movie file will
play on Windows and in PowerPoint. Make sure you have installed the
XviD codecs on the Windows PC or Laptop. See XviD for further
information on installation.

Encoding animations using mencoder and the msmpeg codec

This way works, but requires multiple steps to make the movie play in
PowerPoint?. Choose yourself.

Movies for PowerPoint

Generate movies that Microsoft can read:

::

   mencoder "mf://*.jpg" -mf fps=12.5 -o new.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq

Use the Windows Movie Maker: Import movie and export it again. The
result (\*.wmv) can be used by PowerPoint?.
