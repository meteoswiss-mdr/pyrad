======================================================
SkewT -- Atmospheric Profile Plotting and Diagnostics
======================================================

SkewT provides a few useful tools to help with the plotting and analysis of 
upper atmosphere data. In particular it provides some useful classes to 
handle the awkward skew-x projection.

News
====

**18 June 2015. I have been delaying creating a new release for a long time, 
but here it is! What's new in SkewT version 1.1.0?**

* Fixed some bugs in CAPE calculations: we were mixing ambient 
  temperature-derived significant levels with virtual temperature  
  formulation for CAPE, which was producing some weird results in special 
  cases.
* Fixed plotting and CAPE calculations for soundings that had a 
  tropopause above 100hPa
* A few minor plotting issues, in particular things look better now if you 
  want to plot a different temperature or pressure range (use tmax/tmin 
  pmax/pmin kwargs to ``make_skewt_axes()``).
* Responses to bugs raised by community (Thanks dudes).

**SkewT has undergone a major overhaul in order to implement the new 
features available in Matplotlib 1.4. It also has a bunch of new features 
that we have been meaning to implement since inception. We hope you like it 
more than ever!**

Important Notice
================
* Since version 1.0, SkewT has explicitly included the new SkewX classes 
  that are showcased on the matplotlib website: 
  http://matplotlib.org/mpl_examples/api/skewt.py
  This stuff is completely fundamental to SkewT.py and we are greatful to 
  Ryan May from Unidata for providing this to the Python community.

Sounding Data
=============

The easiest way to get some sounding data is to visit the University of 
Wyoming's website:

http://weather.uwyo.edu/upperair/sounding.html

To get some sounding data, follow the link and find the date, and location 
you are interested in, view the data as a text file and just save the file 
to your system. If you want to get loads of data please be considerate about 
the way you go about doing this! (Lots of wget requests makes the server 
admins unhappy).

You can also pass your own data to SkewT by writing a text file in 
*identical* format to the University of Wyoming files, which are 
constant-width columns. Here's a sample of the first few lines of one of the 
bundled examples::

    94975 YMHB Hobart Airport Observations at 00Z 02 Jul 2013

    -----------------------------------------------------------------------------
       PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
        hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
    -----------------------------------------------------------------------------
     1004.0     27   12.0   10.2     89   7.84    330     14  284.8  306.7  286.2
     1000.0     56   12.4   10.3     87   7.92    325     16  285.6  307.8  286.9
      993.0    115   12.8    9.7     81   7.66    311     22  286.5  308.1  287.9


Alternatively you can create a dictionary with the column headers as keys 
and the data as 1D python arrays (preferably use ``ma.masked_array``). 
There's more about this under the "Running SkewT" section below.

Installing SkewT
================
We recommend that you download the tarball (big green button on this page) 
and run::

    python setup.py install

If you want to put it somewhere different than your system files, you can do::
    
    python setup.py install --prefix=/path/to/local/dir

Just remember, if you use a non-standard location you'll have to tell python 
about where you install it. An easy way to do this is to add the environment 
variable ``PYTHONPATH`` to your ``bashrc`` (you can read about this 
elsewhere).

Because SkewT is written purely in python, you don't even have to install it 
try it out! Just download the tarball and extract it somewhere convenient, 
and navigate to SkewT/skewt, and everything you need is right there.

You can also install using the package manager, but I have had some 
complaints about dependency issues (All you should need is matplotlib and 
numpy).

Running SkewT
=============

There are three basic ways to run SkewT. You can execute it from the command 
line with a text file name as an argument, or you can import it as a module 
and pass it a text file name, or you can pass it data directly.

Running from command line
-------------------------

From the command line (navigate to SkewT/skewt) you can type::

    python SkewT.py /path/to/sounding_filename.txt

What you'll get is all of the default settings. If you do this with the 
bundled example in ``SkewT/skewt/examples/bna_day1.txt``, you'll get this 
`graphical output 
<http://users.monash.edu.au/~tchubb/SkewT_examples/bna_day1_default.png>`_.

If that's what you want, well and good, but if you want to tweak things like 
the colours, read on...

Import SkewT as a module
------------------------

Assuming you have installed the package on your system and your sounding 
file is in your working directory, typical usage of SkewT could look like 
this (I use ``IPython``)::

    In [1]: from skewt import SkewT
    In [2]: S=SkewT.Sounding("bna_day1.txt")
    In [3]: S.plot_skewt(color='r')


The function ``plot_skewt()`` is a wrapper for a bunch of other functions. 
This will give you exactly the same plot as running SkewT from the command 
line, but you have immediate access to all of the ``matplotlib`` plot 
options for the profile traces and the barbs, but you don't get any control 
over anything else.

The full sequence of commands to get what ``plot_skewt`` wraps is this::

    In [1]: S.make_skewt_axes(tmin=-40.,tmax=30.,pmin=100.,pmax=1050.)
    In [2]: S.add_profile(color='r',bloc=0)
    In [3]: parcel=S.get_parcel()
    In [4]: S.lift_parcel(*parcel)

You don't have to put the ``tmin`` and other keyword arguments in to 
``make_skewt_axes()`` unless you want to plot against different values from 
the defaults shown here. The keyword argument ``bloc`` stands for ''barb 
location'' and allows you to shift the wind barbs to the left or right. This 
is handy if you want to plot multiple profiles on the one Skew-T diagram, 
for example, to compare today's and yesterday's soundings::

    In [1]: S=SkewT.Sounding("./skewt/examples/bna_day1.txt")
    In [2]: T=SkewT.Sounding("./skewt/examples/bna_day2.txt")
    In [3]: S.make_skewt_axes()
    In [4]: S.add_profile(color='r',bloc=0)
    In [5]: S.soundingdata=T.soundingdata      # replace the sounding data in S with that from T                      
    In [6]: S.add_profile(color='b',bloc=1)

Import as a module and run with your own data
---------------------------------------------

Got sounding data from another source? Want to make Skew-T diagrams of model 
output? Look no further. All you need to do is define a python dictionary 
like so::

    In [1]: mydata=dict(zip(('hght','pres','temp','dwpt'),(height_m,presssure_pa,temperature_c,dewpoint_c))) 
    In [2]: S=SkewT.Sounding(soundingdata=mydata)

At a minimum we require ``pres``, ``temp`` and ``dwpt`` to generate the 
profile traces, and ``hght`` is required for parcel calculations (although a 
future implementation will use a hydrostatic atmosphere assumption). The other 
keys accepted are those listed in the University of Wyoming sounding data 
header above.

Parcel Ascent
=============

As of version 1.0, SkewT has a full parcel ascent routing including 
automatic parcel definitions and CAPE/CIN and significant level 
calculations.

Automatic Parcel Definition
---------------------------

There are three standard parcel definitions used in predicting severe 
weather (see http://www.spc.noaa.gov/sfctest/help/sfcoa.html):

* Surface Based (``'sb'``): The surface conditions. Found by taking the 
  lowest level where all data is available. This may not represent the 
  convective potential of the sounding very well but is commonly used.
* Mixed Layer (``'ml'``): A parcel representing the mean potential energy in 
  the lowest 100-mb of the atmosphere. Found by averaging potential 
  temperature and the water vapour mixing ratio.
* Most Unstable (``'mu'``): The most unstable parcel of air found within the 
  lowest 300-mb of the atmosphere. Found by calculating CAPE for conditions 
  at all levels in the sounding data, and determining the equivalent surface 
  parcel by adiabatic descent. (Note: if CAPE is 0 for all levels this
  routine defaults to the surface based parcel)

To calculate one of these parcels for your sounding, use the 
``get_parcel()`` routine, which is a wrapper for ``surface_based_parcel()``, 
``mixed_layer_parcel()`` and ``most_unstable_parcel()``. Optionally pass it 
the parcel type you want (default is ``'mu'``)::

    In [1]: S=SkewT.Sounding("./skewt/examples/bna_day1.txt")
    In [2]: parcel=S.get_parcel('mu',depth=300)
    In [3]: parcel
    Out[3]: (1000.0, 23.037, 13.626, 'mu')
    In [4]: S.lift_parcel(*parcel_2)

Or, you can define your own parcel (the fourth item is just some text which 
appears on the Skew-T diagram)::

    In [5]: parcel_2=(1000.0, 25.0, 18, 'user')
    In [6]: S.make_skewt_axes(); S.add_profile(); 
    In [7]: S.lift_parcel(*parcel_2)

CAPE/CIN calculation
--------------------

Definitions in this section are based on Markowsi and Richardson (2010).

The ``lift_parcel()`` routine above is a wrapper for the ``get_cape()`` 
routine, but it also handles the graphics. The ``get_cape()`` routine, by 
itself, will calculate significant levels and CAPE/CIN::

    In [8]: P_lcl,P_lfc,P_el,CAPE,CIN=S.get_cape(*parcel)
    In [9]: print P_lcl,P_lfc,P_el,CAPE,CIN
    870.560154927 859.695806371 382.117602258 427.793216382 -8.64938413185

    In [10]: P_lcl,P_lfc,P_el,CAPE,CIN=S.get_cape(*parcel_2)
    In [11]: print P_lcl,P_lfc,P_el,CAPE,CIN
    902.773891386 902.773891386 178.058628014 2540.55724083 0.0

``get_cape()`` complains a bit if there are any dew point temperatures 
missing in the profile, but its default behaviour is to fill these with the 
minimum dewpoint in the column, and this will have a minimal effect on the 
CAPE calculation. 

The lifted condensation level (LCL) is found by solving for the intersection 
of the temperature for dry adiabatic ascent for the parcel, and a line of 
constant water vapour mixing ratio.

To find the level of free convection (LFC), the parcel is lifted along a 
moist adiabat from the LCL. For details, please see the ``moist_ascent()`` 
routine in ``SkewT.py``. All intersections of the parcel temperature and the 
environmental temperature are identified. Strictly speaking, all such levels 
are `equilibrium levels`. There are basically three possible scenarios:

* Parcel cooler than environment at LCL and no equilibrium levels: There are 
  no unstable levels in the profile above the LCL, so the LFC does not 
  exist.
* Parcel warmer than environment at LCL: This means that LFC=LCL, and there 
  must be at least one stable equilibrium level, which could be as high as 
  the tropopause.
* Parcel cooler than environment at LCL and at least two equilibrium levels: 
  This means that the parcel is initially stable at the LCL, but further 
  lifting will bring it to a condition where it becomes unstable. The LFC is 
  defined as the first point at which this occurs.

The term `Equilibrium Level` (EL) is often used to describe the first 
*stable* equilibrium level above the LFC, if this exists. Once the LCL, LFC 
and EL have been defined, we can calculate the Convective Available 
Potential Energy (CAPE) and Convective Inhibition::

    CAPE=trapz(9.81*(tparcel-tempenv)/tempenv,hght)

This expression only applies to the region where ``T_parcel>T_environment`` 
between the LFC and the EL. ``trapz`` is a basic trapezoidal integration 
routine from ``numpy``.` Similarly for CIN::

    CIN=trapz(9.81*(tparcel-tempenv)/tempenv,hght)

Which applies to the region where ``tparcel<=tempenv`` between the surface 
and the EL.

The example above (``bna_day1.txt``) is a perfect demonstration of why this 
behaviour might not be desirable. Using the `textbook 
<http://users.monash.edu.au/~tchubb/SkewT_examples/bna_day1_textbookcape.png>`_ 
definition (i.e. ``totalcape=False``) of the EL, you get practically no 
CAPE, but it's clear that there is a large layer of instability aloft. 
However, if you define the highest equilibrium level as the EL (i.e. 
``totalcape=True``), you get an answer that is more `representative 
<http://users.monash.edu.au/~tchubb/SkewT_examples/bna_day1_totalcape.png>`_ 
of the conditions of the day.

The keyword argument ``totalcape`` lets you override the default definition 
of the so-called 'Equilibrium Level,' (EL) which I took from Markowsi and 
Richardson (2010, p. 33): "The `equilibrium level` is defined to be the 
height at which a buoyant lifted parcel becomes neutrally buoyant, that is, 
the height above the LFC at which the parcel temperature is equal to the 
environmental temperature."
 
Working Examples
================
We have bundled in a set of example soundings in the ``SkewT/skewt/examples 
directoy``. You can run them like this::

    $ python SkewT.py example1

Substitute digits 1-4 to get the different examples. The code for these is 
right down the end of the SkewT.py file so you can have a look and play 
around with them if you want without affecting how SkewT works on import.

* Example 1: Two soundings from Hobart that I used to develop al ot of the 
  initial code base
* Example 2: Total CAPE vs. Textbook CAPE
* Example 3: Some severe weather events in Australia, with automatic parcel 
  definitions.
* Example 4: Use of custom parcels
* **Example 5 (new in v1.1.0): High tropopause sounding**

The sounding files and output graphics for the examples are all hosted `here 
<http://users.monash.edu.au/~tchubb/SkewT_examples/>`_.


To-Do List
==========
* More column diagnostics.
* Hodographs? Anyone? 

Contributors
============
* Ross Bunn from Monash University is actively developing and finding all my 
  warty bugs.
* Gokhan Sever (North Carolina) is a keen user and has been encouraging me 
  to add more stuff. It's thanks to him that I have finally implemented the 
  CAPE routines.
* Simon Caine.
* Hamish Ramsay (Monash) has promised to at least think about adding some 
  extra diagnostics.
* Holger Wolff as tester


Thanks for your interest in this package and I'd love to hear your feedback: 
thomas.chubb AT monash.edu
