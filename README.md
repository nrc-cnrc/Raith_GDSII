Raith_GDSII MATLAB toolbox
==========================

*Version 1.2*

The four class definitions in this toolbox provide a simple, versatile, and scriptable means of generating patterns for [Raith](http://www.raith.com) electron-beam lithography (EBL) and focused ion beam (FIB) tools using [MATLAB](http://www.mathworks.com/products/matlab/).

This toolbox was developed at the [National Institute for Nanotechnology](http://nint-innt.ca) (NINT), a joint initiative between the Government of Canada, the Government of Alberta, the [National Research Council Canada](http://www.nrc-cnrc.gc.ca) (NRC), and the [University of Alberta](http://www.ualberta.ca).  It is currently maintained by the University of Alberta [nanoFAB](http://nanofab.ualberta.ca) facility.


Features
--------

* Generate Raith-dialect GDSII hierarchy (.csf) and positionlist (.pls) files directly within MATLAB
* Plot patterns in MATLAB using Raith dose factor colouration, from individual GDSII elements to entire positionlists
* Full support for Raith curved elements (circles, ellipses, arcs)  
* Full support for Raith "fixed beam moving stage" exposure elements (paths and circles)  **(new in v. 1.2)**
* Simply-connected font defined to use for text elements
* Export pattern in plain GDSII (Raith curved elements converted to polygons and paths), to use with non-Raith GDSII editors


Requirements and installation
-----------------------------

This toolbox has been tested with MATLAB R2010a and later; use with earlier versions of MATLAB is at your own risk.

To install, simply place `Raith_element.m`, `Raith_structure.m`, `Raith_library.m`, `Raith_positionlist.m` in a folder on your MATLAB path.


Documentation
-------------

Please refer to `Raith_GDSII MATLAB Toolbox - User Guide.pdf` for a quick-start guide, tutorial examples, and a full description of the class definitions.


License
-------

This toolbox is subject to the terms of the [Mozilla Public License, v. 2.0](http://mozilla.org/MPL/2.0/). 