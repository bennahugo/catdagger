**CATDagger**
==============================================================================
A catalog source differential gain tagger based on local noise characteristics

This tool segments regions within residual images that are in need of a differential gain. Preferably the tool is run on stokes V
residuals, which typically contain relatively little real flux and mostly residual calibration errors. In principle it can also be run on Stokes I residuals
if direction independent calibration was successful.

DS9 region maps containing regions and cluster lead information is output by default as shown as example below. Tigger LSM catalogs
can simultaniously be processed and reclustered based on identified dE regions.

.. figure:: https://github.com/bennahugo/catdagger/blob/master/misc/catdagger.png
    :width: 250px
    :height: 250px
    :align: center
