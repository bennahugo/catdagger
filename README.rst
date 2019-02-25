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

Usage
===============================================================================

dagger --help                                                                                              
usage: CATDagger - an automatic differential gain tagger (C) SARAO, Benjamin Hugo 2019
       [-h] [--stokes STOKES] [--min-tiles-region MIN_TILES_REGION]
       [--input-lsm INPUT_LSM] [--ds9-reg-file DS9_REG_FILE]
       [--ds9-tag-reg-file DS9_TAG_REG_FILE] [-s SIGMA]
       [--tile-size TILE_SIZE] [--global-rms-percentile GLOBAL_RMS_PERCENTILE]
       [--de-tag-name DE_TAG_NAME]
       noise_map

positional arguments:
  noise_map             Residual / noise FITS map to use for estimating local
                        RMS

optional arguments:
  -h, --help            show this help message and exit
  --stokes STOKES       Stokes to consider when computing global noise
                        estimates. Ideally this should be 'V', if available
  --min-tiles-region MIN_TILES_REGION
                        Minimum number of tiles per region. Regions with fewer
                        tiles will not be tagged as dE
  --input-lsm INPUT_LSM
                        Tigger LSM to recluster and tag. If this is not
                        specified only DS9 regions will be written out
  --ds9-reg-file DS9_REG_FILE
                        SAODS9 regions filename to write out
  --ds9-tag-reg-file DS9_TAG_REG_FILE
                        SAODS9 regions filename to contain tagged cluster
                        leads as circles
  -s SIGMA, --sigma SIGMA
                        Threshold to use in detecting outlier regions
  --tile-size TILE_SIZE
                        Number of pixels per region tile axis
  --global-rms-percentile GLOBAL_RMS_PERCENTILE
                        Percentile tiles to consider for global rms
                        calculations
  --de-tag-name DE_TAG_NAME
                        Tag name to use for tagged sources in tigger LSM
