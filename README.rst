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
       [--min-distance-from-tracking-centre MIN_DISTANCE_FROM_TRACKING_CENTRE]
       [--add-custom-exclusion-zone ADD_CUSTOM_EXCLUSION_ZONE [ADD_CUSTOM_EXCLUSION_ZONE ...]]
       [--max-region-right-skewness MAX_REGION_RIGHT_SKEWNESS]
       [--psf-image PSF_IMAGE]
       [--remove-tagged-dE-components-from-model-images REMOVE_TAGGED_DE_COMPONENTS_FROM_MODEL_IMAGES]
       [--only-dEs-in-lsm]
       [--max-positive-to-negative-flux MAX_POSITIVE_TO_NEGATIVE_FLUX]
       [--max-region-abs-skewness MAX_REGION_ABS_SKEWNESS]
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
  --min-distance-from-tracking-centre MIN_DISTANCE_FROM_TRACKING_CENTRE
                        Cutoff distance from phase centre in which no tags be
                        raised.This can be used to effectively exclude the
                        FWHM of an parabolic reflector-based interferometer.
  --add-custom-exclusion-zone ADD_CUSTOM_EXCLUSION_ZONE [ADD_CUSTOM_EXCLUSION_ZONE ...]
                        Add manual exclusion zone to which no dE tags shall be
                        added. Expects a tripple of centre X, Y pixel and
                        radius.
  --max-region-right-skewness MAX_REGION_RIGHT_SKEWNESS
                        The maximum tolerance for right skewness of a pixel
                        distribution within a region.A large value (tailed
                        distribution) indicates significant uncleaned flux
                        remaining in the residual. This can be used to
                        effectively control detection sensitivity to uncleaned
                        extended emission, but should be set to 0 if residuals
                        other than stokes I are used
  --psf-image PSF_IMAGE
                        PSF image from which BPA, BMAJ and BMIN may be
                        extracted
  --remove-tagged-dE-components-from-model-images REMOVE_TAGGED_DE_COMPONENTS_FROM_MODEL_IMAGES
                        Blank out model images within resolution of tagged LSM
                        components. Expects list of model FITS files. This
                        option is useful for hybrid DFT-CLEAN component
                        modelling as onlyextended / faint clean components
                        contributes to model.
  --only-dEs-in-lsm     Only store dE tagged sources in lsm. This option is
                        useful for hybrid DFT-CLEAN component modelling, as
                        only bright compact gaussian emission contributes to
                        dE solutions
  --max-positive-to-negative-flux MAX_POSITIVE_TO_NEGATIVE_FLUX
                        The maximum tolerance for the ratio of positive to
                        negative flux. Only to be used with stokes I
  --max-region-abs-skewness MAX_REGION_ABS_SKEWNESS
                        The maximum tolerance for absolute skewness of a pixel
                        distribution within a region.A large value (tailed
                        distribution) indicates significant uncleaned flux
                        remaining in the residual. This can be used to
                        effectively control detection sensitivity to uncleaned
                        extended emission, but should be set to 0 if residuals
                        other than stokes Q,U or V are used

