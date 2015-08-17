otb_module_test()

otb_test_application(NAME MosaicTestLargeFeathering
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestLargeFeathering.tif uint8 
                              -comp.feather large
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestLargeFeathering.tif
                              ${TEMP}/apTvMosaicTestLargeFeathering.tif)


otb_test_application(NAME MosaicTestSlimFeathering
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestSlimFeathering.tif uint8 
                              -comp.feather slim 
                              -comp.feather.slim.lenght 100
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestSlimFeathering.tif
                              ${TEMP}/apTvMosaicTestSlimFeathering.tif)


otb_test_application(NAME MosaicTestSimpleWithHarmoBandRmse
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestSimpleWithHarmoBandRmse.tif uint8 
                              -harmo.method band 
                              -harmo.cost rmse
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestSimpleWithHarmoBandRmse.tif
                              ${TEMP}/apTvMosaicTestSimpleWithHarmoBandRmse.tif)

otb_test_application(NAME MosaicTestSimpleWithHarmoRgbRmse
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestSimpleWithHarmoRgbRmse.tif uint8 
                              -harmo.method rgb 
                              -harmo.cost rmse
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestSimpleWithHarmoRgbRmse.tif
                              ${TEMP}/apTvMosaicTestSimpleWithHarmoRgbRmse.tif)

otb_test_application(NAME MosaicTestSimpleWithCutline
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestSimpleWithCutline.tif uint8 
                              -vdcut ${INPUTDATA}/SP67_FR_subset_1_cutline.shp ${INPUTDATA}/SP67_FR_subset_2_cutline.shp 
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestSimpleWithCutline.tif
                              ${TEMP}/apTvMosaicTestSimpleWithCutline.tif)

otb_test_application(NAME MosaicTestSimpleWithVdstats
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/SP67_FR_subset_1.tif ${INPUTDATA}/SP67_FR_subset_2.tif 
                              -out ${BASELINE}/apTvMosaicTestSimpleWithVdstats.tif uint8 
                              -vdstats ${INPUTDATA}/SP67_FR_subset_1_cutline.shp ${INPUTDATA}/SP67_FR_subset_2_cutline.shp 
                      VALID   --compare-image ${EPSILON_8}
                              ${BASELINE}/apTvMosaicTestSimpleWithVdstats.tif
                              ${TEMP}/apTvMosaicTestSimpleWithVdstats.tif)
