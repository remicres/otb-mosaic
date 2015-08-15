otb_module_test()

otb_test_application(NAME   mosaic1Test1
                      APP  Mosaic
                      OPTIONS -il ${INPUTDATA}/DEM/srtm_directory/N46E008.hgt ${INPUTDATA}/Input/DEM/srtm_directory/N46E007.hgt ${INPUTDATA}/Input/DEM/srtm_directory/N46E006.hgt
                              -out ${TEMP}/mosaic1Test1_result.tif float
                              -output.spacing 0.01
                      VALID   --compare-image ${EPSILON_8}
                              mosaic1Test1_result.tif
                              ${TEMP}/mosaic1Test1_result.tif)
