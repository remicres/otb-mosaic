otb_module_test()

otb_test_application(NAME   mosaic1Test1
                      APP  Mosaic
                      OPTIONS -in1 image1.tif
                              -in2 image2.tif
                              -out ${TEMP}/mosaic1Test1_result.tif uint8
                      VALID   --compare-image ${EPSILON_8}
                              mosaic1Test1_result.tif
                              ${TEMP}/mosaic1Test1_result.tif)

otb_test_application(NAME   mosaic1Test2
                      APP  Mosaic
                      OPTIONS -in1 image1.tif
                              -in2 image2.tif
                              -harmo.method band
                              -out ${TEMP}/mosaic1Test2_result.tif uint8
                      VALID   --compare-image ${EPSILON_8}
                              mosaic1Test2_result.tif
                              ${TEMP}/mosaic1Test2_result.tif)

