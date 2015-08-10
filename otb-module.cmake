set(DOCUMENTATION "Mosaic generation of remote sensing images")

otb_module(Mosaic
  DEPENDS
    OTBBoost
    OTBCommon
	OTBApplicationEngine
	OTBConversion
  TEST_DEPENDS
    OTBTestKernel
    OTBCommandLine
  DESCRIPTION
    "Mosaic images"
)
