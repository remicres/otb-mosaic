cmake_minimum_required (VERSION 2.8.11)

OTB_CREATE_APPLICATION(NAME           Mosaic
                       SOURCES        otbMosaic.cxx
                       LINK_LIBRARIES OTBCommon OTBGdalAdapters)

