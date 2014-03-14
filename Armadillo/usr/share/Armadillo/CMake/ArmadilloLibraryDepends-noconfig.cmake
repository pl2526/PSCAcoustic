#----------------------------------------------------------------
# Generated CMake target import file for configuration "".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "armadillo" for configuration ""
SET_PROPERTY(TARGET armadillo APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
SET_TARGET_PROPERTIES(armadillo PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "/usr/lib/libblas.so;/usr/lib/liblapack.so;/usr/lib/libcblas.so;/usr/lib/liblapack_atlas.so"
  IMPORTED_LOCATION_NOCONFIG "/usr/lib64/libarmadillo.so.2.4.3"
  IMPORTED_SONAME_NOCONFIG "libarmadillo.so.2"
  )

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
