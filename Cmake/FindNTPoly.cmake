# Find the NTPoly Module
# Variables set:
# - NTPoly_FOUND - system found NTPoly.
# - NTPoly_LIBRARIES - the linker line for NTPoly.
# - NTPoly_INCLUDE_DIRS - the path to NTPoly.

# First we search for the libraries
find_library(NTPoly_LIBRARIES NTPoly)
find_path(NTPoly_INCLUDE_DIRS "psmatrixmodule.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NTPoly DEFAULT_MSG
  NTPoly_LIBRARIES NTPoly_INCLUDE_DIRS)
mark_as_advanced(NTPoly_LIBRARIES NTPoly_INCLUDE_DIRS)
