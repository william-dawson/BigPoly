# Find the CheSS Module
# Variables set:
# - CheSS_FOUND - system found CheSS.
# - CheSS_LIBRARIES - the linker line for CheSS.
# - CheSS_INCLUDE_DIRS - the path to CheSS.

# First we search for the libraries
find_library(CheSS_LIBRARIES CheSS-1)
find_path(CheSS_INCLUDE_DIRS "sparsematrix.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CheSS DEFAULT_MSG
  CheSS_LIBRARIES CheSS_INCLUDE_DIRS)
mark_as_advanced(CheSS_LIBRARIES CheSS_INCLUDE_DIRS)
