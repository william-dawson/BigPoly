# Find the Atlab Module
# Variables set:
# - Atlab_FOUND - system found Atlab.
# - Atlab_LIBRARIES - the linker line for Atlab.
# - Atlab_INCLUDE_DIRS - the path to Atlab.

# First we search for the libraries
find_library(Atlab_LIBRARIES atlab-1)
find_path(Atlab_INCLUDE_DIRS "numerics.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Atlab DEFAULT_MSG
  Atlab_LIBRARIES Atlab_INCLUDE_DIRS)
mark_as_advanced(Atlab_LIBRARIES Atlab_INCLUDE_DIRS)
