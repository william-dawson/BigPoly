# Find the Futile Module
# Variables set:
# - Futile_FOUND - system found futile.
# - Futile_LIBRARIES - the linker line for futile.
# - Futile_INCLUDE_DIRS - the path to futile.

# First we search for the libraries
find_library(Futile_LIBRARIES futile-1)
find_path(Futile_INCLUDE_DIRS "dictionaries.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Futile DEFAULT_MSG
  Futile_LIBRARIES Futile_INCLUDE_DIRS)
mark_as_advanced(Futile_LIBRARIES Futile_INCLUDE_DIRS)
