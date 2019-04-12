# Find the YAML Module from BigDFT.
# Variables set:
# - YAML_FOUND - system found YAML.
# - YAML_LIBRARIES - the linker line for YAML.
# - YAML_INCLUDE_DIRS - the path to YAML.

# First we search for the libraries
find_library(YAML_LIBRARIES yaml)
find_path(YAML_INCLUDE_DIRS "yaml_parse.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAML DEFAULT_MSG
  YAML_LIBRARIES YAML_INCLUDE_DIRS)
mark_as_advanced(YAML_LIBRARIES YAML_INCLUDE_DIRS)
