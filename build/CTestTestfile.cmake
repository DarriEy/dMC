# CMake generated Testfile for 
# Source directory: /Users/darrieythorsson/compHydro/code/dmc_design
# Build directory: /Users/darrieythorsson/compHydro/code/dmc_design/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(gradient_verification "/Users/darrieythorsson/compHydro/code/dmc_design/build/test_gradients")
set_tests_properties(gradient_verification PROPERTIES  _BACKTRACE_TRIPLES "/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;155;add_test;/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;0;")
add_test(single_reach_routing "/Users/darrieythorsson/compHydro/code/dmc_design/build/test_single_reach")
set_tests_properties(single_reach_routing PROPERTIES  _BACKTRACE_TRIPLES "/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;160;add_test;/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;0;")
add_test(bmi_interface "/Users/darrieythorsson/compHydro/code/dmc_design/build/test_bmi")
set_tests_properties(bmi_interface PROPERTIES  _BACKTRACE_TRIPLES "/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;165;add_test;/Users/darrieythorsson/compHydro/code/dmc_design/CMakeLists.txt;0;")
subdirs("_deps/codipack-build")
subdirs("_deps/nlohmann_json-build")
