# Grab all *.cpp files in the directory
file(GLOB BENCH_SRC CONFIGURE_DEPENDS "*.cpp")

foreach(BENCH_SRC_FILE ${BENCH_SRC})
  # Strip file extension
  get_filename_component(BENCH_NAME ${BENCH_SRC_FILE} NAME_WLE)
  add_executable(${BENCH_NAME} ${BENCH_SRC_FILE})
  target_link_libraries(${BENCH_NAME} PUBLIC rmat_generator argparse::argparse)
endforeach()