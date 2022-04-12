# TODO add library
add_executable(example main.cpp)
target_link_libraries(example PRIVATE pybind11::embed)
