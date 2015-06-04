project(tan3DExpr)
cmake_minimum_required(VERSION 2.8)
include_directories(../..)
aux_source_directory(. SRC_LIST)
add_definitions("-std=c++0x")
find_package(DGtal 0.9 REQUIRED)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES})
