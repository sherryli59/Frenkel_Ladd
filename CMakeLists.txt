cmake_minimum_required(VERSION 3.0.0)
project(MC VERSION 0.1.0)

include(CTest)
enable_testing()

# add_executable(${PROJECT_NAME} main.cpp)
add_executable(lj_phase1 lj_phase1.cpp)
add_executable(lj_phase2 lj_phase2.cpp)
add_executable(hs hs.cpp)


#target_compile_features(phase1 PUBLIC cxx_std_11)

add_subdirectory(src/)
# target_link_directories(${PROJECT_NAME} PUBLIC src/)
# target_link_libraries(${PROJECT_NAME} simulation)
target_link_directories(lj_phase1 PUBLIC src/)
target_link_libraries(lj_phase1 simulation)
target_link_directories(lj_phase2 PUBLIC src/)
target_link_libraries(lj_phase2 simulation)
target_link_directories(hs PUBLIC src/)
target_link_libraries(hs simulation)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
