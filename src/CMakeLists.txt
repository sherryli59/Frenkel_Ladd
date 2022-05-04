add_library(simulation 
atom.cpp 
box.cpp 
conf.cpp
mc.cpp 
utils.cpp 
lj_potential.cpp 
spring_potential.cpp
hs_potential.cpp
mixed_potential.cpp)
target_include_directories(simulation PUBLIC ../includes)
target_compile_features(simulation PUBLIC cxx_std_11)
