cmake_minimum_required(VERSION 3.16)

project(Takuma_Main_ion_field)

include_directories(/home/hxc1/mfem-4.7/install/include) # 头文件路径
link_directories(/home/hxc1/mfem-4.7/install/lib) # 库路径

set(CMAKE_CXX_STANDARD 17)

# 【新增这一行】：引入第三方单头文件库路径
include_directories(${CMAKE_SOURCE_DIR}/third_party)

# 并行 mfem 依赖文件
# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Hypre
set(HYPRE_DIR /home/hxc1/local/hypre)
include_directories(${HYPRE_DIR}/include)
set(HYPRE_LIBRARIES ${HYPRE_DIR}/lib/libHYPRE.a)

# Metis
set(METIS_DIR /home/hxc1/metis-4.0.3)
include_directories(${METIS_DIR}/Lib)
set(METIS_LIBRARIES ${METIS_DIR}/libmetis.a)

# Define executable name and libraries
set(EXECUTABLE_NAME Takuma_DoubleStack_CPU)
set(TARGET_LIBRARIES
    mfem
    MPI::MPI_CXX
    ${HYPRE_LIBRARIES}
    ${METIS_LIBRARIES}
)

add_executable(${EXECUTABLE_NAME} 
               Main_tubular_TAKUMA_parallel.cpp param_parser.cpp  
                GetNodeE.cpp poisson_FEM_MFEM.cpp ComputeElementVolume.cpp
                ParallelUpwindElementFinder.cpp CVPDE_TAKUMA.cpp ) # 添加可执行程序

target_link_libraries(${EXECUTABLE_NAME} ${TARGET_LIBRARIES}) # 必须指定链接库

## 调试用
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")