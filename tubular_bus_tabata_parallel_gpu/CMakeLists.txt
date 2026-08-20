cmake_minimum_required(VERSION 3.16)

# ============================================================
# GPU 加速版本：CUDA MFEM + Hypre 2.32.0 + 批量 GPU 指针
# ★ 必须在 project() 之前设置 CUDA 编译器，否则 CMake 会找到 CUDA 10.1
set(CMAKE_CUDA_COMPILER /home/hxc1/local/cuda-12.9/bin/nvcc)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_ARCHITECTURES 75)
# ============================================================

project(Tabata_Main_Cylinder2D_monopolar LANGUAGES CXX CUDA)

# --- CUDA 工具包路径 ---
set(CUDA_DIR /home/hxc1/local/cuda-12.9)
set(CUDA_LIB_DIR ${CUDA_DIR}/targets/x86_64-linux/lib)
set(CUDA_INC_DIR ${CUDA_DIR}/targets/x86_64-linux/include)
include_directories(${CUDA_INC_DIR})

# --- MFEM (CUDA 版本, 基于 Hypre 2.32.0, 设备显存模式) ---
set(MFEM_SRC_DIR /home/hxc1/mfem-4.7)
set(MFEM_BUILD_DIR /home/hxc1/mfem-4.7/build_cuda_nouvm)
include_directories(${MFEM_SRC_DIR})                # mfem.hpp, 各模块头文件
include_directories(${MFEM_BUILD_DIR})              # config/_config.hpp
link_directories(${MFEM_BUILD_DIR})                 # libmfem.a

# --- Hypre (CUDA 版本, v2.32.0, 设备显存模式) ---
# GTX 1650 不支持硬件 UM 页迁移, 必须用设备显存模式避免 cudaMemPrefetchAsync 失败
set(HYPRE_DIR /home/hxc1/hypre-fresh/install)
include_directories(${HYPRE_DIR}/include)
set(HYPRE_LIBRARIES ${HYPRE_DIR}/lib/libHYPRE.a)

# --- Metis ---
set(METIS_DIR /home/hxc1/metis-4.0.3)
include_directories(${METIS_DIR}/Lib)
set(METIS_LIBRARIES ${METIS_DIR}/libmetis.a)

# --- MPI ---
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# --- 第三方头文件 ---
include_directories(${CMAKE_SOURCE_DIR}/third_party)

set(CMAKE_CXX_STANDARD 17)

# 使用 MFEM 源码头文件时, 必须手动定义以下宏 (通常由 config.hpp 提供)
add_definitions(-DMFEM_USE_CUDA)
add_definitions(-DMFEM_USE_MPI)
add_definitions(-DMFEM_USE_DOUBLE)
add_definitions(-DMFEM_USE_METIS)

# ============================================================
# 可执行文件
# ============================================================
set(EXECUTABLE_NAME Tabata_GPU)
set(TARGET_LIBRARIES
    mfem
    MPI::MPI_CXX
    ${HYPRE_LIBRARIES}
    ${METIS_LIBRARIES}
    # CUDA 运行时库 (MFEM CUDA 版本需要)
    ${CUDA_LIB_DIR}/libcusparse.so
    ${CUDA_LIB_DIR}/libcurand.so.10
    ${CUDA_LIB_DIR}/libcublas.so
    ${CUDA_LIB_DIR}/libcudart.so
    ${CUDA_LIB_DIR}/libcusolver.so
)

add_executable(${EXECUTABLE_NAME}
               Main_tubular_TABATA_parallel.cpp GetNodeE.cpp poisson_FEM_MFEM.cpp ComputeElementVolume.cpp
                ParallelUpwindElementFinder.cpp CVPDE_TABATA.cpp TabataAssemblyGPU.cu param_parser.cpp)

target_link_libraries(${EXECUTABLE_NAME} ${TARGET_LIBRARIES})

# --- 网格预分区工具 (纯 CPU, 不需要 CUDA kernel) ---
add_executable(partition_mesh
               tools/partition_mesh.cpp)
target_link_libraries(partition_mesh ${TARGET_LIBRARIES})

# GMSH 分区文件 → MFEM 并行格式 (坐标握手, 不需 $PartitionedEntities)
add_executable(gmsh_partition_to_mfem
               tools/gmsh_partition_to_mfem.cpp)
target_link_libraries(gmsh_partition_to_mfem ${TARGET_LIBRARIES})

# 确保运行时能找到 CUDA 共享库
target_link_directories(${EXECUTABLE_NAME} PRIVATE ${CUDA_LIB_DIR})

## 调试用 (需要时取消注释)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")