cmake_minimum_required(VERSION 3.16)
project(IonFlowSimulator VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# 1. 寻找 Qt5
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets OpenGL)
find_package(Qt5 QUIET COMPONENTS Svg)  # 可选: SVG 图标支持

# 2. 寻找 VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# 3. 应对 Ubuntu 的暴力解法
include_directories(/usr/include/opencascade)

# ==========================================================
# 【新增】：让前端 Qt 也能找到 nlohmann/json.hpp 用于导出参数
# ==========================================================
include_directories(${CMAKE_SOURCE_DIR}/third_party)

# 4. 寻找并配置 Gmsh 库
find_library(GMSH_LIBRARY NAMES gmsh)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h)

find_package(Threads REQUIRED)

if(NOT GMSH_LIBRARY OR NOT GMSH_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find Gmsh! Please run: sudo apt install libgmsh-dev")
endif()

include_directories(${GMSH_INCLUDE_DIR})

# 搜集前端源码
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.h" "src/*.ui")

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(IonFlowSimulator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# 链接库
target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL
    ${VTK_LIBRARIES}
    ${GMSH_LIBRARY}
    Threads::Threads
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBO TKBRep
    TKGeomAlgo TKTopAlgo TKShHealing TKMesh TKXMesh
    TKXSBase TKIGES TKSTEP TKSTEPBase TKSTEPAttr TKSTEP209
    TKPrim TKBool TKOffset TKFillet
)
# 可选: SVG 图标支持 (未安装时忽略)
if(TARGET Qt5::Svg)
    target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Svg)
endif()