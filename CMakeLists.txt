PROJECT(PatternViewer)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(src/UI)


INCLUDE_DIRECTORIES(3rdparty/vcglib)
INCLUDE_DIRECTORIES(3rdparty/vcglib/eigenlib)
INCLUDE_DIRECTORIES(3rdparty/LearningQuadrangulationsCore)
INCLUDE_DIRECTORIES(3rdparty/LearningQuadrangulationsCore/SCIP)

#Qt4
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTORCC ON)

FIND_PACKAGE(Qt4 REQUIRED)

# Mac/Linux compiler settings
# IF(UNIX)
# #  SET(CMAKE_BUILD_TYPE Debug)
#   INCLUDE_DIRECTORIES(
#     /usr/local/include
#     /opt/local/include
#   )
#   LINK_DIRECTORIES(
#     /usr/local/lib
#     /opt/local/lib
#   )
  # IF(APPLE)
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -msse4 -stdlib=libc++")
    ADD_DEFINITIONS(-D_MACOSX)
#   ELSE()
#     SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -msse4")
#     ADD_DEFINITIONS(-D_UNIX)
#   ENDIF()
# ENDIF()



FILE(GLOB pattern_viewer_src
	./src/main.cpp
  ./src/UI/*.cpp
)

# linker settings
ADD_EXECUTABLE(PatternViewer
  ${pattern_viewer_src})
FIND_LIBRARY(SQLITE3_LIB sqlite3)
TARGET_LINK_LIBRARIES(PatternViewer
  ${SQLITE3_LIB}
  Qt4::QtGui
  Qt4::QtCore
)
  
