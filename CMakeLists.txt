# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MDM)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# 21/11/15	removed the SiPM module. //
# Load the g4sipm module
# add_subdirectory(externals/g4sipm)
# include_directories(externals/g4sipm/g4sipm/include)
# include_directories(externals/g4sipm/externals/jansson/src)
# include_directories(externals/g4sipm/externals/gtest/include)



#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# 21/11/15	removed the SiPM module. //
# Add persistency sources and headers.
# set(sources ${sources} "${PROJECT_SOURCE_DIR}/src/persistency/PersistencyHandler.cc" 
#   # "${PROJECT_SOURCE_DIR}/src/persistency/PersistencyHandlerMessenger.cc" 
#    "${PROJECT_SOURCE_DIR}/src/persistency/PersistVisitor.cc" 
#    "${PROJECT_SOURCE_DIR}/src/persistency/PersistVisitorFactory.cc" 
#    "${PROJECT_SOURCE_DIR}/src/persistency/JsonPersistVisitor.cc")
# set(headers ${headers} "${PROJECT_SOURCE_DIR}/include/persistency/PersistencyHandler.hh" 
#  # "${PROJECT_SOURCE_DIR}/include/persistency/PersistencyHandlerMessenger.hh" 
#   "${PROJECT_SOURCE_DIR}/include/persistency/PersistVisitor.hh" 
#    "${PROJECT_SOURCE_DIR}/include/persistency/PersistVisitorFactory.hh" 
#   "${PROJECT_SOURCE_DIR}/include/persistency/JsonPersistVisitor.hh")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(MDM_main MDM_main.cc ${sources} ${headers})#
target_link_libraries(MDM_main ${Geant4_LIBRARIES})

# 21/11/15	removed the SiPM module. //
# add the link to G4sipm
# target_link_libraries(MDM_main g4sipm boost_date_time jansson ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MDM_SCRIPTS
  init.mac
  init_vis.mac
  vis_MDM.mac
# 21/11/15	removed the SiPM module. //
# g4sipm.mac

  )

foreach(_script ${MDM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(MDM DEPENDS MDM_main)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS MDM_main DESTINATION bin)


