cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(ODEs_visualizator)

find_package(VTK COMPONENTS 
  vtkInteractionWidgets
  vtkCommonComputationalGeometry
  vtkCommonDataModel
  vtkCommonColor
  vtkCommonTransforms
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingAnnotation
  vtkCommonCore
  vtkIOLegacy
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping ODEs_visualizator: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(ODEs_visualizator MACOSX_BUNDLE ODEs_visualizator.cpp )
  target_link_libraries(ODEs_visualizator PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(ODEs_visualizator MACOSX_BUNDLE ODEs_visualizator.cpp )
  target_link_libraries(ODEs_visualizator PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS ODEs_visualizator
    MODULES ${VTK_LIBRARIES}
    )
endif () 
