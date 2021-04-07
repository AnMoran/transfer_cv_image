cmake_minimum_required(VERSION 2.8)

project(transfer_cv_image)

add_definitions(-std=c++11)
add_definitions("-Wall")
add_definitions(-DUNICODE -D_UNICODE)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)




find_package( OpenCV REQUIRED )

if (NOT OpenCV_FOUND)
    message(FATAL_ERROR "opencv Not Found!")
else()
    message(STATUS "OpenCV_LIBS: ${OpenCV_LIBS}")
    message(STATUS "OpenCV_INCLUDE_DIRS: ${OpenCV_INCLUDE_DIRS}")
endif (NOT OpenCV_FOUND)


set(SRC
        src/main.cpp
        src/norm_4pt.cpp
        )

add_executable(demo ${SRC})

target_link_libraries(demo ${OpenCV_LIBS})

