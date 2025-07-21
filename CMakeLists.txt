cmake_minimum_required(VERSION 3.10)
project(pn25)

FIND_PACKAGE(symmetri REQUIRED)

message(STATUS "Building for ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall -Wextra -pedantic")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -O0")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /fsanitize=address")
endif()

include_directories("include")


add_executable(pn2HDA pn2HDA.cpp)
target_link_libraries(pn2HDA PRIVATE ${symmetri_LIBRARIES})

add_executable(pn2MAXHDA pn2MAXHDA.cpp)
target_link_libraries(pn2MAXHDA PRIVATE ${symmetri_LIBRARIES} z3)
