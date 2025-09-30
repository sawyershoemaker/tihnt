cmake_minimum_required(VERSION 3.20)

set(SRC
    src/main.cpp
    src/ws_server.cpp
    src/board.cpp
    src/json.cpp
    src/solver.cpp
    src/overlay.cpp
)

add_executable(mines_overlay WIN32 ${SRC})

target_include_directories(mines_overlay PRIVATE include)

if(WIN32)
  target_compile_definitions(mines_overlay PRIVATE UNICODE _UNICODE NOMINMAX WIN32_LEAN_AND_MEAN)
  target_link_libraries(mines_overlay PRIVATE user32 gdi32 ws2_32)
endif()


