cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11")

project(miniproject)
add_executable(solve solve.cpp)

# --- Подключаем gmsh для третьего примера, тут делайте так же, как у вас сработало в первой лабе

# - если вам повезло, то сработает простой путь
include_directories(${GMSH_INCLUDE_DIR})

# - если у вас Ubuntu, то возможно, что приехал немного кривой gmsh из репозитория
# - при этом достаточно указать системные пути, куда система положила нужные файлы
# - закомментируйте строку 13 и раскомментируйте эти шесть:
#add_library(gmsh SHARED IMPORTED)
#set_target_properties(gmsh PROPERTIES
#    IMPORTED_LOCATION "/home/amisto/gmsh-4.12.0/lib/libgmsh.so"
#    INTERFACE_INCLUDE_DIRECTORIES "/home/amisto/gmsh-4.12.0/include/"
#)
set(GMSH_LIBRARIES gmsh)

# - если предыдущий вариант не сработал, собирайте gmsh из исходников
# - после этого в предыдущем варианте заменяйте системные пути на те, что вы сделали сами

# --- Закончили с gmsh, остальное править не обязательно (не поможет)


# Аналогично с VTK
find_package(VTK REQUIRED)
# ... подключим её со всеми её компонентами её нативными средствами
include(${VTK_USE_FILE})
# если этот вариант не сработал, собирайтесь из исходников
# можно пользоваться инструкциями от gmsh, все библиотеки собираются по единому принципу

# Укажем библиотеки для линковки
# ... тут только VTK
target_link_libraries(solve ${VTK_LIBRARIES})

# ... а тут ещё и gmsh
target_link_libraries(solve ${GMSH_LIBRARIES})
