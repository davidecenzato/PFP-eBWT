set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/common)

include_directories(${INCLUDE_DIR})

add_executable(circpfp.x circpfp.cpp ${INCLUDE_DIR}/utils.c ${INCLUDE_DIR}/xerrors.c)
target_link_libraries(circpfp.x malloc_count z pthread)

add_executable(circpfpr.x circpfpmw.cpp ${INCLUDE_DIR}/utils.c ${INCLUDE_DIR}/xerrors.c)
target_link_libraries(circpfpr.x malloc_count z pthread)

add_executable(circpfpd.x circpfpd.cpp ${INCLUDE_DIR}/utils.c ${INCLUDE_DIR}/xerrors.c)
target_link_libraries(circpfpd.x malloc_count z pthread)

add_executable(circpfpNT.x circpfp.cpp ${INCLUDE_DIR}/utils.c)
target_link_libraries(circpfpNT.x malloc_count z)
target_compile_options(circpfpNT.x PUBLIC "-DNOTHREADS")

add_executable(circpfprNT.x circpfpmw.cpp ${INCLUDE_DIR}/utils.c)
target_link_libraries(circpfprNT.x malloc_count z)
target_compile_options(circpfprNT.x PUBLIC "-DNOTHREADS")

add_executable(circpfpdNT.x circpfpd.cpp ${INCLUDE_DIR}/utils.c)
target_link_libraries(circpfpdNT.x malloc_count z)
target_compile_options(circpfpdNT.x PUBLIC "-DNOTHREADS")

add_executable(circpfp-gca.x circpfp-gca.cpp ${INCLUDE_DIR}/utils.c ${INCLUDE_DIR}/xerrors.c)
target_link_libraries(circpfp-gca.x malloc_count z pthread)