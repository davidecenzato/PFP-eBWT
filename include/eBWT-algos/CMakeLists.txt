set(INCLUDE_DIR_UT ${PROJECT_SOURCE_DIR}/include/common)
set(INCLUDE_DIR_ALGO ${PROJECT_SOURCE_DIR}/include/eBWT-algos)
set(INCLUDE_DIR_SAIS ${PROJECT_SOURCE_DIR}/include/csais)

include_directories(${INCLUDE_DIR_UT})
include_directories(${INCLUDE_DIR_ALGO})
include_directories(${INCLUDE_DIR_SAIS})


add_executable(parsebwtNT.x ${INCLUDE_DIR_ALGO}/parse.cpp ${INCLUDE_DIR_UT}/utils.c ${INCLUDE_DIR_SAIS}/csais.cpp)
target_link_libraries(parsebwtNT.x malloc_count dl sdsl divsufsort divsufsort64)

add_executable(parsebwtNT64.x ${INCLUDE_DIR_ALGO}/parse.cpp ${INCLUDE_DIR_UT}/utils.c ${INCLUDE_DIR_SAIS}/csais.cpp)
target_link_libraries(parsebwtNT64.x malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(parsebwtNT64.x PUBLIC "-DP64")

add_executable(parsebwtgcaNT.x ${INCLUDE_DIR_ALGO}/parse-gca.cpp ${INCLUDE_DIR_UT}/utils.c ${INCLUDE_DIR_SAIS}/csais.cpp)
target_link_libraries(parsebwtgcaNT.x malloc_count dl sdsl divsufsort divsufsort64)

add_executable(parsebwtgcaNT64.x ${INCLUDE_DIR_ALGO}/parse-gca.cpp ${INCLUDE_DIR_UT}/utils.c ${INCLUDE_DIR_SAIS}/csais.cpp)
target_link_libraries(parsebwtgcaNT64.x malloc_count dl sdsl divsufsort divsufsort64)
target_compile_options(parsebwtgcaNT64.x PUBLIC "-DP64")

add_executable(build-ebwt-gca.x ebwt.cpp ${INCLUDE_DIR_UT}/utils.c)
target_link_libraries(build-ebwt-gca.x malloc_count dl gsacak sdsl divsufsort divsufsort64)

add_executable(build-ebwt-gca_p64.x ebwt.cpp ${INCLUDE_DIR_UT}/utils.c)
target_link_libraries(build-ebwt-gca_p64.x malloc_count dl gsacak sdsl divsufsort divsufsort64)
target_compile_options(build-ebwt-gca_p64.x PUBLIC "-DP64")

add_executable(build-ebwt-gca_d64.x ebwt.cpp ${INCLUDE_DIR_UT}/utils.c)
target_link_libraries(build-ebwt-gca_d64.x malloc_count dl z gsacak64 sdsl divsufsort divsufsort64)
target_compile_options(build-ebwt-gca_d64.x PUBLIC "-DM64")

add_executable(build-ebwt-gca_64.x ebwt.cpp ${INCLUDE_DIR_UT}/utils.c)
target_link_libraries(build-ebwt-gca_64.x malloc_count dl gsacak64 sdsl divsufsort divsufsort64)
target_compile_options(build-ebwt-gca_64.x PUBLIC "-DM64")
target_compile_options(build-ebwt-gca_64.x PUBLIC "-DP64")

add_executable(invertNT.x invertebwt.cpp ${INCLUDE_DIR_UT}/common.hpp)
target_link_libraries(invertNT.x sdsl divsufsort divsufsort64)