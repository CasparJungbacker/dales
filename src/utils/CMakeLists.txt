list(
  APPEND dales_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/modfft_t.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modfft_fftpack_t.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modnvtx.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modtimer.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/submodfft_factory.f90
)

set( dales_srcs "${dales_srcs}" PARENT_SCOPE )