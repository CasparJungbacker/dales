macro( dales_find_cuda )
  # Check if a target architecture has been set. If not, default to 80 (A100 & RTX3090)
  # For H100, use 90
  if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    set( CMAKE_CUDA_ARCHITECTURES 80 )
  endif()
  
  # Look for the CUDA Toolkit, which has cuFFT and NVTX
  find_package( CUDAToolkit )
  if( NOT TARGET CUDA::cufft AND ENABLE_ACC )
    ecbuild_error( "Could not find cuFFT, which is required for GPU builds!" )
  else()
    ecbuild_info( "Found cuFFT: ${CUDA_cufft_LIBRARY}" )
  endif()
  if( NOT TARGET CUDA::nvToolsExt AND ENABLE_NVTX )
    ecbuild_error( "Could not find NVTX library for profiling" )
  else()
    ecbuild_info( "Found NVTX: ${CUDA_nvToolsExt_LIBRARY}" )
  endif() 
endmacro()

