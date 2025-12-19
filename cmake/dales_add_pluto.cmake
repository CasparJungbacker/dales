include(FetchContent)

ecbuild_info( "Fetching Atlas" )

set( ENABLE_FORTRAN ON )

FetchContent_Declare(
    atlas
    GIT_REPOSITORY "https://github.com/ecmwf/atlas.git"
    GIT_TAG 358bfaf6988aa1e90586546cab2c1fcdf70f1ca4 # 0.44.1
    SOURCE_SUBDIR "pluto/"
)

FetchContent_MakeAvailable(atlas)