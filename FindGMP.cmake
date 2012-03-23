# Module for checking for ImaGene
find_path( GMP_INCLUDE_DIR gmp.h $ENV{HOME}/local/include  $ENV{HOME}/include /opt/local/include /usr/local/include /usr/include )
find_library( GMP_LIBRARY gmp PATHS $ENV{HOME}/local/lib $ENV{HOME}/lib /opt/local/lib /usr/local/lib /usr/lib )

if ( GMP_INCLUDE_DIR AND GMP_LIBRARY )
   set ( GMP_FOUND true )
endif ( GMP_INCLUDE_DIR AND GMP_LIBRARY )

if ( GMP_FOUND )
   if ( NOT GMP_FIND_QUIETLY )
       message( STATUS "Found gmp: include=${GMP_INCLUDE_DIR}" )
       include_directories( ${GMP_INCLUDE_DIR} )
   endif ( NOT GMP_FIND_QUIETLY )
else ( GMP_FOUND )
   if ( GMP_FIND_REQUIRED )
       message( FATAL_ERROR "Could not find gmp in $ENV{HOME}/local/inc\
lude  $ENV{HOME}/include /usr/local/include /usr/include" )
   endif ( GMP_FIND_REQUIRED )
endif ( GMP_FOUND )


