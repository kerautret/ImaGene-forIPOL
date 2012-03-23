# Module for checking for ImaGene
find_path( GMPXX_INCLUDE_DIR gmpxx.h $ENV{HOME}/local/include  $ENV{HOME}/include /opt/local/include /usr/local/include /usr/include )
find_library( GMPXX_LIBRARY gmpxx PATHS $ENV{HOME}/local/lib $ENV{HOME}/lib /opt/local/lib  /usr/local/lib /usr/lib )

if ( GMPXX_INCLUDE_DIR AND GMPXX_LIBRARY )
   set ( GMPXX_FOUND true )
endif ( GMPXX_INCLUDE_DIR AND GMPXX_LIBRARY )

if ( GMPXX_FOUND )
   if ( NOT GMPXX_FIND_QUIETLY )
       message( STATUS "Found gmpxx: include=${GMPXX_INCLUDE_DIR}" )
       include_directories( ${GMPXX_INCLUDE_DIR} )
   endif ( NOT GMPXX_FIND_QUIETLY )
else ( GMPXX_FOUND )
   if ( GMPXX_FIND_REQUIRED )
       message( FATAL_ERROR "Could not find gmpxx in $ENV{HOME}/local/inc\
lude  $ENV{HOME}/include /usr/local/include /usr/include" )
   endif ( GMPXX_FIND_REQUIRED )
endif ( GMPXX_FOUND )


