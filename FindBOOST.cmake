# Module for checking for BOOST
find_path( BOOST_INCLUDE_DIR boost/math/distributions.hpp $ENV{HOME}/local/include  $ENV{HOME}/include /usr/local/include /usr/include /opt/local/include )

if ( BOOST_INCLUDE_DIR )
   set ( BOOST_FOUND true )
endif ( BOOST_INCLUDE_DIR )

if ( BOOST_FOUND )
   if ( NOT BOOST_FIND_QUIETLY )
       message( STATUS "Found BOOST: include=${BOOST_INCLUDE_DIR}" )
   endif ( NOT BOOST_FIND_QUIETLY )
else ( BOOST_FOUND )
   if ( BOOST_FIND_REQUIRED )
       message( FATAL_ERROR "Could not find BOOST in $ENV{HOME}/local/include  $ENV{HOME}/include /usr/local/include /usr/include /opt/local/include" )
   endif ( BOOST_FIND_REQUIRED )
endif ( BOOST_FOUND )


