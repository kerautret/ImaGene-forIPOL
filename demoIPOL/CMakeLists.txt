# Make sure the compiler can find include files.
include_directories (${ImaGene_SOURCE_DIR}/include)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
    endif(COMMAND cmake_policy)


set( DEMO_IPOL_EXEC_FILES 
     meaningfulScaleEstim      )

# set( DEMO_IPOL_EXEC_FILES )
foreach( X ${DEMO_IPOL_EXEC_FILES} )
	 set( DEMO_IPOL_EXEC_NAME "${X}${SUFFIXBIN}" )
	 add_executable( ${DEMO_IPOL_EXEC_NAME} ${X}.cxx )
	 target_link_libraries( ${DEMO_IPOL_EXEC_NAME} ${LIBIMAGENE_NAME} )
endforeach( X )





	 











