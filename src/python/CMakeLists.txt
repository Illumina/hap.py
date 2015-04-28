# copy stuff from here to there
# copy_files("${CMAKE_SOURCE_DIR}/src/python/*" "${CMAKE_BINARY_DIR}/python")

file(GLOB COPY_FILES *.py)
add_custom_target(python_bin ALL COMMENT "Copying python scripts...")

foreach(FILENAME ${COPY_FILES})
    get_filename_component(TARGETNAME ${FILENAME} "NAME")
    set(SRC "${FILENAME}")
    set(DST "${CMAKE_BINARY_DIR}/bin/${TARGETNAME}")

    add_custom_command(
      TARGET python_bin
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
    )
endforeach(FILENAME)

SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Haplo")
SET(PYTHON_DST "${CMAKE_BINARY_DIR}/lib/python27/Haplo")
add_custom_target(python_src1 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

add_custom_command(TARGET python_src1
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
)

SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Somatic")
SET(PYTHON_DST "${CMAKE_BINARY_DIR}/lib/python27/Somatic")
add_custom_target(python_src2 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

add_custom_command(TARGET python_src2
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
)

SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Tools")
SET(PYTHON_DST "${CMAKE_BINARY_DIR}/lib/python27/Tools")
add_custom_target(python_src3 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

add_custom_command(TARGET python_src3
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
)

configure_file("Haplo/version.py.in"
               "${CMAKE_BINARY_DIR}/lib/python27/Haplo/version.py")

