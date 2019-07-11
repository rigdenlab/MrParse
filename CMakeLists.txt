cmake_minimum_required(VERSION 2.8)
project(MrParse NONE)

install(
    DIRECTORY
        mrparse
    DESTINATION
        lib/py2
)

install(
    DIRECTORY
        data tests
    DESTINATION
        share/mrparse
)

install(
    FILES
        CHANGELOG.rst CONTRIB.rst LICENSE README.md
    DESTINATION
        share/mrparse
)

if (WIN32)
  install(
      PROGRAMS 
          bin/mrparse 
      DESTINATION 
          bin
  )
else()
  install(
      PROGRAMS
          bin/mrparse
      DESTINATION
          bin
)
endif()
