# from wzdftpd
# via http://fresh.t-systems-sfr.com/unix/src/privat2/wzdftpd-0.8.2.tar.gz:a/wzdftpd-0.8.2/cmake/FindOpenSSL.cmake
# - Try to find the OpenSSL encryption library
# Once done this will define
#
#  CONFIG_FOUND - system has the config++ library
#  CONFIG_INCLUDE_DIR - the config++ include directory
#  CONFIG_LIBRARIES - The libraries needed to use config++

if (CONFIG_INCLUDE_DIR AND CONFIG_LIBRARIES)

  # in cache already
  SET(CONFIG_FOUND TRUE)

else (CONFIG_INCLUDE_DIR AND CONFIG_LIBRARIES)
  FIND_PATH(CONFIG_INCLUDE_DIR libconfig.h++
     ${CONFIG_ROOT}/include
     /usr/include/
     /usr/local/include/
     $ENV{ProgramFiles}/CONFIG/include/
     $ENV{SystemDrive}/CONFIG/include/
  )

  if(WIN32 AND MSVC)
  else(WIN32 AND MSVC)
    FIND_LIBRARY(CONFIG_LIBRARIES NAMES config++
      PATHS
      ${CONFIG_ROOT}/lib
      /usr/lib
      /usr/local/lib
    )
  endif(WIN32 AND MSVC)

  if (CONFIG_INCLUDE_DIR AND CONFIG_LIBRARIES)
     set(CONFIG_FOUND TRUE)
  endif (CONFIG_INCLUDE_DIR AND CONFIG_LIBRARIES)

  if (CONFIG_FOUND)
     if (NOT CONFIG_FIND_QUIETLY)
        message(STATUS "Found config++: ${CONFIG_LIBRARIES}")
     endif (NOT CONFIG_FIND_QUIETLY)
  else (CONFIG_FOUND)
     if (CONFIG_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find config++")
     endif (CONFIG_FIND_REQUIRED)
  endif (CONFIG_FOUND)

  MARK_AS_ADVANCED(CONFIG_INCLUDE_DIR CONFIG_LIBRARIES)

endif (CONFIG_INCLUDE_DIR AND CONFIG_LIBRARIES)