# - Try to find Soap, assumes FindMono has already been run
#
#  SOAP_FOUND
#  SOAP_LIBRARIES
#
# Author: Ewen Cheslack-Postava, 2008
#

# backup PKG_CONFIG_PATH and add Mono's path
SET(SOAP_DEFAULT_LIB System.Runtime.Serialization.Formatters.Soap.dll)
IF(WIN32)

  IF(MONO_ROOT)
    FIND_FILE(SOAP_LIB ${SOAP_DEFAULT_LIB} PATHS ${MONO_ROOT}/lib/mono/assemblies)
    IF(SOAP_LIB)
      SET(SOAP_FOUND TRUE)
      SET(SOAP_LIBRARIES "-r:${SOAP_LIB}")
    ENDIF(SOAP_LIB)
  ELSE(MONO_ROOT)
    SET(SOAP_FOUND TRUE)
    SET(SOAP_LIBRARIES "-r:${SOAP_DEFAULT_LIB}") 
  ENDIF(MONO_ROOT)
ELSE(WIN32)
  SET(SOAP_FOUND TRUE)
  SET(SOAP_LIBRARIES "-r:${SOAP_DEFAULT_LIB}")
ENDIF(WIN32)

# restore PKG_CONFIG_PATH
SET(ENV{PKG_CONFIG_PATH} OLD_PKG_CONFIG_PATH)


# print out what we found
IF(SOAP_FOUND)
  MESSAGE(STATUS "Soap Library Parameters: ${SOAP_LIBRARIES}")
ENDIF(SOAP_FOUND)