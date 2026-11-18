# ExternalProject may rerun PATCH_COMMAND when an sdist extraction path changes.
# Accept an already-applied patch only when its complete reverse check succeeds.
foreach(name IN ITEMS libomp-private.patch libomp-private-windows.patch)
  set(patch "${PATCH_DIR}/${name}")
  execute_process(COMMAND "${GIT_EXECUTABLE}" apply --check "${patch}"
    WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE clean OUTPUT_QUIET ERROR_QUIET)
  if(clean EQUAL 0)
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply "${patch}"
      WORKING_DIRECTORY "${SOURCE_DIR}" COMMAND_ERROR_IS_FATAL ANY)
  else()
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply --reverse --check "${patch}"
      WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE applied OUTPUT_QUIET ERROR_QUIET)
    if(NOT applied EQUAL 0)
      message(FATAL_ERROR "${name} is neither applicable nor already applied; use a clean external build tree")
    endif()
    message(STATUS "${name} already applied (reverse check passed)")
  endif()
endforeach()
