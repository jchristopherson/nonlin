execute_process(
    COMMAND "${TEST_EXECUTABLE}" "${TEST_ARGUMENT}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error_output
)
if (NOT "${result}" STREQUAL "${EXPECTED_ERROR}")
    message(FATAL_ERROR
        "Expected error ${EXPECTED_ERROR}, got ${result}.\n${output}${error_output}")
endif()
