#! /bin/bash


# Usage: grep_asan_failures.sh <ctest_log_file>
#
# Grep for all address sanitizer errors and combine ctest name with the the error that occurred.
# Print statistics about the different failures that occurred.

LOG_FILE=$1

readarray ASAN_DETAILS <<< "$(cat "$LOG_FILE" | grep "SUMMARY: AddressSanitizer" | uniq)"
readarray FAILED_TESTS <<< "$(cat "$LOG_FILE" | grep -A 3000 "The following tests FAILED:")"
unset "FAILED_TESTS[0]"

for f in "${FAILED_TESTS[@]}"; do
    TEST_ID=$(echo $f | awk '{print $1;}')
    echo $f
    echo "${ASAN_DETAILS[*]}" | grep "\s${TEST_ID}: SUMMARY"
done

echo "Summary of failures:"
# Note: "grep ." drops empty lines 
echo "${ASAN_DETAILS[*]}" | awk '{print $4};' | grep . | sort | uniq -c

