#!/bin/sh

GM2CALC=${GM2CALC:-bin/gm2calc.x}

if test ! -x "${GM2CALC}"; then
    echo "Warning: GM2Calc not executable: \"${GM2CALC}\"; skipping test"
    exit
fi

errors=0
passes=0

expect_failure() {
    if [ $1  -eq 0 ] ; then
        printf "%s\n" " [FAIL]"
        errors=$(expr $errors + 1)
    else
        printf "%s\n" " [OK]"
        passes=$(expr $passes + 1)
    fi
}

test_no_input_source() {
    printf "%s" "test_no_input_source $1"
    ${GM2CALC} "--$1-input-file=" 2>/dev/null
    expect_failure "$?"
}

# run tests
test_no_input_source "gm2calc"
test_no_input_source "slha"
test_no_input_source "thdm"

count=$(expr $errors + $passes)

if [ $errors -eq 0 ] ; then
    echo "Test result: OK [${passes}/${count}]"
else
    echo "Test result: FAIL [${passes}/${count}]"
fi

exit ${errors}
