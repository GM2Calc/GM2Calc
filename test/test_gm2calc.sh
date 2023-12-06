#!/bin/sh

BASEDIR="$(dirname $0)"
GM2CALC=${GM2CALC:-bin/gm2calc.x}

if test ! -x "${GM2CALC}"; then
    echo "Warning: GM2Calc not executable: \"${GM2CALC}\"; skipping test"
    exit
fi

errors=0
passes=0

expect_failure() {
    if [ $1 -eq 0 ] ; then
        printf "%s\n" " [FAIL]"
        errors=$(expr $errors + 1)
    else
        printf "%s\n" " [OK]"
        passes=$(expr $passes + 1)
    fi
}

expect_success() {
    if [ $1 -eq 0 ] ; then
        printf "%s\n" " [OK]"
        passes=$(expr $passes + 1)
    else
        printf "%s\n" " [FAIL]"
        errors=$(expr $errors + 1)
    fi
}

test_no_input_source() {
    printf "%s" "test_no_input_source $1"
    ${GM2CALC} "--$1-input-file=" 2>/dev/null
    expect_failure "$?"
}

test_detailed_output() {
    printf "%s" "test_detailed_output"
    { cat $2
      cat <<EOF
Block GM2CalcConfig
     0     1     # detailed output
EOF
    } | ${GM2CALC} "--$1-input-file=-" >/dev/null 2>&1
    expect_success "$?"
}

test_gm2calc_output() {
    printf "%s" "test_gm2calc_output"
    { cat $2
      cat <<EOF
Block GM2CalcConfig
     0     4     # GM2CalcOutput
EOF
    } | ${GM2CALC} "--$1-input-file=-" >/dev/null 2>&1
    expect_success "$?"
}

test_spheno_output() {
    printf "%s" "test_spheno_output"
    { cat $2
      cat <<EOF
Block GM2CalcConfig
     0     3     # SPhenoLowEnergy
EOF
    } | ${GM2CALC} "--$1-input-file=-" >/dev/null 2>&1
    expect_success "$?"
}

test_nmssmtools_output() {
    printf "%s" "test_nmssmtools_output"
    { cat $2
      cat <<EOF
Block GM2CalcConfig
     0     2     # SPhenoLowEnergy
EOF
    } | ${GM2CALC} "--$1-input-file=-" >/dev/null 2>&1
    expect_success "$?"
}

# run tests
test_no_input_source "gm2calc"
test_no_input_source "slha"
test_no_input_source "thdm"
test_detailed_output "gm2calc" "${BASEDIR}/../input/example.gm2"
test_detailed_output "slha" "${BASEDIR}/../input/example.slha"
test_detailed_output "thdm" "${BASEDIR}/../input/example.thdm"
test_gm2calc_output "gm2calc" "${BASEDIR}/../input/example.gm2"
test_gm2calc_output "slha" "${BASEDIR}/../input/example.slha"
test_gm2calc_output "thdm" "${BASEDIR}/../input/example.thdm"
test_nmssmtools_output "gm2calc" "${BASEDIR}/../input/example.gm2"
test_nmssmtools_output "slha" "${BASEDIR}/../input/example.slha"
test_nmssmtools_output "thdm" "${BASEDIR}/../input/example.thdm"
test_spheno_output "gm2calc" "${BASEDIR}/../input/example.gm2"
test_spheno_output "slha" "${BASEDIR}/../input/example.slha"
test_spheno_output "thdm" "${BASEDIR}/../input/example.thdm"

count=$(expr $errors + $passes)

if [ $errors -eq 0 ] ; then
    echo "Test result: OK [${passes}/${count}]"
else
    echo "Test result: FAIL [${passes}/${count}]"
fi

exit ${errors}
