#!/bin/sh

BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
BINDIR="${HOMEDIR}/bin"
GM2CALC="${BINDIR}/gm2calc.x"
frac_diff=0.000000000000001

. "${BASEDIR}/compare.sh"

# points to be tested
# 1: path to parameter point
# 2: format (gm2calc or slha)
# 3: expected exit code
# 4: expected amu (minimal output)
points="\
${BASEDIR}/test_points/problems_funcs_Fa_Fb.in,gm2calc,0,-7.26010839613302e-09
${BASEDIR}/test_points/problems_funcs_I_FC_FN.in,gm2calc,0,-7.38800970725992e-09
"

if test ! -x "${GM2CALC}"; then
    echo "Error: ${GM2CALC} not built"
    exit 1
fi

errors=0

for point in ${points}
do
    error=0
    input="`echo ${point} | tr ',' ' ' | awk '{ print $1 }'`"
    format="`echo ${point} | tr ',' ' ' | awk '{ print $2 }'`"
    expected_exit="`echo ${point} | tr ',' ' ' | awk '{ print $3 }'`"
    expected_amu="`echo ${point} | tr ',' ' ' | awk '{ print $4 }'`"

    amu=$(
        { cat ${input}
          cat <<EOF
Block GM2CalcConfig
     0     0     # minimal output
EOF
        } | ${GM2CALC} "--${format}-input-file=-" 2>/dev/null
    )

    exit_code="$?"

    echo "====== ${input} ======"
    echo "executable     : ${GM2CALC}"
    echo "input          : ${input}"
    echo "format         : ${format}"
    echo "expected exit  : ${expected_exit}"
    echo "expected amu   : ${expected_amu}"
    echo "exit code      : ${exit_code}"
    echo "amu            : ${amu}"

    [ ${expected_exit} -ne ${exit_code} ] && error=1

    if [ $error -eq 0 ] ; then
        CHECK_EQUAL_FRACTION "$expected_amu" "$amu" "$frac_diff"
        [ $? -ne 0 ] && error=1
    fi

    if [ $error -eq 0 ] ; then
        echo "result         : OK"
    else
        echo "result         : FAIL"
        errors=1
    fi

    echo ""
done

if [ $errors -eq 0 ] ; then
    echo "Test result: OK"
else
    echo "Test result: FAIL"
fi

exit ${errors}
