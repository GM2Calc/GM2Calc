#!/bin/sh

BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
BINDIR="${HOMEDIR}/bin"
GM2CALC="${BINDIR}/gm2calc.x"
frac_diff=0.0000000000001

. "${BASEDIR}/compare.sh"

# points to be tested
# 1: path to parameter point
# 2: format (gm2calc or slha)
# 3: expected exit code
# 4: expected amu (minimal output)
points="\
${BASEDIR}/test_points/problems_funcs_Fa_Fb.in,gm2calc,0,-7.26010839613302e-09
${BASEDIR}/test_points/problems_funcs_I_FC_FN.in,gm2calc,0,-7.38800970725992e-09
${BASEDIR}/test_points/P1a_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.80707343687898e-09
${BASEDIR}/test_points/P1a_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.1987194238708e-09
${BASEDIR}/test_points/P1b_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.22472102858175e-09
${BASEDIR}/test_points/P1b_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.11846724283043e-09
${BASEDIR}/test_points/P1c_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.00765380788945e-09
${BASEDIR}/test_points/P1c_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.85987928891313e-09
${BASEDIR}/test_points/P2_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,0,3.1272477101719e-09
${BASEDIR}/test_points/P2_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,0,3.18741953459232e-09
${BASEDIR}/test_points/P3_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.78175706564947e-09
${BASEDIR}/test_points/P3_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.14398237095065e-09
${BASEDIR}/test_points/BM1-1504.05500_1L_resummed.in,gm2calc,1,2.81168064521252e-09
${BASEDIR}/test_points/BM1-1504.05500_2L_resummed.in,gm2calc,1,2.67321328584044e-09
${BASEDIR}/test_points/BM2-1504.05500_1L_resummed.in,gm2calc,1,3.0339292608166e-09
${BASEDIR}/test_points/BM2-1504.05500_2L_resummed.in,gm2calc,1,2.87958582098365e-09
${BASEDIR}/test_points/BM3-1504.05500_1L_resummed.in,gm2calc,0,2.61994966785813e-09
${BASEDIR}/test_points/BM3-1504.05500_2L_resummed.in,gm2calc,0,2.08900966809353e-09
${BASEDIR}/test_points/BM4-1504.05500_1L_resummed.in,gm2calc,1,2.76249843887181e-09
${BASEDIR}/test_points/BM4-1504.05500_2L_resummed.in,gm2calc,1,2.4357222606711e-09
${BASEDIR}/test_points/BM5-1504.05500_1L_resummed.in,gm2calc,0,2.88732972683492e-09
${BASEDIR}/test_points/BM5-1504.05500_2L_resummed.in,gm2calc,0,2.94607528454977e-09
"

if test ! -x "${GM2CALC}"; then
    echo "Error: ${GM2CALC} not built"
    exit 1
fi

errors=0
passes=0
count=0

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
        passes=$(expr $passes + 1)
    else
        echo "result         : FAIL"
        errors=$(expr $errors + 1)
    fi

    count=$(expr $count + 1)

    echo ""
done

if [ $errors -eq 0 ] ; then
    echo "Test result: OK [${passes}/${count}]"
else
    echo "Test result: FAIL [${passes}/${count}]"
fi

exit ${errors}
