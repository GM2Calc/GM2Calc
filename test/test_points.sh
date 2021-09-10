#!/bin/sh

BASEDIR="$(dirname $0)"
GM2CALC=${GM2CALC:-bin/gm2calc.x}
frac_diff=0.000000001

. "${BASEDIR}/compare.sh"

# points to be tested
# 1: path to parameter point
# 2: format (gm2calc or slha)
# 3: expected exit code
# 4: expected amu (minimal output)
points="\
${BASEDIR}/test_points/problems_funcs_Fa_Fb.in,gm2calc,0,-7.26010840E-09
${BASEDIR}/test_points/problems_funcs_I_FC_FN.in,gm2calc,0,-7.38800971E-09
${BASEDIR}/test_points/problems_funcs_Mu_zero_1L.in,gm2calc,0,5.78789509e-04
${BASEDIR}/test_points/problems_funcs_Mu_zero_2L.in,gm2calc,1,
${BASEDIR}/test_points/problems_funcs_M1_zero.in,gm2calc,0,2.57054680e-10
${BASEDIR}/test_points/problems_funcs_xim_zero.in,gm2calc,0,-2.26641081e-10
${BASEDIR}/test_points/problems_hmix_scale.in,slha,0,3.89209714e-10
${BASEDIR}/test_points/problems_throw_me2_convergence.in,slha,0,2.30017838e-10
${BASEDIR}/test_points/problems_negative_soft_mass.in,slha,1,
${BASEDIR}/test_points/problems_infinite_TB.in,slha,1,
${BASEDIR}/test_points/problems_zero_TB.in,slha,1,
${BASEDIR}/test_points/problems_bug_smuon_mixing.in,slha,0,2.33255248e-10
${BASEDIR}/test_points/problems_bug_smuon_mixing-2.in,slha,0,3.17520811E-10
${BASEDIR}/test_points/problems_bino_reordering.in,slha,0,1.67607964E-10
${BASEDIR}/test_points/problems_bino_reordering_pole_running.in,slha,0,1.20357378e-10
${BASEDIR}/test_points/P1a_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.80707344E-09
${BASEDIR}/test_points/P1a_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.16049629E-09
${BASEDIR}/test_points/P1b_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.22472103E-09
${BASEDIR}/test_points/P1b_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.11758474E-09
${BASEDIR}/test_points/P1c_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.00765381E-09
${BASEDIR}/test_points/P1c_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.86038732E-09
${BASEDIR}/test_points/P2_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,0,3.12724771E-09
${BASEDIR}/test_points/P2_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,0,3.19581635E-09
${BASEDIR}/test_points/P3_1L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,2.78175707E-09
${BASEDIR}/test_points/P3_2L_resummed_diploma_thesis_Markus_Bach.in,gm2calc,1,3.09394325E-09
${BASEDIR}/test_points/BM1-1504.05500_1L_resummed.in,gm2calc,1,2.81168065E-09
${BASEDIR}/test_points/BM1-1504.05500_2L_resummed.in,gm2calc,1,2.67369843E-09
${BASEDIR}/test_points/BM2-1504.05500_1L_resummed.in,gm2calc,1,3.03392926E-09
${BASEDIR}/test_points/BM2-1504.05500_2L_resummed.in,gm2calc,1,2.88017459E-09
${BASEDIR}/test_points/BM3-1504.05500_1L_resummed.in,gm2calc,0,2.61994967E-09
${BASEDIR}/test_points/BM3-1504.05500_2L_resummed.in,gm2calc,0,2.09846139E-09
${BASEDIR}/test_points/BM4-1504.05500_1L_resummed.in,gm2calc,1,2.76249844E-09
${BASEDIR}/test_points/BM4-1504.05500_2L_resummed.in,gm2calc,1,2.44726409E-09
${BASEDIR}/test_points/BM5-1504.05500_1L_resummed.in,gm2calc,0,2.88732973E-09
${BASEDIR}/test_points/BM5-1504.05500_2L_resummed.in,gm2calc,0,2.98718794E-09
${BASEDIR}/test_points/thdm_mass-basis.in,thdm,0,1.62362408E-11
${BASEDIR}/test_points/thdm_mass-basis_test_point_1.in,thdm,0,-5.86743344E-12
${BASEDIR}/test_points/thdm_mass-basis_test_point_2.in,thdm,0,1.87875466E-11
${BASEDIR}/test_points/thdm_mass-basis_test_point_3.in,thdm,0,1.87626300E-11
${BASEDIR}/test_points/thdm_mass-basis_test_point_4.in,thdm,0,-5.85790876E-12
${BASEDIR}/test_points/thdm_gauge-basis.in,thdm,0,3.27209495E-11
${BASEDIR}/test_points/thdm_contradictory_input.in,thdm,1,
"

if test ! -x "${GM2CALC}"; then
    echo "Error: GM2Calc not executable: \"${GM2CALC}\""
    exit 1
fi

errors=0
passes=0
count=0

for point in ${points}
do
    error=0
    input="`echo ${point} | awk -F , '{ print $1 }'`"
    format="`echo ${point} | awk -F , '{ print $2 }'`"
    expected_exit="`echo ${point} | awk -F , '{ print $3 }'`"
    expected_amu="`echo ${point} | awk -F , '{ print $4 }'`"

    amu=$(
        { cat ${input}
          cat <<EOF
Block GM2CalcConfig
     0     0     # minimal output
     4     0     # non-verbose
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
