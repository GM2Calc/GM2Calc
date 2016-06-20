#!/bin/sh

BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
BINDIR="${HOMEDIR}/bin"
EXAMPLES="\
${BINDIR}/example-gm2calc.x,${BINDIR}/example-gm2calc_c.x \
${BINDIR}/example-slha.x,${BINDIR}/example-slha_c.x \
"
frac_diff=0.00001

. "${BASEDIR}/compare.sh"

errors=0
passes=0
count=0

for e in ${EXAMPLES}; do
    error=0
    result=

    ex_cpp=$(echo "$e" | tr ',' ' ' | awk '{ print $1 }')
    ex_c=$(echo "$e" | tr ',' ' ' | awk '{ print $2 }')

    [ -x "$ex_cpp" ] || {
        echo "Error: $ex_cpp not built"
        errors=$(expr $errors + 1);
        continue;
    }
    [ -x "$ex_c" ] || {
        echo "Error: $ex_c not built"
        errors=$(expr $errors + 1);
        continue;
    }

    out_cpp=$("$ex_cpp")
    out_c=$("$ex_c")
    amu_cpp=$(echo "$out_cpp" | awk '{ print $3 }')
    amu_c=$(echo "$out_c" | awk '{ print $3 }')
    delta_amu_cpp=$(echo "$out_cpp" | awk '{ print $5 }')
    delta_amu_c=$(echo "$out_c" | awk '{ print $5 }')

    CHECK_EQUAL_FRACTION "$amu_cpp" "$amu_c" "$frac_diff"
    [ $? -ne 0 ] && error=1

    CHECK_EQUAL_FRACTION "$delta_amu_cpp" "$delta_amu_c" "$frac_diff"
    [ $? -ne 0 ] && error=1

    if [ $error -eq 0 ] ; then
        result="OK"
        passes=$(expr $passes + 1)
    else
        result="FAIL"
        errors=$(expr $errors + 1)
    fi

    count=$(expr $count + 1)

    echo "========================="
    echo "running:"
    echo "   $ex_cpp"
    echo "   $ex_c"
    echo "output:"
    echo "   amu       (C++) = $amu_cpp"
    echo "   amu       (C)   = $amu_c"
    echo "   delta amu (C++) = $delta_amu_cpp"
    echo "   delta amu (C)   = $delta_amu_c"
    echo "result: $result"
    echo "========================="
    echo ""
done

if [ $errors -eq 0 ] ; then
    echo "Test result: OK [${passes}/${count}]"
else
    echo "Test result: FAIL [${passes}/${count}]"
fi

exit ${errors}
