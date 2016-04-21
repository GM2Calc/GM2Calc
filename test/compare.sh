# This file contains helper functions for test scripts

# compares two floating point numbers for equality
# with a given maximum relative deviation
#
# Note: scientific notation is allowed
CHECK_EQUAL_FRACTION() {
    if test $# -lt 3 ; then
        echo "Error: CHECK_EQUAL_FRACTION: Too few arguments"
        echo "Usage: CHECK_EQUAL_FRACTION $num1 $num2 $fraction"
        exit 1
    fi

    local num1=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    local num2=$(echo "$2" | sed -e 's/[eE]+*/*10^/')
    local frac=$(echo "$3" | sed -e 's/[eE]+*/*10^/')

    local scale=15

    local error=$(cat <<EOF | bc
define abs(i) {
    if (i < 0) return (-i)
    return i
}

define min(i,j) {
    if (i < j) return i
    return j
}

define max(i,j) {
    if (i > j) return i
    return j
}

# precision of calculation
scale=${scale}

mmin=min($num1,$num2)
mmax=max($num1,$num2)
amax=max(abs($num1),abs($num2))

(mmax - mmin) > $frac * amax
EOF
    )

    if test "x$error" != "x0" ; then
        echo "Test failed: $num1 =r= $num2 with fraction $frac"
    fi

    return $error
}
