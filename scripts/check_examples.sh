#!/bin/sh
#
# Checks the Hecke eigenvalue traces computed by omf5 for the examples
# in the README against the known values.
#
# usage: scripts/check_examples.sh [path to omf5 executable]
#
# Each expected line is "k: traces", where k = 1 for T_p and k = 2 for T_{p^2}.
# Every expected line must appear among the traces printed for that example.
# The order of the forms and any further forms printed are not checked.
#
# The automorphic types of the forms are checked as well. The forms printed
# must have exactly the expected types, up to the order of the forms.

OMF5=${1:-src/omf5}

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

status=0

# check_example name expected-file expected-types omf5-arguments...
# The expected types are given in alphabetical order, such as "F G P".
check_example() {
  name=$1
  expected=$2
  types=$3
  shift 3
  if ! "$OMF5" "$@" > "$tmp/$name.out" 2> "$tmp/$name.err"; then
    echo "FAIL $name: omf5 exited with an error"
    cat "$tmp/$name.err"
    status=1
    return
  fi
  # The traces are printed on the line following the header.
  awk '/traces of hecke eigenvalues of T_\{p\^1\} are:/ {k = 1; next}
       /traces of hecke eigenvalues of T_\{p\^2\} are:/ {k = 2; next}
       k { sub(/[ \t]+$/, ""); print k ": " $0; k = 0 }' \
    "$tmp/$name.err" > "$tmp/$name.got"
  failed=0
  while IFS= read -r line; do
    if ! grep -Fxq -- "$line" "$tmp/$name.got"; then
      echo "FAIL $name: missing traces $line"
      failed=1
    fi
  done < "$expected"
  if [ $failed -ne 0 ]; then
    echo "computed traces were:"
    cat "$tmp/$name.got"
    status=1
  fi
  # The types of the forms, sorted since the order of the forms is not checked.
  got_types=$(grep -o "'aut_rep_type' : '[A-Z]'" "$tmp/$name.out" |
              cut -d "'" -f 4 | sort | tr '\n' ' ' | sed 's/ $//')
  if [ "$got_types" != "$types" ]; then
    echo "FAIL $name: the types of the forms are '$got_types', expected '$types'"
    failed=1
    status=1
  fi
  if [ $failed -eq 0 ]; then
    echo "ok   $name"
  fi
}

# Discriminant 61. Primes dividing the discriminant are skipped.
cat > "$tmp/61.exp" <<'EOF'
1: -7 -3 3 -9 -4 -3 37 -75 10 212 -6 -88 -3 547 -147 -108 -45 -632 -650 859 -978 931 -571 453
2: 7 -9 -9 -42
EOF
check_example disc61 "$tmp/61.exp" "F G P" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=100

# With -nonlifts only the form of type G is computed.
cat > "$tmp/61nl.exp" <<'EOF'
1: -7 -3 3 -9
2: 7 -9
EOF
check_example disc61-nonlifts "$tmp/61nl.exp" "G" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10 -nonlifts

# Discriminant 167, forms of conductor 1 and of conductor 167.
cat > "$tmp/167.exp" <<'EOF'
1: -2 0 -2 2 -14 -34 -15 16 155 40 -152
2: 2 -17 16
1: -3 -9 2 3 92 -41 95 -80 189 -220 -2
2: -3 12 -28
1: -8 -10 -4 -14 -22 -4 -47 -12 41 50 -504
2: 10 11 -44
EOF
check_example disc167 "$tmp/167.exp" "F G G G P" -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -prec=31

# Discriminant 249. The L-polynomial of the second form factors at 2 and at 5,
# but not at 7, so it is only found not to be a lift when 7^2 <= prec.
cat > "$tmp/249.exp" <<'EOF'
1: -1 4 -7 -21 44 167 -86 140 -121 -23 171 -128 -156 -84
2: 0 12 -24
EOF
check_example disc249 "$tmp/249.exp" "F G G G P P P" -quad=1,0,0,-1,0,1,0,0,0,1,1,1,2,-1,11 -format=GG -prec=49

# Discriminant 259, types only. The field of one of the Saito-Kurokawa lifts is
# defined by a polynomial of degree 12 with a repeated factor modulo 3, where
# Z[a] is not maximal. Testing for squares there would find it not to be a lift.
: > "$tmp/259.exp"
check_example disc259 "$tmp/259.exp" "F G G P P P P P Y Y" -quad=1,0,0,0,-1,1,0,1,1,1,-1,0,4,4,6 -format=GG

exit $status
