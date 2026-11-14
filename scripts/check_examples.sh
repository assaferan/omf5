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

# check_error name message omf5-arguments...
# Checks that omf5 exits with an error for these arguments, and reports it
# with a line containing the message. Being killed by a signal, as after
# writing past the end of an array, is not exiting with an error.
check_error() {
  name=$1
  message=$2
  shift 2
  "$OMF5" "$@" > "$tmp/$name.out" 2> "$tmp/$name.err"
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "FAIL $name: omf5 did not exit with an error"
    status=1
  elif [ $rc -gt 128 ] && [ $rc -ne 255 ]; then
    echo "FAIL $name: omf5 was killed by signal $((rc - 128))"
    status=1
  elif ! grep -Fq -- "$message" "$tmp/$name.err"; then
    echo "FAIL $name: the error '$message' was not reported"
    status=1
  else
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

# With -idxs only the forms at the given indices are computed, as type G.
# Here they are the form of type G and the Saito-Kurokawa lift.
cat > "$tmp/61idxs.exp" <<'EOF'
1: -7 -3 3 -9
2: 7 -9
1: 29 59 164 309
2: -3 -4
EOF
check_example disc61-idxs "$tmp/61idxs.exp" "G G" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10 -idxs=0,2

# The indices can be given in any order, and can be repeated. An index at
# which there is no form is ignored.
check_example disc61-idxs-unsorted "$tmp/61idxs.exp" "G G" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10 -idxs=2,0,2,123456789

# An entry which is not an integer is an error. It used to be read as 0.
check_error disc61-idxs-malformed "is not a comma-separated list of integers" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10 -idxs=2,abc
check_error disc61-idxs-empty "is not a comma-separated list of integers" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10 -idxs=

# A matrix which is not given by 15 integers is an error. More than 15 entries
# used to be written past the end of the array, and an entry which is not an
# integer used to be read as 0.
check_error quad-too-long "is not a comma-separated list of 15 integers" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8,1,0,0,1,1,1,0,1,0,1,0,0,1,0,8,1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=10
check_error quad-malformed "is not a comma-separated list of 15 integers" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,x -format=GG -prec=10

# So is a matrix which is not that of a positive definite even lattice. Such a
# matrix used to be accepted, and the computation would not end.
check_error quad-indefinite "the lattice is not positive definite" -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,-8 -format=GG -prec=10
check_error quad-semidefinite "the lattice is not positive definite" -quad=2,0,2,0,0,2,0,0,0,2,0,0,0,0,0 -format=A -prec=10
check_error quad-odd "the lattice is not even" -quad=2,1,2,0,0,2,0,0,0,2,0,0,0,0,3 -format=A -prec=10
# The entries of the Gram matrix have to fit in 32 bits. A larger entry used to
# be truncated, and a coefficient whose double does not fit in 64 bits was read
# as another one, with no error.
check_error quad-large "the entry 4294967296 of the Gram matrix of a lattice does not fit in 32 bits" -quad=1,0,0,0,0,1,0,0,0,1,0,0,1,0,2147483648 -format=GG -prec=10
check_error quad-overflow "the coefficient -9223372036854775807 of a lattice is too large" -quad=1,0,0,0,0,1,0,0,0,1,0,0,1,0,-9223372036854775807 -format=GG -prec=10

# A genus file which cannot be opened is an error, with or without isometries.
check_error genus-missing "could not open the genus file" -genus="$tmp/no-such-file" -disc=61 -prec=10
check_error genus-isom-missing "could not open the genus file" -genus="$tmp/no-such-file" -isom -disc=61 -prec=10

# So is a lattice in the genus file which is not given by 15 integers.
# Here the genus of discriminant 1 has a lattice with 14 entries.
echo '[[],[[1,0,0,1,1,1,0,1,0,1,0,0,1,0,8],[1,0,0,1,1,1,0,1,0,1,0,0,1,0]]]' > "$tmp/bad.db"
check_error genus-malformed "of discriminant 1 is not a list of 15 integers" -genus="$tmp/bad.db" -disc=1 -prec=10
# The same when it is not the last lattice of its genus.
echo '[[],[[1,0,0,1,1,1,0,1,0,1,0,0,1,0],[1,0,0,1,1,1,0,1,0,1,0,0,1,0,8]]]' > "$tmp/bad-first.db"
check_error genus-malformed-first "of discriminant 1 is not a list of 15 integers" -genus="$tmp/bad-first.db" -disc=1 -prec=10

# Discriminant 44, where the genus has two classes, read from a genus file.
# The genus read is the one at the index given by -disc, here after an empty
# genus and a genus which is not complete. The other genera are not checked.
cat > "$tmp/44.exp" <<'EOF'
1: 40 156 400 2380 5220 7240 12720 25260
2: 120 780
1: 5 11 70 110 260 360 445 990
2: -20 -90
1: 7 23 30 234 352 284 579 886
2: -12 -18
EOF
A=[1,0,-1,-1,1,1,-1,1,0,2,1,-1,2,-2,2]
B=[1,-1,0,1,1,1,0,-1,0,1,1,0,2,1,3]
echo "[[],[$B],[$A,$B],[$A]]" > "$tmp/44.db"
check_example genus44 "$tmp/44.exp" "F P P" -genus="$tmp/44.db" -disc=2 -prec=30
# The genus is the same when it is computed from its first lattice.
check_example quad44 "$tmp/44.exp" "F P P" -quad=1,0,-1,-1,1,1,-1,1,0,2,1,-1,2,-2,2 -format=GG -prec=30

# genus_error name message genus
# Checks that the error is reported for a genus file with this genus.
genus_error() {
  echo "[[],$3]" > "$tmp/$1.db"
  check_error "$1" "$2" -genus="$tmp/$1.db" -disc=1 -prec=10
}

# The file has to have a genus at the index given by -disc, and nothing but
# lattices in it. A genus which was not in the file used to be read as an
# empty one.
check_error genus-no-disc "has no genus of discriminant 4" -genus="$tmp/44.db" -disc=4 -prec=10
printf '[[],[%s,[1,-1,0' "$A" > "$tmp/truncated.db"
check_error genus-truncated "ends inside the genus of discriminant 1" -genus="$tmp/truncated.db" -disc=1 -prec=10
genus_error genus-stray "the character 'x' is between the lattices of discriminant 1" "[$A x $B]"
# A lattice longer than the buffer used to be written past its end.
genus_error genus-long "a lattice of discriminant 1 is longer than" "[[$(awk 'BEGIN { for (i = 0; i < 3000; i++) printf "1,"; print "1" }')]]"

# The same genus, read from a file with a Gram matrix and an isometry from the
# first lattice for every lattice, as in qf5db.sage.
MA='[[2, -1, 0, 1, 1], [-1, 2, 0, -1, 0], [0, 0, 2, 1, 0], [1, -1, 1, 4, 1], [1, 0, 0, 1, 6]]'
MB='[[2, 0, -1, -1, 1], [0, 2, -1, 1, 0], [-1, -1, 4, 1, -1], [-1, 1, 1, 4, -2], [1, 0, -1, -2, 4]]'
IA='[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]'
IB='[[1, 0, 0, -1, 1], [2/3, 1, -2/3, 1/3, 1/3], [-1/3, 0, 4/3, 1/3, -2/3], [1/3, 0, -1/3, 2/3, -1/3], [-1/3, 0, 1/3, 1/3, 1/3]]'
echo "[[],[[$MA,$IA],[$MB,$IB]]]" > "$tmp/44isom.db"
check_example genus44-isom "$tmp/44.exp" "F P P" -genus="$tmp/44isom.db" -isom -disc=1 -prec=30
# The white space in the file is not needed. It used to be.
tr -d ' ' < "$tmp/44isom.db" > "$tmp/44isom-compact.db"
check_example genus44-isom-compact "$tmp/44.exp" "F P P" -genus="$tmp/44isom-compact.db" -isom -disc=1 -prec=30
# The fractions do not have to be reduced. The denominators of isometries are
# multiplied when the isometries are composed, so they are reduced when read.
# They used not to be, and the traces of one of the forms were wrong.
echo "[[],[[$MA,$(echo "$IA" | sed 's|1|2147483647/2147483647|')],[$MB,$(echo "$IB" | sed 's|2/3|1000000/1500000|')]]]" > "$tmp/44isom-unreduced.db"
check_example genus44-isom-unreduced "$tmp/44.exp" "F P P" -genus="$tmp/44isom-unreduced.db" -isom -disc=1 -prec=30

# genus_isom_error name message genus
# The same as genus_error, for a genus file with isometries.
genus_isom_error() {
  echo "[[],$3]" > "$tmp/$1.db"
  check_error "$1" "$2" -genus="$tmp/$1.db" -isom -disc=1 -prec=10
}

# A lattice which is not a Gram matrix with an isometry is an error, as is a
# denominator which is not positive. These used to be read past the end of the
# text, or to be accepted.
not_isom="is not a symmetric 5x5 matrix of integers with a 5x5 matrix of rational numbers"
genus_isom_error genus-isom-plain "$not_isom" "[$A,$B]"
genus_isom_error genus-isom-no-isom "$not_isom" "[[$MA],[$MB]]"
genus_isom_error genus-isom-not-symmetric "$not_isom" "[[$MA,$IA],[$(echo "$MB" | sed 's/\[2, 0, -1, -1, 1\]/[2, 1, -1, -1, 1]/'),$IB]]"
genus_isom_error genus-isom-zero-denom "$not_isom" "[[$MA,$IA],[$MB,$(echo "$IB" | sed 's|2/3|2/0|')]]"
genus_isom_error genus-isom-negative-denom "$not_isom" "[[$MA,$IA],[$MB,$(echo "$IB" | sed 's|-2/3|2/-3|')]]"
genus_isom_error genus-isom-four-rows "$not_isom" "[[$MA,$IA],[$MB,[[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0]]]]"
# A file with isometries read without -isom is an error as well.
check_error genus-isom-as-plain "of discriminant 1 is not a list of 15 integers" -genus="$tmp/44isom.db" -disc=1 -prec=10

# The lattices read have to be the classes of a genus. A genus which was not
# complete used to give a smaller space of forms, or a crash.
genus_error genus-incomplete "have mass 1/16, and the mass of the genus of the first lattice is 5/48" "[$A]"
genus_error genus-repeated "at index 2 of the genus is isometric to one of the lattices before it" "[$A,$B,$A]"
genus_error genus-other-disc "at index 1 of the genus has discriminant 61, and the first lattice has discriminant 44" "[$A,[1,0,0,1,1,1,0,1,0,1,0,0,1,0,8]]"
genus_error genus-indefinite "at index 1 of the genus is not positive definite" "[$A,[1,0,0,0,0,1,0,0,0,1,0,0,1,0,-1]]"
# This lattice has the discriminant and the number of automorphisms of the
# second class, so the mass is that of the genus, but it is in another genus.
genus_error genus-other-genus "at index 1 of the genus is not in the genus of the first lattice" "[$A,[1,1,0,0,-1,1,0,0,0,1,1,1,2,2,3]]"
# The spaces of an empty genus are found from its discriminant, which cannot
# be 0. With an empty genus of discriminant 0 the program used not to end.
echo "[[]]" > "$tmp/empty.db"
check_error genus-empty-zero "the genus of discriminant 0 is empty" -genus="$tmp/empty.db" -disc=0 -prec=10

# The same with isometries, where the first lattice is the other class.
genus_isom_error genus-isom-incomplete "have mass 1/24, and the mass of the genus of the first lattice is 5/48" "[[$MA,$IA]]"
genus_isom_error genus-isom-repeated "at index 2 of the genus is isometric to one of the lattices before it" "[[$MA,$IA],[$MB,$IB],[$MA,$IA]]"
# The Gram matrix of the standard lattice, which is not even.
genus_isom_error genus-isom-odd "at index 0 of the genus is not even" "[[$IA,$IA]]"
# An isometry which is not one from the first lattice used to be accepted.
genus_isom_error genus-isom-wrong "the isometry at index 1 of the genus is not an isometry from the first lattice" "[[$MA,$IA],[$MB,$IA]]"
# This one is 14/13 times the identity. The Gram matrix it gives is that of the
# first lattice only after it is rounded.
genus_isom_error genus-isom-rounded "the isometry at index 0 of the genus is not an isometry from the first lattice" "[[$MA,$(echo "$IA" | sed 's|1|14/13|g')],[$MB,$IB]]"
# The spinor norm of an isometry is found modulo the primes which divide the
# discriminant, so the denominator has to be prime to it. This is an isometry
# to the second lattice with the denominator 6. It used to end in an error of
# FLINT.
IH='[[1/2, -1/2, -1, -3/2, 1], [1/6, 1/2, -5/3, -1/6, 1/3], [-5/6, -1/2, 1/3, -1/6, -2/3], [1/3, 0, -1/3, 2/3, -1/3], [-1/3, 0, 1/3, 1/3, 1/3]]'
genus_isom_error genus-isom-denom "the denominator 6 of the isometry at index 1 of the genus is not prime to the discriminant" "[[$MA,$IA],[$MB,$IH]]"

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

# Discriminant 48, where the lattice is not maximal at 2. The traces are those
# of the Eisenstein series, 1 + p + p^2 + p^3 for T_p.
cat > "$tmp/48.exp" <<'EOF'
1: 156 400 1464 2380 5220 7240 12720 25260
2: 780
EOF
check_example disc48 "$tmp/48.exp" "F P P P" -quad=1,0,0,1,1,1,-1,0,-1,1,0,1,1,1,6 -format=GG -prec=30

exit $status
