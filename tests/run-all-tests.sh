#!/bin/sh
countfailed=0
countall=0
# Absolute path to this script, e.g. /home/user/bin/foo.csh
SCRIPT=`readlink -f $0`
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=`dirname $SCRIPT`
TESTPATH="$SCRIPTPATH/unit"
echo "-------------------------------------------"
echo "         start: all tests for F4RT         "
echo "-------------------------------------------"
FILES=$(find $TESTPATH -name \_\* -perm /a+x | sort)
countfailedunit=0
countallunit=0
echo ""
echo "UNIT TESTS"
echo ""
for f in $FILES
do
  fname=`basename $f`
  countallunit=`expr $countallunit + 1`
  countall=`expr $countall + 1`
  $f
  if test $? -eq 0 
  then
    tput setaf 2
    echo -n "PASS: "
    tput sgr0
  else
    tput setaf 3
    echo -n "FAIL: "
    tput sgr0
    countfailedunit=`expr $countfailedunit + 1`
    countfailed=`expr $countfailed + 1`
  fi
  echo "$fname"
done
# Absolute path to this script, e.g. /home/user/bin/foo.csh
SCRIPT=`readlink -f $0`
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=`dirname $SCRIPT`
TESTPATH="$SCRIPTPATH/basic"
FILES=$(find $TESTPATH -name test-\* -perm /a+x | sort)
countfailedbasic=0
countallbasic=0
echo ""
echo ""
echo "BASIC TESTS"
echo ""
for f in $FILES
do
  fname=`basename $f`
  countallbasic=`expr $countallbasic + 1`
  countall=`expr $countall + 1`
  $f
  if test $? -eq 0 
  then
    tput setaf 2
    echo -n "PASS: "
    tput sgr0
  else
    tput setaf 3
    echo -n "FAIL: "
    tput sgr0
    countfailedbasic=`expr $countfailedbasic + 1`
    countfailed=`expr $countfailed + 1`
  fi
  echo "$fname"
done
echo ""
echo ""
printf "%3d of %3d unit  tests failed\n" "$countfailedunit" "$countallunit"
printf "%3d of %3d basic tests failed\n" "$countfailedbasic" "$countallbasic"
echo "-------------------------------------------"
printf "%3d of %3d tests failed\n" "$countfailed" "$countall"
echo "-------------------------------------------"
if test $countfailed -gt 0
then
  echo "Please report to ederc@mathematik.uni-kl.de"
  echo "-------------------------------------------"
  exit 1
fi
