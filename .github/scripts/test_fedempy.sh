#!/bin/bash
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Runs all (or most of) the fedempy regression tests.

set -e

if [ $# -lt 1 ]; then echo "usage: $0 <test-directory>"; exit 1; fi
if [ ! -d $1 ]; then echo "$0: No such directory: $1"; exit 1; fi

FEDEMPY_ROOT=$PWD/PythonAPI
if [ -d $FEDEMPY_ROOT/PythonAPITests ]; then
  export PYTHONPATH=$FEDEMPY_ROOT/src:$FEDEMPY_ROOT/PythonAPITests
  cd $FEDEMPY_ROOT/PythonAPITests
  TEST_DIR=$1/TimeDomain python3 -B test_api.py
  TEST_DIR=$1/TimeDomain python3 -B test_extfunc.py
  TEST_DIR=$1/TimeDomain python3 -B test_Microbatch.py
fi

export PYTHONPATH=$FEDEMPY_ROOT/src:$1/InversePy
for file in `find $1 -name run_API.py | sort`; do
  cd `dirname $file`; pwd
  python3 -B run_API.py --no-print
done
