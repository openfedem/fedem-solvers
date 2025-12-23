#!/bin/bash
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Runs all (or most of) the fedempy regression tests.

if [ $# -lt 1 ]; then exec echo "usage: $0 <test-directory>"; fi
if [ ! -d $1 ]; then exec echo "$0: No such directory: $1"; fi

if [ -d PythonAPI/PythonAPITests ]; then
  cd PythonAPI/PythonAPITests
  TEST_DIR=$1/TimeDomain python3 -B test_api.py
  TEST_DIR=$1/TimeDomain python3 -B test_extfunc.py
  TEST_DIR=$1/TimeDomain python3 -B test_Microbatch.py
fi

for file in `find $1 -name run_API.py`; do
  cd `dirname $file`; pwd
  python3 -B run_API.py --no-print
done
