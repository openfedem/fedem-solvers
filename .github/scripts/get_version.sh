#!/bin/bash
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Small helper to compose the current version string for FEDEM.
# $1 = the file holding the current build number (using the last digit only)
# $2 = the file holding the current Fedem version

if [ $# -lt 2 ]; then exec echo "usage: $0 <VERSION> <version.h>"; fi
if [ ! -e $1 ]; then exec echo "$0: No such file: $1"; fi
if [ ! -e $2 ]; then exec echo "$0: No such file: $2"; fi

buildno=`sed "s/^.*\.//" $1`
sed "1 s/\"//g;1 s/.*$/& (build $buildno)/" $2
