#!/bin/bash
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Checks out previous version of the doc from the gh-pages branch, and
# compares with the newly generated doc while ignoring generated time stamps.
# Since the index.html file now also contains the build date in the heading,
# we need to exclude this file from the diff (assuming it is rarely changed).
# Therefore moving this file around.

git checkout gh-pages
mv docs/solver/index.html index.old
mv doc/solver_html/index.html index.new
git rm -rf docs/solver
mv doc/solver_html docs/solver
mv index.old docs/solver/index.html
git add docs/solver
if git diff -I "^Generated on .* for FEDEM" HEAD --quiet; then
  echo No changes in source code documentation
  git reset --hard HEAD
  exit 0
else
  mv index.new docs/solver/index.html
  git add docs/solver/index.html
  exit 1
fi
