#!/usr/bin/env bash
coverage run --source=blimpy -m pytest
EXITCODE=$?
if [ $EXITCODE -ne 0 ]; then
    echo
    echo '*** Oops, coverage pytest failed, exit code = '$EXITCODE' ***'
    echo
    exit $EXITCODE
fi
coverage html -d coverage_html
EXITCODE=$?
if [ $EXITCODE -ne 0 ]; then
    echo
    echo '*** Oops, coverage report failed, exit code = '$EXITCODE' ***'
    echo
    exit $EXITCODE
fi

