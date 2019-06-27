#!/usr/bin/env bash
if [[ $DIST == *"py3"* ]]; then
    python3 setup.py install; cd tests
    pip3 install coverage coveralls
else
    python setup.py install; cd tests
    pip install coverage coveralls
fi
coverage run --source=blimpy -m py.test
coverage report
ls -a