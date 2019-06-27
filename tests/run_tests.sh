#!/usr/bin/env bash
if [[ $DIST == *"py3"* ]]; then
    python3 setup.py install; cd tests
    pip3 install coverage python-coveralls pyyaml
else
    python setup.py install; cd tests
    pip install coverage python-coveralls pyyaml
fi
apt-get install git
cd ..
coverage run --source=blimpy -m pytest
coverage report
coveralls --nogit