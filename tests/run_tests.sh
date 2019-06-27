#!/usr/bin/env bash
python setup.py install; cd tests
pip install coverage coveralls
coverage run --source=blimpy -m py.test
coverage report
ls