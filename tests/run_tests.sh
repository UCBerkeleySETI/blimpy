#!/usr/bin/env bash
python setup.py install; cd tests
echo 'running'
coverage run --source=blimpy -m py.test
coverage report
