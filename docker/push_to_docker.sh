#!/usr/bin/env bash
docker login -u $DOCKER_USER -p $DOCKER_PASS
export REPO=ucberkeleyseti/blimpy
export TAG=`if [ "$TRAVIS_BRANCH" == "master" ]; then echo "latest"; else echo "$TRAVIS_BRANCH" ; fi`
docker tag ${DIST}:latest $REPO:${DIST%.docker}_$TAG
docker push $REPO
