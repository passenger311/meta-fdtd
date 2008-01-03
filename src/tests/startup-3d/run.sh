#!/bin/bash

../run-meta-engine ../build/meta3
[ $? != 0 ] && exit 1

exit 0
