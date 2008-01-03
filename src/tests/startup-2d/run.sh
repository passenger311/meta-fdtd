#!/bin/bash

../run-meta-engine ../build/meta2
[ $? != 0 ] && exit 1

exit 0
