#!/bin/bash

../run-meta-engine ../build/meta3-cf
[ $? != 0 ] && exit 1

exit 0
