#!/bin/bash

../run-meta-engine ../build/meta2-cf
[ $? != 0 ] && exit 1

exit 0
