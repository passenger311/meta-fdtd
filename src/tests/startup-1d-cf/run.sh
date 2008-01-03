#!/bin/bash

../run-meta-engine ../build/meta1-cf
[ $? != 0 ] && exit 1

exit 0
