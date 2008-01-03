#!/bin/bash

../run-meta-engine ../build/meta1
[ $? != 0 ] && exit 1

exit 0
