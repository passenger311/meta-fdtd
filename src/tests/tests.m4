define(`M4_RUN_ENGINE',`../run-meta-engine $1; ERR=$?; [ $ERR != 0 ] && exit $ERR')
define(`M4_EXIT_CHECK',`ERR=$?; [ $ERR != 0 ] && exit $ERR')
define(`M4_EXIT_OK',`exit 0')
define(`M4_COMP_DATA',`../compdata.pl; [ $? != 0 ] && exit 1')
