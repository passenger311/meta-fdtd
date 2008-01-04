define(`M4_RUN_ENGINE',`../run-meta-engine $1')
define(`M4_EXIT_CHECK',`ERR=$?; [ $ERR != 0 ] && exit $ERR')
define(`M4_EXIT_OK',`exit 0')
