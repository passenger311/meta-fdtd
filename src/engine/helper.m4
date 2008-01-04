define({M4_QUOTE}, {ifelse({$#}, {0}, {}, {{$*}})})
define({M4_DQUOTE}, {{$@}})
define({M4_DQUOTE_ELT}, {ifelse({$#}, {0}, {}, {$#}, {1}, {{{$1}}}, {{{$1}},$0(shift($@))})})
define({M4_FOREACHQ}, {pushdef({$1})_M4_FOREACHQ($@)popdef({$1})})
define({_M4_ARG1}, {$1})
define({_M4_FOREACHQ}, {ifelse(M4_QUOTE($2), {}, {},
{define({$1}, {_M4_ARG1($2)})$3{}$0({$1}, {shift($2)}, {$3})})})