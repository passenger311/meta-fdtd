define({M4_VERSION},{0.0.2008.1.6})
define({M4_DATE},{2007-2008})
define({M4_AUTHORS},{J.Hamm A.Klaedtke C.Hermann S.Scholz})
define({M4_BUILD},{M4_ARCH() M4_IFELSE_DBG({DBG })M4_IFELSE_OMP({OMP })M4_IFELSE_MPI({MPI })M4_IFELSE_MPELOG({MPELOG })})
define({M4_FLAVOUR},{M4_SDIM({D })M4_IFELSE_CF({CF })M4_IFELSE_NG({NG })M4_IFELSE_TE({TE })})
