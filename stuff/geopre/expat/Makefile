ifeq ($(ARCH),)
include ../../Make.default
else 
include ../../Make.$(ARCH)
endif 

# If you know what your system's byte order is, define XML_BYTE_ORDER:
# use -DXML_BYTE_ORDER=12 for little-endian byte order;
# use -DXML_BYTE_ORDER=21 for big-endian (network) byte order.
# -DXML_NS adds support for checking of lexical aspects of XML namespaces spec
# -DXML_MIN_SIZE makes a smaller but slower parser
# -DXML_DTD adds full support for parsing DTDs

INC=-Ixmltok -Ixmlparse -DXML_NS -DXML_DTD
AR=ar
# Use one of the next two lines; unixfilemap is better if it works.
FILEMAP_OBJ=xmlwf/unixfilemap.o
#FILEMAP_OBJ=xmlwf/readfilemap.o
LIBOBJS=xmltok/xmltok.o \
  xmltok/xmlrole.o \
  xmlparse/xmlparse.o

OBJS=xmlwf/xmlwf.o \
  xmlwf/xmlfile.o \
  xmlwf/codepage.o \
  $(FILEMAP_OBJ)
LIB=xmlparse/libexpat.a
EXE=
XMLWF=xmlwf/xmlwf$(EXE)

all: $(XMLWF)

$(XMLWF): $(OBJS) $(LIB)
	@echo [LD] $@
	@$(CC) $(CFLAGS) $(INC) -o $@ $(OBJS) $(LIB)

$(LIB): $(LIBOBJS)
	@echo [AR] $@
	@$(AR) rc $(LIB) $(LIBOBJS)

clean:
	@rm -f $(OBJS) $(LIBOBJS) $(LIB) $(XMLWF)

xmltok/nametab.h: gennmtab/gennmtab$(EXE)
	@rm -f $@
	@gennmtab/gennmtab$(EXE) >$@

gennmtab/gennmtab$(EXE): gennmtab/gennmtab.c
	@echo [CC] $@
	@$(CC) $(CFLAGS) $(INC)  -o $@ gennmtab/gennmtab.c

xmltok/xmltok.o: xmltok/nametab.h

.c.o:
	@echo [CC] $@
	@$(CC) $(CFLAGS) $(INC)  -c -o $@ $<
