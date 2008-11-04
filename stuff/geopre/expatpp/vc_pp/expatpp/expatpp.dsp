# Microsoft Developer Studio Project File - Name="expatpp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=expatpp - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "expatpp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "expatpp.mak" CFG="expatpp - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "expatpp - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "expatpp - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "\oofile\expatpp\expat\xmlparse" /I "\oofile\expatpp\expat\xmltok" /I "\oofile\expatpp\src_pp" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D XMLTOKAPI="" /D XMLPARSEAPI="" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_Release"
# PROP BASE Intermediate_Dir "expatpp___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /I "\oofile\expatpp\expat\xmlparse" /I "\oofile\expatpp\expat\xmltok" /I "\oofile\expatpp\src_pp" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D XMLTOKAPI="" /D XMLPARSEAPI="" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /ML /W3 /GX /O2 /I "\oofile\expatpp\expat\xmlparse" /I "\oofile\expatpp\expat\xmltok" /I "\oofile\expatpp\src_pp" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D XMLTOKAPI="" /D XMLPARSEAPI="" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "expatpp - Win32 Debug"
# Name "expatpp - Win32 Release"
# Begin Group "xmlparse files"

# PROP Default_Filter ""
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Source File

SOURCE=..\..\expat\xmlparse\xmlparse.c
# End Source File
# End Group
# Begin Group "xmltok files"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\expat\xmltok\xmlrole.c
# End Source File
# Begin Source File

SOURCE=..\..\expat\xmltok\xmltok.c
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\src_pp\expatpp.cpp
# End Source File
# End Target
# End Project
