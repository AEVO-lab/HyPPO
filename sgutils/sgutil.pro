#-------------------------------------------------
#
# Project created by QtCreator 2016-05-25T20:39:09
#
#-------------------------------------------------

#QT       -= core
QT       -= gui

TARGET = sgutils
CONFIG   += console
CONFIG   -= app_bundle

QMAKE_CXXFLAGS += -std=c++0x

TEMPLATE = app


SOURCES += main.cpp \
    trees/newicklex.cpp \
    trees/node.cpp \
    trees/treeinfo.cpp \
    trees/treeiterator.cpp \
    trees/genespeciestreeutil.cpp

HEADERS += \
    trees/autester.h \
    trees/newicklex.h \
    trees/node.h \
    trees/treeinfo.h \
    trees/treeiterator.h \
    div/tinydir.h \
    div/util.h \
    trees/genespeciestreeutil.h


DEFINES -= UNICODE
DEFINES += _MBCS
