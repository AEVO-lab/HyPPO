#-------------------------------------------------
#
# Project created by QtCreator 2017-04-01T09:59:44
#
#-------------------------------------------------

QT       -= core
QT       -= gui

TARGET = OCR
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    trees/genespeciestreeutil.cpp \
    trees/newicklex.cpp \
    trees/node.cpp \
    trees/treeinfo.cpp \
    trees/treeiterator.cpp \
    main_archive.cpp \
    supergenetreemaker.cpp \
    cluster.cpp \
    speciestreemaker.cpp \
    div/longcounter.cpp \
    clusterfinder.cpp


QMAKE_CXXFLAGS += -std=c++0x

HEADERS += \
    trees/genespeciestreeutil.h \
    trees/newicklex.h \
    trees/node.h \
    trees/treeinfo.h \
    trees/treeiterator.h \
    div/util.h \
    div/tinydir.h \
    supergenetreemaker.h \
    cluster.h \
    speciestreemaker.h \
    div/longcounter.h \
    clusterfinder.h


DEFINES -= UNICODE
DEFINES += _MBCS
