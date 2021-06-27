#-------------------------------------------------
#
# Project created by QtCreator 2021-05-04T13:31:34
#
#-------------------------------------------------

QT       += core gui
QT       += serialport multimedia network
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
#踩过的坑之，必须包含printsupport，否侧应用第三方的绘图库qcustomplot报错
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport
TARGET = GNSSDateProcessing
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    rinexdata.cpp \
    stellarmap.cpp \
    qcustomplot.cpp \
    satepositiontable.cpp

HEADERS += \
        mainwindow.h \
    rinexdata.h \
    stellarmap.h \
    qcustomplot.h \
    satepositiontable.h

FORMS += \
        mainwindow.ui
QT += concurrent testlib #使用多线程
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
