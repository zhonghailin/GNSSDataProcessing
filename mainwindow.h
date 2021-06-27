#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <rinexdata.h>   //包含星历数据
#include <QGridLayout>
#include <QHBoxLayout>
#include <QStackedLayout>
#include <stellarmap.h>  //包含星空图
#include <satepositiontable.h>  //包含卫星位置表格
namespace Ui {
class MainWindow;
}
class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
private slots:
    void on_open_file_action_triggered();
private:
    Ui::MainWindow *ui;
    NavigationObservationData *data;//this object is a singleton pattern
    SatePositionTable *tableWid;
};

#endif // MAINWINDOW_H
