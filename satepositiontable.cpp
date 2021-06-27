#include "satepositiontable.h"

SatePositionTable::SatePositionTable(QTableWidget *parent) : QTableWidget(parent)
{
         initialization();
}
void SatePositionTable::initialization()
{   
    //创建一个tablewidget
    this->setColumnCount(7);//设置列数
    //设置表头内容
    QStringList header;
    header<<"PRN"<<"X(m)"<<"Y(m)"<<"Z(m)"<<"L"<<"B"<<"H(m)";
    this->setHorizontalHeaderLabels(header);
    //去除选中虚线框
    this->setFocusPolicy(Qt::NoFocus);
    //设置垂直头不可见
    this->verticalHeader()->setVisible(false);
    //根据内容调整列宽
    this->resizeColumnsToContents();
    //修改表格编辑状态权限
    this->setEditTriggers(QAbstractItemView::DoubleClicked);
    //行和列的大小设为与内容相匹配
    this->resizeColumnsToContents();
    this->resizeRowsToContents();
    //设置表头自动填充
    horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
      //设置表格选择方式
      this->setSelectionBehavior(QAbstractItemView::SelectRows);
    //单个选中和多个选中的设置
    this->setSelectionMode(QAbstractItemView::SingleSelection);
}
void SatePositionTable::updateTable(QMap<QString,QList<double>> table)
{

    //先清除表格中的数据
     clearALLData();
    //每个QList含6值，其顺序为：X,Y,Z,L,B,H
     QMap<QString, QList<double>>::iterator iter=table.begin();
     while (iter!= table.end())
     {
         int num=this->rowCount();
         this->setRowCount(num+1); //增加一行
         this->setItem(num,0,new QTableWidgetItem(iter.key()));  //输入PRN号
         this->setItem(num,1,new QTableWidgetItem(QString::number(table[iter.key()][0],'f',4)));  //输入X
         this->setItem(num,2,new QTableWidgetItem(QString::number(table[iter.key()][1],'f',4)));  //输入Y
         this->setItem(num,3,new QTableWidgetItem(QString::number(table[iter.key()][2],'f',4)));  //输入Z
         this->setItem(num,4,new QTableWidgetItem(QString::number(table[iter.key()][3],'f',4)));  //输入L
         this->setItem(num,5,new QTableWidgetItem(QString::number(table[iter.key()][4],'f',4)));  //输入B
         this->setItem(num,6,new QTableWidgetItem(QString::number(table[iter.key()][5],'f',4)));  //输入H
         iter++;
     }
}
void SatePositionTable::clearALLData()
{
    for(int i=this->rowCount()-1;i>=0;i--)
        this->removeRow(i);
}
SatePositionTable:: ~SatePositionTable()
{
    qDebug()<<"表格类销毁！";
}
