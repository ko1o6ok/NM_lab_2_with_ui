import numpy as np
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QTableWidgetItem, QLabel
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, \
    NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

from UI.infoWindow import UI_infoWindow

import ctypes
import os
import time

from table_columns import columns, extra_info_rows

ui_file = './UI/MainWindow.ui'

MAX_SIZE = 40.0 # Максимальное выводимое число разбиений


def create_plot(parent):
    parent.fig = Figure(figsize=(parent.width() / 100, parent.height() / 100))
    parent.canvas = FigureCanvas(parent.fig)
    parent.plot = parent.fig.add_subplot(projection='3d')
    return parent.plot


class UI_mainWindow(QMainWindow):
    def __init__(self):
        super(UI_mainWindow, self).__init__()
        uic.loadUi(ui_file, self)

        # создание окон для графиков
        self.plt = create_plot(self.plot_widget_1)  # Функция и сплайн на одном графике
        self.plt_PS = create_plot(self.plot_widget_2)  # График первых производных функции и сплайна

        # присвоение мест для окон
        self.plot_widget_1.canvas.setParent(self.plot_widget_1)
        self.plot_widget_2.canvas.setParent(self.plot_widget_2)

        self.tabWidget.currentChanged.connect(
            self.toolBar_changing)  # задание функционала. В данной строке: Меняет тулбар при переходе на другую вклвдку
        self.plot_toolBar = NavigationToolbar(self.plot_widget_1.canvas, self)

        self.addToolBar(self.plot_toolBar)  # создание тулбара

        #функционал кнопок
        self.plot_button.clicked.connect(
            self.plotting)  # задание функционала. В данной строке: построение графика при нажатии на кнопку "Построить"
        self.delete_plot.clicked.connect(
            self.clear_plots)  # задание функционала. В данной строке: очистка окон от ВСЕХ графиков (чистит все окна(графики и таблицу))
        #функционал списков
        #self.task_type.currentTextChanged.connect(self.standart_params)
        #self.task_num.currentTextChanged.connect(self.standart_params)
        #self.standart_params()
        self.task_index = 0
        # Названия осей
        self.plot_widget_1.plot.set_xlabel("x")
        self.plot_widget_1.plot.set_ylabel("y")
        self.plot_widget_1.plot.set_zlabel("z")

        self.plot_widget_2.plot.set_xlabel("x")
        self.plot_widget_2.plot.set_ylabel("y")
        self.plot_widget_2.plot.set_zlabel("z")

        # self.plot_widget_3.plot.set_xlabel("x")
        # self.plot_widget_3.plot.set_ylabel("Значение")

        # настройка включения второго окна
        self.info_button.triggered.connect(lambda: self.info_window("m_i.pdf"))

        # Чтобы сразу запустить
        self.input_n.setText(str("6"))
        self.input_m.setText(str("6"))
        self.Nmax.setText(str("1000"))
        self.eps_met.setText(str("0.01"))
        self.Nmax2.setText(str("1000"))
        self.eps_met2.setText(str("0.01"))

    # def standart_params(self):
    #     task_type = self.task_type.currentIndex()
    #     task_num = self.task_num.currentIndex()
    #     self.task_index =0 # Индекс задачи в массиве
    #     if task_type == 0:
    #         self.task_index = 0
    #     elif task_type == 1:
    #         self.task_index = task_num + 1
    #     elif task_type == 2:
    #         self.task_index = task_num + 4
    #     a, b, M1, M2 = task_standart_params[self.task_index]
    #     self.input_a.setText(str(a))
    #     self.input_b.setText(str(b))
    #     self.input_mu_0.setText(str(M1))
    #     self.input_mu_1.setText(str(M2))

    def info_window(self, file_name):
        self.i_window = QMainWindow()
        self.i_window.ui = UI_infoWindow(file_name)
        self.i_window.ui.show()

    def clear_plots(self):
        self.clear_plot(self.plot_widget_1)
        self.clear_plot(self.plot_widget_2)

        self.clear_table(self.info_table)
        self.clear_table(self.info_table_2)
        self.clear_table(self.info_table_3)
        self.clear_table(self.info_table_4)
        self.clear_table(self.info_table_5)
    def clear_plot(self, cur_plot_widget):
        cur_plot_widget.plot.cla()
        cur_plot_widget.canvas.draw()  # обновление окна

        # Названия осей
        cur_plot_widget.plot.set_xlabel("x")
        cur_plot_widget.plot.set_ylabel("y")
        cur_plot_widget.plot.set_zlabel("z")
        cur_plot_widget.canvas.draw()

    def toolBar_changing(self, index):  # изменение привязки тулбара
        self.removeToolBar(self.plot_toolBar)
        if index == 0:  # тулбал для вкладки # Функция и сплайн на одном графике
            self.plot_toolBar = NavigationToolbar(self.plot_widget_1.canvas, self)
        elif index == 1:  # тулбар для вкладки # График первых производных функции и сплайна
            self.plot_toolBar = NavigationToolbar(self.plot_widget_2.canvas, self)
        self.addToolBar(self.plot_toolBar)

    def file_to_table(self, file_name,n,m):  # из str делает list(list(str))
        if len(file_name.split('.')) == 1:
            file_name += '.txt'
        table = []
        if(min(n,m)<=MAX_SIZE):
            with open(file_name, 'r') as f:
                for line in f:
                    table.append(line.replace(" \n",'').split(' '))
        else:
            # Если узлов много, то выводим равномерно с шагом + обязательно границу
            step_x = int(np.floor(n/MAX_SIZE))
            step_y = int(np.floor(m/MAX_SIZE))
            y_counter = 0
            with open(file_name, 'r') as f:
                for line in f:
                    if y_counter%step_y == 0 or y_counter == (m+1) or y_counter == 1:
                        ln = line.replace(" \n",'').split(' ')
                        a = [ln[0]]
                        a_ = [ln[1]]
                        b = ln[2:n:step_x]
                        c = [ln[n+1]]

                        r = a
                        r.extend(a_)
                        r.extend(b)
                        r.extend(c)
                        table.append(r)
                    y_counter +=1
        return table

    def clear_exrta_info_table(self):
        while self.extra_info_layout.count():
            item = self.extra_info_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def update_extra_info_table(self, task_index, table):
        self.clear_exrta_info_table()

        table = table[0]
        i = 0
        cur_table = extra_info_rows[task_index]
        for elem in table:
            cur_text = f"{cur_table[i]} {elem}"
            self.extra_info_layout.addWidget(QLabel(cur_text, self))
            i += 1
    def prepare_for_3d_drawing(self,table):
        # 0 x0 x1 x2 ... xn
        # y1
        # y2
        # ...
        x_coords = [float(x) for x in table[0][1:]]

        y_coords = []
        z_coords = []
        for el in table[1:]:
            y_coords.append(float(el[0]))
            z_coords.append([float(val) for val in el[1:]])

        x_coords = np.asarray(x_coords)
        y_coords = np.asarray(y_coords)
        z_coords = np.asarray(z_coords)
        return x_coords, y_coords, z_coords
    def plotting(self):
        lib = ctypes.windll.LoadLibrary("libNM1_lib.dll")

        task_type = int(self.task_type.currentIndex())

        if(task_type == 0):
            my_func = lib.solve_test
            my_func.argtypes = [ctypes.c_int, ctypes.c_int,
                                ctypes.c_int, ctypes.c_double,
                                ]
            my_func.restype = ctypes.c_void_p

            n = int(self.input_n.text())
            m = int(self.input_m.text())
            Nmax = int(self.Nmax.text())
            eps_met = float(self.eps_met.text())

            my_func(n,m,Nmax,eps_met)
            #time.sleep(0.01)

            file_name_1 = "nach_priblizh_test"
            file_name_2 = "tochn_resh_test"

            file_name_3 = "chisl_resh_test"
            file_name_4 = "diff_resh_test"
            file_name_extra_info = 'spravka_test'

            self.toolBox.setItemText(0, "Начальное приближение")
            self.toolBox.setItemText(1, "Точное решение")
            self.toolBox.setItemText(2, "Численное решение")
            self.toolBox.setItemText(3, "Разность решений")

            self.clear_table(self.info_table)
            table_1 = self.file_to_table(file_name_1,n,m)
            self.set_table(self.info_table, table_1, file_name_1,n)

            x1, y1, z1 = self.prepare_for_3d_drawing(table_1)
            X1, Y1 = np.meshgrid(x1, y1)
            surface3 = self.plt_PS.plot_surface(X1, Y1, z1, cmap="plasma", label="Начальное приближение")

            self.clear_table(self.info_table_2)
            table_2 = self.file_to_table(file_name_2,n,m)
            self.set_table(self.info_table_2, table_2, file_name_2,n)

            x2,y2,z2 = self.prepare_for_3d_drawing(table_2)
            X2,Y2 = np.meshgrid(x2,y2)
            surface1 = self.plt.plot_surface(X2,Y2,z2,color = "r",label = "Точное решение")

            self.clear_table(self.info_table_3)
            table_3 = self.file_to_table(file_name_3,n,m)
            self.set_table(self.info_table_3, table_3, file_name_3,n)

            x3, y3, z3 = self.prepare_for_3d_drawing(table_3)
            X3, Y3 = np.meshgrid(x3, y3)
            surface2 = self.plt.plot_surface(X3, Y3, z3,cmap="plasma", label = "Численное решение")

            # Добавление цветовой шкалы для каждой поверхности
            # self.plt.colorbar(surface1, shrink=0.5, aspect=5, label='Первая поверхность')
            # self.plt.colorbar(surface2, shrink=0.5, aspect=5, label='Вторая поверхность')

            self.clear_table(self.info_table_4)
            table_4 = self.file_to_table(file_name_4,n,m)
            self.set_table(self.info_table_4, table_4, file_name_4,n)

            x4, y4, z4 = self.prepare_for_3d_drawing(table_4)
            X4, Y4 = np.meshgrid(x4, y4)
            surface3 = self.plt.plot_surface(X4, Y4, z4,cmap = "plasma",label="Разность решений")


        else:
            my_func = lib.solve_main
            my_func.argtypes = [ctypes.c_int, ctypes.c_int,
                                ctypes.c_int, ctypes.c_double,
                                ctypes.c_int, ctypes.c_double,
                                ]
            my_func.restype = ctypes.c_void_p

            n = int(self.input_n.text())
            m = int(self.input_m.text())
            Nmax = int(self.Nmax.text())
            eps_met = float(self.eps_met.text())
            Nmax2 = int(self.Nmax2.text())
            eps_met2 = float(self.eps_met2.text())

            my_func(n, m, Nmax, eps_met,Nmax2,eps_met2)
            #time.sleep(0.01)

            file_name_1 = "nach_priblizh1_main"
            file_name_2 = "nach_priblizh2_main"
            file_name_3 = "chisl_resh1_main"
            file_name_4 = "chisl_resh2_main"
            file_name_5 = "diff_resh_main"
            file_name_extra_info = 'spravka_main'

            self.toolBox.setItemText(0, "Начальное приближение на осн. сетке")
            self.toolBox.setItemText(1, "Начальное приближение на контр. сетке")
            self.toolBox.setItemText(2, "Численное решение на осн. сетке")
            self.toolBox.setItemText(3, "Численное решение на контр. сетке")
            self.toolBox.setItemText(4, "Разность решений")

            self.clear_table(self.info_table)
            table_1 = self.file_to_table(file_name_1,n,m)
            self.set_table(self.info_table, table_1, file_name_1,n)

            x1, y1, z1 = self.prepare_for_3d_drawing(table_1)
            X1, Y1 = np.meshgrid(x1, y1)
            surface1 = self.plt_PS.plot_surface(X1, Y1, z1, color = 'r',label="Начальное приближение на осн. сетке")

            self.clear_table(self.info_table_2)
            table_2 = self.file_to_table(file_name_2,n,m)
            self.set_table(self.info_table_2, table_2, file_name_2,n)

            x2, y2, z2 = self.prepare_for_3d_drawing(table_2)
            X2, Y2 = np.meshgrid(x2, y2)
            surface2 = self.plt_PS.plot_surface(X2, Y2, z2, cmap="plasma",label="Начальное приближение на контр. сетке")

            self.clear_table(self.info_table_3)
            table_3 = self.file_to_table(file_name_3,n,m)
            self.set_table(self.info_table_3, table_3, file_name_3,n)

            x3, y3, z3 = self.prepare_for_3d_drawing(table_3)
            X3, Y3 = np.meshgrid(x3, y3)
            surface3 = self.plt.plot_surface(X3, Y3, z3, color = 'r',label="Числ. реш. на осн. сетке")

            self.clear_table(self.info_table_4)
            table_4 = self.file_to_table(file_name_4,n,m)
            self.set_table(self.info_table_4, table_4, file_name_4,n)

            x4, y4, z4 = self.prepare_for_3d_drawing(table_4)
            X4, Y4 = np.meshgrid(x4, y4)
            surface4 = self.plt.plot_surface(X4, Y4, z4, cmap="plasma",label="Числ. реш. на контр. сетке")

            self.clear_table(self.info_table_5)
            table_5 = self.file_to_table(file_name_5,n,m)
            self.set_table(self.info_table_5, table_5, file_name_5,n)

            x5, y5, z5 = self.prepare_for_3d_drawing(table_5)
            X5, Y5 = np.meshgrid(x5, y5)
            surface5 = self.plt.plot_surface(X5, Y5, z5, cmap = "plasma", label="Разность решений")



        table_extra_info = self.file_to_table(file_name_extra_info,0,0)

        self.update_extra_info_table(file_name_extra_info,
                                     table_extra_info)  # заполнение вспомогательной информации(правый нижний угол)
        self.plot_widget_1.canvas.draw()
        self.plot_widget_2.canvas.draw()

    def set_row(self, table, row):
        max_row_index = table.rowCount()
        table.insertRow(max_row_index)  # создание строки
        for i in range(len(row)):
            table.setItem(max_row_index, i, QTableWidgetItem(str(row[i])))  # заполнение элементами

    def set_columns(self, table, task_index,n):
        cols = [str(k) for k in range(1,n+3)]
        table.setColumnCount(n+2)  # создание пустых колонок, в количестве len(cols) штук
        table.setHorizontalHeaderLabels(cols)  # присвоение имен для колонок

    def set_table(self, table, data, task_index,n):
        self.set_columns(table, task_index,len(data[0])-2)
        for row in data:
            self.set_row(table, row)

    def clear_table(self, table):
        while (table.rowCount() > 0):
            table.removeRow(0)
