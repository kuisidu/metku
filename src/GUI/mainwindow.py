import sys
import numpy as np

from PySide2.QtCore import Qt, Slot
from PySide2.QtGui import QPainter, QMouseEvent, QDoubleValidator, QPen
from PySide2.QtWidgets import *
from PySide2.QtCharts import QtCharts

from frame2d.frame2d import Frame2D





class MainWindow(QMainWindow):

    def __init__(self, widget):
        super().__init__()
        self.setWindowTitle("Frame2D")

        # Menu
        self.menu = self.menuBar()
        self.file_menu = self.menu.addMenu("File")

        # Exit QAction
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.exit_app)

        # Add actions to menu
        self.file_menu.addAction(exit_action)

        # Set Central Widget
        self.setCentralWidget(widget)

    @Slot()
    def exit_app(self):
        QApplication.quit()


class TabWidget(QTabWidget):

    def __init__(self, widgets):
        super().__init__()

        for widget, tab_name in widgets.items():
            self.addTab(widget, tab_name)


class LoadTab(QWidget):

    def __init__(self):
        super().__init__()
        # Buttons
        self.add_pointload_button = QPushButton("Add Pointload")
        self.add_lineload_button = QPushButton("Add Lineload")

        # Labels

        # Line Edits
        validator = QDoubleValidator(1, 1e5, 2)
        self.bay_line_edit = QLineEdit("2000")
        self.bay_line_edit.setMaxLength(5)
        self.bay_line_edit.setValidator(validator)
        self.bay_line_edit.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.storey_line_edit = QLineEdit("2000")
        self.storey_line_edit.setMaxLength(5)
        self.storey_line_edit.setValidator(validator)
        self.storey_line_edit.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Spinboxes
        self.bays_spinbox = QSpinBox()
        self.bays_spinbox.setRange(1, 10)
        self.storeys_spinbox = QSpinBox()
        self.storeys_spinbox.setRange(1, 10)
        # Checkboxes
        self.beams_checkbox = QCheckBox()
        self.beams_checkbox.setChecked(True)
        # Load tables
        self.pointload_table = QTableWidget()
        self.pointload_table.setColumnCount(5)
        self.pointload_table.setHorizontalHeaderLabels(["Load Id", "Coordinate", "Fx [kN]", "Fy [kN]", "M [kNm]"])

        self.lineload_table = QTableWidget()
        self.lineload_table.setColumnCount(5)
        self.lineload_table.setHorizontalHeaderLabels(["Load Id", "Member", "q0 [kN/m]", "q1 [kN/m]", "Direction"])

        # Test item
        self.pointload_table.insertRow(0)
        self.pointload_table.setItem(0, 0, QTableWidgetItem("1"))
        self.pointload_table.setItem(0, 1, QTableWidgetItem("[0, 1000]"))
        self.pointload_table.setItem(0, 2, QTableWidgetItem("-10"))
        self.pointload_table.setItem(0, 3, QTableWidgetItem("0"))
        self.pointload_table.setItem(0, 4, QTableWidgetItem("0"))

        self.lineload_table.insertRow(0)
        self.lineload_table.setItem(0, 0, QTableWidgetItem("1"))
        self.lineload_table.setItem(0, 1, QTableWidgetItem("0"))
        self.lineload_table.setItem(0, 2, QTableWidgetItem("-10"))
        self.lineload_table.setItem(0, 3, QTableWidgetItem("-10"))
        self.lineload_table.setItem(0, 4, QTableWidgetItem("y"))



        # Connect buttons

        # Add widgets to grid layout
        self.layout = QGridLayout()
        self.layout.setSpacing(5)
        self.layout.addWidget(self.lineload_table)
        self.layout.addWidget(self.add_lineload_button)
        self.layout.addWidget(self.pointload_table)
        self.layout.addWidget(self.add_pointload_button)



        self.setLayout(self.layout)



class Widget(QWidget):

    def __init__(self):
        super().__init__()

        # Buttons
        self.generate_button = QPushButton("Generate")
        self.calculate_button = QPushButton("Calculate")

        # Labels
        self.bays_label = QLabel("Bays: ")
        self.bay_length_label = QLabel("Bay Lenght: ")
        self.storeys_label = QLabel("Storeys: ")
        self.storey_height_label = QLabel("Storey Height: ")
        self.beams_label = QLabel("Insert beams")
        # Line Edits
        validator = QDoubleValidator(1, 1e5, 2)
        self.bay_line_edit = QLineEdit("2000")
        self.bay_line_edit.setMaxLength(5)
        self.bay_line_edit.setValidator(validator)
        self.bay_line_edit.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.storey_line_edit = QLineEdit("2000")
        self.storey_line_edit.setMaxLength(5)
        self.storey_line_edit.setValidator(validator)
        self.storey_line_edit.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Spinboxes
        self.bays_spinbox = QSpinBox()
        self.bays_spinbox.setRange(1, 10)
        self.storeys_spinbox = QSpinBox()
        self.storeys_spinbox.setRange(1, 10)

        # Checkboxes
        self.beams_checkbox = QCheckBox()
        self.beams_checkbox.setChecked(True)

        # Connect buttons
        self.generate_button.clicked.connect(self.generate)
        self.calculate_button.clicked.connect(self.calculate)

        # Add widgets to grid layout
        self.left = QGridLayout()
        self.left.setSpacing(5)
        self.left.addWidget(self.bays_label, 0, 0)
        self.left.addWidget(self.bays_spinbox, 0, 1)
        self.left.addWidget(self.bay_length_label, 1, 0)
        self.left.addWidget(self.bay_line_edit, 1, 1)
        self.left.addWidget(self.storeys_label, 2, 0)
        self.left.addWidget(self.storeys_spinbox, 2, 1)
        self.left.addWidget(self.storey_height_label, 3, 0)
        self.left.addWidget(self.storey_line_edit, 3, 1)
        self.left.addWidget(self.beams_label, 4, 0)
        self.left.addWidget(self.beams_checkbox, 4, 1)
        self.left.addWidget(self.generate_button, self.left.rowCount() + 1, 0)
        self.left.addWidget(self.calculate_button, self.left.rowCount() - 1, 1)

        self.left.setRowStretch(self.left.rowCount() +1, 1)


        # Graphics
        self.pen = QPen(Qt.black, 1)
        self.scene = QGraphicsScene()
        self.view = QGraphicsView(self.scene)
        self.view.setMinimumSize(400, 400)


        self.right = QVBoxLayout()
        self.right.addWidget(self.view)

        # Set layout
        self.layout = QHBoxLayout()
        self.layout.addLayout(self.left)
        self.layout.addLayout(self.right)

        self.setLayout(self.layout)




    def draw_frame(self):
        self.scene.clear()
        H = self.view.height()
        L = self.view.width()
        scl = min(H / self.frame.H, L /  self.frame.L) * 0.75
        for mem in self.frame.members.values():
            size = mem.A * scl / 20
            self.pen.setWidth(size)
            [x1, y1], [x2, y2] = np.asarray(mem.coordinates) * scl

            # Upper left corner is (0, 0) so Y-values need to be subtracted from
            # QGraphicsView height
            y1 = H - y1
            y2 = H - y2


            item = QGraphicsLineItem(x1, y1, x2, y2)
            item.setPen(self.pen)
            self.scene.addItem(item)





    @Slot()
    def generate(self):
        bays = self.bays_spinbox.value()
        storeys = self.storeys_spinbox.value()
        bay_lenght = float(self.bay_line_edit.text())
        storey_height = float(self.storey_line_edit.text())
        self.frame = Frame2D(simple=[storeys, bays, storey_height, bay_lenght],
                        beams=self.beams_checkbox.checkState())
        self.frame.generate()
        self.draw_frame()




    @Slot()
    def calculate(self):
        print("Calculate pressed!")
        self.frame.calculate()



    def sceneMouse(self, event:QGraphicsSceneMouseEvent):
        print("Mouse Pressed!")


if __name__ == "__main__":
    # Qt Application
    app = QApplication(sys.argv)
    # QWidget
    widget = Widget()
    # LoadTab
    load_tab = LoadTab()
    # All tabs
    widgets = {
        widget: 'Main',
        load_tab: 'Loads'
    }
    # TabWidget
    tab_widget = TabWidget(widgets)
    # TabWidget as the central widget
    window = MainWindow(tab_widget)
    window.show()

    # Execute application
    sys.exit(app.exec_())

