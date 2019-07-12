from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import QRectF

import numpy as np

class MemberItem(QGraphicsLineItem):

    def __init__(self, member, scl, H):
        super().__init__()

        self.member = member
        self.scl = scl
        self.H = H
        self.pen = QPen("black")
        self.setAcceptHoverEvents(True)
        self.hover = False
        self.setFlag(self.acceptHoverEvents(True))





    def hoverEnterEvent(self, event:QGraphicsSceneHoverEvent):
        print("HOVER")
        self.hover = True
        self.pen.setColor("yellow")
        self.pen.setWidth(self.pen.width() + 1)
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event:QGraphicsSceneHoverEvent):
        self.hover = False
        self.pen.setColor("black")
        self.pen.setWidth(self.pen.width() - 1)
        self.update()


    def boundingRect(self):
        penWidth = self.pen.width()
        [x1, y1], [x2, y2] = np.asarray(self.member.coordinates) * self.scl

        return QRectF(x1, y1, 50, 50)


    def paint(self, painter:QPainter, option:QStyleOptionGraphicsItem, widget:QWidget=...):
        [x1, y1], [x2, y2] = np.asarray(self.member.coordinates) * self.scl
        print("PAINTING")

        # Upper left corner is (0, 0) so Y-values need to be subtracted from
        # graphics view height

        self.setPos(0, self.H/2)

        painter.drawLine(x1, y1, x2, y2)
        painter.drawRect(x1, y1, 50, 50)
        painter.fillRect(x1, y1, 50, 50, self.pen.color())



if __name__ == '__main__':

    scene = QGraphicsScene()
    view = QGraphicsView(scene)


