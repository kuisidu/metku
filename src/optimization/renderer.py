import pygame
import numpy as np
from src.frame2d.frame2d import *
from src.truss2d import Truss2D
from src.sections.steel.catalogue import *

PROFILES = list(rhs_profiles.keys())

GREEN = (0, 255, 0)
RED = (255, 0, 0)



frame = Frame2D(simple=[2,2, 1e3, 1e3], supports='fixed', num_elements=10)

for mem in frame.members.values():
    if mem.mtype == 'beam':
        frame.add(LineLoad(mem, [-1e3, -1e3], 'y'))

frame.add(PointLoad([0, 1e3], [2e5, 0, 0]))
frame.generate()
#frame.hinge_joints()
frame.calculate()


WIDTH = 500
HEIGHT = 500
pygame.init()

SCALE = frame.L / WIDTH * 2

HEIGHT = max(HEIGHT, int(frame.H / SCALE) + 100)

screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Frame2D")

pressed = 0
RUN = True

def draw_frame():

    # Members
    for mem in frame.members.values():
        start, end = np.asarray(mem.coordinates) / SCALE
        start = np.array([3/2*frame.L/SCALE , HEIGHT - 30]) - start
        end = np.array([3/2*frame.L/SCALE , HEIGHT - 30]) - end
        pygame.draw.line(screen, (0, 0, 0), start, end,3)

        # Hinges
        if  mem.Sj1 <= 1e-4:
            size = 50 / SCALE
            pygame.draw.ellipse(screen,  (0,0,0), (start[0]-size/2, start[1] - size / 2, size, size))
            pygame.draw.ellipse(screen, (255, 255, 255), (start[0] - size / 4, start[1] - size / 4, size/2, size/2))

        if  mem.Sj2 <= 1e-4:
            size = 50 / SCALE
            pygame.draw.ellipse(screen, (0, 0, 0), (end[0] - size / 2, end[1] - size / 2, size, size))
            pygame.draw.ellipse(screen, (255, 255, 255), (end[0] - size / 4, end[1] - size / 4, size / 2, size / 2))


    # Supports
    for sup in frame.supports.values():
        rect_size = 60 / SCALE
        x, y = np.asarray(sup.coordinate) / SCALE
        x  = 3/2*frame.L/SCALE - x - rect_size/2
        y = HEIGHT - 30 - y
        rect = pygame.Rect(x, y, rect_size, rect_size)
        pygame.draw.rect(screen , (0,0,0), rect)


def draw_deflected(scl=1):

    for mem in frame.members.values():

        if max(mem.r[:-1]) < 1:
            color = GREEN
        else:
            color = RED

        x_vals = np.asarray([coord[0] for coord in mem.nodal_coordinates])
        y_vals = np.asarray([coord[1] for coord in mem.nodal_coordinates])


        x_disp_vals = np.asarray([disp[0] for disp in mem.nodal_displacements.values()])
        y_disp_vals = np.asarray([disp[1] for disp in mem.nodal_displacements.values()])


        x_vals -= x_disp_vals
        y_vals += y_disp_vals

        x_vals /= SCALE
        y_vals /= SCALE

        x_vals = 3/2*frame.L/SCALE - x_vals
        y_vals = HEIGHT - 30 - y_vals

        pointlist = []

        for x, y in zip(x_vals, y_vals):
            pointlist.append([int(x), int(y)])
        pygame.draw.lines(screen, color, False,  pointlist, 2)



def update_frame(i=[]):

    idx = len(i)

    if idx < len(PROFILES):
        for mem in frame.members.values():
            mem.profile = PROFILES[idx]
        frame.calculate()
        i.append(i)






while RUN:
    screen.fill((255, 255, 255))
    pygame.time.delay(300)
    if pressed:
        pressed = 0

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            RUN = False

    draw_frame()
    draw_deflected()
    update_frame()

    pygame.display.flip()

pygame.quit()
