import time
import pygame
import numpy as np
import sys

BACKGROUND = (255, 255, 255)
CIRCLE = (0, 0, 0)

def get_data():
    data = open('data.txt', 'r').read()
    lines = np.array(data.split('\n'))
    Positions = []
    for line in lines:
        if len(line) > 1:
            Positions.append((float(line.split(' ')[0]), float(line.split(' ')[1])))
    return Positions

def visualize(boxSize, radius):
    pygame.init()
    Display = pygame.display.set_mode()
    windowSize = (pygame.display.Info().current_w, pygame.display.Info().current_h)

    while(1):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
        data = get_data()
        if len(data)==0:
            continue
        drawFrame(data, radius, boxSize, windowSize, Display)


def drawFrame(positions, radius, boxSize, windowSize, Display):
    Display.fill(BACKGROUND)
    offset = ((windowSize[0]-boxSize[0]-4)/2, (windowSize[1]-boxSize[1]-4)/2)
    pygame.draw.rect(Display, CIRCLE, pygame.Rect(offset[0], offset[1], boxSize[0]+5, boxSize[1]+5), 3)
    for pos in positions:
        pygame.draw.circle(Display, CIRCLE, (pos[0]+offset[0], pos[1]+offset[1]), radius)
    pygame.display.update()