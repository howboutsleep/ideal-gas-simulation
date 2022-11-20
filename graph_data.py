import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
from distfit import distfit
import numpy as np
from random import random

HIST_BINS = 50

def get_data():
    data = open('data.txt', 'r').read()
    lines = np.array(data.split('\n'))
    if len(lines) < 2:
        return ([], 0)
    time = float(lines[0])
    Vels = []
    for line in lines[1:]:
        if len(line) > 1:
            Vels.append(float(line.split(' ')[0]))
    return (np.array(Vels), time)

def graph(vel_lim):
    def update(i):
        GD = get_data()
        data, time = GD[0],GD[1]
        if not time:
            return updatable
        len_data = len(data)
        n, _ = np.histogram(data, HIST_BINS)
        for count, rect in zip(n, bar_container):
            rect.set_height(count/len_data*100)
        time_text.set_text(f'time = {int(time)}ms')
        return updatable

    style.use("fivethirtyeight")
    data = get_data()[0]
    n, _ = np.histogram(data, HIST_BINS)
    fig, ax = plt.subplots()

    n, bins, bar_container = ax.hist(data, HIST_BINS,
                                  ec="yellow", fc="green", alpha=0.5)

    plt.xlabel('Velocity')
    plt.ylabel('Probability (%)')
    plt.title('Histogram of velocities')
    plt.grid(True)
    ax.set_ylim(top=10)
    ax.set_xlim(right=vel_lim)
    line = ax.plot([], [], color='red', alpha=0.75, linewidth=1)[0]
    time_text = ax.text(0.05, 0.95, f'time = 0ms', transform=ax.transAxes)

    updatable = list(bar_container) + [time_text]
    anim = animation.FuncAnimation(fig, update, repeat=False, blit=True)
    plt.show()
    dist = distfit()
    GD = get_data()
    #dist.fit_transform(GD[0])
    dist.fit_transform(np.array([vel**2 for vel in GD[0]]))
    dist.plot()
    #time_text = plt.text(2, 0.95, f'time = {int(GD[1])}ms', transform=ax.transAxes)
    plt.show()