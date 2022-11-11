import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib import style

HIST_BINS = 50

def get_data():
    data = open('data.txt', 'r').read()
    lines = np.array(data.split('\n'))
    Vels = []
    for line in lines:
        if len(line) > 1:
            Vels.append(float(line.split(' ')[2]))
    #print(Vels)
    return Vels

def ret_update(bar_container):
    def update(i):
        data = get_data()
        #print(data)
        n, _ = np.histogram(data, HIST_BINS)
        for count, rect in zip(n, bar_container.patches):
            rect.set_height(count)
        return bar_container.patches
    return update



def graph(vel_lim):
    style.use("fivethirtyeight")
    data = get_data()
    #data = np.random.randn(1000)
    n, _ = np.histogram(data, HIST_BINS)
    fig, ax = plt.subplots()


    _, _, bar_container = ax.hist(data, HIST_BINS, lw=1,
                                  ec="yellow", fc="green", alpha=0.5)

    plt.xlabel('Velocity')
    plt.ylabel('Number')
    plt.title('Histogram of velocities')
    plt.grid(True)
    ax.set_ylim(top=100)
    ax.set_xlim(right=vel_lim)
    ani = animation.FuncAnimation(fig, ret_update(bar_container), repeat=False, blit=True)
    plt.show()