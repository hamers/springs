import argparse
import numpy as np

from matplotlib import pyplot

import csv

def add_bool_arg(parser, name, default=False,help=None):
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true',help="Enable %s"%help)
    group.add_argument('--no-' + name, dest=name, action='store_false',help="Disable %s"%help)
    parser.set_defaults(**{name:default})


def parse_arguments():
    

    parser = argparse.ArgumentParser()

    parser.add_argument("--file",                 type=str,       dest="filename",        default="node_data.txt",       help=" ")

    args = parser.parse_args()

    return args

def plot_function(args):
    """
    Plot node output data (positions) in the form of animations
    """
    

    filename = args.filename

    with open(filename, 'r') as read_obj:
        
        csv_reader = csv.reader(read_obj, delimiter = '\t')
        for index,row in enumerate(csv_reader):
            if index==0:
                index_nodes = [row.index(x) for x in row if x == 'N_nodes'][0] + 1
                index_m = [row.index(x) for x in row if x == 'm'][0] + 1
                index_k = [row.index(x) for x in row if x == 'k'][0] + 1
                index_b = [row.index(x) for x in row if x == 'b'][0] + 1
                N_nodes = int(row[index_nodes])
                m = float(row[index_m])
                k = float(row[index_k])
                b = float(row[index_b])

                break

        ts = []
        snap_pos_x = [[] for x in range(N_nodes)]
        snap_pos_y = [[] for x in range(N_nodes)]
        snap_pos_z = [[] for x in range(N_nodes)]
 
        for index,row in enumerate(csv_reader):
            t = float(row[0])
            ts.append(t)
            for i in range(N_nodes):
                data = []
                for j in range(3):
                    data.append(float(row[1+3*i+j]))
                snap_pos_x[i].append(data[0])
                snap_pos_y[i].append(data[1])
                snap_pos_z[i].append(data[2])

    N = N_nodes
    Nt = len(ts)

    snap_pos_x_rel = [[0 for x in range(Nt)] for x in range(N)]
    snap_pos_y_rel = [[0 for x in range(Nt)] for x in range(N)]
    
    for i_t in range(Nt):
        
        x_av = np.mean(np.array([snap_pos_x[j][i_t] for j in range(N)]))
        y_av = np.mean(np.array([snap_pos_y[j][i_t] for j in range(N)]))
        
        for j in range(N):
            snap_pos_x_rel[j][i_t] = snap_pos_x[j][i_t] - x_av
            snap_pos_y_rel[j][i_t] = snap_pos_y[j][i_t] - y_av

    fig = pyplot.figure()
    plot = fig.add_subplot(1,1,1)
    plot.set_xlabel("x")
    plot.set_ylabel("y")
    plot.set_title("m = " + str(m) + "; k = " + str(k) + "; b = " + str(b))
    for i in range(N):
        plot.scatter(snap_pos_x[i][0],snap_pos_y[i][0],s=20)
        plot.plot(snap_pos_x[i],snap_pos_y[i])
        
    
    from matplotlib.animation import FuncAnimation
    pyplot.style.use('seaborn-pastel')


    fig = pyplot.figure()
    r = 10
    ax = pyplot.axes(xlim=(-r, r), ylim=(-r, r))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("m = " + str(m) + "; k = " + str(k) + "; b = " + str(b))

    line = ax.scatter([],[],s=20)

    def init():
        line.set_offsets(np.c_[0,0])
        return line,
    def animate(i):
        X = [snap_pos_x[j][i] for j in range(N)]
        Y = [snap_pos_y[j][i] for j in range(N)]

        line.set_offsets(np.c_[X,Y])

        return line,

    anim = FuncAnimation(fig, animate, init_func = init,
                                   frames=Nt, interval=10, blit=True)

    import matplotlib.animation as animation
    FFwriter = animation.PillowWriter(fps=30) 
    anim.save('animation.gif', writer = FFwriter)
    
    fig2 = pyplot.figure()
    ax2 = pyplot.axes(xlim=(-5, 5), ylim=(-5, 5))
    ax2.set_title("Relative to center; m = " + str(m) + "; k = " + str(k) + "; b = " + str(b))
    
    line2 = ax2.scatter([],[],s=20)

    def init2():
        line2.set_offsets(np.c_[0,0])
        return line,
    def animate2(i):

        X = [snap_pos_x_rel[j][i] for j in range(N)]
        Y = [snap_pos_y_rel[j][i] for j in range(N)]
        
        line2.set_offsets(np.c_[X,Y])

        return line2,

    anim2 = FuncAnimation(fig2, animate2, init_func = init2,
                                   frames=Nt, interval=10, blit=True)

    import matplotlib.animation as animation

    FFwriter = animation.PillowWriter(fps=30) 
    anim2.save('animation_rel.gif', writer = FFwriter)
    
    pyplot.show()
    
if __name__ == '__main__':
    args = parse_arguments()
    plot_function(args)

