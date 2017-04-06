#! encoding:utf8 #

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time, glob, os

############################# Set Constants ###################################

H = 137100.               # Energy difference [J/mol]
R = 8.314               # Gas constant [J/(mol K)]

b = 250. * 10**(-9)    # nearest neighbour spacing [m]

############################## Functions ######################################

def remove_pictures():
    # Clean up old frames
    for name in glob.glob('tmp*.png'):
        os.remove(name)

def create_environment(n):
    # create environment in which the cells move
    # on the left hand side (first half) 10 % of the cells shall be
    # of different species (Aluminum cells)
    matrix = np.zeros((n))
    matrix[:n/2] = 40
    matrix = matrix.reshape((1, n))
    # matrix = 40*matrix
    print matrix
    return matrix

def choose_cell(matrix):
    '''Choose one of the cells where the diffusion event will happen. Also
    consider the amount of aluminum cells inside each cell for calculating
    the probabilities'''

    probability_matrix = matrix / np.sum(matrix)
    acumulated_probability = np.cumsum(probability_matrix)

    r = np.random.uniform(0.,1.)

    cell = 0
    while acumulated_probability[cell] < r:
        cell += 1

    return cell

def move_atom(matrix, cell):
    '''Move one atom in either the left hand side or right hand side direction;
    take also into account the number of available sites in the neighbor cells'''

    global move_attempts
    move_attempts += 1
    print "Atom no. {} moved".format(move_attempts)

    left_cell = matrix[0][cell-1]
    right_cell = matrix[0][cell+1]
    middle_cell = matrix[0][cell]

    # number of all lattice sites per cell is 10 times more the number of Al-
    # occupied ones
    all_lattice_sites = 10 * z

    # number of available lattice sites is the sum of free sites on the left
    # and right
    left_free_lattice_sites = all_lattice_sites - left_cell
    right_free_lattice_sites = all_lattice_sites - right_cell
    available_lattice_sites = left_free_lattice_sites + right_free_lattice_sites

    # calculate probabilities
    prob_left = left_free_lattice_sites / available_lattice_sites

    # move the cell back or forth
    # check if cell at the left border is chosen
    if cell == 0:
        if right_cell == 0:
            middle_cell += -1
            right_cell += 1
        else:
            pass
    # check if cell at the right border is chosen
    elif cell == n-1:
        if left_cell == 0:
            middle_cell += -1
            left_cell += 1
        else:
            pass
    # check if it collides with another cell which is already full of aluminum
    elif right_cell <= 140 and left_cell <= 140:
        middle_cell += -1
        if np.random.uniform(0., 1.) <= prob_left:
            left_cell += 1
        else:
            right_cell += 1
    elif right_cell >= 143 and left_cell <= 143:
        if np.random.uniform(0., 1.) <= prob_left:
            middle_cell += -1
            left_cell += 1
        else:
            pass
    elif left_cell >= 143 and right_cell <= 143:
        if np.random.uniform(0., 1.) > prob_left:
            middle_cell += -1
            right_cell += 1
        else:
            pass
    else:
        pass

    return matrix

def calc_time_increment(n, z):
    '''calculate the average time per jump using the total amount of
    aluminum atoms n/2 * z'''
    tau = b**2 / (2 * D)
    random_R = np.log(np.random.uniform())
    B = n/2. * z * tau                        # after Lesar
    t_step = - B * random_R
    return t_step

def make_plots(matrix, vis_matrix, temp, t, n, intervall):

    # years = t_k // 31536000                             # number of years
    # weeks = (t_k % 31536000.) // (86400.*7)             # number of weeks
    # hours = (t_k % (86400*7)) // 3600.                   # number of hours
    # seconds = t_k % 3600.                               # number of seconds

    # create a figure
    fig1 = plt.figure(figsize=(12,1), facecolor='white')
    # add plot to the figure
    fig1.subplots_adjust(left=0.06, right=0.98, top=0.7, bottom=0.3, hspace=0.3, wspace=0.3)
    ax1 = fig1.add_subplot(211)
    im = ax1.imshow(vis_matrix, aspect = 'auto')
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    # plt.title("Temperature = {:.0f},\t $\Delta$t = {:.0f} years \t\t {:.0f} weeks \t\
    #  {:.0f} hours  \t\t{:.0f} sec".format(temp, years, weeks, hours, seconds))
    plt.title("Temperature = {:.0f}, \t\t{:.1f} sec".format(
     temp, t))

    x = np.arange(n)
    y = np.zeros(n)
    p = n/10
    q = n/10
    ratios = np.zeros(n)
    for k in range(p, len(matrix[0][p:-p]) + p):
        intervall1 = 0.
        for j in matrix[0][k-q:k+q]:
            if abs(j - 1) < 1e-6:
                intervall1+= 1
        ratios[k] = intervall1 / p

    ax2 = fig1.add_subplot(212)
    im2 = ax2.plot(ratios)
    ax2.set_xlim([0, n])
    ax2.set_ylim([0, 0.35])
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])
    plt.savefig('tmp_{:04d}.png'.format(intervall))
    plt.clf()       # to delete the figure from the cache
    plt.close()     # to close the window

def make_plot2(average_over_position, positions):
    # create a figure
    plt.figure(facecolor='white')
    # add plot to the figure
    plt.plot(np.arange(n)*0.25 - (n*0.25 / 2), average_over_position)
    plt.title("Average Number of Atoms at particular Site")
    plt.xlabel('Distance from Interface [$\mu m$]')
    plt.ylabel('Aluminum Concentration')
    plt.xlim(- (n*0.25 / 2), (n*0.25 / 2))
    print 'Anzahl der Elemente in positions[:, 0]:', len(positions[:, 0])
    plt.savefig('Average_Concentration_{}.png'.format(
     len(positions[:, 0])))

def make_movie():
    from scitools.std import movie
    movie('tmp_*.png', encoder='convert', fps=4,
          output_file='tmp_diffusion_task1.gif')

############################### Precommands ###################################

remove_pictures()
time0 = time.time()

############################### main code #####################################

# set number of cells in diffusion direction
n = 80

# set number of Al cells per cell
z = 14 # 143 atom sites per cell --> see calculations

t_end = 10.**6
temperature_steps = 31             # number of steps in temperature intervall

temperature = np.linspace(300, 600, temperature_steps)

intervall = 0
move_attempts = 0
move_attempts_list = []

msg = ""

round_list = [102]#, 5002]#, 10000]

# calculate the time intervall
t_k = t_end / temperature_steps

for rounds in round_list:
    '''loop over every possible number of rounds'''
    temperature_averages = np.zeros((rounds, n))
    for j in range(rounds):
        t = 0.                  # initial time

        if j % 100 == 0:
            print j
            time.sleep(0.3)

        position_average = np.zeros((temperature_steps, n))
        matrix = create_environment(n)

        for k, temp in enumerate(temperature):
            '''loop over all temperatures'''
            t_temp = 0.
            # set diffusion coefficient
            D = 1.49 * 10**(-7) * np.exp(-H/(R*temp)) # in [m**2/s]

            # positions is the array of the model copied many
            # times on below another
            position_average[k] = matrix

            temperature_counter = 0

            while t_temp < t_k:
                '''loop inside of every temp intervall'''
                # choose cell randomly
                cell = choose_cell(matrix)
                matrix = move_atom(matrix, cell)

                # calculate the increment of time per jump
                delta_t = calc_time_increment(n, z)
                print delta_t
                time.sleep(0.1)
                t_temp += delta_t

                temperature_counter += 1

            # make big time steps
            else:
                if temperature_counter == 1:
                    t += t_k
                else:
                    t += t_temp

            if t > t_end:
                break

            msg += "{} attempts to move an atom in temperature intervall {}\
              \n".format(temperature_counter, intervall)

            intervall += 1

        temp_average = np.average(position_average, axis = 0) #* temperature_steps
        temperature_averages[j] = temp_average

        # make_plot2(temperature_averages, temp_average)

    round_average = np.average(temperature_averages, axis = 0) #* temperature_steps
    make_plot2(round_average, temperature_averages)

    move_attempts_list.append(move_attempts)
    move_attempts = 0

    ############################### End Commanfs #################################

    # print out time, make movie, remove picture files, print time to convert
    # pictures into gif file


    time1 = time.time()

    print 'Simulation time: {} minutes and {:2d} seconds'.format(
     int((time1 - time0) // 60), int((time1 - time0) % 60))

    # if rounds == round_list[-1]:
    #     make_movie()

    remove_pictures()

    print 'Conversion time to GIF: {} minutes and {:2d} seconds'.format(
     int((time.time() - time1) // 60), int((time.time() - time1) % 60))


print 'Atoms moved:', move_attempts_list

print msg
