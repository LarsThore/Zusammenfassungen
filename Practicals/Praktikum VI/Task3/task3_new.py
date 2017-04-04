#! encoding:utf8 #

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time, glob, os

############################### Set Constants #################################

H = 137100.               # Energy difference [J/mol]
R = 8.314               # Gas constant [J/(mol K)]

b = 250. * 10**(-9)    # nearest neighbour spacing [m]
T = (1050 + 273.15)     # Temperature [K]

kickouts = ["Au_in", "Au_out"]
diffusion = ["Au_diffusion"]
events = ["Au_in", "Au_out", "Au_diffusion"] # kick_out_1, kick_out_2 or diffusion

Au_out_rate = 10.       # Si-atom kicking out an Au-atom [s^(-1)]
Au_in_rate = 600.      # Au-atom kicking out an Si-atom [s^(-1)]
Au_diffusion_rate = 2700./2    # Au interstitial diffusion rate [s^(-1)]
                               # divided by 2 because of double amount
                               # of interstitial sites
################################ Functions ####################################

def remove_pictures():
    # Clean up old frames
    for name in glob.glob('tmp*.png'):
        os.remove(name)

def create_environment(n):
    # create environment in which there is a constant equilibrium
    # concentration at the left hand side boundary
    matrix = np.zeros((n, n))
    matrix[::2, 0] = 1
    return matrix

def choose_atom(matrix):
    '''Choose an Au atom randomly'''
    # find out all places where an Au atom sits
    Au_atoms = np.where(matrix[::2] == 1)
    print Au_atoms
    Au_x = Au_atoms[1]
    Au_y = Au_atoms[0]
    # choose one of these atoms randomly
    k = np.random.randint(len(Au_x))
    print 'k = ', k
    atom = (Au_x[k], Au_y[k])
    print atom
    return atom

def choose_atom_kickout(matrix):
    '''Choose an Au atom randomly which sits on a Si site'''
    # find out the particular places
    Au_atoms = np.where(matrix[1::2] == 1)
    print Au_atoms
    Au_x = Au_atoms[1]
    Au_y = Au_atoms[0]
    # choose one of these atoms randomly
    k = np.random.randint(len(Au_x))
    print 'k = ', k
    atom = (Au_x[k], Au_y[k])
    print atom
    return atom

def choose_event(events):
    '''Choose one of the possible three events'''

    global Au_out_rate, Au_in_rate, Au_diffusion_rate

    rate1 = Au_in_rate
    rate2 = Au_out_rate
    rate3 = Au_diffusion_rate

    rate_sum = sum([rate1, rate2, rate3])

    a = np.random.rand()

    try:
        if 0 <= a < rate1/rate_sum:
            return events[0]
        elif rate1/rate_sum <= a < rate2/rate_sum:
            return events[1]
        elif (rate1+rate2)/rate_sum <= a < 1:
            return events[2]
    except:
        raise ValueError("No happening could be chosen. Maybe the rates are\
                         calculated wrong.")

def choose_kickout(kickouts):
    rate1 = Au_out_rate
    rate2 = Au_in_rate
    rate_sum = sum([rate1, rate2])

    a = np.random.rand()

    try:
        if 0 <= a < rate1/rate_sum:
            return kickouts[0]
        elif rate1/rate_sum <= a < 1:
            return kickouts[1]
    except:
        raise ValueError("No happening could be chosen. Maybe the rates are\
                         calculated wrong.")

def make_diffusion(matrix):
    '''Carries out a diffusion event of an Au atom'''

    atom = choose_atom(matrix)

    i = atom[1]
    j = atom[0]

    global diffusionCounter, frustDiffusionCounter

    # move the atom back or forth

    # check if atom at the left border is chosen
    if j == 0:
        if matrix[i, j+1] == 0:
            matrix[i, j] = 1            # atom at left border is kept
            matrix[i, j+1] = 1
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
        else:
            print "Didn't move atom at this step"
            frustDiffusionCounter += 1
    # check if atom at the right border is chosen
    elif j == n-1:
        if matrix[i, j-1] == 0:
            matrix[i, j] = 0            # atom diffuses into the body
            matrix[i, j-1] = 1
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
        else:
            frustDiffusionCounter += 1
            print "Didn't move atom at this step"
    # check if it collides with another Au atom
    elif matrix[i, j+1] == 0 and matrix[i, j-1] == 0:
        matrix[i, j] = 0
        if np.random.randint(2) == 1:
            matrix[i, j+1] = 1              # atom moves to right hand side
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
        else:
            matrix[i, j-1] = 1              # atom moves to left hand side
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
    elif matrix[i, j+1] == 1 and matrix[i, j-1] == 0:
        if np.random.randint(2) == 1:
            matrix[i, j] = 0
            matrix[i, j-1] = 1              # atom moves to left hand side
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
        else:
            frustDiffusionCounter += 1
            print "Didn't move atom at this step"
    elif matrix[i, j-1] == 1 and matrix[i, j+1] == 0:
        if np.random.randint(2) == 1:
            matrix[i, j] = 0
            matrix[i, j+1] = 1              # atom moves to right hand side
            diffusionCounter += 1
            print "Atom no. {} moved".format(diffusionCounter)
        else:
            frustDiffusionCounter += 1
            print "Didn't move atom at this step"
    # check if atom sits on Si site
    elif j % 2 == 1:
        frustDiffusionCounter += 1
        print "Didn't move atom at this step"
    else:
        frustDiffusionCounter += 1
        print "Didn't move atom at this step"

    return matrix

def make_Au_kick_in(matrix):
    '''Makes the event of Au atom kick in mechanism happen'''

    atom = choose_atom(matrix)

    i = atom[1]
    j = atom[0]

    # check if atom at the very top is chosen
    if i == 0:
        matrix[i, j] = 0
        matrix[i-1, j] = 1              # atom moves to bottom side
    elif matrix[i+1, j] == 0 and matrix[i-1, j] == 0:
        matrix[i, j] = 0
        if np.random.randint(2) == 1:
            matrix[i+1, j] = 1              # atom moves to top side
        else:
            matrix[i-1, j] = 1              # atom moves to bottom side

    return matrix

def make_Au_kick_out(atom_site):
    matrix[i, j] = 0
    matrix[i-1, j] = 1

    return matrix

def make_kickout_move(atom_site, kickouts):
    kick_out = choose_kickout(kickouts)

    if kick_out == "Au_in":
        matrix = make_Au_kick_in(atom_site)
    elif kick_out == "Au_out":
        matrix = make_Au_kick_out(atom_site)

    return matrix

def calc_time_increment(jump_art):
    # average time per jump
    N_diffusors = 1         # TODO find a reasonable value for number of diffusors

    tau = b**2 / (2 * D)
    random_R = np.log(np.random.uniform())
    tau_ = N_diffusors * tau            # after Lesar
    t_step = - tau_ * random_R
    return t_step

def make_plots(matrix, t, n, counter):

    # years = t_k // 31536000                             # number of years
    # weeks = (t_k % 31536000.) // (86400.*7)             # number of weeks
    # hours = (t_k % (86400*7)) // 3600.                   # number of hours
    # seconds = t_k % 3600.                               # number of seconds

    # create a figure
    fig1 = plt.figure(figsize=(12,1), facecolor='white')
    # add plot to the figure
    fig1.subplots_adjust(left=0.06, right=0.98, top=0.7, bottom=0.3, hspace=0.3, wspace=0.3)
    ax1 = fig1.add_subplot(211)
    im = ax1.imshow(matrix, aspect = 'auto')
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    # plt.title("Temperature = {:.0f},\t $\Delta$t = {:.0f} years \t\t {:.0f} weeks \t\
    #  {:.0f} hours  \t\t{:.0f} sec".format(temp, years, weeks, hours, seconds))
    plt.title("Time = {:.1f} sec".format(t))

    plt.savefig('tmp_{:04d}.png'.format(counter))
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

def make_movie(n):
    from scitools.std import movie
    movie('tmp_*.png', encoder='convert', fps=4,
          output_file='tmp_diffusion_task3_system_size{}x{}.gif'.format(n, n))

############################### Precommands ###################################

remove_pictures()
time0 = time.time()

############################### main code #####################################

# set number of atoms in diffusion direction
n = 10                         # system size, length of one axis
t_end = 1000.                  # 1000 seconds annealing time

pictureCounter = 0                    # counter for pictures
diffusionCounter = 0                  # counter for number of diffused atoms
frustDiffusionCounter = 0 # counter for number of attempted but not-diffused atoms
diffusionCounterList = []

round_list = [20]#2, 5002]#, 10000]

for rounds in round_list:
    '''Loop über unterschiedlich viele Runden für Vergleich der Performance'''

    positions_rounds = np.zeros((rounds, n))
    for j in range(rounds):
        '''Loop über alle Runden für statistische Auswertung'''

        t = 0.                  # initial time

        if j % 100 == 0:
            '''printing current progress of simulation'''
            print j
            time.sleep(0.3)

        # make numpy array for statistical averaging
        # positions_temp = np.zeros((temp_steps, n))

        # create the diffusion environment
        matrix = create_environment(n)
        print matrix

        # set diffusion coefficient
        D = 1.49 * 10**(-7) * np.exp(-H/(R*T)) # in [m**2/s]

        # choose randomly an event
        event = choose_event(events)

        for i in range(20):
            matrix = make_diffusion(matrix)
            print matrix

        for i in range(10):
            matrix = make_Au_kick_in(matrix)
            print matrix

        time.sleep(10)

        if event == 'Au_in':
            '''kick in mechanism'''
            matrix = make_Au_kick_in(matrix)
        elif event == 'Au_out':
            '''kick out mechanism'''
            matrix = make_Au_kick_out(matrix)
        else:
            '''diffusion event'''
            matrix = make_diffusion(matrix)

        if matrix[i, j] == 1:

            # happening1 = choose_happening1(events)

            for k in range(10):  # small random value
                if i % 2 == 0:      # interstitial sites
                    matrix = make_kickout_move((i, j), kickouts)
                for l in range(1000): # big random value
                    if i % 2 == 1:      # lattice sites
                        matrix = make_diffusion((i, j))

                print matrix

            # calculate the increment of time per jump
            delta_t = 1 #calc_time_increment(jump_art)
            t += delta_t

            make_plots(matrix, delta_t, n, pictureCounter)

        # breche ab wenn maximale Zeitdauer erreicht ist
        if t > t_end:
            break

        # if j % 100 == 0:
            # print "T: {} \t D: {:.2g} \t time: {:.2g} sec".format(
            #  temp, D, t)

        pictureCounter += 1

        # average_position_temp = np.average(positions_temp, axis = 0) #* temp_steps
        # positions_rounds[j] = average_position_temp

    # average_position_rounds = np.average(positions_rounds, axis = 0) #* temp_steps

    # make_plot2(average_position_rounds, positions_rounds)

    diffusionCounterList.append(diffusionCounter)
    diffusionCounter = 0

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


print 'Atoms moved:', diffusionCounterList
print 'Atoms not moved:', frustDiffusionCounter
