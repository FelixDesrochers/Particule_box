#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.animation as animation
from matplotlib.patches import Rectangle

class particule:
    def __init__(self,x,y,vx,vy):
        self.r = np.array([x,y])
        self.v = np.array([vx,vy])

    def dist(self,part2):
        return np.linalg.norm(self.r-part2.r)

    def actualiser(self,dt):
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt


class Box:
    def __init__(self,x_lim,y_lim,masse,rayon,T,N,dt):
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.N = N
        self.dt = dt
        self.masse = masse
        self.rayon = rayon
        self.T = T

        self.liste_particule= []
        liste_particule = np.random.random((N,6))
        liste_particule[:,0] *= x_lim
        liste_particule[:,1] *= y_lim
        liste_particule[:,2] *= 2*np.pi
        a = np.sqrt((1.38*10**-23)*T/masse)
        liste_particule[:,3] = stats.maxwell.rvs(loc=0, scale=a, size=N)

        for param in liste_particule:
            self.liste_particule.append(particule(param[0],param[1],param[3]*np.cos(param[2]),param[3]*np.sin(param[3])))


    def step(self):
        #Vérifier pour collision avec le mur
        for part in self.liste_particule:
            if (part.r[0] + self.rayon) > self.x_lim or (part.r[0] - self.rayon) < 0:
                part.v[0] *= -1
            if (part.r[1] + self.rayon) > self.y_lim or (part.r[1] - self.rayon) < 0:
                part.v[1] *= -1

        #Vérifier pour collision entre les particules
        for i,part1 in enumerate(self.liste_particule):
            for j,part2 in enumerate(self.liste_particule):
                if i >= j:
                    pass

                elif part1.dist(part2) < 2*self.rayon:
                    part1.v,part2.v = (part1.v - np.dot(part1.v-part2.v, part1.r-part2.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part1.r-part2.r),
                                       part2.v - np.dot(part2.v-part1.v, part2.r-part1.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part2.r-part1.r))

        #Actualiser le système
        for part in self.liste_particule:
            part.actualiser(self.dt)


#Création du système
dt = 0.0001
sys = Box(10,10,1,0.1,300,10,dt)

#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                xlim=(-3.2, 3.2), ylim=(-2.4, 2.4))

# particles holds the locations of the particles
particles, = ax.plot([], [], 'ro', ms=6)

# rect is the box edge
ax.add_patch(Rectangle((10 - 0.1, 10 - 0.1), 0.2, 0.2, alpha=1, facecolor='none'))


planetes_espace = [ plt.plot(part.r[0], part.r[1], 'o', color='r', markersize=2) for part in sys.liste_particule ]

def animate():
    # update pieces of the animation
    sys.step()
    for part,plot_part in zip(sys.liste_particule,planetes_espace):
        plot_part[0].set_data(part.r[0], part.r[1])
    return plot_part

anim = animation.FuncAnimation(fig, animate, interval=5, blit=False, repeat=False, save_count=300,)

plt.show()
