#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

class particule(Object):
    def __init__(self,x,y,z,vx,vy,vz):
        self.r = np.array[x,y,z]
        self.v = np.array[vx,vy,vz]

    def dist(self,part2):
        return np.linalg.norm(self.r-part2.r)

    def actualiser(self,dt):
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt
        self.z = self.z + self.vz*dt


class Box(Object):
    def __init__(self,x_lim,y_lim,z_lim,masse,rayon,T,N,dt):
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.z_lim = z_lim
        self.N = N
        self.dt = dt
        self.masse = masse
        self.rayon = rayon
        self.T = T

        self.liste_particule= []
        liste_particule = np.random((N,6))
        liste_particule[:,0] *= x_lim
        liste_particule[:,1] *= y_lim
        liste_particule[:,2] *= z_lim
        liste_particule[:,3] *= 2*np.pi
        liste_particule[:,4] *= np.pi
        a = np.sqrt((1.38*10**-23)*T/masse)
        liste_particule[:,5] *= data = stats.maxwell.rvs(loc=0, scale=a, size=N)

        for param in liste_particule:
            self.liste_particule.append(particule(param[0],param[1],param[2],param[5]*np.sin(param[4])*np.cos(param[3]),param[5]*np.sin(param[4])*np.sin(param[3]),param[5]*np.cos(param[4])))


    def step(self):
        #Vérifier pour collision avec le mur
        for part in self.liste_particule:
            if (part.x + self.rayon) > self.x_lim or (part.x - self.rayon) < 0:
                part.v[0] *= -1
            if (part.y + self.rayon) > self.y_lim or (part.y - self.rayon) < 0:
                part.v[1] *= -1
            if (part.z + self.rayon) > self.z_lim or (part.z - self.rayon) < 0:
                part.v[2] *= -1

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
sys = Box(10,10,10,1,0.1,300,10,0.0001)

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

lines = [ax.plot(part.r[0], part.r[1], part.r[2])[0] for part in sys.liste_particule]

def update_lines(num, dataLines, lines) :
    for line, data in zip(lines, dataLines) :
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2,:num])
    return lines


# Setting the axes properties
ax.set_xlim3d([0.0, sys.x_lim])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, sys.y_lim])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, sys.z_lim])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),
                              interval=50, blit=False)

plt.show()


