#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.animation as animation


class particule:
    def __init__(self,x,y,vx,vy):
        self.r = np.array([x,y])
        self.v = np.array([vx,vy])

    def dist(self,part2):
        return np.linalg.norm(self.r-part2.r)

    def actualiser(self,dt):
        self.r[0] = self.r[0] + self.v[0]*dt
        self.r[1] = self.r[1] + self.v[1]*dt


class Box:
    def __init__(self,x_lim,y_lim,masse,rayon,T,N,dt,pression=0):
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.N = N
        self.dt = dt
        self.masse = masse
        self.rayon = rayon
        self.T = T
        self.pression = pression

        self.liste_particule= []
        liste_particule = np.random.random((N,6))
        liste_particule[:,0] = (liste_particule[:,0]*(x_lim-2*self.rayon)) + self.rayon
        liste_particule[:,1] = (liste_particule[:,1]*(y_lim-2*self.rayon)) + self.rayon
        liste_particule[:,2] *= 2*np.pi
        a = np.sqrt((1.38*10**-23)*T/masse)
        liste_particule[:,3] = stats.maxwell.rvs(loc=0, scale=a, size=N)

        for param in liste_particule:
            self.liste_particule.append(particule(param[0],param[1],param[3]*np.cos(param[2]),param[3]*np.sin(param[3])))


    def step(self):

        #self.pression = 0

        #Vérifier pour collision avec le mur
        for part in self.liste_particule:
            if (part.r[0] + self.rayon) > self.x_lim :
                part.v[0] = -abs(part.v[0])
                #self.pression += 2*self.masse*abs(part.v[0])/(self.dt*self.y_lim*1)
            elif (part.r[0] - self.rayon) < 0:
                part.v[0] = abs(part.v[0])
                #self.pression += 2*self.masse*abs(part.v[0])/(self.dt*self.y_lim*1)

            if (part.r[1] + self.rayon) > self.y_lim:
                part.v[1] = -abs(part.v[1])
                #self.pression += 2*self.masse*abs(part.v[1])/(self.dt*self.x_lim*1)
            elif (part.r[1] - self.rayon) < 0:
                part.v[1] = abs(part.v[1])
                #self.pression += 2*self.masse*abs(part.v[1])/(self.dt*self.x_lim*1)

        #Vérifier pour collision entre les particules
        for i,part1 in enumerate(self.liste_particule):
            for j,part2 in enumerate(self.liste_particule):
                if i >= j:
                    pass
                elif part1.dist(part2) < 2*self.rayon:
                    part1.v,part2.v = (part1.v - np.dot(part1.v-part2.v, part1.r-part2.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part1.r-part2.r),
                                       part2.v - np.dot(part2.v-part1.v, part2.r-part1.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part2.r-part1.r))

        #Actualiser le système
        for i,part in enumerate(self.liste_particule):
            part.actualiser(self.dt)

        v2_average = np.array([np.linalg.norm(i.v)**2 for i in self.liste_particule]).mean()
        self.pression = self.N * self.masse * v2_average / (self.x_lim * self.y_lim)


#Création du système
dt = 1
sys = Box(10,10,10**-15,0.1,300,10,dt)


#Définition d'une fonction pour pour actualiser la position de plusieurs planètes à la fois
def actualiser_systeme(systeme):
    while True:
        systeme.step()
        yield systeme

t = 0
#Programme principal
def main():
    x_lim = 15
    y_lim = 15
    masse = 2*1.008/(6.022*10**23)
    rayon = 0.12
    dt = 0.005
    N = 50
    T = 300

    sys = Box(x_lim,y_lim,masse,rayon,T,N,dt)

    #Initialisation de la figure
    #fig, (ax1,ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
    fig = plt.figure(figsize=(10,5))
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    ax2 = plt.subplot2grid((2, 3), (0, 2))
    ax3 = plt.subplot2grid((2, 3), (1, 2))

    #Paramètres esthétiques
    ax1.set_xlim([0,x_lim])
    ax1.set_ylim([0,y_lim])
    ax1.set_xlabel('x (m)', fontsize = 16)
    ax1.set_ylabel('y (m)', fontsize = 16)
    ax1.set_aspect(1)
    plt.axis('equal')
    points_espace = [ax1.plot([], [], 'ro',  markersize=(rayon*450/x_lim)) for i in sys.liste_particule]
    T_texte = ax1.text(0.02, 0.95, '0:.2f'.format(''), color='k', transform=ax1.transAxes)

    ax2.set_xlabel('Vitesse (m/s)',fontsize = 14)
    ax2.set_ylabel('Nombre de particules', fontsize = 12)
    ax2.hist([], bins=round(N/10), facecolor='blue', alpha=1, edgecolor='black', linewidth=1.2)

    ax3.set_xlabel('Temps (s)', fontsize = 14)
    ax3.set_ylabel('Pression (Pa)', fontsize = 14)
    p_line = ax3.plot([], [])
    #ax3.plot([], [])

    #Définition de la fonction d'animation du système
    p = []
    time = []
    def run(data):
        #actualisation du graphique
        for i,points in enumerate(points_espace):
            points[0].set_data(data.liste_particule[i].r[0], data.liste_particule[i].r[1])

        #print([i.v[0] for i in data.liste_particule])
        #print([i.v[1] for i in data.liste_particule])
        #print([i.r[0] for i in data.liste_particule])
        #print([i.r[1] for i in data.liste_particule])

        global t
        t += dt
        time.append(t)
        T_texte.set_text('t = {0:.2f} s'.format(t))

        ax2.cla()
        ax2.hist([np.linalg.norm(i.v) for i in data.liste_particule], bins=round(N/10), facecolor='blue', alpha=1, edgecolor='black', linewidth=1.2)
        ax2.set_xlabel('Vitesse (m/s)', fontsize = 14)
        ax2.set_ylabel('Nombre de particules', fontsize = 12)

        p.append(data.pression*10**20)
        print(p)
        print(time)
        #ax3.plot(time,p)
        p_line[0].set_data(time, p)
        #ax3.set_xlim([0,max(time)])
        #ax3.set_ylim([0,max(p)])

        return points

    #Animation
    anim = animation.FuncAnimation(fig, run, actualiser_systeme(sys), interval=10, blit=False, repeat=True)

    #Traçage de l'animation
    #plt.figure(figsize=(10,10))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
