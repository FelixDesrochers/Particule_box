#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.animation as animation


class particule:
    def __init__(self,x,y,z,vx,vy,vz):
        self.r = np.array([x,y,z])
        self.v = np.array([vx,vy,vz])

    def dist(self,part2):
        return np.linalg.norm(self.r-part2.r)

    def actualiser(self,dt):
        self.r[0] = self.r[0] + self.v[0]*dt
        self.r[1] = self.r[1] + self.v[1]*dt
        self.r[2] = self.r[2] + self.v[2]*dt


class Box:
    def __init__(self,x_lim,y_lim,z_lim,masse,rayon,T,N,dt,pression=0):
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.z_lim = z_lim
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
        liste_particule[:,2] = (liste_particule[:,2]*(y_lim-2*self.rayon)) + self.rayon
        liste_particule[:,3] *= 2*np.pi
        liste_particule[:,4] *= np.pi
        a = np.sqrt((1.38*10**-23)*T/masse)
        liste_particule[:,5] = stats.maxwell.rvs(loc=0, scale=a, size=N)

        for param in liste_particule:
            self.liste_particule.append(particule(param[0],param[1],param[2],param[5]*np.sin(param[4])*np.cos(param[3]),param[5]*np.sin(param[4])*np.sin(param[3]),param[5]*np.cos(param[4])))


    def step(self):

        self.pression = 0

        #Vérifier pour collision avec le mur
        for part in self.liste_particule:
            if (part.r[0] + self.rayon) > self.x_lim :
                part.v[0] = -abs(part.v[0])
                self.pression += 2*self.masse*abs(part.v[0])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))
            elif (part.r[0] - self.rayon) < 0:
                part.v[0] = abs(part.v[0])
                self.pression += 2*self.masse*abs(part.v[0])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))

            if (part.r[1] + self.rayon) > self.y_lim:
                part.v[1] = -abs(part.v[1])
                self.pression += 2*self.masse*abs(part.v[1])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))
            elif (part.r[1] - self.rayon) < 0:
                part.v[1] = abs(part.v[1])
                self.pression += 2*self.masse*abs(part.v[1])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))
                
            if (part.r[2] + self.rayon) > self.z_lim:
                part.v[2] = -abs(part.v[2])
                self.pression += 2*self.masse*abs(part.v[2])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))
            elif (part.r[2] - self.rayon) < 0:
                part.v[2] = abs(part.v[2])
                self.pression += 2*self.masse*abs(part.v[2])/(self.dt*(2*self.x_lim*self.y_lim + 2*self.x_lim*self.z_lim + 2*self.y_lim*self.z_lim))

        #Vérifier pour collision entre les particules
        for i,part1 in enumerate(self.liste_particule):
            for j,part2 in enumerate(self.liste_particule):
                if i >= j:
                    pass
                elif part1.dist(part2) < 2*self.rayon:
                    part1.v,part2.v = (part1.v - np.dot(part1.v-part2.v, part1.r-part2.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part1.r-part2.r),
                                       part2.v - np.dot(part2.v-part1.v, part2.r-part1.r)/(np.linalg.norm(part1.r-part2.r)**2)*(part2.r-part1.r))
                    part1.r,part2.r = ((part1.r + part2.r)/2) + self.rayon * (part1.r - part2.r)/np.linalg.norm(part1.r-part2.r), ((part1.r + part2.r)/2) - self.rayon * (part1.r - part2.r)/np.linalg.norm(part1.r-part2.r)
        #Actualiser le système
        for i,part in enumerate(self.liste_particule):
            part.actualiser(self.dt)


#Définition d'une fonction pour pour actualiser la position de plusieurs planètes à la fois
def actualiser_systeme(systeme):
    while True:
        systeme.step()
        yield systeme

t = 0
#Programme principal
def main():
    x_lim = 15*10**-9
    y_lim = 15*10**-9
    z_lim = 15*10**-9
    masse = 2*1.008/(6.022*10**23)
    rayon = 0.22*10**-9
    dt = 0.000000000002
    N = 50
    T = 300

    sys = Box(x_lim,y_lim,z_lim,masse,rayon,T,N,dt)

    #Initialisation de la figure
    #fig, (ax1,ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
    fig = plt.figure(figsize=(10,5))
    plt.suptitle('Gaz Parfait en 3D', fontsize=20)
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    ax2 = plt.subplot2grid((2, 3), (0, 2))
    ax3 = plt.subplot2grid((2, 3), (1, 2))
    fig.text(.99, .005, 'Par Félix Desrochers, Gregory Giard et Olivier Gamache', ha='right', fontsize=6)

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
    ax2.hist([], bins=round(N/3), facecolor='blue', alpha=1, edgecolor='black', linewidth=1.2)

    ax3.set_xlabel('Temps (ns)', fontsize = 14)
    ax3.set_ylabel('Pression (kPa)', fontsize = 14)
    ax3.set_ylim(0,300)
    ax3.set_xlim(left=0)
    P_texte = ax3.text(0.02, 0.90, '0:.2f'.format(''), color='k', transform=ax3.transAxes)
    p_line = ax3.plot([], [], 'g-')

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
        time.append(t*10**9)
        T_texte.set_text('t = {0:.3e} s'.format(t))

        ax2.cla()
        ax2.hist([np.linalg.norm(i.v) for i in data.liste_particule], bins=round(N/3), facecolor='blue', alpha=1, edgecolor='black', linewidth=1.2)
        ax2.set_xlabel('Vitesse (m/s)', fontsize = 14)
        ax2.set_ylabel('Nombre de particules', fontsize = 12)
        
        p.append(data.pression/1000)
        ax3.set_ylim([0,300])
        ax3.set_ylim(bottom = 0)
        ax3.set_xlim(0,t*10**9)
        p_line[0].set_data(time, p)
        P_texte.set_text('P = {0:.2f} kPa'.format(np.mean(p)))
        
        return points

    #Animation
    anim = animation.FuncAnimation(fig, run, actualiser_systeme(sys), interval=30, blit=False, repeat=True)

    #Traçage de l'animation
    #plt.figure(figsize=(10,10))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    #print(np.mean([p[i] for i in np.nonzero(p)[0]]))
    print('La pression moyenne est de '+str(np.mean(p))+ ' kPa')

if __name__ == "__main__":
    main()
    