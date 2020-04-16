import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import animation

dict = {}
dict['bh'] = [1,5.91824521E-7]#Solarmass and AU
dict['earth'] = [0.000003003,4.26354*1E-5]
print('Please choose from the following:')
print('bh --- earth ----  ')
object1 = input('Object 1:')
object2 = input('Object 2:')
sunmass = dict[object1][0]
mass_plan1 = dict[object2][0]

radius_1 = dict[object1][1]#0.00005267828
radius_2 = dict[object2][1]#0.00005267828

Gr = 39.47841760435743
c = 63239.7263 #Au/yr
r1= 6.32397263E12 #AU


#Integration Loop

def integrator(planets,sun,planets_index):

    def gravity_on_planet(x_sun,y_sun,x_planet,y_planet,index):

        r = np.asarray([x_sun-x_planet, y_sun-y_planet])
        grav = (Gr*sunmass*masses[index])/(np.linalg.norm(r))**2
        ax = grav*r[0] /(np.linalg.norm(r))
        ay = grav*r[1]/(np.linalg.norm(r))

        return np.asarray([ax,ay],float)

    print("--------------")


    time = float(input("Choose simulation time in yrs: "))
    dt = 0.0001
    N = int(time/dt)

    t = dt
    number_of_planets = int(len(planets))
    cm = np.array([0,0])
    masses = [mass_plan1]
    x0xsun,x0ysun,v0xsun,v0ysun = sun
    h = np.zeros(N)
    x_planets = np.asarray([np.zeros((N,2),float) for i in range(number_of_planets)])
    v_planets = [np.asarray([planets[i,2],planets[i,3]]) for i in range(number_of_planets)]

    for i in range(number_of_planets):
        x_planets[i,0,:] = [planets[i,0],planets[i,1]]
    a_i_planets = np.zeros((number_of_planets,2))
    a_iplus1_planets = np.zeros((number_of_planets,2))

    for i in range(number_of_planets):
        a_i_planets[i] = gravity_on_planet(x0xsun,x0ysun,planets[i,0],planets[i,1],i)
    x_sun = np.zeros((N,2),float)
    x_sun[0,:] = [x0xsun,x0ysun]
    v_sun = np.asarray([v0xsun,v0ysun],float)
    a_i_sun = -sum(a_i_planets)/sunmass
    count = 0
    for i in range(N-1):
        t+=dt
        cm_x = [x_planets[g,i,0]*masses[g] for g in range(number_of_planets)]
        cm_y = [x_planets[g,i,1]*masses[g] for g in range(number_of_planets)]
        cm = [(sum(cm_x)+sunmass*x_sun[i,0])/(sum(masses)+sunmass), \
        (sum(cm_y)+sunmass*x_sun[i,1])/(sum(masses)+sunmass)]

        for j in range(number_of_planets):
            x_planets[j,i+1,:] = x_planets[j,i,:]+(v_planets[j]*dt + 0.5*a_i_planets[j]/masses[j]*dt**2)-cm

        x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm

        for k in range(number_of_planets):
            a_iplus1_planets[k] = gravity_on_planet(x_sun[i+1,0],x_sun[i+1,1],x_planets[k,i+1,0],\
            x_planets[k,i+1,1],k)

        a_iplus1_sun = -sum(a_i_planets)/sunmass

        for l in range(number_of_planets):
            v_planets[l] += 0.5*(a_i_planets[l] + a_iplus1_planets[l])/masses[l]*dt

        a_i_planets = a_iplus1_planets
        v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
        a_i_sun = a_iplus1_sun


        R = np.linalg.norm(x_sun[i+1,:] - x_planets[0,i+1,:])
        omega = np.sqrt(Gr*sunmass/(4*R**3))

        h[i+1] = h_stretch(r1,omega,R,t)


        if R < (radius_1+radius_2):
            break

        count +=1
        #i += 1
    return x_planets,x_sun, count, h, t,N


def h_stretch(r,omega,R,t):
    return 8*Gr*sunmass/(r*c**4)*omega**2*R**2*np.cos(2*omega*(t-r/c))

#Input index of planets wanted in simulation
planets_index = [0]
r = 5


answer= input("If you want to collide press 'c'. If you want to orbit press 'o'")
hyp = np.asarray([-r,0])
tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
hyp[0]/np.linalg.norm(hyp)])
vel = np.sqrt(Gr*mass_plan1*sunmass/np.linalg.norm(hyp))
v_orbit = vel*tan_vec
v_colide = v_orbit*0.5 + vel*np.asarray([-0.5,0])
if answer == 'c':
    planet_1 = np.asarray([[r,0,v_colide[0],v_colide[1]]])
    planet_2 = np.asarray([0,0,0,0])
elif answer == 'o':
    planet_1 = np.asarray([[r,0,v_orbit[0],v_orbit[1]]])
    planet_2 = np.asarray([0,0,0,0])

#Final results and plotting

fig, (ax1,ax2) = plt.subplots(1,2)

planet_orbit, sun_orbit, count, h, t,N = integrator(planet_1,planet_2,planets_index)

planet_orbit = planet_orbit[0,:count+1,:]

sun_orbit = sun_orbit[:count+1,:]



r = np.linalg.norm(planet_orbit[-1,:] - sun_orbit[-1,:])/2
print("Distanse is = {} AU".format(r))
P = -32/5*Gr**4/c**5*((sunmass*mass_plan1)**2*(sunmass + mass_plan1))/r**5
print("Energi radiert = {}".format(P))

#for i in range(len(planets_index)):

ax1.plot(planet_orbit[:-1,0], planet_orbit[:-1,1], "b--",label=object1)
ax1.plot(sun_orbit[:-1,0], sun_orbit[:-1,1],"--", label=object2)
ax1.set(xlabel=("x-position (AU)"),ylabel= ("y-position (AU)"))
ax1.legend(loc='lower right')



#plt.savefig("garvity-sun-pos.jpeg")
ax2.plot(np.linspace(0,t,len(h)),h)
ax2.set(xlabel=('Time(y)'),ylabel=('Distortion(Au)'))
plt.show() #viser banen til objektene
"""
Animation kode
"""


intr = int(1e3)
# First set up the figure, the axis, and the plot element we want to animate
fig, (ax1,ax2) = plt.subplots(1,2)
ax1.set(xlim=(0,10),ylim=(np.min(h)*1.5,np.max(h)*1.5), ylabel=('Distortion'))
ax2.set(xlim=(-2,2),ylim=(-2,2))
line1, = ax1.plot(np.linspace(0,10,intr), h[0:int(intr)], lw=2)
line2, = ax2.plot(planet_orbit[0:int(intr),0],planet_orbit[0:int(intr),1])
line3, = ax2.plot(sun_orbit[0:int(intr),0],sun_orbit[0:int(intr),1])

# initialization function: plot the background of each frame
def init1():
    line1.set_data([], [])
    return line1,

# animation function.  This is called sequentially
def animate1(i):
    #x = np.linspace(int(i*t/(1e2)),int((i+1)*t/(1e2)), 1e2)
    x = np.linspace(0,10,intr)
    y = h[int((i*intr)):int((i+1)*intr)]
    line1.set_data(x, y)
    return line1,
def init2():
    line2.set_data([], [])
    return line2,

# animation function.  This is called sequentially
def animate2(i):
    #x = np.linspace(int(i*t/(1e2)),int((i+1)*t/(1e2)), 1e2)
    x = planet_orbit[int(i*intr):int((i+1)*intr),0]
    y = planet_orbit[int(i*intr):int((i+1)*intr),1]
    line2.set_data(x, y)
    return line2,
def init3():
    line3.set_data([], [])
    return line3,

# animation function.  This is called sequentially
def animate3(i):
    #x = np.linspace(int(i*t/(1e2)),int((i+1)*t/(1e2)), 1e2)
    x = sun_orbit[int(i*intr):int((i+1)*intr),0]
    y = sun_orbit[int(i*intr):int((i+1)*intr),1]
    line3.set_data(x, y)
    return line3,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim1 = animation.FuncAnimation(fig, animate1, init_func=init1,
                               frames=int(N/intr), interval=100, blit=True)
anim2 = animation.FuncAnimation(fig, animate2, init_func=init2,
                               frames=int(N/intr), interval=100, blit=True)
anim3 = animation.FuncAnimation(fig, animate3, init_func=init3,
                               frames=int(N/intr), interval=100, blit=True)
plt.show()
