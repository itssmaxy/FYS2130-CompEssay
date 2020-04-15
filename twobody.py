import numpy as np
import matplotlib.pyplot as plt
import sys 


sunmass = 1#6.38437672E-06
mass_plan1 = 1#6.38437672E-06

radius = 9.02419262E-8#0.00005267828
Gr = 39.47841760435743
c = 63239.7263 #Au/yr
r = 6.32397263E12 #AU


#Integration Loop

def integrator(planets,sun,planets_index):

    def gravity_on_planet(x_sun,y_sun,x_planet,y_planet,index):

        r = np.asarray([x_sun-x_planet, y_sun-y_planet])
        grav = (Gr*sunmass*masses[index])/(np.linalg.norm(r))**2
        ax = grav*r[0] /(np.linalg.norm(r))
        ay = grav*r[1]/(np.linalg.norm(r))

        return np.asarray([ax,ay],float)

    print("--------------")
    N = int(1E6)
    
    time = float(input("Choose simulation time in yrs: "))
    dt = time/N
    
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

        h[i+1] = h_stretch(r,omega,R,t)


        if R < 2*radius:
            break
        
        count +=1
        
        #i += 1



    return x_planets,x_sun, count, h, t


def h_stretch(r,omega,R,t):
    return 8*Gr*sunmass/(r*c**4)*omega**2*R**2*np.cos(2*omega*(t-r/c))

#Input index of planets wanted in simulation
planets_index = [0]

print("Please enter initial conditions for Object 1: ")
x_planet = float(input("X position for Object 1: "))
y_planet = float(input("Y position for Object 1: "))
vx_planet = float(input("VX velocity for Object 1: "))
vy_planet = float(input("VY velocity for Object 1: "))

planets = [[x_planet,y_planet,vx_planet,vy_planet]]

print("--------------")

print("Please enter initial conditions for Object 2: ")
x_sun = float(input("X position for Object 2: "))
y_sun = float(input("Y position for Object 2: "))
vx_sun = float(input("VX velocity for Object 2: "))
vy_sun = float(input("VY velocity for Object 2: "))

planets= np.asarray(planets)
Sun = [x_sun,y_sun,vx_sun,vy_sun]


#Final results and plotting

fig, ax = plt.subplots()

planet_orbit, sun_orbit, count, h, t = integrator(planets,Sun,planets_index)

planet_orbit = planet_orbit[0,:count+1,:]
print(np.shape(planet_orbit))
sun_orbit = sun_orbit[:count+1,:]
print(np.shape(sun_orbit))


r = np.linalg.norm(planet_orbit[-1,:] - sun_orbit[-1,:])/2
print("Distanse is = {} AU".format(r))
P = -32/5*Gr**4/c**5*((sunmass*mass_plan1)**2*(sunmass + mass_plan1))/r**5
print("Energi radiert = {}".format(P))

#for i in range(len(planets_index)):
ax.plot(planet_orbit[:,0], planet_orbit[:,1], "b--",label="Object 0")
ax.plot(sun_orbit[:,0], sun_orbit[:,1],"--", label="Object 1")

circle1 = plt.Circle((sun_orbit[-1,0], sun_orbit[-1,1]), radius, color="Green")
circle2 = plt.Circle((planet_orbit[-1,0], planet_orbit[-1,1]), radius, color="Black")

ax.add_artist(circle1)
ax.add_artist(circle2)

plt.xlabel("x-position (AU)")
plt.ylabel("y-position (AU)")
plt.legend(loc='lower right')
#plt.savefig("garvity-sun-pos.jpeg")
plt.show() #viser banen til objektene

fig, ax = plt.subplots()
ax.plot(np.linspace(0, t, len(h)), h)
plt.show()


"""
plt.plot(sun_orbit[:,0], sun_orbit[:,1], label="Sun")
plt.xlabel("x-position (AU)")
plt.ylabel("y-position (AU)")
plt.legend(loc='lower left')
plt.savefig("sun-all-grav.jpeg")
plt.show()
"""