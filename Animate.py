import numpy as np
import matplotlib.pyplot as plt
#from progress.bar import ChargingBar
from matplotlib import animation
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.image as mpimg


dict = {}
"""
Define objects first by weight and then radius, First Solarmass and then radius
"""
dict['bh'] = [1,5.91824521E-7,'k']#BlackHole
dict['bbh'] = [30,5.924*1e-7,'k']
dict['ns'] = [1.4,1.00268807*1E-7,'grey'] #Neutral Star
dict['sun'] = [1,0.0046524726,'y'] # Sun
dict['rg'] = [0.8,0.46524726,'r'] #Red Gianta
print('Please choose from the following:')
print('bh (Black Hole)---bbh (Big Black Hole) --- sun --- ns (Neutron Star)---rg (Red Giant)')
object1 = input('Object 1:')
object2 = input('Object 2:')
sunmass = dict[object2][0]
mass_plan1 = dict[object1][0]

radius_1 = dict[object1][1]#0.00005267828
radius_2 = dict[object2][1]#0.00005267828

Gr = 39.47841760435743
c = 63239.7263 #Au/yr
r1= 6.32397263E12 #AU Distance from the objects where we detect the distortion/wave.




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
    dt = 0.00001 #år
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
    answer2='o'
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

        if np.linalg.norm(v_planets[0])>100:
            print('THE OBJECTS HAVE COLLIDED!')
            answer2='c'
            break
        if R<(radius_1+radius_2):
            print('THE OBJECTS HAVE COLLIDED!')
            answer2='c'
            break
        count +=1

    return x_planets,x_sun, count, h, t,N ,dt,answer2


def h_stretch(r,omega,R,t):
    return 8*Gr*sunmass/(r*c**4)*omega**2*R**2*np.cos(2*omega*(t-r/c))

#Input index of planets wanted in simulation
planets_index = [0]
dist = 1
if 5*(radius_1 +radius_2)>1:
    dist= 1e1*(radius_1 +radius_2)


"""
Collision or orbit computing
"""
answer= input("If you want to collide press 'c'. If you want to orbit press 'o'")
hyp = np.asarray([-dist,0])
tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
hyp[0]/np.linalg.norm(hyp)])
vel = np.sqrt(Gr*mass_plan1*sunmass/np.linalg.norm(hyp))
v_orbit = vel*tan_vec
v_colide = v_orbit*0.3 + vel*np.asarray([-0.7,0])

if answer == 'c':
    planet_1 = np.asarray([[dist,0,v_colide[0],v_colide[1]]])
    planet_2 = np.asarray([0,0,0,0])
    action = "collision"
elif answer == 'o':
    planet_1 = np.asarray([[dist,0,v_orbit[0],v_orbit[1]]])
    planet_2 = np.asarray([0,0,0,0])
    action = "orbit"


planet_orbit, sun_orbit, count, h, t,N, dt,answer2 = integrator(planet_1,planet_2,planets_index)
Sampling = int(2*len(h[:count+1]))
#Fourier Transform
Hdef = int(31556926)
def MonsieurFourier(h,dt,t):

    dt = dt
    N = h.size
    Yr = N*dt
    Hdef = int(31556926)
    frate = N/(Yr*Hdef)
    Gfrate = frate
    print("Yr: ", Yr)
    print(frate, "Samples pr. second")
    FourierTransform = np.fft.fft(h)*frate
    Sample_Frequency = np.linspace(0,Sampling/2, Sampling//2)*frate

    plt.plot(Sample_Frequency, 2/Sampling*np.abs(FourierTransform[:Sampling//2]),label="Fourier Analysis")
    plt.xlabel('Frekvens [Hz]')
    plt.ylabel('Utslag')
    plt.legend()
    plt.show()

def ask():
    print("Do you wish to do a Fourier analysis? (y,n)")
    ans = input("Answer: ")
    if ans == "y":
        MonsieurFourier(h,dt,t)
    elif ans == "n":
        pass
    else:
        ask()
ask()

# Wavelet Analasys
def Wavelet_Transform(sp, w, h, fs, N, Func, w_a, K):
    """
    runs a specific wavelet through the signal
    """

    #sp = np.fft.fft(h)
    #w = np.fft.fftfreq(Func.size, 1/fs)
    #w = np.linspace(0, fs, N)*2*np.pi
    wavelet = 2*(np.exp(-(K*(w-w_a)/w_a)**2) - np.exp(-K**2)*np.exp(-(K*w/w_a)**2))

    return np.fft.ifft(sp*wavelet)

def Wavelet_diagram(h, t, Sampling, Fmax, Fmin, dt, K=8):
    """
    Runs Wavelet_Transform for all wavelets across the signal and compiles them into one diagram
    """
    dt = dt
    global action
    N = h.size
    Yr = dt*N
    Gfrate = N/(Yr*Hdef)
    fs = Sampling*Gfrate
    print(fs)
    Func = np.linspace(0, N/fs, N)
    omega_a = np.arange(Fmin, Fmax)*2*np.pi
    sp = np.fft.fft(h)*Gfrate
    w = np.linspace(0, fs, N)*2*np.pi
    wavelet_stuff = np.zeros((len(omega_a), len(Func)))
    print("Running Wavelet analasys!")

    #bar = ChargingBar('Processing', max = len(omega_a))
    import time
    StartT = time.perf_counter()
    for i in range(len(omega_a)):
        wavelet_stuff[i,:] = np.abs(Wavelet_Transform(sp, w, h, fs, N, Func, omega_a[i], K))
        print(f" {i/len(omega_a)*100:.1f}%", end='\r', flush=True)
        #bar.next()
    #bar.finish()
    print("")
    print("DONE!")
    print("")
    FinCompT = time.perf_counter()
    print(f"Elapsed time for computation: {(FinCompT - StartT):.1f} seconds")
    print("Starting matplotlib plotting..")
    X, Y = np.meshgrid(Func, omega_a/(2*np.pi))
    plt.title("Wavelet Analasys")
    plt.contourf(X, Y, wavelet_stuff,)
    plt.colorbar()
    plt.ylabel("Frekvens ['Hz']")
    plt.xlabel("tid")
    plt.title("{} with {} and {}, with k = {}".format(action,object1, object2, K))
    #plt.savefig("{}_{}_{}.png".format(object1,object2,action))
    EndT = time.perf_counter()
    print(f"Total elapsed time for wavelet process: {(EndT - StartT):.1f} seconds")
    print("")
    plt.show()


def ask2():
    print("Do you wish to do Wavelet analysis? (y,n)")
    ans = input("Answer: ")
    if ans == "y":
        print("Choose upper freq: ")
        Fmax = int(input("Fmax: "))
        print("Choose start freq:  ")
        Fmin = int(input("Fmin: "))
        print("Choose K-value")
        K = int(input("K: "))
        Wavelet_diagram(h[:(count+int((count*0.1)))], t, Sampling, Fmax, Fmin, dt, K)
    elif ans == "n":
        pass
    else:
        ask2()
ask2()





#Final results and plotting
fig, (ax1,ax2) = plt.subplots(1,2)

planet_orbit = planet_orbit[0,1:count+1,:]
sun_orbit = sun_orbit[1:count+1,:]


r = np.linalg.norm(planet_orbit[-1,:] - sun_orbit[-1,:])/2
print("Distanse is = {} AU".format(r))
P = -32/5*Gr**4/c**5*((sunmass*mass_plan1)**2*(sunmass + mass_plan1))/r**5
print("Energi radiert = {}".format(P))



ax1.plot(planet_orbit[:-1,0], planet_orbit[:-1,1], "--",label=object1)
ax1.plot(sun_orbit[:-1,0], sun_orbit[:-1,1],"--", label=object2)
ax1.set(xlabel=("x-position (AU)"),ylabel= ("y-position (AU)"))
ax1.legend(loc='lower right')


#If you wish to savefig, un-comment the line below
#plt.savefig("garvity-sun-pos.jpeg")
ax2.plot(np.linspace(0,t,count),h[:count])
ax2.set(xlabel=('Time(y)'),ylabel=('Distortion(Au)'))
plt.show() #viser banen til objektene

"""
Animation kode
"""
rat = 0.4
if object1 == 'bbh' or object2 == 'bbh':
    rat = 2

intr = int(1e3)
# First set up the figure, the axis, and the plot element we want to animate
fig, (ax1,ax2) = plt.subplots(1,2)
ax1.set(xlim=(0,10),ylim=(np.min(h)*1.5,np.max(h)*1.5), ylabel=('Distortion'))
ax2.set(xlim=(-rat*dist,rat*dist),ylim=(-rat*dist,rat*dist),xlabel='AU',ylabel='AU',title='Objects trajectory')
line1, = ax1.plot(np.linspace(0,10,intr), h[0:int(intr)], lw=2)
line2, = ax2.plot(planet_orbit[0:int(intr),0],planet_orbit[0:int(intr),1],label=object1,color=dict[object1][2])
line3, = ax2.plot(sun_orbit[0:int(intr),0],sun_orbit[0:int(intr),1],label=object2,color=dict[object2][2])
patch1 = plt.Circle((planet_orbit[0:int(intr),0],planet_orbit[0:int(intr),1]), radius_1,color=dict[object1][2])
patch2 = plt.Circle((sun_orbit[0:int(intr),0],sun_orbit[0:int(intr),1]), radius_2,color=dict[object2][2])
line = [line1,line2,line3]
centr =sun_orbit[-1,:] +  radius_2*(planet_orbit[-1,:]-sun_orbit[-1,:])/np.linalg.norm(planet_orbit[-1,:]-sun_orbit[-1,:])

def init3():
    line[0].set_data([], [])
    line[1].set_data([],[])
    line[2].set_data([], [])
    return line

def update(i,line,patch1,patch2):
    x = np.linspace(0,10,intr)
    y = h[int((i*intr)):int((i+1)*intr)]
    line[0].set_data(x, y)
    
    x = planet_orbit[int(i*intr):int((i+1)*intr),0]
    y = planet_orbit[int(i*intr):int((i+1)*intr),1]
    line[1].set_data(x, y)
    patch1.center= (x[-1],y[-1])
    ax2.add_patch(patch1)
    x = sun_orbit[int(i*intr):int((i+1)*intr),0]
    y = sun_orbit[int(i*intr):int((i+1)*intr),1]
    line[2].set_data(x, y)
    patch2.center= (x[-1],y[-1])
    ax2.add_patch(patch2)
    if (i >96*int(count/intr)/100 and answer2=='c'):
        imagebox = OffsetImage(mpimg.imread('boom.png'),zoom=(rat*3))
        ab = AnnotationBbox(imagebox, (centr),frameon=False)
        ax2.add_artist(ab)
        plt.draw()
    return line

# call the animator.  blit=True means only re-draw the parts that have changed.
anim1 = animation.FuncAnimation(fig, update, init_func=init3,
                               frames=int(count/intr), fargs=[line,patch1,patch2], interval=150, blit=True)

#anim1.save('anim1.gif', writer='imagemagick', fps = 60)
#anim2.save('anim2.gif', writer='imagemagick',fps = 60)
plt.legend()
plt.show()
