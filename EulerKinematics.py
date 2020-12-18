# matplotlib 3.3.1
import matplotlib.pyplot as plt     #Plotting Library
import numpy as np      #Math Library
import math     #python's math module

def main():
    """
    This script takes any initial values and constants and calculates and plots the trajectory
    of an object assuming constant gravitational force and typical drag force equation.

    Uses the Euler method, assuming that acceleration is constant between time intervals.

    Units: kg, m, sec, rad
    """
    #Drag Force Equation: 1/2 * rho * Cd * A * v^2

    #User-Defined Constants
    global m
    global v0
    global theta
    global rho  #Fluid Density
    global A    #Cross-sectional Area
    global Cd   #Drag coefficient
    global tStep
    global g

    m = 1
    v0 = 30
    theta = math.radians(45)
    rho = 1.225
    A = 0.05
    Cd = 0.5    #A ball is approx. 0.5
    tStep = 0.005
    g = 9.8


    #Data Structures
    global tHist
    global xHist
    global yHist
    global vHist
    global vXHist
    global vYHist
    global aXHist
    global aYHist
    tHist = []      #list for all time steps
    xHist = []      #list for all x position steps
    yHist = []      #list for all y position steps
    vHist = []      #list for all velocities at every time step
    vXHist = []     #list for all x-axis velocities at every time step
    vYHist = []     #list for all y-axis velocities at every time step
    aXHist = []     #list for all x-axis accelerations at every time step
    aYHist = []     #list for all y-axis accelerations at every time step

    #Initialize initial values
    tHist.append(0.0)
    xHist.append(0.0)
    yHist.append(0.0)
    vHist.append(v0)
    vXHist.append(v0 * math.cos(theta))
    vYHist.append(v0 * math.sin(theta))
    vTheta = math.atan(vYHist[0] / vXHist[0])
    aXHist.append(-0.5*rho*Cd*A*vHist[0]**2*math.cos(vTheta) / m)
    aYHist.append(-g - 0.5*rho*Cd*A*vHist[0]**2*math.sin(vTheta) / m)
    # print("t: " + str(tHist[0]))
    # print("x: " + str(xHist[0]))
    # print("y: " + str(yHist[0]))
    # print("Vx: " + str(vXHist[0]))
    # print("Vy: " + str(vYHist[0]))
    # print("Ax: " + str(aXHist[0]))
    # print("Ay: " + str(aYHist[0]))

    counter = 1
    #Loop until the y-displacement becomes negative (projectile reaches ground again)
    while True:
        tHist.append(counter * tStep)   #increment time
        vXHist.append(vXHist[counter-1] + aXHist[counter-1]*tStep)  #calculate new x velocity using v = v0 + at
        vYHist.append(vYHist[counter-1] + aYHist[counter-1]*tStep)  #calculate new x velocity using v = v0 + at
        vHist.append(math.hypot(vXHist[counter], vYHist[counter]))   #calculate new net velocity from x and y velocities

        vTheta = math.atan(vYHist[counter] / vXHist[counter])   #calculates current angle of motion using arctan
        vSquared = vHist[counter]**2    #temporary, convenience variable
        xDragAccel = -0.5*rho*Cd*A*vSquared*math.cos(vTheta) / m                                    #calculation for acceleration due to drag in the x-axis
        yDragAccel = -math.copysign(0.5*rho*Cd*A*vSquared*math.sin(vTheta) / m, vYHist[counter])    #calculation for net acceleration in the y-axis

        #This velocity calculation is incorrect; double counts cos and sin
        # vXSquared = vXHist[counter]**2
        # vYSquared = vYHist[counter]**2
        # xDragAccel = -0.5*rho*Cd*A*vXSquared / m
        # yDragAccel = -math.copysign(0.5*rho*Cd*A*vYSquared / m, vYHist[counter])

        #Add the accelerations to their respective lists
        aXHist.append(xDragAccel)
        aYHist.append(-g + yDragAccel)

        tSquared = tStep**2     #temporary, convenience variable
        xHist.append(xHist[counter-1] + vXHist[counter]*tStep + 0.5*aXHist[counter]*tSquared)   #x = x0 + v0t + 1/2 * at^2
        yHist.append(yHist[counter-1] + vYHist[counter]*tStep + 0.5*aYHist[counter]*tSquared)   #y = y0 + v0t + 1/2 * at^2

        # xHist.append(xHist[counter-1] + vXHist[counter]*tStep)
        # yHist.append(yHist[counter-1] + vYHist[counter]*tStep)

        if vYHist[counter-1] > 0.0 and vYHist[counter] < 0.0:   #Find when the max height has been reached by finding a critical point
            print("max height reached at time=" + str(tHist[counter]))

        print("t:  " + str(tHist[counter]))
        print("x:  " + str(xHist[counter]))
        print("y:  " + str(yHist[counter]))
        print("0:  " + str(vTheta))
        print("Vx: " + str(vXHist[counter]))
        print("Vy: " + str(vYHist[counter]))
        print("Ax: " + str(aXHist[counter]))
        print("Ay: " + str(aYHist[counter]))
        print("")

        if yHist[counter] < 0 or counter > 99999:   #End the loop if the projectile has reached the ground (or limit the number of iterations to avoid computer death)
            break

        counter += 1

    plotData()


def plotData():
    print("t: " + str(tHist[len(tHist)-1]))
    plt.plot(tHist, xHist, marker='o', markersize=2, label='x')
    plt.plot(tHist, yHist, marker='o', markersize=2, label='y')
    plt.grid()
    plt.xlabel('time')
    plt.legend()
    plt.show()

    plt.plot(xHist, yHist, marker='o', markersize=2, label='position')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    plt.plot(tHist, vXHist, marker='o', markersize=2, label='x velocity')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('x velocity')
    plt.show()

    plt.plot(tHist, aXHist, marker='o', markersize=2, label='x acceleration')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('x acceleration')
    plt.show()

    plt.plot(tHist, vYHist, marker='o', markersize=2, label='y velocity')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('y velocity')
    plt.show()

    plt.plot(tHist, aYHist, marker='o', markersize=2, label='y acceleration')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('y acceleration')
    plt.show()


if __name__ == "__main__":  #Run the main function
    main()
