# matplotlib 3.3.1
import matplotlib.pyplot as plt #Plotting Library
import numpy as np  #Math library
import scipy.integrate as integrate     #Integration library
import scipy.special as special         #Integration library
from sympy import Symbol, sin, cos, tan, sec, csc, cot  #Trig Library
from sympy.integrals.trigonometry import trigintegrate  #Integration library
from sympy.abc import x     #Integration library
import math     #python's math module

def main():
    """
    This script takes any initial values and constants and calculates and plots the trajectory
    of an object assuming constant gravitational force and typical drag force equation.

    This is another method that I invented. It uses separation of variables to solve
    for the velocity at every time step, although it assumes a constant angle between
    time steps.

    Note: I made some use of online integral calculators.

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
    global thetaHist
    global vHist
    global vXHist
    global vYHist
    tHist = []      #list for all time steps
    xHist = []      #list for all x position steps
    yHist = []      #list for all y position steps
    thetaHist = []  #List for all theta at every time step
    vHist = []      #list for all velocities at every time step
    vXHist = []     #list for all x-axis velocities at every time step
    vYHist = []     #list for all y-axis velocities at every time step

    #Initialize intial values
    tHist.append(0.0)
    xHist.append(0.0)
    yHist.append(0.0)
    thetaHist.append(theta)
    vHist.append(v0)
    vXHist.append(v0 * math.cos(theta))
    vYHist.append(v0 * math.sin(theta))
    vTheta = math.atan(vYHist[0] / vXHist[0])
    # print("t: " + str(tHist[0]))
    # print("x: " + str(xHist[0]))
    # print("y: " + str(yHist[0]))
    # print("v: " + str(vHist[0]))
    # print("Vx: " + str(vXHist[0]))
    # print("Vy: " + str(vYHist[0]))

    #Convenience variables
    global k

    counter = 1
    #Loop until the y-displacement becomes negative (projectile reaches ground again)
    while True:
        tHist.append(counter * tStep)   #increment time
        print("t:  " + str(tHist[counter]))

        #This large hunk is the solution to the net force differential equation in the x-axis
        # oneOverVX = (1/vXHist[counter-1]) + (((rho*A*Cd*math.cos(thetaHist[counter-1]))/(2*m))*(tStep))   #STABLE
        # oneOverVX = (1/vXHist[counter-1]) + (((rho*A*Cd)/(2*m))*(tStep))
        # oneOverVX = (1/vHist[counter-1]) + (((rho*A*Cd*math.cos(thetaHist[counter-1]))/(2*m))*(tStep))
        oneOverVX = (1/vXHist[counter-1]) + ((rho*A*Cd)/(2*m*math.cos(thetaHist[counter-1]))*(tStep))   #This is one over the solution for velocity in the x-axis net force differential equation
        vXHist.append(1 / oneOverVX)    #Adding the velocity to the list of velocities

        vY0 = vYHist[counter-1]     #Convenience variable
        # k = 0.5 * rho * A * Cd * math.sin(abs(thetaHist[counter-1]))  #STABLE
        # k = 0.5 * rho * A * Cd
        k = (rho * A * Cd) / (2 * math.sin(abs(thetaHist[counter-1])))  #Convenience variable
        print("k:  " + str(k))
        print("vX: " + str(vXHist[counter]))
        rootGMK = math.sqrt(g*m*k)  #Convenience variable
        if vYHist[counter-1] > 0.0:     #If the projectile is going upwards
            #Solving the y-axis differential equation for velocity
            equationRight = -rootGMK * ((tStep/m) - (math.atan((k*vY0)/(rootGMK))/rootGMK))
            vYHist.append((math.tan(equationRight) * rootGMK) / k)
        elif vYHist[counter-1] < 0.0:   #If the projectile is going downwards
            #Solving the y-axis differential equation for velocity

            # Hand-solved integral
            # exponent = -(2*tStep*rootGMK)/m
            # numerator = g*m*math.exp(exponent) - math.exp(exponent)*vY0*rootGMK - vY0*rootGMK - g*m
            # denominator = math.exp(exponent)*(vY0-rootGMK) - vY0*k - rootGMK
            # vYHist.append(numerator / denominator)

            #Wolfram Alpha arctanh integral
            arctanh =(vY0*math.sqrt(k))/(math.sqrt(g*m))
            print("arctanh: " + str(arctanh))
            equationRight = (np.arctanh(arctanh))/(rootGMK) - (tStep/m)
            vYHist.append(np.tanh(rootGMK * equationRight) * ((math.sqrt(g*m))/(math.sqrt(k))))
        else:   #If current y velocity is 0
            vYHist.append(vY0 - g*tStep)
        print("vY: " + str(vYHist[counter]))

        vHist.append(math.hypot(vXHist[counter], vYHist[counter]))  #Calculate the net velocity and add it to the velocities list
        print("v:  " + str(vHist[counter]))
        thetaHist.append(math.atan(vYHist[counter]/vXHist[counter]))    #Calculate the current angle based on the velocities and add it to the theta list
        print("0:  " + str(math.degrees(thetaHist[counter])))

        x0 = xHist[counter-1]
        y0 = yHist[counter-1]

        # yIntegral = trigintegrate()

        """
        Note: What I wanted to do here was to integrate the velocity functions over the time interval to find the exact
        changes in position. Unfortunately, I was running short of time and decided it was not worth it to move forward with
        this final step, and instead worked on the presentation and testing different cases.
        """
        xHist.append(x0 + vXHist[counter]*tStep)    #Calculate new x position using x = x0 + vt
        yHist.append(y0 + vYHist[counter]*tStep)    #Calculate new y position using y = y0 + vt
        print("x:  " + str(xHist[counter]))
        print("y:  " + str(yHist[counter]))
        print()

        # xHist.append(xHist[counter-1] + vXHist[counter-1]*tStep + 0.5*aXHist[counter-1]*tStep**2)
        # yHist.append(yHist[counter-1] + vYHist[counter-1]*tStep + 0.5*aYHist[counter-1]*tStep**2)
        # vXHist.append(vXHist[counter-1] + aXHist[counter-1]*tStep)
        # vYHist.append(vYHist[counter-1] + aYHist[counter-1]*tStep)
        # vHist.append(math.hypot(vXHist[counter], vYHist[counter]))
        #
        # vTheta = math.atan(vYHist[counter] / vXHist[counter])
        # xDragAccel = -0.5*rho*Cd*A*vHist[counter]**2*math.cos(vTheta) / m
        # yDragAccel = -math.copysign(0.5*rho*Cd*A*vHist[counter]**2*math.sin(vTheta) / m, vYHist[counter])
        #
        # aXHist.append(xDragAccel)
        # aYHist.append(-g*tStep + yDragAccel)

        if vYHist[counter-1] > 0.0 and vYHist[counter] < 0.0:   #Check if the projectile has reached it's peak by checking for a critical point
            print("max height reached at time=" + str(tHist[counter]))
            # break

        # print("t: " + str(tHist[counter]))
        # print("x: " + str(xHist[counter]))
        # print("y: " + str(yHist[counter]))
        # print("Vx: " + str(vXHist[counter]))
        # print("Vy: " + str(vYHist[counter]))
        # print("Ax: " + str(aXHist[counter]))
        # print("Ay: " + str(aYHist[counter]))

        if yHist[counter] < 0 or counter > 99999:   #End the loop if the projectile has reached the ground (or limit the number of iterations to avoid computer death)
            break

        counter += 1

    plotData()


def plotData():
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

    for i in range(len(thetaHist)):
        thetaHist[i] = math.degrees(thetaHist[i])
    plt.plot(tHist, thetaHist, marker='o', markersize=2, label='theta')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('theta')
    plt.legend()
    plt.show()


    plt.plot(tHist, vXHist, marker='o', markersize=2, label='x velocity')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('x velocity')
    plt.show()

    plt.plot(tHist, vYHist, marker='o', markersize=2, label='y velocity')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('y velocity')
    plt.show()


if __name__ == "__main__":  #Run the main function
    main()
