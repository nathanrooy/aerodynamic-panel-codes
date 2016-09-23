#------------------------------------------------------------------------------+
#
#   Nathan A. Rooy
#   Inviscid Constant Strength Doublet Numerical Panel Code
#   Low-Speed Aerodynamics by Joseph Katz and Allen Plotkin (2nd Edition)
#   2016-SEPT-07
#
#------------------------------------------------------------------------------+

#--- IMPORT DEPENDENCIES ------------------------------------------------------+

from __future__ import division
from collections import defaultdict
from math import atan2
from math import cos
from math import pi
from math import sin
from math import sqrt

import numpy as np
import re

#--- FUNCTIONS ----------------------------------------------------------------+

# DETERMINE VELOCITY POTENTIAL
def phi(panel_i,panel_j):

    # calculate x2
    x2=cos(panel_j.aoa)*(panel_j.x2-panel_j.x1)-sin(panel_j.aoa)*(panel_j.y2-panel_j.y1)

    # translate from global coordinates to panel coordinates (eqn. 11.23a, pg.277)
    x=cos(panel_j.aoa)*(panel_i.xc-panel_j.x1)-sin(panel_j.aoa)*(panel_i.yc-panel_j.y1)   
    y=sin(panel_j.aoa)*(panel_i.xc-panel_j.x1)+cos(panel_j.aoa)*(panel_i.yc-panel_j.y1)

    # calculate velocity potential (eqn. 10.28, pg. 235)
    phi=-(panel_j.mu/(2*pi))*((atan2(y,x-x2))-(atan2(y,x)));

    return phi

# PANEL CLASS
class panel:
    def __init__(self,x1,y1,x2,y2):
        self.chord=sqrt((x2-x1)**2+(y2-y1)**2)
        self.x1=x1                          # panel left x-coordinate
        self.y1=y1                          # panel left y-coordinate
        self.x2=x2                          # panel right x-coordinate
        self.y2=y2                          # panel right y-coordinate
        self.xc=(x2+x1)/2                   # panel x-coordinate collocation point
        self.yc=(y2+y1)/2                   # panel y-coordinate collocation point
        self.aoa=atan2(y1-y2,x2-x1)         # panel angle of attack
        self.mu=1                           # panel doublet strength
        self.v_t=0                          # panel tangential velocity
        self.cp=0                           # panel Cp
        
# DETERMINE PANELS
def determinePanels(x,y):
    panels=[]
    M=len(x)
    # populate panel list with panel objects
    for i in range(M):
        # add additional wake panel
        if i==M-1:
            panels.append(panel(x[-1],y[-1],100000,y[-1]))
        # normal panel
        else:
            panels.append(panel(x[i],y[i],x[i+1],y[i+1]))
    return panels

# READ IN THE AIRFOIL COORDINATES FROM A TEXT FILE
def getCoords(fileName):
    x_list=[]
    y_list=[]
    with file(settingsDict['FileName'],'r') as csvFile:
        for line in csvFile:
            x,y=re.findall(r'(^-?\d*.?\d*)\s*(-?\d*.?\d*)',line)[0]
            try:
                x=float(x)
                y=float(y)
                x_list.append(x)
                y_list.append(y)
            except:
                pass
    x_list.reverse()
    y_list.reverse()
    return x_list,y_list

#--- MAIN ---------------------------------------------------------------------+

def main(settingsDict):
    #--- FREE STREAM INITIALIZATION -----------------------+
    settingsDict['Alpha radians']=settingsDict['Alpha degrees']*pi/180                              # Angle of Attack
    settingsDict['u-inf']=settingsDict['Freestream Velocity']*cos(settingsDict['Alpha radians'])    # x-velocity    
    settingsDict['v-inf']=settingsDict['Freestream Velocity']*sin(settingsDict['Alpha radians'])    # y-velocity

    #--- READ IN AIRFOILL COORDINATES FROM TEXT FILE ------+
    x,y=getCoords(settingsDict['FileName'])
    M=len(x)
    print 'NUMBER OF POINTS:',M

    #--- DETERMINE PANELS AND GEOMETETRIC PROPERTIES ------+
    panels=determinePanels(x,y)
    print 'NUMBER OF PANELS:',len(panels)

    #--- CREATE LHS ---------------------------------------+
    A=np.zeros((M,M),dtype=float)
    for i in range(M-1):
        for j in range(M):
            if i==j:
                A[i][j]=0.5
            else:
                A[i][j]=phi(panels[i],panels[j])

    #--- EXPLICIT KUTTA CONDITION -------------------------+
    A[-1][0]=1
    A[-1][-2]=-1
    A[-1][-1]=1

    #--- CREATE RHS ---------------------------------------+
    RHS=np.zeros((M,1),dtype=float)
    for i in range(M-1):
        RHS[i]=-settingsDict['u-inf']*panels[i].xc-settingsDict['v-inf']*panels[i].yc

    #--- SOLVE LINEAR SET OF EQUATIONS --------------------+
    mu=np.linalg.solve(A,RHS)   # eqn. 11.7, pg. 266

    #--- CALCULATE Cp AND Cl ------------------------------+
    Cl=0
    for j in range(M-1):
        # first panel on airfoil
        if i==0:
            delta_l=sqrt((panels[0].xc-panels[1].xc)**2+(panels[0].yc-panels[1].yc)**2);
            panels[j].v_t=(mu[0]-mu[1])/delta_l;
            
        # last panel on airfoil (excluding wake panel)
        elif j==M-2:
            delta_l=sqrt((panels[-3].xc-panels[-2].xc)**2+(panels[-3].yc-panels[-2].yc)**2);
            panels[j].v_t=(mu[-3]-mu[-2])/delta_l;
            
        # standard panel
        else:
            delta_l=sqrt((panels[j].xc-panels[j+1].xc)**2+(panels[j].yc-panels[j+1].yc)**2);
            panels[j].v_t=(mu[j]-mu[j+1])/delta_l;           

        panels[j].cp=1-(panels[j].v_t**2/settingsDict['Freestream Velocity']**2)    # eq. 11.77, pg. 293
        Cl=Cl-panels[j].v_t*delta_l

    #--- PLOT Cp ------------------------------------------+
    if settingsDict['Plot Cp?']==True:
        import matplotlib.pyplot as plt
        pxc=[]
        pcp=[]
        for i in range(M-1):
            pxc.append(panels[i].xc)
            pcp.append(panels[i].cp)
        plt.figure()
        plt.plot(pxc,pcp)
        plt.gca().invert_yaxis()
        #plt.title('Coeffcient of Pressure Distribution over the surface')
        plt.xlabel('x')
        plt.ylabel("C_p")
        plt.xlim(0,1)
        plt.show()

    #--- PRINT RESULTS ------------------------------------+
    print 'Cl: %0.5f'%Cl
    return

#--- RUN ----------------------------------------------------------------------+

settingsDict=defaultdict(int)

settingsDict['FileName']='naca2412.dat'
settingsDict['ReynoldsNumber']=2000
settingsDict['Alpha degrees']=0
settingsDict['Freestream Velocity']=1.0
settingsDict['Plot']=False
settingsDict['Close Trailing Edge']=True
settingsDict['Plot Cp?']=True

main(settingsDict)

#--- END ----------------------------------------------------------------------+
