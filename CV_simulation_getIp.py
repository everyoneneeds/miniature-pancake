##/*-------------------------------------------------------------------
##  Copyright 2020-03-11 SUN, Sheng
##//Last Change:  2017-06-09 18:56:12
##
##The code is used to simulate cyclic voltametry curves
##
##The code is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU General Public License for more details.
##
##You should have received a copy of the GNU General Public License
##along with the code.  If not, see <http://www.gnu.org/licenses/>.
##----------------------------------------------------------------*/

#!/bin/python
import numpy as np
import matplotlib.pyplot as mpl
from scipy.interpolate import spline

class cycVol:
    def __init__(self, vol=1, Einit=0.5, Eend=-0.5):
        self.vol = vol
        self.Einit = Einit
        self.Eend = Eend
        self.htime = np.abs((Einit-Eend)/vol)
        print "One half scan time: ", self.htime
        self.ksi = []
        self.time = []
        self.SLambda = []
        self.V = []

    def initParams(self, E0=0, k0=1, alpha=0.5, DO=1e-5, DR=1e-5, nume=1, COinit=1e-6, area=1, Psi=0.5):
        self.E0 = E0   # default is Fe3+/Fe2+
        self.k0 = k0
        self.DO = DO
        self.DR = DR
        self.alpha = alpha
        self.gamma = np.sqrt(self.DO/self.DR)
        self.nume = nume
        self.COinit = COinit
        self.area = area
        Farad_const = 96485 #J/mol V
        Gas_const = 8.314 #J/mol K
        Temp = 300 #K
        self.a = self.nume*self.vol*Farad_const/Gas_const/Temp  #a=nFv/RT
        self.theta = np.exp(self.a/self.vol*(self.Einit-self.E0))
        # self.theta=self.COinit/(1e-100)
        self.Psi = self.gamma**(self.alpha)*self.k0/np.sqrt(np.pi*self.a*self.DO)
        # self.Psi = Psi

    def getIV(self, tnumber=2):
        self.delta = self.a*self.htime/tnumber
        # Single direction
        for currNu in np.arange(1, 2*tnumber):
            # print "Processing: ", currNu
            Slamb = self.get_V(currNu)
            self.get_ksi(currNu, Slamb)
        #print self.ksi

    def get_V(self, currN):
        currat = currN*self.delta
        self.time.append(currat/self.a)
        if currat <= self.a*self.htime:
            self.V.append(self.Einit-self.vol*currat/self.a)
            Slamb = np.exp(-1*currat)
            self.SLambda.append(Slamb)
        else:
            self.V.append(self.Einit-2*self.vol*self.htime+self.vol*currat/self.a)
            Slamb = np.exp(currat-2*self.a*self.htime)
            self.SLambda.append(Slamb)
        return Slamb

    def get_ksi(self, currN, Slamb):
        coeff1 = (self.gamma*self.theta*Slamb)**(self.alpha)/self.Psi
        coeff2 = 1-Slamb
        coeff3 = 1+self.gamma*self.theta*Slamb

        '''
        print "C1: ", coeff1
        print "C2: ", coeff2
        print "C3: ", coeff3
        '''

        if currN == 1:
            self.ksi.append(coeff2/(coeff1+2*np.sqrt(self.delta)*coeff3))
        else:
            coeffKsi = self.accumKsi(currN)
            ksiN = (coeff2-2*np.sqrt(self.delta)*coeffKsi*coeff3)/(coeff1+2*np.sqrt(self.delta)*coeff3)
            self.ksi.append(ksiN)

    def accumKsi(self, currN):
        n = len(self.ksi)
        assert n+1 == currN
        ksi1 = self.ksi[0]*np.sqrt(currN)
        accum = 0
        if currN == 2:
            return (np.sqrt(currN)-1)*self.ksi[0]
        else:
            aKsi1 = ksi1-self.ksi[-1]
            for k in np.arange(1, currN-1):  # Note that the last value is currN-2 by np.arange()
                accum += np.sqrt(currN-k)*(self.ksi[k]-self.ksi[k-1])
            return aKsi1+accum

    def printParams(self):

        print "\n Print parameters begin: "
        Farad_const = 96500 #J/mol V
        Gas_const = 8.314 #J/mol K
        Temp = 300 #K

        print "Theta: ", self.theta
        print "Gamma: ", self.gamma
        print "Psi: ", self.Psi
        print "a:", self.a
        Varray = np.array(self.V)
        print "time:", np.array(self.time)
        print "V:", Varray
        Iarray = np.array(self.ksi)*self.nume*Farad_const*self.area*self.COinit*np.sqrt(np.pi*self.a*self.DO)
        print  "Ksi: ", np.array(self.ksi)
        print  "Iarray:", Iarray
        print   "Salambda: ", np.array(self.SLambda)
        print "a=", self.a, "TimeLambda=", self.htime, "Delta is: ", self.delta
    
    def getEppIpf(self):
        assert len(self.V) == len(self.ksi)
        Farad_const = 96485 #J/mol V
        Gas_const = 8.314 #J/mol K
        Temp = 300 #K
        Varray =-1* np.array(self.V)
        ''' Pls note that when flux is positive, which means that oxides move toward the electrode and a cathodic current, the current is negative.''' 
        Iarray = -1* np.array(self.ksi)*-1.0*self.nume*Farad_const*self.area*self.COinit*np.sqrt(np.pi*self.a*self.DO)
        assert len(Varray) == len(Iarray)
        Imax = Iarray.max()*1000
        ixImax = Iarray.argmax()
        Imin = Iarray.min()*1000
        ixImin = Iarray.argmin()
        Vmax = Varray[ixImax]
        Vmin = Varray[ixImin]
        DeltaEpp = np.abs(Vmax-Vmin)
        return Imax, Imin, DeltaEpp, Vmax, Vmin

    def plotIV(self, ax, myls='-', mymarker=None, mylegend=r"example") :
        
        assert len(self.V) == len(self.ksi)
        Farad_const = 96485 #J/mol V
        Gas_const = 8.314 #J/mol K
        Temp = 300 #K
        Varray =-1* np.array(self.V)
        ''' Pls note that when flux is positive, which means that oxides move toward the electrode and a cathodic current, the current is negative.''' 
        Iarray = -1* np.array(self.ksi)*-1.0*self.nume*Farad_const*self.area*self.COinit*np.sqrt(np.pi*self.a*self.DO)
        ax.plot(Varray, Iarray*1000, 'b', ls=myls, marker=mymarker, mfc='none', label=mylegend)

def main():

    # from matplotlib import rcParams
    # rcParams.update({'figure.autolayout': True})
    # mpl.rc('text',usetex=True)
    # font={'family':"sans-serif"}
    # mpl.rc('font',**font)
    # mpl.rc('xtick',labelsize=20)
    # mpl.rc('ytick',labelsize=20)

    # fig=mpl.figure()
    # ax=mpl.gca()
    # ax.set_xlabel(r"-E/V", fontsize=20)
    # ax.set_ylabel(r"-I/mA", fontsize=20)
    # minor_ticks = np.arange(-1.0, 1.0, 0.2)
    # major_ticks = np.arange(-1.0, 1.0, 0.4)
    # ax.set_yticks(major_ticks)
    # ax.set_yticks(minor_ticks, minor=True)

    myDl=1e-6
    myDs=1e-6
    myVol=1e-2
    myk0=1e-0
    myCmax = 5 
    myCinit =  1e-6

    Epp = []
    for myDl in np.arange(1e-7, 1e-6, (1e-6-1e-7)/20):
        for myVol in np.logspace(-4, 4, num=20):
            for myk0 in np.logspace(-18.42, 4.7, num=20):
                for myCinit in np.arange(5e-7, 5e-6, (5e-6-5e-7)/20):
                    myCV = cycVol(vol=myVol)
                    myCV.initParams(k0=myk0, DR=myDs, DO=myDl, COinit=myCinit)
                    myCV.getIV(100)
                    Imax, Imin, DVpp, Vmax, Vmin = myCV.getEppIpf()
                    Epp.append([Imax, Imin, DVpp, myk0, myVol, myDl, Vmax, Vmin, myCinit])
    npEpp = np.array(Epp)
    np.savetxt("Ip.dat", npEpp)
   
# import subprocess

if __name__ == "__main__":
    main()
    # subprocess.call('evince ./CV_curve.pdf', shell = True)
