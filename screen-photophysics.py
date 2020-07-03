import os
import sys
import matplotlib.pyplot as plt
import glob
import pandas
import numpy as np
from scipy.interpolate import UnivariateSpline

#define cis and trans
cislog=sys.argv[1]
translog=sys.argv[2]
cis=cislog.split('.log')[0]
trans=translog.split('.log')[0]

#how to name output files/plot
name=cis

#get the energies
os.system("""cat {0} |grep Excited|grep -v singles|sed 's/=/\t/g'|gawk '{{print $7,$10}}' > {1}.impulses""".format(cislog,cis))
os.system("""cat {0} |grep Excited|grep -v singles|sed 's/=/\t/g'|gawk '{{print $7,$10}}' > {1}.impulses""".format(translog,trans))

path=os.getcwd()

names = ['{0}.impulses'.format(trans), '{0}.impulses'.format(cis)]

#constants, broadening parameter, and number of states to expect
nmol = 1
nstates=10
sigma=0.2
sigmain=sigma
h=6.62607015e-34
cmeters=299792458
c=cmeters*100
eV2J=1.602176634e-19
sigmacm=((sigma*eV2J)/(h*c))

#setting up spectra
start=380
stop=800
wavelengths=np.arange(start,stop,0.5)
extra=len(names)+1
linearcomb=np.zeros((len(wavelengths),extra))

tableau100 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),   
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),   
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),   
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),   
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),(31, 119, 180),(31, 119, 180)]

for i in range(len(tableau100)):
    r, g, b = tableau100[i]
    tableau100[i] = (r / 255., g / 255., b / 255.)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
fig =  plt.figure()
ax1 = fig.add_subplot(111)

#actual calculation starts here
for x,n in enumerate(names):

    inputfile = n
    lines=file_len(inputfile)
    data=np.loadtxt(inputfile)
    bound1 = int((lines/nmol)-1)
    bound2 = int((lines)-1)
    data1=data[0:bound1,:]
    nmol2=data.shape[0]
    nsampled=(nmol2/nstates)
  
#Generate gaussian curves

    gauss=np.zeros((bound1,len(wavelengths)))

    for i in range(0,len(wavelengths)):

        gauss[:,i]=1.3062974e8 * (data1[:,1]/sigmacm) * np.exp(-(((1/wavelengths[i]-1/data1[:,0])/(sigmacm*10e-8))**2))

    linearcomb[:,x] = np.sum(gauss,0)

 
#linearcomb=linearcomb/100
max=np.amax(linearcomb)
negmax=max*(-1)
linearcomb[:,1]=linearcomb[:,1]*(-1)  

#PLOTTING
for x,n in enumerate(names):

    ev=1240/wavelengths

    colors=[1,7]
    ax1.scatter(wavelengths,linearcomb[:,x],s=10,alpha='0.8', c=tableau100[colors[x]] )

difference=linearcomb[:,0] + linearcomb[:,1]

#difference=difference/(np.amax(difference))


ax1.scatter(wavelengths,difference,s=2, c='k', label = 'difference' )
ax1.scatter(wavelengths,np.zeros(len(wavelengths)), s=5, c='k')
y=difference
x=wavelengths
y_spl = UnivariateSpline(x,y,s=0,k=4)
#ax1.scatter(wavelengths,y_spl(wavelengths),s=1,label = 'fit')
y_spl_d =  y_spl.derivative(n=1)
droots=y_spl_d.roots()
droots=np.round(droots)
ax1.scatter(wavelengths,y_spl_d(wavelengths),s=1,c=tableau100[22],label = 'dy/dx')

y_spl_d2 =  y_spl.derivative(n=2)
d2rootval=y_spl_d2(droots)

ax1.scatter(wavelengths,y_spl_d2(wavelengths),s=1,label = 'd2y/dx2',c=tableau100[40])

#Determine whether passes or fails
transbool=False
cisbool=False


cispeakvalue=np.zeros(droots.size)
transpeakvalue=np.zeros(droots.size)
ciscount=0
transcount=0
for x,i in enumerate(d2rootval):
    index=np.where(wavelengths == droots[x])
    if i > 0:
        if (difference[index]/linearcomb[index,1]) > 0.6:
            cisbool=True
            cispeakvalue[ciscount] = wavelengths[index]
            ciscount+=1
                
    elif i < 0: 
        if difference[index]/linearcomb[index,0] > 0.6:
            transbool=True
            transpeakvalue[transcount]= wavelengths[index]
            transcount+=1
    else:
        pass

transpeakvalue=transpeakvalue[transpeakvalue != 0]
cispeakvalue=cispeakvalue[cispeakvalue != 0]

with open('{0}.txt'.format(name),'w') as screen:
    screen.write("""trans: {0} at {2}
cis: {1} at {3}""".format(trans,cis, transpeakvalue,cispeakvalue))

for i in transpeakvalue:
    i=int(i)
    index=np.where(wavelengths == i)
    ax1.scatter(i,difference[index],s=100, c='y', label = 'HIT-trans',marker='*')
for i in cispeakvalue:
    i=int(i)
    index=np.where(wavelengths == i)
    ax1.scatter(i,difference[index],s=100, c='y', label = 'HIT-cis', marker='*')

#make plot png
plt.legend()

ax1.set_xlabel('Wavelength (nm)',fontsize=15)
ax1.set_ylabel('Normalized Intensity',fontsize=15)
plt.title('{0}'.format(name))
plt.tight_layout()
plt.savefig('{0}.png'.format(name))
plt.show()

