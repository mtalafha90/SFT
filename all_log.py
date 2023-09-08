#!/usr/bin/env/ python
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from gekko import GEKKO 
from lmfit.models import PolynomialModel,ExpressionModel
from scipy.optimize import curve_fit
'''
def fitting_KK(ampl,findip):
	model = PolynomialModel(4, method='least_squares')
	params = model.guess(findip, x=ampl)
	result = model.fit(findip, params, x=ampl)
	return result.best_fit

'''
def lin(x,y0,x0):
	return y0/x0*x

def pol(x,c1,c2,c3):
	return c1*(x)**2+c2*(x)+(c3)

def fitting_KK(ampl,findip,c0):
	init_vals = [1,1,1]
	best_vals, covar = curve_fit(pol, ampl, findip, p0=init_vals)
	return best_vals

def fitting_lin(ampl,findip):
	init_vals = [1,1]
	best_vals, covar = curve_fit(lin, ampl, findip, p0=init_vals)
	return best_vals

def fitting_lm(ampl,findip):
	model = PolynomialModel(1, method='least_squares')
	params = model.guess(findip, x=ampl)
	result = model.fit(findip, params, x=ampl)
	#print(result.fit_report())
	c0 = result.params['c0'].value
	return result.best_fit,c0
'''
def fit_exp(ampl,findip,c0):
	#gmod = ExpressionModel("c1*(1+b)*x-c1*(b/x0)*x**2")
	gmod = ExpressionModel("c1*x-c2*(x-x0)**2+(c0*x0)")
	params = gmod.guess(findip, x=ampl)
	result = gmod.fit(final_dip, params,x=ampl)
	c1 = result.params['c1'].value
	c2 = result.params['c2'].value
	#print("b=",b)
	print("c1=",c1)
	print("c2=",c2)
	return result.best_fit
'''
##########################################################################
maxamp=[]
final_dip=[]
u0=[]
eta=[]
tau=[]
v00=[]
blat=[]

datafile='params2.dat'
f2=open(datafile,"r")
lines1=f2.readlines()
for x in lines1:
	maxamp.append(float(x.split('	')[0]))
	final_dip.append(float(x.split('	')[13]))
	u0.append(float(x.split('	')[14]))
	eta.append(float(x.split('	')[15]))
	tau.append(float(x.split('	')[16]))
	blat.append(float(x.split('	')[17]))
	v00.append(float(x.split('	')[19]))	
f2.close()

#maxamp=np.log(np.array((maxamp)))#*100
#final_dip=np.log(np.array((final_dip)))#*0.15

maxamp_aft=[]
final_dip_aft=[]
u0_sel = float(sys.argv[1])
eta_sel=int(sys.argv[2])
tau_sel=int(sys.argv[3])


for i in range (len(u0)):
	if (v00[i]==0 and blat[i]==0 and tau[i]==tau_sel and eta[i]==eta_sel):
		maxamp_aft.append(maxamp[i])
		final_dip_aft.append(final_dip[i])


fig = plt.figure()
fig.suptitle('$u_{0} =$'+str(u0_sel)+ ',$\eta=$'+str(eta_sel)+', $ \\tau=$'+str(tau_sel)+', $v_{0}=5$,'+' $\lambda_{c}=15$' )
plt.ylabel('$|\\Delta D(t)|$', fontsize=14)
plt.xlabel('$S_{max},n$', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

maxamp=(np.array((maxamp_aft))*400)#
final_dip=(np.array((final_dip_aft)))#

f=fitting_lin(maxamp,final_dip)
c0=f[0]/f[1]
y0=f[0]
x0=f[1]
#print('y0 =',y0,'x0 =',x0)
#x=np.linspace(4,6.5,100)
standard_dev = np.std(final_dip, axis=0)
fitted=lin(maxamp,f[0],f[1])
xs, ys = zip(*sorted(zip(maxamp, fitted)))

#plt.fill_between(x, fitted-standard_dev, fitted+standard_dev,alpha=0.2,color='#1f77b4')
plt.plot(xs,ys,'b',markersize=2,label='NoQ + No inflow')
plt.scatter(maxamp,final_dip,s=1.9,alpha=0.5,c='#1f77b4')
#maxamp*=0
#final_dip*=0
###############################################################################

maxamp1=[]
final_dip1=[]
u01=[]
eta1=[]
tau1=[]
v001=[]
blat1=[]
datafile='params2.dat'
f3=open(datafile,"r")
lines2=f3.readlines()
for x in lines2:
	maxamp1.append(float(x.split('	')[0]))
	final_dip1.append(float(x.split('	')[13]))
	u01.append(float(x.split('	')[14]))
	eta1.append(float(x.split('	')[15]))	
	tau1.append(float(x.split('	')[16]))
	blat1.append(float(x.split('	')[17]))
	v001.append(float(x.split('	')[19]))	
f3.close()


maxamp_aft1=[]
final_dip_aft1=[]

for i in range (len(u01)):
	if (v001[i]==5 and blat1[i]==0 and tau1[i]==tau_sel and eta1[i]==eta_sel):
		maxamp_aft1.append(maxamp1[i])
		final_dip_aft1.append(final_dip1[i])

maxamp1=(np.array((maxamp_aft1))*400)#
final_dip1=(np.array((final_dip_aft1)))#



s=fitting_KK(maxamp1,final_dip1,c0)
#print ('c1=',s[0],'c2=',s[1],'x0=',s[2],'c0=',s[3])
'''

y=np.linspace(-1.5,0.9,100)
standard_dev = np.std(final_dip1, axis=0)
fitted=pol(y,s[0],s[1])
#plt.fill_between(y, fitted-standard_dev, fitted+standard_dev,alpha=0.2,color='#3cff33')
'''
fitted1=pol(maxamp1,s[0],s[1],s[2])
xs1, ys1 = zip(*sorted(zip(maxamp1, fitted1)))

plt.plot(xs1,ys1,'g',markersize=2,label='Inflow')
plt.scatter(maxamp1,final_dip1,s=1.9,alpha=0.5,c='#3cff33')

##############################################################################
maxamp2=[]
final_dip2=[]
u02=[]
eta2=[]
tau2=[]
v002=[]
blat2=[]
datafile='params2.dat'
f4=open(datafile,"r")
lines3=f4.readlines()
for x in lines3:
	maxamp2.append(float(x.split('	')[0]))
	final_dip2.append(float(x.split('	')[13]))
	u02.append(float(x.split('	')[14]))
	eta2.append(float(x.split('	')[15]))	
	tau2.append(float(x.split('	')[16]))
	blat2.append(float(x.split('	')[17]))
	v002.append(float(x.split('	')[19]))
f4.close()

maxamp_aft2=[]
final_dip_aft2=[]

for i in range (len(u02)):
	if (v002[i]==0 and blat2[i]==2.4 and tau2[i]==tau_sel and eta2[i]==eta_sel):
		maxamp_aft2.append(maxamp2[i])
		final_dip_aft2.append(final_dip2[i])

maxamp2=(np.array((maxamp_aft2))*400)
final_dip2=(np.array((final_dip_aft2)))

#print ('c1=',t[0],'c2=',t[1],'x0=',t[2],'c0=',t[3])

'''
x=np.linspace(4,6.5,100)
standard_dev = np.std(final_dip, axis=0)
fitted=pol(x,t[0],t[1],t[2],t[3],t[4])
plt.fill_between(x, fitted-standard_dev, fitted+standard_dev,alpha=0.2,color='#e933ff')
'''
#print(final_dip2)
t=fitting_KK(maxamp2,final_dip2,c0)
fitted2=pol(maxamp2,t[0],t[1],t[2])
xs2, ys2 = zip(*sorted(zip(maxamp2, fitted2)))

plt.plot(xs2,ys2,'m',markersize=2,label='LQ')
plt.scatter(maxamp2,final_dip2,s=1.9,alpha=0.5,c='#e933ff')

##############################################################################

maxamp3=[]
final_dip3=[]
u03=[]
eta3=[]
tau3=[]
v003=[]
blat3=[]
datafile='params2.dat'
f5=open(datafile,"r")
lines4=f5.readlines()
for x in lines4:
	maxamp3.append(float(x.split('	')[0]))
	final_dip3.append(float(x.split('	')[13]))
	u03.append(float(x.split('	')[14]))
	eta3.append(float(x.split('	')[15]))	
	tau3.append(float(x.split('	')[16]))
	blat3.append(float(x.split('	')[17]))
	v003.append(float(x.split('	')[19]))
f5.close()

maxamp_aft3=[]
final_dip_aft3=[]

for i in range (len(u03)):
	if (v003[i]==5 and blat3[i]==2.4 and tau3[i]==tau_sel and eta3[i]==eta_sel):
		maxamp_aft3.append(maxamp3[i])
		final_dip_aft3.append(final_dip3[i])

maxamp3=(np.array((maxamp_aft3))*400)
final_dip3=(np.array((final_dip_aft3)))
print (len(final_dip3))

n=fitting_KK(maxamp3,final_dip3,c0)
#print ('c1=',n[0],'c2=',n[1],'x0=',n[2],'c0=',n[3])

'''
x=np.linspace(4,6.5,100)
standard_dev = np.std(final_dip, axis=0)
fitted=pol(x,n[0],n[1],n[2],n[3],n[4])
plt.fill_between(x, fitted-standard_dev, fitted+standard_dev,alpha=0.2,color='#000000')
'''
fitted3=pol(maxamp3,n[0],n[1],n[2])
xs3, ys3 = zip(*sorted(zip(maxamp3, fitted3)))
plt.plot(xs3, ys3,'k',markersize=2,label='LQ + Inflow')
plt.scatter(maxamp3,final_dip3,s=1.9,alpha=0.5,c='#000000')


datafilename='c1_c2' + '.csv'
printout=str(u0_sel)+'\t' + str(eta_sel)+'\t' + str(tau_sel)+'\t' + str(n[0])+'\t' + str(n[1])+'\t' + str(s[0])+'\t' + str(s[1])+'\t' + str(t[0])+'\t' + str(t[1])+'\t'+'\n'
with open(datafilename, 'a') as outtable:
	outtable.write(printout)
outtable.closed

print(len(ys),len(ys1),len(ys2),len(ys3))
line1=np.array(ys)
line2=np.array(ys1)
line3=np.array(ys2)
line4=np.array(ys3)

d=line1-line4
d=np.absolute(d)

#plt.plot(maxamp[np.argmin(d)],fitted_no[np.argmin(d)],'ro')
x_inter=xs[np.argmin(d)]
doub_x0=2*x_inter
sig_pls=doub_x0+standard_dev*5
sig_min=doub_x0-standard_dev*5
#plt.axvline(x=sig_pls)
plt.axvline(x=x_inter, linestyle= '--')
plt.axvline(x=doub_x0)
print(doub_x0)
#plt.axvspan(sig_min, sig_pls, alpha=0.3, color='red')


dev=[]

def find_horizantal(maxamp,fitted):
	avrgfit=[]
	for count,value in enumerate(maxamp):
		if (value>sig_min and value<sig_pls):
			avrgfit.append(fitted[count])
		#elif (value>sig_min-7 and value<sig_pls+7):
		#	avrgfit.append(fitted[count])
	avr=sum(avrgfit)/len(avrgfit)
	dev.append(avr)
	plt.axhline(y=avr, linestyle='--')
	avr = round(avr, 4)
	plt.text(40,avr,'%s'%(avr))
	avrgfit*=0


find_horizantal(xs,ys)

find_horizantal(xs1,ys1)

find_horizantal(xs2,ys2)

find_horizantal(xs3,ys3)

TQ_dev=dev[0]-dev[1]
LQ_dev=dev[0]-dev[2]
LQTQ_dev=dev[0]-dev[3]
dev*=0

datafilename='dev' + '.csv'
printout=str(u0_sel)+'\t' + str(eta_sel)+'\t' + str(tau_sel)+'\t'+str(TQ_dev)+'\t' + str(LQ_dev)+'\t' + str(LQTQ_dev)+'\n'
with open(datafilename, 'a') as outtable:
	outtable.write(printout)
outtable.closed

plotname='Scatter'+'/case' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.png'

plt.legend()
plt.savefig(plotname, dpi=300)
#plt.show()
