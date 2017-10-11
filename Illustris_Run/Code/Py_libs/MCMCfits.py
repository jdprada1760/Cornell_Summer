
# coding: utf-8

# In[1]:
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from IPython.display import display, Math, Latex
import emcee
import matplotlib.gridspec as gridspec
#import corner
from scipy import stats

# In[2]:

#mass units
M = (10.0**10)
#vol of illustris run
vol = (75.0/(0.7))**3
print(75./0.7)
# Folder to save
fold = "Filter19/F19_"

# Gives the x,y of the histogram for logmasses from the given array of mass
def histog( mf, nbins ):
    hist, bin_edges = np.histogram(mf, bins = nbins)
    # Obtains the center point for each bin
    xcenter = (bin_edges[:-1] + bin_edges[1:])/2
    return np.array(xcenter), np.array(hist/vol)

# Truncates float to n decimal parts
def truncate(d,n):
    before_dec, after_dec = str(d).split('.')
    d = float('.'.join((before_dec, after_dec[0:n-1])))
    return d

# Schechter function
def Schechter(m,phi_s,alpha,m_s):
    MHI = 10.**m
    M_s = 10.**m_s
    return np.log(10.)*phi_s*((MHI/M_s)**(alpha+1))*np.exp(-MHI/M_s)

# loads the file with halo masses and makes the quartiles
# file formated env,mass
# j defines if we treat gass or dm j = 0,1
# j = 0 -> gas
# j = 1 -> dm
# j = 2 -> stars and dust
# j = 3 -> black holes
# j = 4 -> central black hole
# minMass minimum mass allowed to have into account a halo
def loadHaloes(name, massname, neigh, j, minMass):
    table = np.loadtxt(name, skiprows = 0, delimiter = ',').T
    mtable = np.loadtxt(massname, skiprows = 0, delimiter  = ',').T
    #table[0] = neigh/(np.pi*table[0]+1e-9)
    print("number of haloes before: " + str(len(table.T)))
    
    # Only take the haloes with not null mass
    if j == 0:
    	lsigma = table[0][np.where(((M*mtable[j]) > minMass)& ((M*mtable[j]) < 10**(11.5)))]
    	mass = M*mtable[j][np.where((M*mtable[j] > minMass) & ((M*mtable[j]) < 10**(11.5)))]
    else:
        lsigma = table[0][np.where(((M*mtable[j]) > minMass))]
        mass = M*mtable[j][np.where(M*mtable[j] > minMass)]  

    print("number of haloes after: " + str(len(mass)))

    # Calculates the first quartile of logsigma
    global q1
    global q2
    global q3
    global q4
    global q5
    q5 = M*mtable[j][np.where(M*mtable[j]> minMass)]

    
    q = np.percentile(lsigma,25)
    # Keeps first quartile
    q1 = mass[np.where(lsigma <= q)]
    # Removes fist quartile from lsigma
    mass = mass[np.where(lsigma > q)]
    lsigma = lsigma[np.where(lsigma > q)]
    # Second quartile
    q = np.percentile(lsigma,33.333333)
    q2 = mass[np.where(lsigma <= q)]
    mass = mass[np.where(lsigma > q)]
    lsigma = lsigma[np.where(lsigma > q)]
    # Third and Fourth quartile
    q = np.percentile(lsigma,50)
    q3 = mass[np.where(lsigma <= q)]
    q4 = mass[np.where(lsigma >= q)]
    
    # verify the length of each quartile
    print(len(q1),len(q2),len(q3),len(q4), sum([len(q1),len(q2),len(q3),len(q4)]))

# Loads the file with environment classification T-Web
# file formated env,massgas,masdm
# j defines if we treat gass or dm j = 0,1
# minMass minimum mass allowed to have into account a halo
def loadTweb(name, massname, neigh, j, minMass):
    table = np.loadtxt(name, skiprows = 0, delimiter = ',').T
    mtable = np.loadtxt(massname, skiprows = 0, delimiter  = ',').T
    #table[0] = neigh/(np.pi*table[0]+1e-9)
    print("number of haloes before: " + str(len(table.T)))
    
    # Only take the haloes with not null mass
    if j == 0:
    	lsigma = table[0][np.where(((M*mtable[j]) > minMass)& ((M*mtable[j]) < 10**(11.5)))]
    	mass = M*mtable[j][np.where((M*mtable[j] > minMass) & ((M*mtable[j]) < 10**(11.5)))]
    else:
        lsigma = table[0][np.where(((M*mtable[j]) > minMass))]
        mass = M*mtable[j][np.where(M*mtable[j] > minMass)]
    
    print("number of haloes after: " + str(len(mass)))

    # Calculates the first quartile of logsigma
    global q1 
    global q2 
    global q3 
    global q4 
    global q5 
    q1 = mass[np.where(lsigma == 3)]
    q2 = mass[np.where(lsigma == 2)]
    q3 = mass[np.where(lsigma == 1)]
    q4 = mass[np.where(lsigma == 0)]
    print("Number of Voids: " ,len(q4))
    print("######################################################################")
    q4 = -1
    q5 = mass
    # verify the length of each quartile
    print(len(q1),len(q2),len(q3), sum([len(q1),len(q2),len(q3)]))
    
def SchtrGraph2(bins,filename,title, boolean, indi):
    # indi is 0 if graphing quartiles, 1 if graphing Tweb
    # Graphics of the mass function for each quartile
    fig, axes = plt.subplots( nrows=4, ncols=2, figsize=(10,15) )
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.3)
    gs = gridspec.GridSpec(4,2)
    ax = plt.subplot(gs[0:2,0:])
    ax1 = plt.subplot(gs[2:,0:])
    #ax2 = plt.subplot(gs[1,1])
    #ax0 = plt.subplot(gs[0, 0])
    
    # Plots Total mass function
    '''
    fig, axes = plt.subplots( nrows=2, ncols=1, figsize=(15,20) )
    gs = gridspec.GridSpec(2,1)
    ax0 = plt.subplot(gs[0, 0])
    ax = plt.subplot(gs[1,0])
    tmp = histog(log10(q5),bins+5)
    tmp = [tmp[0][where(tmp[1]!=0)],tmp[1][where(tmp[1]!=0)]]
    ax0.errorbar( tmp[0], log10(tmp[1]),  1./(sqrt(tmp[1]*vol)), fmt='o', ecolor = 'b',
                        c = 'g', markersize=7,elinewidth=1.5)
    '''
    
    # Callable lists
    if indi:
        qs = [q1,q2,q3,q4]
        indx = ['1st quartile','2nd quartile','3rd quartile','4th quartile']
        cs = ['b','g','r','y']
        fmts = ['o','>','<','s']
    else: 
        qs = [q1,q2,q3]
        indx = ['Cluster','Sheet','Fillament']
        cs = ['b','g','r']
        fmts = ['o','>','<']
    alph = []
    dalph = []
    ms = []
    dms = []
    minx = 40
    maxx = -1
    miny = 100
    maxy = -100
    # Chooses normalization factor as first phi_s
    tmp = histog(np.log10(q2),15)
    #print(tmp)
    tmp = [tmp[0][np.where(tmp[1]!=0)],tmp[1][np.where(tmp[1]!=0)]]
    phi = 1
    if ( boolean ):
        fiiit,sampler = fitMCMC(tmp)
        #print(fiiit)
        ft = fiiit.T[0]
        errors = fiiit.T[1:]
        phi = ft[0]
    ind = 0
    for q,i,co,fm in zip(qs,indx,cs,fmts):
        # Histo for first quartile
        tmp = histog(np.log10(q),bins)
        tmp = [tmp[0][np.where(tmp[1]!=0)],tmp[1][np.where(tmp[1]!=0)]]
        #print tmp
        if ( boolean ):
            logm = np.linspace(min(tmp[0])-0.3,max(tmp[0])+0.1, 200)
            fiiit,sampler = fitMCMC(tmp)
            #print(fiiit)
            ft = fiiit.T[0]
            errors = fiiit.T[1:].T
            print("......................................................")
            print(ft[0],errors[0])
            print(ft[1],errors[1])
            print(ft[2],errors[2])
            ax.errorbar( tmp[0], np.log10(tmp[1]*(phi/ft[0]))-ind,  1./(np.sqrt(tmp[1]*vol)), fmt=fm, ecolor = co,
                        c = co, markersize=7,elinewidth=1.5)
        else: 
            ax.errorbar( tmp[0], np.log10(tmp[1]),  1./(np.sqrt(tmp[1]*vol))-ind, fmt=fm, ecolor = co,
                        c = co, markersize=7,elinewidth=1.5)
        # Decimals to show
        if ( boolean ):
            m = 5
            ax.plot(logm, np.log10(Schechter(logm,ft[0],ft[1],ft[2])*(phi/ft[0])) -ind,co
                 #,label = r"$\alpha$ = " + str(truncate(ft[1],m))+"\n"
                 #+ r"$m_\ast$ = " + str(truncate(ft[2],m))+"\n"  )
                 )
            alph.append(ft[1])
            dalph.append(errors[1])
            ms.append(ft[2])
            dms.append(errors[2])
            #ax.bar( tmp[0] , log10(tmp[1]), width = tmp[0][1]-tmp[0][2], alpha = 0.4, align = 'center' )
            #ax.legend(prop={'size':13})
        #savetxt(name + "q2.csv",(array([tmp[0], tmp[1], sqrt(vol*tmp[1])/(vol*tmp[1])]).T), delimiter = ',')
        #ax.set_title(i, fontsize = 18)
        ax.set_ylabel('$Log_{10}(n)$', fontsize = 20)
        ax.set_xlabel('$Log_{10}(M/M_{\odot})$', fontsize = 20)
        minx = min([minx,min(tmp[0])])
        maxx = max([maxx,max(tmp[0])])
        if ( boolean ):
            miny = min([miny,min(np.log10(tmp[1]*(phi/ft[0])))-ind])
            maxy = max([maxy,max(np.log10(tmp[1]*(phi/ft[0])))-ind])
        else:
            miny = min([miny,min(np.log10(tmp[1]))-ind])
            maxy = max([maxy,max(np.log10(tmp[1]))-ind])
        ################################################################
        ind += 0.5
    ax.set_title(title,fontsize = 20)
    tmp = histog(np.log10(q5),bins)
    ax.set_xlim(minx-0.3,maxx+0.3)
    if ( boolean ):
        ax.set_xlim(minx-0.3,maxx+1.5)
    ax.set_ylim(miny-0.3,maxy+0.3)
    #ax.set_ylim(-6.,0.)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ###############################################################################3
    # Plots knee mass and alpha
    ###################################################################################
    labs = []
    if indi:
        labs = ['4th', '3rd', '2nd','1st']
        cs = ['b','g','r','y']
        fmts = ['o','>','<','s']
    else: 
        labs = ['Cluster','Sheet','Fillament']
        cs = ['b','g','r']
        fmts = ['o','>','<']
    for al,dal,ma,dma,fmtt,col,lab in zip(alph,dalph,ms,dms,fmts,cs,labs):
        print(al,dal,ma,dma,fmtt,col,lab)
        ax1.errorbar( np.array([al]), np.array([ma]),xerr=np.array([dal]).T, yerr = np.array([dma]).T
                     , fmt=fmtt, ecolor = 'r', 
                     c = col, markersize=10,elinewidth=1.5,label =lab)
        ax1.legend(prop={'size':13})
    #ax.set_xlim(0.5,4.5)
    d = 0.1
    d2 = 0.1
    ax1.set_xlim(min(alph)-d,max(alph)+d)
    ax1.set_ylim(min(ms)-d2,max(ms)+d2)
    ax1.set_ylabel(r"$Log(m_\ast)$", fontsize = 30)
    ax1.set_xlabel(r"$\alpha$", fontsize = 30)
    ax1.tick_params(axis='x', labelsize=20)
    ax1.tick_params(axis='y', labelsize=20)
    #ax.set_title(title, fontsize = 30)
    #plt.savefig("../../data/pics/Filter/"+filename)
    plt.savefig("../../data/pics/" + fold+filename)
    return np.array(alph),np.array(dalph).T,np.array(ms),np.array(dms).T

    
def graph(x,dx,d,name,title,filename):
    fig, ax = plt.subplots( nrows=1, ncols=1, figsize=(10,7) )
    q = [1,2,3,4]
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.3)
    ax.errorbar( q, x, yerr = dx, fmt='o', ecolor = 'r', c = 'b', markersize=7,elinewidth=1.5)
    ax.set_xlim(0.5,4.5)
    ax.set_ylim(min(x)-d,max(x)+d)
    ax.set_ylabel(name, fontsize = 30)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_xlabel(r'$Environment$', fontsize = 30)
    ax.set_title(title, fontsize = 30)
    #plt.savefig("../../data/pics/Filter/"+filename)
    plt.savefig("../../data/pics/"+fold+filename)


# In[3]:

# Defines functions for emcee
def lnprior(p):
    # The parameters are stored as a vector of values, so unpack them
    norm,alpha,m = p
    # We're using only uniform priors, and only eps has a lower bound
    if norm > 0 and   8 < m < 20 and -3 < alpha < 1:
        return 0
    return -np.inf

# Defines the likelihood
def lnlike(p, x, y, yerr):
    model = SchechterMC(x,p)
    # the likelihood is sum of the lot of normal distributions
    denom = np.power(yerr,2.)
    lp = -0.5*sum(np.power((y - model),2.)/denom + np.log(denom)) 
    return lp

# The probability
def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, x, y, yerr)

# MCMC fit
p0 = [0.003, -1., 10.]
def fitSci(tmp,yerr,p):
    return optimize.curve_fit(Schechter,tmp[0] ,tmp[1],p0=[0,-2.0,12], sigma = np.sqrt(tmp[1]/vol),
    #bounds = (np.array([0,-np.inf,8.]),np.array([np.inf,np.inf,20.])),method = 'trf',
    maxfev = 100000)


# Schechter function
def SchechterMC(m,params):
    MHI = 10.**m
    M_s = 10.**params[2]

    return np.log(10.)*params[0]*((MHI/M_s)**(params[1]+1))*np.exp(-MHI/M_s)

# Gets the error taking into account the confidence interval
# Arr -> the array
# center -> central value to calculate the interval
# perc -> percent of data within the interval
def getErr(arr, center, perc):
    #conf_int  = stats.norm.interval(perc, loc = center, scale = 1 )
    conf_int  = stats.norm.interval(perc, loc = center, scale = np.std(arr))
    return conf_int

# Uses MCMC emcee to improve the SchechterMC fit and error estimation
# return the best fit values and their respective errors within confidence interval given by perc = (0.68,0.95)
# temp -> data to fit
# Nsteps -> steps for MCMC
# Nwalkers -> number of MCMC walkers
# suggestion: fitMCMC(temp,0.68,1000,50)
def fitMCMC(tmp):
    Nsteps = 1000
    Nwalker = 50
    tmp = [tmp[0][np.where(tmp[1]!=0)],tmp[1][np.where(tmp[1]!=0)]]
    # The error in y assuming poisson dist
    yerr = np.sqrt(tmp[1]/vol)
    # Uses scipy fit as initial guess
    p00 = [0.003,-1.,10.]
    ft,errors = fitSci(tmp,yerr,p00)
    #print ft
    # Creates walkers around initial guess
    Ndim = 3
    walkers = [ft+[1.e-7*np.random.randn(),1.e-4*np.random.randn(),1.e-3*np.random.randn()] for i in range(Nwalker)]
    # Initialises sampler
    sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob,args=(tmp[0],tmp[1],yerr))
    pos,prob,state = sampler.run_mcmc(walkers, 500)
    sampler.reset()
    #res=plot(sampler.chain[:,:,1].T, '-', color='k', alpha=0.3)
    #axhline(ft[1], color='blue')
    # Runs the walks
    pos,prob,state = sampler.run_mcmc(pos, Nsteps)
    # fit and Errors in format (fit,upper,lower)
    samples = sampler.chain[:, 50:, :].reshape((-1, Ndim))
    #corner.corner(samples, labels=['normV','alpha','m*'], 
     #           truths=[ft[0], ft[1], ft[2]])
    #print shape(samples)
    fit = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
    return np.array(list(fit)),samples



# Number of bins
nbins = 10


# # GAS NEAREST NEIGHBOR

# In[4]:

# Parameters of loadHaloes
# 1) Environment data (distance to nn, gas mass, dm)
# 2) Masses (mgas,mdm,x,x,mstellar,mbh) as stated in Illustris doc
# 3) # nearest neighbor to calculate environment number desity
# 4) index of mass related to Masses file
# 5) Mass cut in solar masses
#loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,0,10**(8.7))
loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,0,10**(10.3))
#loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,0,0)
#histo(15)

a = SchtrGraph2(nbins,"quartilesGas","Nearest Neighbor Gas Mass Functions",True,1)


# In[5]:

alph,dalph,ms,dms = a
print(dms.T)

# # GAS TWEB

# In[8]:

#loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,0,10**(8.7))
loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,0,10**(10.3))
#loadTweb("../data/twEnv/Tweb1.csv","../data/Masses1.csv",3.0,0,0)

a = SchtrGraph2(nbins,"T-Web_Gas","T-Web Gas Mass Functions",True,0)

# In[7]:

#loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,1,10**(10))
loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,1,10**(11.2))
#loadHaloes("../data/nnEnv/env1.data","../data/Masses1.csv",3.0,1,0)

#histo(15)
a = SchtrGraph2(nbins,"quartilesDM","Nearest Neighbor DM Mass Functions",True,1)

# In[16]:

#loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,1,10**(10))
loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,1,10**(11.2))
#loadTweb("../data/twEnv/Tweb1.csv","../data/Masses1.csv",3.0,1,0)

a = SchtrGraph2(nbins,"T-Web_DM","T-Web DM Mass Functions",True,0)

# In[10]:
#loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,4,10**(7.5))
loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,2,10**(9.4))
#histo(15)
a = SchtrGraph2(nbins,"quartilesSellar","Nearest Neighbor Stellar Mass Functions",True,1)

# In[24]:

#loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,4,10**(7.5))
loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,2,10**(9.4))

a = SchtrGraph2(nbins,"T-Web_Stellar","T-Web Stellar Mass Functions",True,0)


# # NEAREST NEIGHBOR BH MF

# In[9]:

#loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,5,10**6.8)
loadHaloes("../../data/nnEnv/env1.data","../../data/Galaxy1.csv",3.0,3,10**(6.8))
#6.4
#histo(15)
a = SchtrGraph2(15,"quartilesBH","Nearest Neighbor BH Mass Functions",True,1)

# # TWEB BH MF

# In[28]:

#loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,5,10**6.8)
loadTweb("../../data/twEnv/Tweb1.csv","../../data/Galaxy1.csv",3.0,3,10**(6.8))
a = SchtrGraph2(15,"T-Web_BH","T-Web BH Mass Functions",True,0)




