# =============================================================================
# Module for well logs
# =============================================================================
import numpy as np
import scipy as sc
from scipy import io
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
import scipy.stats as st
import matplotlib.colors as colors
import scipy.stats as stats
import seaborn as sns

# Plotting definitions
font_ = 24
f_size=(15,12)


# Correlation between x and y
def corr(x,y):
    results = sm.OLS(y,sm.add_constant(x)).fit()
    plt.scatter(x,y,s=3.0)
    plt.plot(x, results.fittedvalues, 'k')
#    plt.title(np.str(results.params[0]b ) + np.str(results.params[1]))
    plt.title(r'%3.3f + %3.3f x'%(results.params[0] , results.params[1]), fontsize=18)
    return results.params[0] , results.params[1]

# =============================================================================
# Rock physics modelling
# =============================================================================

# Lithology classification
def lith_class(log, depth, vshale, phi, sw):
    lith = np.zeros(log.shape[0])
    for i in range(log.shape[0]):
        if log[i,vshale] < 0.1 and log[i,phi] >= 0.2 and log[i,sw] <= 0.6:
            lith[i] = 0                                                                                 # 'High porosity oil limestone'; 
        elif log[i,vshale] < 0.1 and log[i,phi] < 0.2 and log[i,sw] <= 0.6:
            lith[i] = 1                                                                                 # 'Low porosity oil limestone'; 
        elif log[i,vshale] < 0.1 and log[i,phi] >= 0.2 and log[i,sw] > 0.6:
            lith[i] = 2                                                                                 # 'High porosity water limestone';
        elif log[i,vshale] < 0.1 and log[i,phi] < 0.2 and log[i,sw] > 0.6:  
            lith[i] = 3                                                                                 # 'Low porosity water limestone';
        elif log[i,vshale] >= 0.1:
            lith[i] = 4                                                                                 # 'Marly Chalk'
        
    # LFC classification names and colours for plotting
    ccc = np.array(['red','black','green','blue','brown',])
    cmap_facies = colors.ListedColormap(ccc[0:len(ccc)], 'indexed')
    lith_var = ['HP oil LS', 'LP oil LS', 'HP water LS', 'LP water LS', 'Marly Chalk']

    lith_depth = [0,10, log[-1,depth],log[0,depth]] # (X,Y ) adjusting the litholog to the same depth as the rest of the logs
        
    #size of second dimension for the lith-array (only plotting related)
    y_size = 100    
    lith = np.repeat(np.expand_dims(lith,1), y_size, 1)    
    X = np.zeros([log[:,phi].shape[0],3])
    X[:,0] = log[:,phi]
    X[:,1] = log[:,vshale] #SD.lit
    X[:,2] = log[:,sw]    
    y = lith[:,0]    
    # if-condition to replace np.single arrays of litho classification --> cannot be used in discriminant analysis
    for i in range(np.unique(y).size):
        if np.where(lith[:,0]==np.unique(lith)[i])[0].size == 1:
            y[np.where(lith[:,0] == np.unique(lith)[i])[0]]   =    y[np.where(lith[:,0]==np.unique(lith)[i])[0]-1]
    clf = QuadraticDiscriminantAnalysis()
    clf.fit(X, y)
    newlith = clf.predict(X)
    newlith = np.repeat(newlith, y_size, axis = 0)
    newlith = np.reshape(newlith, lith.shape)  
    return lith, newlith, lith_depth, ccc, cmap_facies, lith_var


# Bulk modulus
def bulk_shear(vp, vs, rho):
    vp  = vp / 1000.        #[km/s]
    vs  = vs / 1000.        #[km/s]
    rho = rho / 1000.       #[g/cm^3]
    mu  = rho * vs**2.      #shear modulus [GPa]
    k_s = rho * vp**2. - (4./3.)*mu #bulk modulus [GPa]
    return k_s, mu

# Voigt and Reuss bounds + Hill averaging
def bounds_vrh(f,k,u):
#% [K_U,K_L,U_U,U_L,KA,UA]=BOUND(IB,F,K,U)
#% calculate the elastic bounds (upper & lower) of an aggregate.
#%
#% input:    IB - type of bound (ib=0: Voigt-Reuss; ib=1: Hashin-Shtrikman);
#%           F - volume fractions (<=1); K - bulk moduli; U - shear moduli.
#% output:   K_U, K_L, U_U, U_L - elastic bounds of the aggregate.
#%           KA, UA - arithmetic average of the upper and lower bounds
#%                    (equals the Hill average for Hashin-Shtrikman bounds)
#% 
#% note:     1. VR bounds are the simplest;
#%           2. HS bounds are the narrowest possible;
#%           3. assumption: rock is isotropic.
#
#% source:   Berryman, J.G., 1993, Mixture theories for rock properties
#%           Mavko, G., 1993, Rock physics formulas
#%  

#% Written by Xingzhou 'Frank' Liu
#% Modified by Isao Takahashi 4/27/99
    
    k = np.ones(f.size)*k
    u = np.ones(f.size)*u
#    lf=length(f);lk=length(k);lu=length(u);
#if lf~=lk|lf~=lu,
#	error('Input f, k, and u must have the same length')
##end
#if sum(f)~=1,
#	error('F must sum up to 1')
#end

#if ib==0		% use Voigt-Reuss bounds

    k_u=np.sum(f*k)			#% Voigt bound
    k_l=1./np.sum(f/k)		#% Reuss bound
    
    u_u=np.sum(f*u)			#% Voigt bound
    u_l=1./np.sum(f/u)		#% Reuss bound
    
    ka=(k_u+k_l)/2.			#% Hill average
    ua=(u_u+u_l)/2.
        
    
#    
#    
#    f = np.array(volumes).T
#    k = np.resize(np.array(k),np.shape(f))
#    mu = np.resize(np.array(mu),np.shape(f))
#    k_u = np.sum(f*k, axis=1)           #upper bound
#    k_l = 1. / np.sum(f/k, axis=1)      #lower bound
#    mu_u = np.sum(f*mu, axis=1)         #upper bound
#    mu_l = 1. / np.sum(f/mu, axis=1)    #lower bound
#    k0 = (k_u+k_l) / 2.                 #average
#    mu0 = (mu_u+mu_l) / 2.
#    return k_u, k_l, mu_u, mu_l, k0, mu0

    return k_u,k_l,u_u,u_l,ka,ua



# Kernel density estimation (and scatterplot)
    
def kde(par1, par2, par3, cmap_facies, lith_var,ccc, kde_contour, size_marker):
    fig = plt.figure(figsize=(15,12)) 

    plt.scatter(par1, par2,  s = size_marker, c=par3[:,0], marker='o', edgecolors='k',alpha=0.2, cmap=cmap_facies, vmin=0, vmax=4)
    plt.xlabel(r'AI $(km/s * g/cm^3)$'); plt.ylabel('Vp/Vs')
#    plt.grid()
    extra_gridx = par1.min()/5 # Defined for extra space for the density estimations (avoiding discontinuities)
    extra_gridy = par2.min()/5 # Defined for extra space for the density estimations (avoiding discontinuities)
    cbar = plt.colorbar()
    cbar.set_label((15*' ').join(lith_var))
    cbar.set_ticks(range(0,1)); cbar.set_ticklabels('')
    xx, yy = np.mgrid[par1.min()-extra_gridx:par1.max()+extra_gridx:100j, par2.min()-extra_gridy:par2.max()+extra_gridy:100j]
    
    for i in range(5):
        ind_ = np.where(par3[:,0]==i)
        if ind_[0].size == 0:
            continue
        x = par1[ind_]
        y = par2[ind_]
#        xmin, xmax = np.min(x), np.max(x)
#        ymin, ymax = np.min(y), np.max(y)        
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values)
#        kernel.set_bandwidth(bw_method='silverman')
        f = np.reshape(kernel(positions).T, xx.shape)
        ax = fig.gca()
        if kde_contour == True: #If false, only scatterplot is shown
            cset = ax.contour(xx, yy, f, colors=ccc[i], alpha = 0.8)
    
    ax.set_xlim(par1.min()-extra_gridx, par1.max()+extra_gridx)
    ax.set_ylim(par2.min()-extra_gridy, par2.max()+extra_gridy)
    return



def HS(k1, mu1, k2, mu2):
    dx = 104 #resolution of phi, sat, and vshale
    res_phi, res_vshale, res_sat = dx, dx, dx
    phi_max = 0.45 # Critical porosity
    #rock.com.solidMix.model = 'hswl';
    #rock.com.fluidMix.model = 'voigt';
    
    # Por, Lit, Sat arrays
    phi = np.linspace(0,phi_max, num = res_phi)
    vshale = np.linspace(0,1, num = res_vshale)
    sat = np.linspace(0,1, num = res_sat)
    
    f1, f2 = sat, 1.-sat
    
    
    # Lower - and upper Hashkin-Strikman bounds for K
    k_HS_lower = k1 + (f1 / ((k2-k1)**-1. + f2 * (k1+ 4./3.  * mu2)**-1. ))
    k_HS_upper = k1 + (f1 / ((k2-k1)**-1. + f2 * (k1+ 4./3.  * mu1)**-1. ))
    
    # Lower - and upper Hashkin-Strikman bounds for Mu
    
#    mu_HS_lower = mu2 + (f2/( (mu1-mu2)**-1. + 2.*f1*(k2 + 2.*mu2) / (5.*mu2 *(k2 + 4.*mu2/3.)) ) ) 
#    mu_HS_upper  = mu1 + (f2/( (mu2-mu1)**-1. + 2.*f1*(k1 + 2.*mu1) / (5.*mu1 *(k1 + 4.*mu1/3.)) ) ) 

    fgu = mu1*(9.*k1+8.*mu1)/(6.*(k1+2.*mu1))
    fgl = mu2*(9.*k2+8.*mu2)/(6.*(k2+2.*mu2))
    mu_HS_upper = mu2+(mu1-mu2)*f2*(mu2+fgu) / (mu2+fgu+f1*(mu1-mu2))
    mu_HS_lower = mu2+(mu1-mu2)*f2*(mu2+fgl) / (mu2+fgl+f1*(mu1-mu2))
    
    return k_HS_lower, k_HS_upper, mu_HS_lower, mu_HS_upper

def voigt(k1,f1,k2):
    k = np.array([k1,k2])
    f = np.array([f1,1-f1]).T
    k_u=np.sum(f*k,axis=1)			#% Voigt bound
    return k_u


def reuss(k1,f1,k2):
    k = np.array([k1,k2])
    f = np.array([f1,1-f1]).T
    k_l=1./np.sum(f/k , axis=1)		#% Reuss bound
    return k_l





def avopp(vp1,vs1,d1,vp2,vs2,d2,ang,approx):
    #%Rpp=AVOPP(vp1,vs1,d1,vp2,vs2,d2,ang,approx);
    #%
    #%Calculates P-to-P reflectivity (Rpp) as a function of
    #%the angle of incidence (ang).
    #%input parameters:
    #%  layer 1 (top): vp1, vs1, density1 (d1)
    #%  layer 2 (bottom): vp2, vs2, density2 (d2)
    #% ang: vector with angles(DEG)
    #% approx: 1)Full Zoeppritz(A&R)
    #%	  2)Aki&Richards
    #%         3)Shuey's paper
    #%         4)Castagna's paper->Shuey (slightly different formulation of Shuey)
    #%
    #% With no output arguments, plots Rpp vs. angle.
    #%
    #% See also AVOPS, AVO_ABE
    
    t=ang*np.pi/180.;	p=np.sin(t)/vp1;	ct=np.cos(t);
    da=(d1+d2)/2.;     Dd=(d2-d1);
    vpa=(vp1+vp2)/2.;  Dvp=(vp2-vp1);
    vsa=(vs1+vs2)/2.;  Dvs=(vs2-vs1);
    
    if approx == 'zoeppritz':
#       case 1,		%FULL Zoeppritz (A&K)
    	ct2=np.sqrt(1.-(np.sin(t)**2.*(vp2**2./vp1**2.))) 
    	cj1=np.sqrt(1.-(np.sin(t)**2.*(vs1**2./vp1**2.))) 
    	cj2=np.sqrt(1.-(np.sin(t)**2.*(vs2**2./vp1**2.))) 
    	a=(d2 *(1.-(2. *vs2**2. *p**2.)))-(d1 *(1.-(2 *vs1**2. *p**2.))) 
    	b=(d2 *(1.-(2. *vs2**2. *p**2.)))+(2. *d1 *vs1**2. *p**2.) 
    	c=(d1 *(1.-(2. *vs1**2. *p**2.)))+(2. *d2 *vs2**2. *p**2.) 
    	d=2. *((d2 *vs2**2.)-(d1 *vs1**2.)) 
    	E=(b *ct/vp1)+(c *ct2/vp2) 
    	F=(b *cj1/vs1)+(c *cj2/vs2) 
    	G=a-(d *ct *cj2/(vp1 *vs2)) 
    	H=a-(d *ct2 *cj1/(vp2 *vs1)) 
    	D=(E *F)+(G *H *p**2.) 
    	Rpp = ( (((b *ct/vp1)-(c *ct2/vp2)) *F) - ((a+(d *ct *cj2/(vp1 *vs2))) *H *p**2.) ) / D 
        
    elif approx == 'aki_richard': #%Aki & Richard (aprox) , assuming (angles) i=i1
    	Rpp=(0.5 *(1.-(4. *p**2. *vsa**2.)) *Dd/da) + (Dvp/(2. *ct**2. *vpa)) - (4. *p**2. *vsa *Dvs) 
    elif approx == 'shuey':
    	poi1=((0.5 *(vp1/vs1)**2.)-1.)/((vp1/vs1)**2.-1.) 
    	poi2=((0.5 *(vp2/vs2)**2.)-1.)/((vp2/vs2)**2.-1.) 
    	poia=(poi1+poi2)/2.; Dpoi=(poi2-poi1)
    	Ro=0.5 *((Dvp/vpa)+(Dd/da)) 
    	Bx=(Dvp/vpa)/((Dvp/vpa)+(Dd/da)) 
    	Ax=Bx-(2. *(1.+Bx) *(1.-2. *poia)/(1.-poia)) 
    	Rpp= Ro + (((Ax *Ro)+(Dpoi/(1.-poia)**2.)) *np.sin(t)**2.) + (0.5 *Dvp *(np.tan(t)**2.-np.sin(t)**2.)/vpa) 
    elif approx == 'shuey_lin':		#%Shuey linear
    	A=0.5 *((Dvp/vpa)+(Dd/da)) 
    	B=(-2. *vsa**2. *Dd/(vpa**2. *da)) + (0.5 *Dvp/vpa) - (4. *vsa *Dvs/(vpa**2.)) 
    	Rpp=A+(B*np.sin(t)**2.) 	
    

    return Rpp

def sec_mean(vp,vs,rho):
    vp_m = np.mean(vp)
    vs_m = np.mean(vs)
    rho_m = np.mean(rho)
    return vp_m, vs_m, rho_m







# =============================================================================
# PANDAS SCRIPTS
# =============================================================================

def multivariateGrid(col_x, col_y, col_k, df, k_is_color=False, scatter_alpha=.5):
    def colored_scatter(x, y, c=None):
        def scatter(*args, **kwargs):
            args = (x, y)
            if c is not None:
                kwargs['c'] = c
            kwargs['alpha'] = scatter_alpha
            plt.scatter(*args, **kwargs)

        return scatter
    g = sns.JointGrid(
        x=col_x,
        y=col_y,
        data=df
    )
    color = None
    legends=[]
    for name, df_group in df.groupby(col_k):
        legends.append(name)
        if k_is_color:
            color=name
        g.plot_joint(
            colored_scatter(df_group[col_x],df_group[col_y],color),
        )
        sns.distplot(
            df_group[col_x].values,
            ax=g.ax_marg_x,
            color=color,
        )
        sns.distplot(
            df_group[col_y].values,
            ax=g.ax_marg_y,
            color=color,            
            vertical=True
        )
    # Do also global Hist:
    sns.distplot(
        df[col_x].values,
        ax=g.ax_marg_x,
        color='grey'
    )
    sns.distplot(
        df[col_y].values.ravel(),
        ax=g.ax_marg_y,
        color='grey',
        vertical=True
    )
    plt.legend(legends)
    plt.gcf().set_size_inches(11.7, 8.27) ## Change figsize
    
    return




    
# =============================================================================
#     Jupyter Script Definitions:
# =============================================================================
    

def plot_vawig(axhdl, data, t, excursion):

    import numpy as np
    import matplotlib.pyplot as plt

    [ntrc, nsamp] = data.shape
    

    
    
    t = np.hstack([0, t, t.max()])
    
    for i in range(0, ntrc):
        tbuf = excursion * data[i,:] / np.max(np.abs(data)) + i
        
        tbuf = np.hstack([i, tbuf, i])
            
        axhdl.plot(tbuf, t, color='black', linewidth=0.5)
        plt.fill_betweenx(t, tbuf, i, where=tbuf>i, facecolor='k', linewidth=0, alpha=0.5)
        #plt.fill_betweenx(t, tbuf, i, where=tbuf<i, facecolor=[0.6,0.6,1.0], linewidth=0)
    
    axhdl.set_xlim((-excursion, ntrc+excursion))
    axhdl.xaxis.tick_top()
    axhdl.xaxis.set_label_position('top')
    axhdl.invert_yaxis()
    
  
    
def ricker(cfreq, phase, dt, wvlt_length):
    '''
    Calculate a zero-phase ricker wavelet
    
    Usage:
    ------
    t, wvlt = wvlt_ricker(cfreq, dt, wvlt_length)
    
    cfreq: central frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    '''
    
    import numpy as np
    import scipy.signal as signal
    
    nsamp = int(wvlt_length/dt + 1)
    t_max = wvlt_length*0.5
    t_min = -t_max
    
    t = np.arange(t_min, t_max, dt)
    
    t = np.linspace(-wvlt_length/2, (wvlt_length-dt)/2, wvlt_length/dt)
    wvlt = (1.0 - 2.0*(np.pi**2)*(cfreq**2)*(t**2)) * np.exp(-(np.pi**2)*(cfreq**2)*(t**2))
    
    if phase != 0:
        phase = phase*np.pi/180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase)*wvlt - np.sin(phase)*wvlth
    
    return t, wvlt



def wvlt_bpass(f1, f2, f3, f4, phase, dt, wvlt_length):
    '''
    Calculate a trapezoidal bandpass wavelet
    
    Usage:
    ------
    t, wvlt = wvlt_ricker(f1, f2, f3, f4, phase, dt, wvlt_length)
    
    f1: Low truncation frequency of wavelet in Hz
    f2: Low cut frequency of wavelet in Hz
    f3: High cut frequency of wavelet in Hz
    f4: High truncation frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    '''
    
    from numpy.fft import fft, ifft, fftfreq, fftshift, ifftshift
    
    nsamp = int(wvlt_length/dt + 1)
    
    
    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    aspec = freq*0.0
    pspec = freq*0.0
    
    # Calculate slope and y-int for low frequency ramp
    M1 = 1/(f2-f1)
    b1 = -M1*f1
    
    # Calculate slop and y-int for high frequency ramp
    M2 = -1/(f4-f3)
    b2 = -M2*f4
    
    # Build initial frequency and filter arrays
    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    filt = np.zeros(nsamp)
    
    # Build LF ramp
    idx = np.nonzero((np.abs(freq)>=f1) & (np.abs(freq)<f2))
    filt[idx] = M1*np.abs(freq)[idx]+b1
    
    # Build central filter flat
    idx = np.nonzero((np.abs(freq)>=f2) & (np.abs(freq)<=f3))
    filt[idx] = 1.0
    
    # Build HF ramp
    idx = np.nonzero((np.abs(freq)>f3) & (np.abs(freq)<=f4))
    filt[idx] = M2*np.abs(freq)[idx]+b2
    
    # Unshift the frequencies and convert filter to fourier coefficients
    filt2 = ifftshift(filt)
    Af = filt2*np.exp(np.zeros(filt2.shape)*1j)
    
    # Convert filter to time-domain wavelet
    wvlt = fftshift(ifft(Af))
    wvlt = np.real(wvlt)
    wvlt = wvlt/np.max(np.abs(wvlt)) # normalize wavelet by peak amplitude

    # Generate array of wavelet times
    t = np.linspace(-wvlt_length*0.5, wvlt_length*0.5, nsamp)
    
    
    # Apply phase rotation if desired
    if phase != 0:
        phase = phase*np.pi/180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase)*wvlt - np.sin(phase)*wvlth
    
    return t, wvlt
    
    

def calc_times(z_int, vp_mod):
    '''
    Calculate two-way travel time through a layered model
    
    Usage:
    -----
    t_int = calc_times(z_int, vp_mod)
    
    '''
    
    nlayers = len(vp_mod)
    nint = nlayers - 1

    t_int = []
    for i in range(0, nint):
        if i == 0:
            tbuf = z_int[i]/vp_mod[i]
            t_int.append(tbuf)
        else:
            zdiff = z_int[i]-z_int[i-1]
            zdiff = zdiff*2.0   # multiply by 2 for two-way traveltimes
            tbuf = zdiff/vp_mod[i] + t_int[i-1]
            tbuf = tbuf
            t_int.append(tbuf)
    
    return t_int


def digitize_model(rc_int, t_int, t):
    '''
    Sample a simple layered reflectivity model
    
    Usage:
    ------
    rc = digitize_model(rc, t_int, t)
    
    rc = reflection coefficients corresponding to interface times
    t_int = interface times
    t = regularly sampled time series defining model sampling
    '''
    
    import numpy as np
    
    nlayers = len(rc_int)
    nint = nlayers - 1
    nsamp = len(t)
    
    rc = list(np.zeros(nsamp,dtype='float'))
    lyr = 0
    
    for i in range(0, nsamp):

        if t[i] >= t_int[lyr]:
            rc[i] = rc_int[lyr]
            lyr = lyr + 1    

        if lyr > nint:
            break
            
    return rc
    

def rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
    Reflection & Transmission coefficients calculated using full Zoeppritz
    equations.
    
    Usage:
    ------
    R = rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta1)
    
    Reference:
    ----------
    The Rock Physics Handbook, Dvorkin et al.
    '''
    
    import math
    
    # Cast inputs to floats
    vp1  = float(vp1)
    vp2  = float(vp2)
    vs1  = float(vs1)
    vs2  = float(vs2)
    rho1 = float(rho1)
    rho2 = float(rho2)
    theta1 = float(theta1)
    
    # Calculate reflection & transmission angles
    theta1 = math.radians(theta1)   # Convert theta1 to radians
    p      = ray_param(vp1, math.degrees(theta1)) # Ray parameter
    theta2 = math.asin(p*vp2);      # Transmission angle of P-wave
    phi1   = math.asin(p*vs1);      # Reflection angle of converted S-wave
    phi2   = math.asin(p*vs2);      # Transmission angle of converted S-wave
    
    # Matrix form of Zoeppritz Equations... M & N are two of the matricies
    M = np.array([ \
        [-math.sin(theta1), -math.cos(phi1), math.sin(theta2), math.cos(phi2)],\
        [math.cos(theta1), -math.sin(phi1), math.cos(theta2), -math.sin(phi2)],\
        [2*rho1*vs1*math.sin(phi1)*math.cos(theta1), rho1*vs1*(1-2*math.sin(phi1)**2),\
            2*rho2*vs2*math.sin(phi2)*math.cos(theta2), rho2*vs2*(1-2*math.sin(phi2)**2)],\
        [-rho1*vp1*(1-2*math.sin(phi1)**2), rho1*vs1*math.sin(2*phi1), \
            rho2*vp2*(1-2*math.sin(phi2)**2), -rho2*vs2*math.sin(2*phi2)]
        ], dtype='float')
    
    N = np.array([ \
        [math.sin(theta1), math.cos(phi1), -math.sin(theta2), -math.cos(phi2)],\
        [math.cos(theta1), -math.sin(phi1), math.cos(theta2), -math.sin(phi2)],\
        [2*rho1*vs1*math.sin(phi1)*math.cos(theta1), rho1*vs1*(1-2*math.sin(phi1)**2),\
            2*rho2*vs2*math.sin(phi2)*math.cos(theta2), rho2*vs2*(1-2*math.sin(phi2)**2)],\
        [rho1*vp1*(1-2*math.sin(phi1)**2), -rho1*vs1*math.sin(2*phi1),\
            -rho2*vp2*(1-2*math.sin(phi2)**2), rho2*vs2*math.sin(2*phi2)]\
        ], dtype='float')
    
    # This is the important step, calculating coefficients for all modes and rays
    R = np.dot(np.linalg.inv(M), N);
    
    return R


def ray_param(v, theta):
    '''
    Calculates the ray parameter p
    
    Usage:
    ------
        p = ray_param(v, theta)
    
    Inputs:
    -------
            v = interval velocity
        theta = incidence angle of ray (degrees)
    
    Output:
    -------
        p = ray parameter (i.e. sin(theta)/v )
    '''
    
    import math
    
    # Cast inputs to floats
    theta = float(theta)
    v = float(v)
    
    p = math.sin(math.radians(theta))/v # ray parameter calculation
    
    return p
    
    

def hill(k1, f1, k2):
    k = np.array([k1,k2])
    f = np.array([f1,1-f1])
    
    k_u = np.sum(f*k)			#% Voigt bound
    k_l = 1./np.sum(f/k)		#% Reuss bound
    
    ka = (k_u+k_l)/2.			#% Hill average
    
    return ka

def IF():
    
    K_a=0.000131 ## bulk modulus of air

    K_ca=71  
    G_ca=30       ## chalk
    Rho_ca=2.71

    K_w=2.2
    G_w=0         ## water
    Rho_w=1.0


    Por= np.arange(0,.41,0.01)
    Rho=Rho_ca*(1-Por)+Por*Rho_w

    #WT = input('wet or dry sample? 0-dry, 1-wet:\n')
    WT = 1

    fig1 = plt.figure(figsize=(15,12))
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(222)

    ax1.set_xlabel('Porosity')
    ax1.set_ylabel('K [Gpa]')

    ax2.set_xlabel('Porosity')
    ax2.set_ylabel('G [Gpa]')


    IF = np.arange(0,1.1,0.1)

    K_out = np.zeros([IF.size,Por.size])
    G_out = np.zeros([IF.size,Por.size])

    for i in range(IF.size):
        #print(i)
        f1=IF[i]*(1-Por)
        f2=Por+(1-IF[i])*(1-Por)

        if WT == 0:
            K_sus=K_a                                     ### for dry sample
        elif WT==1:
            K_sus= ( (Por/K_w) + (((1-Por)+Por*(1-IF[i])) / K_ca) )**(-1)   ###  for wet sample
        #K_sus_start = ( (Por/K_w) + (((1-Por)*(1-IF[i])) / K_ca) )**(-1)   ###  for wet sample
    
    
        
        #K=K_ca+f2/((K_sus-K_ca)**(-1)+f1*(K_ca+4/3*G_ca)**(-1))
        #G=G_ca+ f2/(2*f1*(K_ca+2*G_ca)/(5*G_ca*(K_ca+4/3*G_ca))-1/G_ca)

        zeta=G_ca/6*((9*K_ca+8*G_ca)/(K_ca+2*G_ca))
        K=((Por+(1-IF[i])*(1-Por))/(K_sus+4/3*G_ca)+(IF[i]*(1-Por)/(K_ca+4/3*G_ca)))**(-1)-40; ## calibrated Ida -perfect, 40 replaced as 4/3*G_ca?
        G=((Por+(1-IF[i])*(1-Por))/zeta+IF[i]*(1-Por)/(G_ca+zeta))**(-1)-zeta;                 ## calibrated Ida -perfect, zeta=33.55
    
        M=K+4/3*G
        Vp=(M/Rho)**0.5 * 1000.
        Vs=(G/Rho)**0.5 * 1000.

        K_out[i,:]=K 
        G_out[i,:]=G 


        ax1.plot(Por,K_out[i])
        ax2.plot(Por,G_out[i])

    return


def rolling_window(a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        rolled = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
        
        return rolled
    
    
    # Despike function
def despike(curve, curve_sm, max_clip): 
    spikes = np.where(curve - curve_sm > max_clip)[0]
    spukes = np.where(curve_sm - curve > max_clip)[0]
    out = np.copy(curve)
    out[spikes] = curve_sm[spikes] + max_clip  # Clip at the max allowed diff
    out[spukes] = curve_sm[spukes] - max_clip  # Clip at the min allowed diff
    return out