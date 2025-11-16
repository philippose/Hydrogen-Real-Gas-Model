import numpy as np; 
import matplotlib.pyplot as plt; plt.close("all")
import matplotlib.axes as plax;


################################ -- Ideal toroidal CFVNs -- ###############################################################
def Tor_ideal(l, al, nz):
    zc = np.sqrt(2.0**2 - (2.0-(2.5-1.0)/2)**2);
    zal = 2.0*np.sin(al*np.pi/180.0);
    ze = l;
    ntest = round((nz-5)/4);
    nz = 4*ntest + 5;
    BSz = np.linspace(0.0, ze, nz);
    BSy = np.linspace(0.0, 0.0, nz);
    for i in range(0,nz):
        if BSz[i] < zc + zal:
            BSy[i] = 5/2 - np.sqrt(2.0**2 - (BSz[i]-zc)**2);
        else:
            BSy[i] = 1/2 + (1.0-np.cos(al*np.pi/180.0))*2.0 + np.tan(al*np.pi/180.0)*(BSz[i]-(zc+zal)) ;
    BSd = int(np.round(np.median(np.where(BSy == BSy.min()))));
    return BSz, BSy, BSd
###########################################################################################################################


################################ -- Ideal cylindrical CFVNs -- ############################################################
def Cyl_ideal(l, al, nz):
    ze = l;
    ntest = round((nz-5)/4);
    nz = 4*ntest + 5;
    BSz = np.linspace(0.0, ze, nz);
    BSy = np.linspace(0.0, 0.0, nz);
    for i in range(0,nz):
        if BSz[i] < 1.0:
            BSy[i] = 3/2 - np.sqrt(1.0**2 - (BSz[i]-1.0)**2);
        elif BSz[i] > 2.0:
            BSy[i] = 1/2 + np.tan(al*np.pi/180.0)*(BSz[i]-2.0) ;
        else:
            BSy[i] = 1/2;
    BSd = int(np.round(np.median(np.where(BSy == BSy.min()))));        
    return BSz, BSy, BSd
###########################################################################################################################


################################ -- Measured CFVNs -- #####################################################################
def Meas_CFVN(meas_data, NType, l, n_in, n_out):
    ## Measurement input
    BSSMeas = np.loadtxt(meas_data, delimiter='\t', skiprows = 1);
    Shift_cyl = 1.5;
    Shift_tor = np.sqrt(39)/4;
    # NType = 0 --> Cyl; NType = 1 --> Tor;
    BSzMeas = BSSMeas[:,0] + (1 - NType) * Shift_cyl + NType * Shift_tor;
    BSyMeas = BSSMeas[:,1];
    BSdMeas = int(np.round(np.median(np.where(BSyMeas == BSyMeas.min()))));
    ze = l;
    nz = len(BSzMeas)*l/(BSzMeas[len(BSzMeas)-1]-BSzMeas[0]); 
    ntest = round((nz-5)/4);
    nz = 4*ntest + 5;

    ## Inlet extension
    nz1 = round(nz*(BSzMeas[0] - 0.0)/l) + 1;
    BSz1 = np.linspace(0.0, BSzMeas[0], nz1);
    BSy1 = np.linspace(0.0, 0.0, nz1);    
    
    zn_in = BSzMeas[0:n_in];
    yn_in = BSyMeas[0:n_in];
    
    ## Solving Linear System (A^T*Ax = A^Tb) on n points
    A = np.ones((3, 3));
    b = np.ones(3);
    
    A[0,0] = sum(zn_in**2); A[0,1] = sum(zn_in*yn_in); A[0,2] = sum(zn_in);
    A[1,0] = sum(zn_in*yn_in); A[1,1] = sum(yn_in**2); A[1,2] = sum(yn_in);
    A[2,0] = sum(zn_in); A[2,1] = sum(yn_in); A[2,2] = n_in;
    
    b[0] = sum(zn_in*(zn_in**2 + yn_in**2));
    b[1] = sum(yn_in*(zn_in**2 + yn_in**2));
    b[2] = sum(zn_in**2 + yn_in**2);
    
    A_T = A.transpose();
    A_new = np.matmul(A_T,A);
    b_new = np.matmul(A_T,b);
    n_opt = np.linalg.solve(A_new,b_new);

    pz = n_opt[0]/2;
    py = n_opt[1]/2;
    R = np.sqrt(n_opt[2] + (n_opt[0]**2 + n_opt[1]**2)/4);
    
    BSy1 = py - np.sqrt(R**2 - (BSz1 - pz)**2);
    n_nan = np.count_nonzero(np.isnan(BSy1));
    if (n_nan > 0):
        BSzmax = BSz1[-1];
        BSy1 = BSy1[n_nan:];
        BSz1 = BSz1[0:-n_nan];
        delta = BSzmax - BSz1[-1];
        BSzMeas = BSzMeas - delta;
        nz1 = nz1 - n_nan;

    ## Outlet extension
    nz2 = nz - len(BSzMeas) - nz1 + 2;
    BSz2 = np.linspace(BSzMeas[-1], ze, nz2);
    BSy2 = np.linspace(0.0, 0.0, nz2);
    
    m = (BSyMeas[-(n_out+1)] - BSyMeas[-(n_out+2)])/(BSzMeas[-(n_out+1)] - BSzMeas[-(n_out+2)]);
    n = BSyMeas[-1] - m * BSzMeas[-1];
    BSy2 = m * BSz2 + n;

    ## Creation of full BS vectors
    BSz = np.concatenate((BSz1[0:-1], BSzMeas, BSz2[1:]), axis=None);
    BSy = np.concatenate((BSy1[0:-1], BSyMeas, BSy2[1:]), axis=None);
    BSd = int(np.round(np.median(np.where(BSy == BSy.min()))));
    return BSz, BSy, BSd, BSzMeas, BSyMeas
###########################################################################################################################


################################ -- Creation of geometry file -- ##########################################################
def create_geometry_file(name, BSz, BSy, BSd): 
    f = open(name, 'w');
    f.write('BSz ( ');
    np.savetxt(f, BSz, fmt='%.6f', newline = ' ');
    f.write(' ); \n \nBSy ( ');
    np.savetxt(f, BSy, fmt='%.6f', newline = ' ');
    f.write(' ); \n \nBSd ' + str(BSd) + ';');
    f.close();
########################################################################################################################   


### Input parameters ###
meas_data = 'Cyl_D1_meas.txt'; # measured contour data
NType = 0;  # nozzle type (Cyl. type = 0, Tor. type = 1)
l = 9.0;    # nozzle length
n_in = 40;  # point index of inlet circle
n_out = 25; # point index of outlet slope

### Nozzle contour (Meas. cylindrical nozzle) ###
zMeas, yMeas, dMeas, zMeasRaw, yMeasRaw  = Meas_CFVN(meas_data, NType, l, n_in, n_out); 


### Create geometry file ###
create_geometry_file('Meas_Cyl.out', zMeas, yMeas, dMeas);








