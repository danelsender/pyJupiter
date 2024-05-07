import matplotlib.ticker
from numpy import *
import matplotlib
#matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class output: #name of the class

    def __init__(self, N):
        '''
        init

        IN: N       = number of output file to analse
        '''
        self.N = N
        self.folder = 'output%.5d'%N

        #read the descriptor file
        with open(self.folder+'/Descriptor%.1d.dat'%N,'r') as f:
             self.lines = f.readlines()

        #highest refinement level
        self.maxreflvl=int(self.lines[-3][-2:-1]) #number of grid refinement levels
        
        # initialise np arrays for grid info
        self.xdim = empty(self.maxreflvl+1,dtype=int) 
        self.ydim = empty(self.maxreflvl+1,dtype=int) 
        self.zdim = empty(self.maxreflvl+1,dtype=int) 
        self.xlim = empty((self.maxreflvl+1,2)) 
        self.ylim = empty((self.maxreflvl+1,2)) 
        self.zlim = empty((self.maxreflvl+1,2)) 

        for lvl in range(self.maxreflvl+1): #here we loop over all the grid refinement levels.
            
            #current refinement level
            self.dim_line = lvl*11+6

            dim = self.lines[self.dim_line]

            self.xdim[lvl] = int(dim.split(" ")[0])
            # azimuth upper and lower bounds
            self.xlim[lvl,0] = float(self.lines[self.dim_line+2].split(" ")[2])
            self.xlim[lvl,1] = float(self.lines[self.dim_line+2].split(" ")[-4])

            self.ydim[lvl] = int(dim.split(" ")[1])
            # radial upper and lower bounds
            self.ylim[lvl,0] = float(self.lines[self.dim_line+3].split(" ")[2])
            self.ylim[lvl,1] = float(self.lines[self.dim_line+3].split(" ")[-4])

            self.zdim[lvl] = int(dim.split(" ")[2])
            # latitude upper and lower bounds
            self.zlim[lvl,0] = float(self.lines[self.dim_line+4].split(" ")[2])
            self.zlim[lvl,1] = float(self.lines[self.dim_line+4].split(" ")[-4])

        # check whether the simulation is 2D or 3D
        if self.zdim[lvl] == 1:
            self.dim=2
        else:
            self.dim=3

        #create some epmty instance variables
        self.field = "" #name of the field
        self.azimuth = 0.0
        self.polar = False
        self.velocityfield = False
        self.verticalplot = False
        self.midplaneplot = False
        self.planetmass = 0.001
        self.physical = False
        self.cmap=plt.get_cmap('hot')

    def get_physical_units(self):
        '''Caclulates the physical units based on code
        untits
        '''
        ###########################################
        self.AU     = 1.496e+13                         # Astronomical unit in cm 
        self.XMSOL  = 1.989e+33                         # Solar mass in grams
        self.XMSTAR = 1.0                               # Stellar mass
        self.RGAS   = 8.314e+07                         # Gas constant
        self.GRAVC  = 6.67e-08                          # Gravitational constant
        self.R0     = 30.0 * self.AU                         # Planet semi major axis
        self.VOL0   = (self.R0*self.R0*self.R0)                        # Unit volume
        self.XM0    = (self.XMSOL * self.XMSTAR)                  # Unit mass
        self.RHO0   = (self.XM0 / self.VOL0)                      # Unit density
        self.TIME0  = sqrt(self.VOL0/self.GRAVC/(self.XMSTAR*self.XMSOL))   # Unit time
        self.V0     = (self.R0 / self.TIME0)                      # Unit velocity
        self.TEMP0  = (self.V0*self.V0 / self.RGAS)                    # Unit temperature
        self.P0     = self.XM0 / self.R0 / self.TIME0*self.TIME0            # Unit pressure
        ###########################################

    def convert_to_physical_units(self,data,field):
        '''Convert from code units to physical cgs units
        
        IN: data                = JUPITER data

        CALLS: self.get_physical_units()
        
        Assumptions:
        Star mass XMSTAR        = 1 M_sun
        Planet orbital radius   = 30 au
        
        OUT: phys_data'''

        self.field = field
        
        self.get_physical_units()

        # phys_data = data
        
        if (self.field == 'gasdensity'): #density field
            phys_data = data*self.RHO0
        
        if (self.field == 'gastemperature'): #temperature field
            phys_data = data*self.TEMP0

        if (self.field == 'gasvelocity'): #temperature field
            phys_data = data*self.V0
        
        return phys_data


    def readdata(self,field,reflvl,physical=False):
        '''
        reads in .dat file and returns z*y*x dimensional array with simulation data

        IN: field       = quantity to be analysed
            reflvl      = current level of mesh refinement

        OUT: fielddata

        '''
        #read .dat file
        data = fromfile(self.folder+'/%s%.1d_%.1d_%.1d.dat'%(field,self.N,reflvl,reflvl))
        # #check the data dimensions
        # print(f"dimensions of the data: {data.shape}")
        #rearange
        xdim = self.xdim[reflvl]
        ydim = self.ydim[reflvl]
        zdim = self.zdim[reflvl]
        if self.velocityfield==True:
            fielddata = reshape(data,(3,zdim,ydim,xdim))
        else:
            fielddata = reshape(data,(zdim,ydim,xdim))

        if (field=='dustdensity'): # I chose a different colormap for dust
            self.cmap=plt.get_cmap('inferno')
        else:
            self.cmap=plt.get_cmap('hot')

        self.physical = physical
        print(self.physical)
        if self.physical:
            fielddata = self.convert_to_physical_units(fielddata,field)

        return fielddata
    
    def get_minmax(self,field,reflvls=-1,velocity=False,physical=False):
        '''
        gets the max and min of field data in the midplane
        
        IN: field       = quantity of which max and min are recorded
            reflvls     = levels of refinement
            
        OUT: min        = miniumum of field data
             max        = maximum of field data'''
        
        self.velocityfield = velocity
        self.physical = physical
        # print(self.physical)
        
        if reflvls==-1: #plot all the refinement levels
            reflvls = arange(0,self.maxreflvl+1)
        
        self.upddate_maxmin(field,reflvls,vel=0,z=0,azigrid=0,avg=False,vert=False,
                            physical=self.physical)

        return self.midmin, self.midmax

    def upddate_maxmin(self,field,reflvls,vel=0,z=0,azigrid=0,avg=True,vert=False,
                       physical=False):
        '''
        computes the maximum and minimum values over all the refinement levels for the colormap

        IN: field       = quantity of which min and max will be found
            reflvls     = current mesh refinement level
            vel         = ?
            z           = height in the disc
            azigrid     = ?
            avg         = azimuthal average of field data
            vert        = take a vertical slice

        OUT: None
        '''
        self.midmax = -1.0e10
        self.vertmax = -1.0e10
        self.midmin = 1.0e10
        self.vertmin = 1.0e10
        self.physical = physical
        # print(self.velocityfield)
        for lvl in reflvls:
            if self.velocityfield==True:
                    middata = self.readdata(field,lvl,self.physical)[vel,-1-z,:,:]
                    vertdata = self.readdata(field,lvl,self.physical)[vel,:,:,int(self.xdim[lvl]/2)]
                    data = self.readdata(field,lvl,self.physical)[vel,:,:]
            else:
                    middata = self.readdata(field,lvl,self.physical)[-1-z,:,:]
                    vertdata = self.readdata(field,lvl,self.physical)[:,:,int(self.xdim[lvl]/2)]
                    data = self.readdata(field,lvl,self.physical)
            midmaximum = amax(middata)
            midminimum = amin(middata)
            vertmaximum = amax(vertdata)
            vertminimum = amin(vertdata)

            if midmaximum > self.midmax:
                self.midmax = midmaximum
            if midminimum < self.midmin:
                self.midmin = midminimum
            if vertmaximum > self.vertmax:
                self.vertmax = vertmaximum
            if vertminimum < self.vertmin:
                self.vertmin = vertminimum
# Commented out this lines as they make no sense really
            # if (self.midmin<1e-50):
            #     self.midmin = 1e-50

            if vert==True:
                if avg==True:
                    self.vertmax = amax(data.mean(-1))
                    self.vertmin = amin(data.mean(-1))
                else:
                    self.vertmax = amax(data[:,:,azigrid])
                    self.vertmin = amin(data[:,:,azigrid])
# Commented out this lines as they make no sense really
                # if (self.vertmin<1e-50):
                #     self.vertmin = 1e-50

            return None


    def plot_setlimits(self,polar=False,vertical=False):
        if self.polar:
            if vertical:
                self.ax.set_thetamin(80)
                self.ax.set_thetamax(100)
                self.ax.set_rorigin(-0.2)
                self.ax.set_rmin(0.2)
                self.ax.set_rmax(2.6)
                self.ax.set_theta_zero_location('S')
                self.ax.grid(False)
                self.ax.set_xticklabels([])
            else:
                self.ax.set_theta_zero_location('E')

                # self.ax.set_xlabel('r [au]')
                # self.ax.set_ylabel('colatidute [deg]')
                # remove the polar grid
                self.ax.grid(False)
                #remove radial ticks
                self.ax.set_yticklabels([])
                #remove azi ticks
                self.ax.set_xticklabels([])
                self.ax.axis("off")
            
        else:
            if self.verticalplot:
                self.ax.set_ylabel('polar [rad]')
                self.ax.set_xlabel(r'r [$r_p$]')
            if self.midplaneplot:
                self.ax.set_ylabel(r'r [$r_p$]')
                self.ax.set_xlabel('azimuth [rad]')
        if vertical:
            None
            #self.ax.set_ylim(pi/2-0.01,pi/2+0.01)
            #self.ax.set_xlim(0.99,1.01)
            #self.ax.set_xlim(self.zlim[0][0],(self.zlim[0][1]-self.zlim[0][0])+self.zlim[0][1]+0.01)
            #self.ax.set_ylim(self.zlim[0][0],(self.zlim[0][1]-self.zlim[0][0])+self.zlim[0][1]+0.01)

        # if self.azimuth != 0.0:
        #     self.ax.set_title(self.folder+' '+self.field+', azimuth=%.2f'%self.azimuth)
        # else:
        #     self.ax.set_title(self.folder+' '+self.field+' ')

        return None


    def plotmidlvl(self,lvl,polar,grid):
        slcdata = self.slice
        xcsize = (self.xlim[lvl,1]-self.xlim[lvl,0])/self.xdim[lvl]
        ycsize = (self.ylim[lvl,1]-self.ylim[lvl,0])/self.ydim[lvl]
        xplot = linspace(self.xlim[lvl,0],self.xlim[lvl,1],self.xdim[lvl]+1)
        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl]+1)
        xplot_2d, yplot_2d = meshgrid(xplot,yplot)
        #print x and y lengths
        print(f"x length is: {len(xplot)}")
        print(f"y length is: {len(yplot)}")
        if polar:
            if self.velocityfield==True:
                cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.PowerNorm(gamma=0.4325,vmin=-400000,vmax=1.6e6),
                                      shading='flat',cmap='seismic')
                cx.set_edgecolor('face')
            else:
                cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.LogNorm(vmin=self.midmin, vmax=self.midmax),shading='flat',cmap=self.cmap)
                # cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.LogNorm(vmin=0.0000003,  vmax=3),shading='flat',cmap=self.cmap)
            if grid:
                    grid_corners = [[xplot[0],xplot[-1],xplot[-1],xplot[0],xplot[0]],[yplot[0],yplot[0],yplot[-1],yplot[-1],yplot[0]]]
                    print(grid_corners)
                    self.ax.plot(grid_corners[0],grid_corners[1],c='black')
        else:
            if self.velocityfield==True:
                cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.PowerNorm(gamma=0.5,vmin=-1,vmax=3),
                                      shading='flat',cmap='seismic')
            else:
                cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.LogNorm(vmin=self.midmin, vmax=self.midmax),shading='flat',cmap=self.cmap)
                if grid:
                    grid_corners = [[xplot[0],xplot[-1],xplot[-1],xplot[0],xplot[0]],[yplot[0],yplot[0],yplot[-1],yplot[-1],yplot[0]]]
                    print(grid_corners)
                    self.ax.plot(grid_corners[0],grid_corners[1],c='black')

        

        if lvl==0:
            if polar:
                cb = self.fig.colorbar(cx,ax=self.ax,pad=0.15,format='%4.1e',
                                  label='Gas velocity [cm/s]',
                                  ticks=[1.6e6,1e6,3e5,0,-3e5,-4e5])
                # minorticks = cx.norm(arange(-4e5, 2e6, 1e5))
                # cb.ax.xaxis.set_ticks(minorticks, minor=True)
                cb.ax.minorticks_on()
            else:
                self.fig.colorbar(cx,ax=self.ax)
        return None

    def plotmidslice(self,field='gasdensity',reflvls=-1,polar=False,vel=0,z=0,
                     save=False,filename='',grid=False,physical=False):
        self.field=field
        self.polar = polar
        self.midplaneplot = True
        self.save = save
        self.physical = physical
        if ((field=='gasvelocity') or (field=='dustvelocity')):
            self.velocityfield=True
        if reflvls==-1: #plot all the refinement levels
            reflvls = arange(0,self.maxreflvl+1)
        self.upddate_maxmin(field,reflvls,vel,z,avg=False,vert=False,
                            physical=self.physical)

        if polar:
            figsz = (9,9)
        else:
            figsz=(9,5)
        # self.fig, self.ax  = plt.subplots()
        # self.ax.set_aspect('equal')
        
        if polar:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(polar=polar)
            # self.ax1.set_aspect('equal')
            # self.axp = self.fig.add_axes(self.ax1.get_position().bounds, polar=polar,frameon = False)
            # self.fig.subplots_adjust(right=0.8)
            # self.cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
        else:    
            self.fig, self.ax  = plt.subplots()
            # self.ax = self.fig.add_subplot(111)
            # self.ax = self.fig.add_subplot(111,polar=polar)
        for lvl in reflvls:
            if self.velocityfield==True:
                self.slice = self.readdata(field,lvl,self.physical)[vel,-1-z,:,:]
            else:
                self.slice = self.readdata(field,lvl,self.physical)[-1-z,:,:] #midplane data of reflevel 'lvl'
            # self.ax.set_xlim(-self.ylim[lvl,1]+0.2,self.ylim[lvl,1]+0.2)
            # self.ax.set_ylim(-self.ylim[lvl,1]+0.2,self.ylim[lvl,1]+0.2)
            # _len = self.ylim[lvl,1]
            # if lvl == reflvls[-1]:
            #     self.ax1.set_xlim(-self.ylim[lvl,1],self.ylim[lvl,1])
            #     self.ax1.set_ylim(-self.ylim[lvl,1],self.ylim[lvl,1])
            #     self.cax = self.ax1.inset_axes([1.1,-0.9,0.5,2])
            self.plotmidlvl(lvl,polar,grid) #plot midplane data of reflevel 'lvl'
        # self.ax1.set_xlim(-self.ylim[lvl,1],self.ylim[lvl,1])
        # self.ax1.set_ylim(-self.ylim[lvl,1],self.ylim[lvl,1])
        if polar:
        #     self.ax1.set_aspect('equal')
        #     self.axp.set_ylim(0,self.ylim[0,1])
        #     self.ax1.axis("off")
            self.ax.text(x=14*pi/48,y=2.5,s=f"Orbit {self.N}")
        # else:
        #     self.ax.axis("off")
        #print the shape of self.slice
        print(f"shape of self.slice: {self.slice.shape}")
        print(f"values of self.slice: {self.slice[0,0]}")
        self.plot_setlimits(polar,vertical=False)
        plt.tight_layout()
        if save:
            plt.savefig(f"/home/delsender/figures/intro/{filename}_{self.N:03d}.pdf")
        plt.show()
        return None


    def plotvertlvl(self,lvl):
        slcdata = self.slice
        ycsize = (self.ylim[lvl,1]-self.ylim[lvl,0])/self.ydim[lvl]
        zcsize = (self.zlim[lvl,1]-self.zlim[lvl,0])/self.zdim[lvl]
        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl]+1)
        zplot = linspace(self.zlim[lvl,0],self.zlim[lvl,1]+(self.zlim[lvl,1]-self.zlim[lvl,0]),self.zdim[lvl]*2+1)
        yplot_2d, zplot_2d = meshgrid(yplot, zplot)

        if self.polar:
            cont = transpose(zplot_2d)
            zplot_2d = transpose(yplot_2d)
            yplot_2d = cont
            slcdata = transpose(slcdata)
            if self.velocityfield==True:
                cx=self.ax.pcolormesh(yplot_2d,zplot_2d,slcdata,norm=colors.PowerNorm(gamma=0.44),
                                      shading='flat',cmap='seismic')
            else:
                cx=self.ax.pcolormesh(yplot_2d,zplot_2d,slcdata,norm=colors.LogNorm(vmin=self.midmin, vmax=self.midmax),shading='flat',cmap=self.cmap)
                # cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.LogNorm(vmin=0.0000003,  vmax=3),shading='flat',cmap=self.cmap)
        else:
            cx = self.ax.pcolormesh(yplot_2d,zplot_2d,slcdata,norm=colors.LogNorm(vmin=0.0001, vmax=0.3),
                          shading='flat',cmap=self.cmap,alpha=1.0)

        self.fig.colorbar(cx,ax=self.ax,fraction=0.026, pad=0.04)

        return cx

    def plotvertslice(self,field='gasdensity',reflvls=-1,polar=False,azi=0.0,vel=0,
                      filename='',save=False):
        if self.dim!=3:
            print('no vertical slice for 2 dimensional array')
            return None
        self.field=field
        self.polar = polar
        self.azimuth = azi
        self.verticalplot = True
        if ((field=='gasvelocity') or (field=='dustvelocity')):
            self.velocityfield=True
        if reflvls==-1:
            reflvls = arange(0,self.maxreflvl+1)
        self.upddate_maxmin(field,reflvls,vel,z=0,avg=False,vert=True)
        if polar:
            figsz = (9,9)
        else:
            figsz=(12,4)
        self.fig=plt.figure(figsize=figsz)
        self.ax = self.fig.add_subplot(111,polar=polar)
        if azi != 0.0:
            reflvls = arange(0,1);
        for lvl in reflvls:
            azigrid = int(self.xdim[lvl]/2)
            if (lvl==0):
                self.upddate_maxmin(field,[lvl],vel,z=0,azigrid=azigrid,avg=False,vert=True)
            if azi != 0.0:
                azigrid = int(self.xdim[lvl]*(azi + pi)/(2.*pi))
            if self.velocityfield==True:
                self.slice = self.readdata(field,lvl)[vel,:,:,azigrid]
            else:
                self.slice = self.readdata(field,lvl)[:,:,azigrid]
            self.slice = concatenate((self.slice,self.slice[::-1,:]),axis=0)
            cx = self.plotvertlvl(lvl) #plot midplane data of reflevel 'lvl'
        self.plot_setlimits(polar,vertical=True)
        self.ax.text(x=47*pi/80,y=2,s=f"Orbit {self.N}")
        if save:
            plt.savefig(f"/home/delsender/figures/intro/{filename}_{self.N:03d}.png")
        plt.show(block=False)
        return None

    def plotvertavg(self,field='gasdensity',polar=False):
        self.verticalplot = True
        if self.dim!=3:
            print('no vertical average for 2 dimensional array')
            return None
        self.field=field
        self.polar = polar
        lvl = ([0])
        self.upddate_maxmin(field,lvl,avg=True,vert=True)

        if polar:
            figsz = (9,9)
        else:
            figsz=(12,4)
        fig=plt.figure(figsize=figsz)
        self.ax = fig.add_subplot(111,polar=polar)

        self.slice = average(self.readdata(field,0),axis=2)
        self.slice = concatenate((self.slice,self.slice[::-1,:]),axis=0)
        cx = self.plotvertlvl(0) #plot midplane data of reflevel 'lvl'
        fig.colorbar(cx)
        self.plot_setlimits(polar,vertical=True)
        plt.show(block=False)
        return None

    def surfdens(self,field):
        lvl = 0
        zcsize = (self.zlim[lvl,1]-self.zlim[lvl,0])/self.zdim[lvl]
        y = zcsize*linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl])
        dz = reshape(tile(y,self.xdim[0]),(self.xdim[0],self.ydim[0]))
        dz = transpose(dz)
        fielddata = 2.*sum(self.readdata(field,0),axis=0)*dz
        self.surfdens = fielddata
        return fielddata

    def plotgapprofile(self,field='gasdensity',ylog=True,avg=False,azi=0.0):
        lvl=0
        self.azimuth=azi
        self.field=field
        if avg:
            if self.dim == 2:
                data = average(self.readdata(field,0)[0,:,:],axis=1)
            if self.dim == 3:
                data = average(self.surfdens(field),axis=1)
        else:
            azigrid = int(self.xdim[lvl]*(pi+ self.azimuth)/(2.*pi))



            if self.dim == 2:
                data = self.readdata(field,0)[0,:,azigrid]
            if self.dim == 3:
                data = self.surfdens(field)[:,azigrid]

        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl])
        plt.figure(figsize=(8,4))
        plt.plot(yplot,data,color = 'k')
        if ylog:
            plt.yscale('log')
        plt.ylabel(r'$\Sigma$ [$M_{*}/au^2$]')
        plt.xlabel('r [au]')
        if avg:
            plt.title('surafce density profile (azim avg), '+self.folder+' '+self.field+'')
        else:
            plt.title('surafce density profile, azimuth=%.2f'%self.azimuth)
        plt.show(block=False)
        return None
    
    def hill_sphere_mass(self,field='gasdensity',planetmass=0.001,physical=False):
        '''Calculates the mass of gas contained within the Hill
        of the planet.
        
        IN: field       = field qunatity, set to gas density
            planetmass  = mass of planet in code units
            physical    = convert mass into grams
        
        Planet location is assumed to be at:
        r               = 1
        azi             = 0
        co-lat          = pi/2

        Planet mass assumed to be: 
        planetmass = 0.001

        Star mass assumed to be:
        XMSTAR = 1
        
        OUT: rHmass     = mass of gas contained with the Hill sphere
                        in code units unless physical=True'''
        
        # load in variables
        self.field = field 
        self.planetmass = planetmass
        # initialise the Hill mass
        self.rHmass = 0
        # initialise the volume summation
        volsum = 0
        # hardcoding to work for the base mesh
        lvl = 0

        # calculate the Hill radius
        # UPDATE THE 1 TO RP
        # AND THE STAR MASS XMSTAR
        rH = 1 * (planetmass/(3*(1+planetmass)))**(1/3) 
        # calculate the volume of the Hill sphere
        rHvolume = 4*pi*rH**3 / 3

        # load in the field data (gas density by default)
        gasdata  = self.readdata(field,lvl)[:,:,:]
        # get the delta r, theta, and phi
        deltax = float((self.xlim[lvl,1] - self.xlim[lvl,0])/self.xdim)
        deltay = float((self.ylim[lvl,1] - self.ylim[lvl,0])/self.ydim)
        deltaz = float((self.zlim[lvl,1] - self.zlim[lvl,0])/self.zdim)
        #
        # loop over all cells
        # could be done more efficiently with prior knowledge to the data structure
        #
        for i, xi in enumerate(self.lines[self.dim_line+2].split(" ")[2:-4]):
            for j, yi in enumerate(self.lines[self.dim_line+3].split(" ")[2:-4]):
                for k, zi in enumerate(self.lines[self.dim_line+4].split(" ")[2:-4]):
                    # we know where the planet is, so just ignore cells far away
                    # no such condition on z as the disc is very thin!
                    # x, y, z = float(xi), float(yi), float(zi)
                    if abs(float(xi)) < pi/20 and 1-rH < float(yi) < 1+rH:
                        dist = sqrt(1*1 + float(yi)*float(yi) - 2*1*float(yi)*
                                    (sin(pi/2)*sin(float(zi))*cos(0-float(xi)) + 
                                     cos(pi/2)*cos(float(zi))))
                        if dist < rH:
                            # multiple of two comes for the fact we are running a half disc
                            self.rHmass += 2*gasdata[k,j,i]
                            # sum the volumes of each cell
                            # multiply by 2 as we are in a hemisphere
                            x, y, z = float(xi), float(yi), float(zi)
                            volsum += 2 * y*y * sin(z) * deltay * deltax * deltaz
                            # print(zi)
        # we have calculated the density within the Hill sphere
        # multilpy this by Hill sphere volume                    
        self.rHmass = self.rHmass
        if physical:
            self.get_physical_units()
            self.rHmass *= XM0
        print(f"The Hill volume is: {rHvolume}")
        print(f"The cell sum volume is: {volsum}")
        print(f"The Hill mass sphere is: {self.rHmass * rHvolume}")
        print(f"Hill mass cell volume is: {self.rHmass * volsum}")
        return self.rHmass * volsum

class get:
    '''A class to get values out of the simulation data'''
    
    def __init__(self,N):
        '''Set the output up to which to search
        
        IN:     N       = Final output, i.e. search all outputs up to N'''

        self.N = N

    def global_minmax(self,field='gasdensity',physical=False):
        '''Gets the global minimum and maximum for field data
        
        IN:     field           = field of which global minimum and maximum are to be found
            
        OUT:    globmin         = minimum of the global field data after N outputs
                globmax         = maximum of the global field data after N outputs'''
        
        self.globmin = 1e+10
        self.globmax = -1e+10
        self.field = field
        self.velocityfield = False
        self.physical = physical

        if ((field=='gasvelocity') or (field=='dustvelocity')):
            self.velocityfield = True

        print(self.physical)
        # print(self.velocityfield)

        for i in range(self.N):
            tempmin, tempmax = output(i+1).get_minmax(field=self.field,
                                                      velocity=self.velocityfield,
                                                      physical=self.physical)

            if tempmin < self.globmin:
                self.globmin = tempmin
            if tempmax > self.globmax:
                self.globmax = tempmax
        
        return self.globmin, self.globmax


