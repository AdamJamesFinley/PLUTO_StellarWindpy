# -*- coding: utf-8 -*-

"""Main module."""
#================ CLASS DECLARATION ====================
class StellarWindData:

	def __init__(self,wdir,dimensions,dataset=-1,STRENGTH=1,Vphi_Vkep=4.46e-3,datatype='double',RF=True):
		'''
		Reads in the PLUTO dataset using pyPLUTO, calculates coordinates and
		simple derived quantities.
		'''
		plutodir = os.environ['PLUTO_DIR']
		time_data = pp.nlast_info(w_dir=wdir,datatype=datatype)
		self.lastStep = time_data['nlast']
		if dataset == -1:
			self.currentStep = self.lastStep
		else:
			self.currentStep = dataset
		if int(self.currentStep) > int(self.lastStep):
			sys.exit('Required dataset must be less than : '+str(self.lastStep))
		self.dimensions = dimensions
		self.omega = Vphi_Vkep
		self.bstar = STRENGTH #or use D.B1_tot[0,0]
		if RF == True:
			self.RF = True
		else:
			self.RF = False

		#Output Variables
		D = pp.pload(self.currentStep,w_dir=wdir,datatype=datatype)
		self.rCoord = D.x1
		self.thCoord = D.x2
		self.thDelta = D.dx2
		if self.dimensions == 3:
			self.phiCoord = D.x3

		if self.dimensions == 2:
			self.thGrid, self.rGrid = np.meshgrid(D.x2, D.x1)
			self.xGrid = np.outer(D.x1, np.sin(D.x2))
			self.yGrid = np.outer(D.x1, np.cos(D.x2))
			self.rcylGrid = np.abs(self.xGrid)
		elif self.dimensions == 3:
			self.thGrid, self.rGrid, self.phiGrid = np.meshgrid(D.x2, D.x1, D.x3)
			self.dthGrid, self.drGrid, self.dphiGrid = np.meshgrid(D.dx2, D.dx1, D.dx3)
			self.xGrid = np.outer(D.x1,np.outer(np.sin(D.x2),np.cos(D.x3))).reshape(D.x1.shape[0],D.x2.shape[0],D.x3.shape[0])
			self.yGrid = np.outer(D.x1,np.outer(np.cos(D.x2),np.ones_like(D.x3))).reshape(D.x1.shape[0],D.x2.shape[0],D.x3.shape[0])
			self.zGrid = np.outer(D.x1,np.outer(np.sin(D.x2),np.sin(D.x3))).reshape(D.x1.shape[0],D.x2.shape[0],D.x3.shape[0])
			self.rcylGrid = np.sqrt(self.xGrid**2+self.zGrid**2)

		self.rho = D.rho
		self.logrho = np.log10(D.rho)
		self.prs = D.prs
		self.BrGrid = D.B1_tot
		self.BthGrid = D.B2_tot
		self.BphiGrid = D.B3_tot
		self.VrGrid = D.vx1
		self.VthGrid = D.vx2
		self.VphiGrid = D.vx3
		if self.RF == True:
			self.VphiRF = self.VphiGrid+(self.omega)*self.rcylGrid
		else:
			self.VphiRF = self.VphiGrid

		if self.dimensions == 2:
			self.BxGrid = self.BrGrid*np.sin(self.thGrid)*np.cos(0.) + self.BthGrid*np.cos(self.thGrid)*np.cos(0.) - self.BphiGrid*np.sin(0.)
			self.BzGrid = self.BrGrid*np.sin(self.thGrid)*np.sin(0.) + self.BthGrid*np.cos(self.thGrid)*np.sin(0.) + self.BphiGrid*np.cos(0.)
			self.ByGrid = self.BrGrid*np.cos(self.thGrid) - self.BthGrid*np.sin(self.thGrid)
		elif self.dimensions == 3:
			self.BxGrid = self.BrGrid*np.sin(self.thGrid)*np.cos(self.phiGrid) + self.BthGrid*np.cos(self.thGrid)*np.cos(self.phiGrid) - self.BphiGrid*np.sin(self.phiGrid)
			self.BzGrid = self.BrGrid*np.sin(self.thGrid)*np.sin(self.phiGrid) + self.BthGrid*np.cos(self.thGrid)*np.sin(self.phiGrid) + self.BphiGrid*np.cos(self.phiGrid)
			self.ByGrid = self.BrGrid*np.cos(self.thGrid) - self.BthGrid*np.sin(self.thGrid)

		#Derived Variables:
		self.Vpol = np.sqrt(self.VrGrid*self.VrGrid+self.VthGrid*self.VthGrid)
		self.Bpol = np.sqrt(self.BrGrid*self.BrGrid + self.BthGrid*self.BthGrid)
		self.Vmag =np.sqrt(self.VrGrid*self.VrGrid+self.VthGrid*self.VthGrid+self.VphiGrid*self.VphiGrid)
		self.Bmag = np.sqrt(self.BrGrid*self.BrGrid + self.BthGrid*self.BthGrid + self.BphiGrid*self.BphiGrid)
		self.Cs = np.sqrt(1.05*self.prs/self.rho)
		self.VApol = self.Bpol/np.sqrt(self.rho)
		self.VAphi = self.BphiGrid/np.sqrt(self.rho)
		self.Mach = self.Vpol/self.Cs
		self.MachB = self.Vpol/self.VApol
		self.MachBf = 2.*self.Vpol*self.Vpol/(self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi-np.sqrt((self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi)*(self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi)-4*self.Cs*self.Cs*self.VApol*self.VApol))
		self.MachBs = 2.*self.Vpol*self.Vpol/(self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi+np.sqrt((self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi)*(self.Cs*self.Cs+self.VApol*self.VApol+self.VAphi*self.VAphi)-4*self.Cs*self.Cs*self.VApol*self.VApol))
		self.Temp = self.prs/self.rho

		#Arrow Stuff (for 2D atm)
		I = pp.Image()
		self.arrowR,self.arrowZ,SphData = I.getSphData(D, w_dir=wdir)
		self.arrowVr = SphData['v1c']
		self.arrowVth = SphData['v2c']

	def calcWindProperties(self):
		'''
		Lambda - Quantity is related to the angular momentum flux. The total flux of angular
		momentum per unit of poloidal magnetic flux is constant along a fieldline.
		Mdot - Mass loss rate.
		Jdot - Angular momentum loss rate.
		Bflux - Magnetic flux.
		'''

		#Functions To Integrate
		if self.dimensions == 2:
			self.VB = (self.VrGrid*self.BrGrid+self.VthGrid*self.BthGrid)/self.Bpol
			dA = 2*np.pi*np.outer(self.rCoord*self.rCoord, np.sin(self.thCoord))
		if self.dimensions == 3:
			self.VB = (self.VrGrid*self.BrGrid+self.VthGrid*self.BthGrid)/self.Bpol
			dA = np.outer(self.rCoord*self.rCoord,np.outer(np.sin(self.thCoord),np.ones_like(self.phiCoord))).reshape(self.rCoord.shape[0],self.thCoord.shape[0],self.phiCoord.shape[0])
		self.Lambda = self.rcylGrid*(self.VphiRF-self.BphiGrid*self.Bmag/self.rho/self.VB)
		self.LambdaPlasma = self.rcylGrid*(self.VphiRF)
		self.LambdaField = self.rcylGrid*(-self.BphiGrid*self.Bmag/self.rho/self.VB)
		F1 = self.Lambda*self.rho*self.VrGrid*dA #Integrand F(theta)
		F2 = self.rho*self.VrGrid*dA
		F3 = np.abs(self.BrGrid*dA)
		if self.dimensions == 2:
			#Trapezium rule for Jdot, Mdot etc
			def integrateFunction(self,F):
				'''
				Uses Trapezium rule to integrate the function F(theta).
				'''
				Result = 0.5*(self.thDelta[0]*F[0:len(self.rCoord),0]+ self.thDelta[len(self.thCoord)-1]*F[0:len(self.rCoord),len(self.thCoord)-1]) #Initial and final components
				for j in range(1,len(self.thCoord)-2):                 #Sum central terms
						tmp = self.thDelta[j]*F[0:len(self.rCoord),j]
						Result += tmp
				return Result
			self.jdotVsRadius = integrateFunction(self,F1)
			self.mdotVsRadius = integrateFunction(self,F2)
			self.bfluxVsRadius = integrateFunction(self,F3)
		if self.dimensions == 3:
			#Integrate over a sphercal shell
			def sphericalSum(self,F):
				Result = [np.zeros(len(self.rCoord))]
				tmp=[]
				for i in range(0,len(self.thCoord)-1):
					for j in range(0,len(self.phiCoord)-1):
						tmp = F[:,i,j]*self.dthGrid[:,i,j]*self.dphiGrid[:,i,j]
						Result+=tmp
				return Result
			self.jdotVsRadius = sphericalSum(self,F1)[0]
			self.mdotVsRadius = sphericalSum(self,F2)[0]
			self.bfluxVsRadius = sphericalSum(self,F3)[0]

		Res = len(self.rCoord)
		UPBOUND = int(0.9*Res)
		LOWBOUND = int(0.6*Res)
		self.Jdot = np.median(self.jdotVsRadius[LOWBOUND:UPBOUND])
		self.Mdot = np.median(self.mdotVsRadius[LOWBOUND:UPBOUND])
		self.Bopen = np.median(self.bfluxVsRadius[LOWBOUND:UPBOUND])
		self.surfaceFlux = self.bfluxVsRadius[0]
		#Output Variables
		self.AlfvRad = np.sqrt(self.Jdot/self.Mdot/self.omega)
		self.Upsilon = 4.*np.pi*self.bstar**2/(self.Mdot*np.sqrt(2.))
		self.Upsilon_open = 4.*np.pi*self.Bopen**2/(self.Mdot*np.sqrt(2.))
		#Upsilon_star = 2*pi*B_FLUX**2/(Mdot*sqrt(2))
		print("Alfven Radius = "+str(self.AlfvRad))
		print("Upsilon = "+str(self.Upsilon))
		print("Upsilon_open = "+str(self.Upsilon_open))
		if self.Upsilon > 4e4 :
			print('You may need to concider time variance.')

	def plotColormap(self,variable,label='',slice=0,cmap='jet',vmin=False,vmax=False,zoom=20.):
		'''
		Simple colourmap plotting tool. For quick looks at the data.
		'''
		variable = getattr(self, variable)
		if (vmin == False and vmax == False):
			vmin = np.min(variable)
			vmax = np.max(variable)
		plt.clf()
		fig, ax = plt.subplots(figsize=[10,10])
		#Plot Colormap
		if self.dimensions == 2:
			plt.pcolormesh(self.xGrid,self.yGrid,variable, cmap=cmap)
		if self.dimensions == 3:
			plt.pcolormesh(self.xGrid[:,:,slice],self.yGrid[:,:,slice],variable[:,:,slice], cmap=cmap)
		plt.colorbar(label=label)
		plt.clim(vmin,vmax)
		ax.set_ylabel('Stellar Radii')
		ax.set_xlabel('Stellar Radii')
		ax.set_ylim(-zoom,zoom)
		ax.set_xlim(0,zoom)
		ax.set_aspect('equal')

		#Draw The Star
		circle1 = plt.Circle((0, 0), 1.0, color='orange')
		ax.add_artist(circle1)
		plt.show()

	def calcStreamFunction(self):
		'''
		Calculates the magnetic vector potential in the toroidal direction, which
		can be used to produce the scalar phi, whose contours represent fieldlines
		in axsymmetric MHD simulations.
		'''
		if self.dimensions == 3:
			sys.exit('StreamFunction does not work for 3D fields')
		Ap = np.zeros_like(self.BrGrid) #Create blank array for Aphi (and sets Ap[0,0]=0)
		for i in range(0,len(self.rCoord)-1): #Create initial array from which the rest of the potential can be calculated
			Ap[i+1,0]=-self.rCoord[i]*self.BthGrid[i,0]*(self.rCoord[i+1]-self.rCoord[i])/self.rCoord[i+1] + self.rCoord[i]*Ap[i,0]/self.rCoord[i+1]
		for i in range(0,len(self.rCoord)-1):
			for j in range(0,len(self.thCoord)-1):
				Ap[i,j+1]=self.BrGrid[i,j]*self.rCoord[i]*np.sin(self.thCoord[j])*(self.thCoord[j+1]-self.thCoord[j])/np.sin(self.thCoord[j+1]) + np.sin(self.thCoord[j])*Ap[i,j]/np.sin(self.thCoord[j+1])
				Ap[len(self.rCoord)-1,j] = Ap[len(self.rCoord)-2,j]*self.rCoord[len(self.rCoord)-2]/self.rCoord[len(self.rCoord)-1]
		self.streamFunction = Ap*self.xGrid

	def plotStreamFunction(self):
		'''
		Plots the magnetic field lines only. calcStreamFunction() must be called
		prior to use.
		'''
		if self.dimensions == 3:
			sys.exit('Stream Function does not work for 3D fields')
		plt.clf()
		fig, ax = plt.subplots(figsize=[10,10])
		ax.set_ylabel('Stellar Radii')
		ax.set_xlabel('Stellar Radii')
		ax.set_ylim(-20,20)
		ax.set_xlim(0,20)
		ax.set_aspect('equal')
		#Plot Stream  Function
		plt.contour(self.xGrid, self.yGrid, self.streamFunction, 30, colors='k', linestyles='-',linedwidth=0.1)
		#Draw The Star
		circle1 = plt.Circle((0, 0), 1.0, color='orange')
		ax.add_artist(circle1)
		plt.show()

	def calcOmegaEff(self):
		'''
		Effective rotation of the plasma, conserved along streamlines (fieldlines).
		'''
		if self.dimensions == 2:
			self.omegaEff = 1./self.rcylGrid*(self.VphiRF-self.VB*self.BphiGrid/self.Bpol)/(self.omega)
			self.omegaEffPlasma = 1./self.rcylGrid*(self.VphiRF)/(self.omega)
			self.omegaEffField = 1./self.rcylGrid*(-self.VB*self.BphiGrid/self.Bpol)/(self.omega)
		if self.dimensions == 3:
			self.omegaEff = 1./self.rcylGrid*(self.VphiRF-self.VB*self.BphiGrid/self.Bpol)/(self.omega)
			self.omegaEffPlasma = 1./self.rcylGrid*(self.VphiRF)/(self.omega)
			self.omegaEffField = 1./self.rcylGrid*(-self.VB*self.BphiGrid/self.Bpol)/(self.omega)

	def plotOmegaEffective(self,slice=0,all=False,pltlimit=False,averaged=False,smoothed=False,max=False):
		'''
		Plot the effective rotation rate versus theta angle. calcOmegaEff() must be called
		prior to use.
		'''
		plt.clf()
		if self.dimensions == 2:
			plt.scatter(self.thGrid,self.omegaEff,c=self.omegaEff, cmap=Balance_20.mpl_colormap,vmin=0.5, vmax=1.5,s=0.1)
		if self.dimensions == 3:
			if averaged==True:
				plt.scatter(np.mean(self.thGrid,axis=2),np.mean(self.omegaEff,axis=2),c=np.mean(self.omegaEff,axis=2), cmap=Balance_20.mpl_colormap,vmin=0.5, vmax=1.5,s=0.1)
			else:
				if all==True:
					if smoothed == True:
						plt.scatter(np.mean(self.thGrid[:,:,:],axis=0),np.mean(self.omegaEff[:,:,:],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,:],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*1./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*1./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*1./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*2./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*2./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*2./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*3./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*3./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*3./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*4./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*4./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*4./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*5./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*5./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*5./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*6./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*6./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*6./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						# plt.scatter(np.mean(self.thGrid[:,:,int(len(self.phiCoord)*7./8.+slice)],axis=0),np.mean(self.omegaEff[:,:,int(len(self.phiCoord)*7./8.+slice)],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,int(len(self.phiCoord)*7./8.+slice)],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
					else:
						if max == True:
							plt.scatter(np.mean(self.thGrid[:,:,:],axis=0),np.max(self.omegaEff[:,:,:],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,:],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
							plt.scatter(np.mean(self.thGrid[:,:,:],axis=0),np.min(self.omegaEff[:,:,:],axis=0),cmap='jet',c=np.mean(self.phiGrid[:,:,:],axis=0),s=1,alpha=1,vmin=0,vmax=2.*np.pi)
						else:
							plt.scatter(self.thGrid[:,:,slice],self.omegaEff[:,:,slice],cmap='jet',c=self.phiGrid[:,:,slice],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*1./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*1./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*1./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*2./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*2./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*2./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*3./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*3./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*3./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*4./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*4./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*4./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*5./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*5./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*5./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*6./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*6./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*6./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
							plt.scatter(self.thGrid[:,:,int(len(self.phiCoord)*7./8.+slice)],self.omegaEff[:,:,int(len(self.phiCoord)*7./8.+slice)],cmap='jet',c=self.phiGrid[:,:,int(len(self.phiCoord)*7./8.+slice)],s=0.1,alpha=0.6,vmin=0,vmax=2.*np.pi)
					plt.colorbar(label=r'$\phi Slice$')
				else:
					plt.scatter(self.thGrid[:,:,slice],self.omegaEff[:,:,slice],c=self.omegaEff[:,:,slice], cmap=Balance_20.mpl_colormap,vmin=0.5, vmax=1.5,s=0.1)
		plt.ylabel(r'$\Omega_{eff}/\Omega_*$')
		if averaged==True:
			plt.ylabel(r'$\langle\Omega_{eff}/\Omega_*\rangle_{\phi}$')
		plt.xlabel(r'$\theta$')
		if pltlimit==False:
			plt.ylim(0.95*np.min(self.omegaEff),1.05*np.max(self.omegaEff))
		else:
			plt.ylim(1.-pltlimit,1.+pltlimit)
		plt.show()

	def plotLambda(self,averaged=True):
		'''
		'''
		plt.clf()
		fig, ax = plt.subplots()
		if self.dimensions == 2:
			ax.plot(self.thCoord/np.pi*180.-90., self.Lambda[int(len(self.rCoord)*0.85),:]/np.max(self.Lambda[int(len(self.rCoord)*0.85),:]), c='k')
			ax.plot(self.thCoord/np.pi*180.-90., 1*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
			ax.plot(self.thCoord/np.pi*180.-90., 2*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
			axb = ax.twinx()
			axb.plot(self.thCoord/np.pi*180.-90., self.rho[int(len(self.rCoord)*0.85),:]*self.VrGrid[int(len(self.rCoord)*0.85),:]/np.max(self.rho[int(len(self.rCoord)*0.85),:]*self.VrGrid[int(len(self.rCoord)*0.85),:]), c='red')
		if self.dimensions == 3:
			if averaged == True:
				ax.plot(self.thCoord/np.pi*180.-90., np.mean(self.Lambda[int(len(self.rCoord)*0.85),:,:]/np.max(np.mean(self.Lambda[int(len(self.rCoord)*0.85),:,:],axis=1)),axis=1), c='k')
				ax.plot(self.thCoord/np.pi*180.-90., 1*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
				ax.plot(self.thCoord/np.pi*180.-90., 2*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
				axb = ax.twinx()
				axb.plot(self.thCoord/np.pi*180.-90., np.mean(self.rho[int(len(self.rCoord)*0.85),:,:]*self.VrGrid[int(len(self.rCoord)*0.85),:,:]/np.max(np.mean(self.rho[int(len(self.rCoord)*0.85),:,:]*self.VrGrid[int(len(self.rCoord)*0.85),:,:],axis=1)),axis=1), c='red')
			else:
				ax.plot(self.thCoord/np.pi*180.-90., self.Lambda[int(len(self.rCoord)*0.85),:,:]/np.max(self.Lambda[int(len(self.rCoord)*0.85),:,:]), c='k')
				ax.plot(self.thCoord/np.pi*180.-90., 1*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
				ax.plot(self.thCoord/np.pi*180.-90., 2*np.sin(self.thCoord)**2,c='grey',ls='--',zorder=0)
				axb = ax.twinx()
				axb.plot(self.thCoord/np.pi*180.-90., self.rho[int(len(self.rCoord)*0.85),:,:]*self.VrGrid[int(len(self.rCoord)*0.85),:,:]/np.max(self.rho[int(len(self.rCoord)*0.85),:,:]*self.VrGrid[int(len(self.rCoord)*0.85),:,:]), c='red')
		ax.set_ylim(0,2.1)
		axb.set_ylim(0,1.1)
		ax.set_ylabel(r'$\langle\Lambda\rangle_{\phi}$')
		axb.set_ylabel(r'$\langle\rho v_r\rangle_{\phi}$',color='red')
		axb.tick_params(axis='y', labelcolor='red')
		ax.set_xlabel(r'Latitude')
		plt.show()

	def plotOutputFig(self,slice=0,zoom=20.):
		'''
		Creates a "standard" output.
		'''
		plt.clf()
		fig, ax = plt.subplots(figsize=[10,10])
		#Plot Colormap
		if self.dimensions == 2:
			plt.pcolormesh(self.xGrid,self.yGrid,self.omegaEff, cmap=Balance_20.mpl_colormap)
		if self.dimensions == 3:
			plt.pcolormesh(self.xGrid[:,:,slice],self.yGrid[:,:,slice],self.omegaEff[:,:,slice], cmap=Balance_20.mpl_colormap)
		plt.colorbar(label=r'$\Omega_{eff}/\Omega_*$')
		plt.clim(0.5,1.5)
		ax.set_ylabel('Stellar Radii')
		ax.set_xlabel('Stellar Radii')
		ax.set_ylim(-zoom,zoom)
		ax.set_xlim(0,zoom)
		ax.set_aspect('equal')
		#Plot Stream  Function
		if self.dimensions == 2:
			plt.contour(self.xGrid, self.yGrid, self.streamFunction, 30, colors='grey', linestyles='-',linedwidth=0.1)
		if self.dimensions == 3:
			pass
		#Plot Flow Arrows
		newdims = 2*(20,)
		ppt=pp.Tools()
		xcong = ppt.congrid(self.arrowR,newdims,method='linear')
		ycong = ppt.congrid(self.arrowZ,newdims,method='linear')
		velxcong = ppt.congrid(self.arrowVr,newdims,method='linear')
		velycong = ppt.congrid(self.arrowVth,newdims,method='linear')
		plt.gca().quiver(xcong, ycong, velxcong, velycong,color='k')
		#Plot Sonic and Alfv\'en Surfaces
		if self.dimensions == 2:
			Sound = plt.contour(self.xGrid, self.yGrid, self.Mach, [1], linewidths=2, colors='black')
			Alfven = plt.contour(self.xGrid, self.yGrid, self.MachB, [1], linewidths=2, colors='blue')
			AlFast = plt.contour(self.xGrid, self.yGrid, self.MachBf, [1], linewidths=1, linestyles='--', colors='white')
			AlSlow = plt.contour(self.xGrid, self.yGrid, self.MachBs, [1], linewidths=1, linestyles='-.', colors='white')
		if self.dimensions == 3:
			Sound = plt.contour(self.xGrid[:,:,slice], self.yGrid[:,:,slice], self.Mach[:,:,slice], [1], linewidths=2, colors='black')
			Alfven = plt.contour(self.xGrid[:,:,slice], self.yGrid[:,:,slice], self.MachB[:,:,slice], [1], linewidths=2, colors='blue')
			AlFast = plt.contour(self.xGrid[:,:,slice], self.yGrid[:,:,slice], self.MachBf[:,:,slice], [1], linewidths=1, linestyles='--', colors='white')
			AlSlow = plt.contour(self.xGrid[:,:,slice], self.yGrid[:,:,slice], self.MachBs[:,:,slice], [1], linewidths=1, linestyles='-.', colors='white')
		#Indicate Average Alfv\'en Radius
		ax.axvline(x=self.AlfvRad, color='grey', ls='--', lw=5)
		#Draw The Star
		circle1 = plt.Circle((0, 0), 1.0, color='orange')
		ax.add_artist(circle1)
		plt.show()

	def plotRadialFunction(self,variable,ylabel=''):
		plt.clf()
		plt.plot(self.rCoord,getattr(self, variable), c='k')
		plt.xlabel('Stellar Radii')
		plt.ylabel(ylabel)
		if (variable == 'bfluxVsRadius'):
			plt.loglog()
		plt.show()

	def plotSurface(self,variable,vmin=False,vmax=False,slice=0,label='',cmap=Balance_20.mpl_colormap):
		'''
		Plots a quantity on a radial surface.
		'''
		if self.dimensions == 2:
			sys.exit('Surface cannot be shown as a colourmap in 2D')
		plt.clf()
		if vmin == False and vmax == False:
			vmin = np.min(getattr(self, variable)[slice,:,:])
			vmax = np.max(getattr(self, variable)[slice,:,:])
		fig, ax = plt.subplots(figsize=[10,10])
		ax.set_ylabel('Latitude')
		ax.set_xlabel('Longitude')
		ax.set_aspect('equal')
		C=plt.pcolormesh(self.phiGrid[slice,:,:]/np.pi*180., self.thGrid[slice,:,:]/np.pi*180.-90., getattr(self, variable)[slice,:,:], cmap=cmap,vmin=vmin,vmax=vmax)
		ax.set_ylim(-90,90)
		ax.set_xlim(0,360)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="2%", pad=0.05)
		plt.colorbar(C,label=label,cax=cax)
		plt.clim()
		plt.show()

	def calcAvRA(self,angle=0.):
		'''
		'''
		if self.dimensions == 2:
			# Make 3D Grid  ###########
			TH = self.thGrid[0,:]
			PHI = np.linspace(0.,2.*np.pi,2*len(self.thGrid[0,:]))
			PHI_grid, TH_grid = np.meshgrid(PHI,TH)

			# Rotate Grid  ###########
			TH_prime = (angle)*np.pi/180.
			TH_rot=np.empty_like(TH_grid)
			PHI_rot=np.empty_like(PHI_grid)

			for i in range(0,len(self.thGrid[0,:]) ):
			    for j in range(0,2*len(self.thGrid[0,:]) ):
			        old_points = np.array([np.sin(TH_grid[i,j])*np.cos(PHI_grid[i,j]),np.sin(TH_grid[i,j])*np.sin(PHI_grid[i,j]),np.cos(TH_grid[i,j])])
			        old_axis = np.array([np.sin(-np.pi/2)*np.cos(0),np.sin(-np.pi/2)*np.sin(0),np.cos(-np.pi/2)])
			        transform = np.cos(TH_prime)*old_points+np.sin(TH_prime)*np.cross(old_axis,old_points)+np.dot(old_axis,old_points)*(1-np.cos(TH_prime))*old_axis
			        TH_rot[i,j] = np.pi-np.arctan2(np.sqrt(transform[0:1]**2+transform[1:2]**2),transform[2:3])
			        PHI_rot[i,j] = np.arctan2(transform[1:2],transform[0:1])

			RAsum = 0
			Mdotsum = 0
			for k in np.arange(0, 2*len(self.thGrid[0,:]) ):
				for j in np.arange(0, len(self.thGrid[0,:]) ):
					for i in np.arange(0, len(self.rCoord) ):
						if (self.MachB[i,j] > 1.):
							if (j>1 and j<len(self.thGrid[0,:])-1 and k>1 and k<2*len(self.thGrid[0,:])-1):
								dVa_r = (self.MachB[i+1,j]-self.MachB[i-1,j])/(self.rCoord[i+1]-self.rCoord[i-1])
								dVa_th = (self.MachB[i,j+1]-self.MachB[i,j-1])/(self.thCoord[j+1]-self.thCoord[j-1])
								dVa_phi = (self.MachB[i,j]-self.MachB[i,j])/(PHI[k+1]-PHI[k-1])
								weightR = dVa_r/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								weightTh = dVa_th/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								weightPhi = dVa_phi/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								AlfvenRadius = self.rCoord[i-1] + (self.rCoord[i-1]-self.rCoord[i])*(1-self.MachB[i-1,j])/(self.MachB[i-1,j]-self.MachB[i,j])
								AlfvenRho = self.rho[i-1,j] + (self.rho[i-1,j]-self.rho[i,j])*(1-self.MachB[i-1,j])/(self.MachB[i-1,j]-self.MachB[i,j])
								AlfvenVr = self.VrGrid[i-1,j] + (self.VrGrid[i-1,j]-self.VrGrid[i,j])*(1-self.MachB[i-1,j])/(self.MachB[i-1,j]-self.MachB[i,j])
								AlfvenVth = self.VthGrid[i-1,j] + (self.VthGrid[i-1,j]-self.VthGrid[i,j])*(1-self.MachB[i-1,j])/(self.MachB[i-1,j]-self.MachB[i,j])
								AlfvenVphi = self.VphiRF[i-1,j] + (self.VphiRF[i-1,j]-self.VphiRF[i,j])*(1-self.MachB[i-1,j])/(self.MachB[i-1,j]-self.MachB[i,j])
								RAsum = RAsum + (AlfvenRadius)**2.*np.sin(TH_rot[j,k])**2.*AlfvenRho*( weightR*AlfvenVr+weightTh*AlfvenVth-weightPhi*AlfvenVphi )*2.*(AlfvenRadius)**2.*np.sin(self.thCoord[j])*np.sin(np.pi/2./len(TH))
								Mdotsum = Mdotsum + AlfvenRho*( weightR*AlfvenVr+weightTh*AlfvenVth-weightPhi*AlfvenVphi )*2.*(AlfvenRadius)**2.*np.sin(self.thCoord[j])*np.sin(np.pi/2./len(TH))
							else:
								RAsum = RAsum + self.rCoord[i]**2.*np.sin(TH_rot[j,k])**2.*self.rho[i,j]*self.VrGrid[i,j]*2.*self.rCoord[i]**2.*np.sin(self.thCoord[j])*np.sin(np.pi/2./len(TH))
								Mdotsum = Mdotsum + self.rho[i,j]*self.VrGrid[i,j]*2.*self.rCoord[i]**2.*np.sin(self.thCoord[j])*np.sin(np.pi/2./len(TH))
							break
		elif self.dimensions == 3:
			# Rotate Grid  ###########
			TH_prime = (angle)*np.pi/180.
			TH_rot=np.empty_like(self.thGrid[0,:,:])
			PHI_rot=np.empty_like(self.phiGrid[0,:,:])

			for i in range(0,len(self.thCoord) ):
			    for j in range(0,len(self.phiCoord) ):
			        old_points = np.array([np.sin(self.thGrid[0,i,j])*np.cos(self.phiGrid[0,i,j]),np.sin(self.thGrid[0,i,j])*np.sin(self.phiGrid[0,i,j]),np.cos(self.thGrid[0,i,j])])
			        old_axis = np.array([np.sin(-np.pi/2)*np.cos(0),np.sin(-np.pi/2)*np.sin(0),np.cos(-np.pi/2)])
			        transform = np.cos(TH_prime)*old_points+np.sin(TH_prime)*np.cross(old_axis,old_points)+np.dot(old_axis,old_points)*(1-np.cos(TH_prime))*old_axis
			        TH_rot[i,j] = np.pi-np.arctan2(np.sqrt(transform[0:1]**2+transform[1:2]**2),transform[2:3])
			        PHI_rot[i,j] = np.arctan2(transform[1:2],transform[0:1])

			RAsum = 0
			Mdotsum = 0
			for k in np.arange(0, len(self.phiCoord) ):
				for j in np.arange(0, len(self.thCoord) ):
					for i in np.arange(0, len(self.rCoord) ):
						if (self.MachB[i,j,k] > 1.):
							if (j>1 and j<127 and k>1 and k<255):
								dVa_r = (self.MachB[i+1,j,k]-self.MachB[i-1,j,k])/(self.rCoord[i+1]-self.rCoord[i-1])
								dVa_th = (self.MachB[i,j+1,k]-self.MachB[i,j-1,k])/(self.thCoord[j+1]-self.thCoord[j-1])
								dVa_phi = (self.MachB[i,j,k+1]-self.MachB[i,j,k-1])/(self.phiCoord[k+1]-self.phiCoord[k-1])
								weightR = dVa_r/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								weightTh = dVa_th/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								weightPhi = dVa_phi/np.sqrt(dVa_r**2.+dVa_th**2.+dVa_phi**2.)
								AlfvenRadius = self.rCoord[i-1] + (self.rCoord[i-1]-self.rCoord[i])*(1-self.MachB[i-1,j,k])/(self.MachB[i-1,j,k]-self.MachB[i,j,k])
								AlfvenRho = self.rho[i-1,j,k] + (self.rho[i-1,j,k]-self.rho[i,j,k])*(1-self.MachB[i-1,j,k])/(self.MachB[i-1,j,k]-self.MachB[i,j,k])
								AlfvenVr = self.VrGrid[i-1,j,k] + (self.VrGrid[i-1,j,k]-self.VrGrid[i,j,k])*(1-self.MachB[i-1,j,k])/(self.MachB[i-1,j,k]-self.MachB[i,j,k])
								AlfvenVth = self.VthGrid[i-1,j,k] + (self.VthGrid[i-1,j,k]-self.VthGrid[i,j,k])*(1-self.MachB[i-1,j,k])/(self.MachB[i-1,j,k]-self.MachB[i,j,k])
								AlfvenVphi = self.VphiRF[i-1,j,k] + (self.VphiRF[i-1,j,k]-self.VphiRF[i,j,k])*(1-self.MachB[i-1,j,k])/(self.MachB[i-1,j,k]-self.MachB[i,j,k])
								RAsum = RAsum + (AlfvenRadius)**2.*np.sin(TH_rot[j,k])**2.*AlfvenRho*( weightR*AlfvenVr+weightTh*AlfvenVth-weightPhi*AlfvenVphi )*2.*(AlfvenRadius)**2.*np.sin(self.thCoord[j])*np.sin(np.pi/256.)
								Mdotsum = Mdotsum + AlfvenRho*( weightR*AlfvenVr+weightTh*AlfvenVth-weightPhi*AlfvenVphi )*2.*(AlfvenRadius)**2.*np.sin(self.thCoord[j])*np.sin(np.pi/256.)
							else:
								RAsum = RAsum + self.rCoord[i]**2.*np.sin(TH_rot[j,k])**2.*self.rho[i,j,k]*self.VrGrid[i,j,k]*2.*self.rCoord[i]**2.*np.sin(self.thCoord[j])*np.sin(np.pi/256.)
								Mdotsum = Mdotsum + self.rho[i,j,k]*self.VrGrid[i,j,k]*2.*self.rCoord[i]**2.*np.sin(self.thCoord[j])*np.sin(np.pi/256.)
							break
		return np.sqrt(RAsum/Mdotsum)
#================ EXAMPLES =============================
#Load Data And Derive Quantites
# DIPOLE = StellarWindData(wdir=wdir3Dcone,dimensions=3,STRENGTH=1.,Vphi_Vkep=4.46e-3)
# DIPOLE.calcWindProperties()
# DIPOLE.calcStreamFunction()
# DIPOLE.calcOmegaEff()

#Plotting Functions
# DIPOLE.plotColormap(variable='zGrid',label='Radial Magnetic Field',cmap=Balance_20.mpl_colormap,coneAngle=0.1,slice=0)
# DIPOLE.plotOmegaEffective()
# DIPOLE.plotStreamFunction()
# DIPOLE.plotOutputFig()
# DIPOLE.plotRadialFunction('jdotVsRadius',ylabel=r'$\int_{A}(\Lambda\rho v_r)dA$')
