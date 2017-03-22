#import nova.blanketsegments.utilities as ut
import numpy as np
#import nova.streamfunction as SF
#import nova.setup as Setup

def curvature(x,y,sind,eind,calltype='normal'):
	if len(x)!=len(y):
		raise('Curves are not the same length!')
	sind=int(round(sind))
	eind=int(round(eind))
	dl=np.sqrt((x[eind]-x[sind])**2+(y[eind]-y[sind])**2)
	grad=np.gradient(y)/np.gradient(x)
	acc=np.gradient(grad)/np.gradient(x)
	e=len(x)-1
	grad[0]=(y[1]-y[0])/(x[1]-x[0])
	grad[1]=(y[2]-y[1])/(x[2]-x[1])
	grad[e-1]=(y[e-1]-y[e-2])/(x[e-1]-x[e-2])
	grad[e]=(y[e]-y[e-1])/(x[e]-x[e-1])
	acc[0]=(grad[1]-grad[0])/(x[1]-x[0])
	acc[1]=(grad[2]-grad[1])/(x[2]-x[1])
	acc[e-1]=acc[e-2]#(grad[e-1]-grad[e-2])/(x[e-1]-x[e-2])
	acc[e]=acc[e-2]#(grad[e]-grad[e-1])/(x[e]-x[e-1])
	r=[]
	r.append(abs((1+grad[sind]**2)**1.5)/acc[sind])
	r.append(abs((1+grad[eind]**2)**1.5)/acc[eind])
	theta=[2*np.arcsin(dl/(2*ri)) for ri in r]
	inftycheck=[True if grad[ri]==float('inf') else False for ri in [sind,eind]]
	if all(inftycheck):
		avtheta=0
	else:
		avtheta=np.mean(theta)
	if calltype=='normal':
		return avtheta
	elif calltype=='acc':
		return acc
	elif calltype=='grad':
		return grad
	elif calltype=='check':
		return inftycheck
	elif calltype=='r':
		return r

def blanketloop(x,y,ml=1.5,mt=20,subdiv=5,accuracy=0.5,**kwargs):
	'''
	x and y are the loop coordinates of the blanket segment
	ml is the maximum length of an individual segment
	mt is the maximum angle covered by a blanket segment, default to 20degrees
	acc is accuracy limit for curve averaging; if points are within acc% curvature, they are
		considered to be part of a constant curve and can have segments averaged
	possible kwarg is a thickness array, same length as x and y
	'''
	limit=accuracy/100 #converting accuracy
	
	if 'dr' in kwargs:
		dr=kwargs.get('dr')
		if type(dr)==float:
			dr=[dr for i in x]
	else:
		dr=[0.2 for i in x]
	
	
	mt=mt*np.pi/180 #convert mt to radians
	sind=[] # start index list
	eind=[] # end index list
	dl=[] # segment length
	ac=[] # average curvature
	mp=[] # midpoint array
	nml=[] # array of normals
	xpath=[] # x coordinate objects
	ypath=[] # y coordinate objects
	hoz=[] # determine if line is horizontal or vertical
	ver=[] # determine if line is horizontal or vertical
	mp=[] # midpoint of segment
	xii=[] # x, inner, start coordinate
	xie=[] # x, inner, final coordinate
	xoi=[] # x, outer, start coordinate
	xoe=[] # x, outer, final coordinate
	yii=[] # y, inner, start coordinate
	yie=[] # y, inner, final coordinate
	yoi=[] # y, outer, start coordinate
	yoe=[] # y, outer, final coordinate
	theta0=[] # angle of incline to chord
	
	laveraged=[] #length of averaged section
	avstartpoint=[] #start index of an averaged section
	avendpoint=[] # end index of an averaged section
	navgd=[] #number of segments to be averaged
	
	'''
	This is a fudged call, don't use the 'grad' calltype normally
	This calltype still calculates curvature as an angle, but it doesn't affect the output
	So the 'start' and 'end' indices don't matter
	It simply returns d2y/dx2 for the entire curves of x and y
	
	This section determines if the blanket is on a straight or vertical, and 
	'''
	d2ydx2=curvature(x,y,0,len(x)-1,calltype='acc')
	dydx=curvature(x,y,0,len(x)-1,calltype='grad')
	nody=[True if y[i+1]==y[i] else False for i,z in enumerate(y[:-1])]
	nodx=[True if x[i+1]==x[i] else False for i,z in enumerate(x[:-1])]
	
	#checking for constant gradient
	consgrad=[]
	#checking for constant curvature
	conscurve=[]
	fullcurve=[curvature(x,y,i-1,i) for i,a in enumerate(x[:])]
	
	
	pv=1 # positional variable, or could use 'i' for index
	si=0 # segment position index
	sn=0 # segment number
	while si+pv<len(x)-1:
		sind.append(si)
		dl.append(np.sqrt((x[pv+si]-x[si])**2+(y[pv+si]-y[si])**2))
		while dl[sn]<ml:
			pv=pv+1
			if pv+si<len(x)-1:
				dl[sn]=np.sqrt((x[pv+si]-x[si])**2+(y[pv+si]-y[si])**2)
			else:
				break
		# find first position more than ml away from the initial sgment	
		pv=pv-1
		#calculate a preliminary end point
		eind.append(pv+si)
		#subdivide it into checkpoints
		checkpoints=[round(a) for a in np.linspace(sind[sn],eind[sn],subdiv)]
		cc=[curvature(x,y,int(a)+1,int(a)) for a in checkpoints[:-1]]
		#calculate curvature at checkpoints, which must be below threshold
		while max(cc)>ml:
			'''
			if curvature is not below threshold, then move back a node and try again
			This can probably be sped up, by increasing the rate checks move back
			'''
			pv=pv-1
			eind[sn]=pv+si
			checkpoints=[round(a) for a in np.linspace(sind[sn],eind[sn])]
			cc=[curvature(x,y,a+1,a) for a in checkpoints[:-1]]
		# calculate average curvature
		ac.append(curvature(x,y,sind[sn],eind[sn]))
		'''
		now we have start and end positions, length, thickness and average curvature: 
		the equivalent rhombus can be defined from these
		
		first, we start from the initial starting point of the segment
		- generate the close outer point
		- generate the far outer point
		
		These are calculated from the point locations and angles between lines
		
		equation is as follows:
		th_tot = th_0 + pi/2 + th_av/2
		[x,y]_o,i = [x,y]_i,i + dr*[cos(th_tot),sin(th_tot)]/cos(th_av/2)
		[x,y]_o,e = [x,y]_o,i + (dl + 2dr * tan(th_av/2)) * [cos(th_0/2),sin(th_0/2)]
		'''
		mp.append(int(round(np.mean([sind[sn],eind[sn]]))))
		xii.append(x[sind[sn]])
		yii.append(y[sind[sn]])
		xie.append(x[eind[sn]])
		yie.append(y[eind[sn]])
		# finally commit segment length
		dl[sn]=np.sqrt((xie[sn]-xii[sn])**2+(yie[sn]-yii[sn])**2)
		theta0.append(np.arctan((yie[sn]-yii[sn])/(xie[sn]-xii[sn])))
		if sn==0:
			xoi.append(xii[sn]+dr[int(mp[sn])]*np.cos(theta0[sn]+np.pi/2+ac[sn]/2)/np.cos(ac[sn]/2))
			yoi.append(yii[sn]+dr[int(mp[sn])]*np.sin(theta0[sn]+np.pi/2+ac[sn]/2)/np.cos(ac[sn]/2))
		else:
			xoi.append(xoe[sn-1])
			yoi.append(yoe[sn-1])
		xoe.append(xoi[sn]+(dl[sn]+2*dr[int(mp[sn])]*np.tan(ac[sn]/2))*np.cos(theta0[sn]))
		yoe.append(yoi[sn]+(dl[sn]+2*dr[int(mp[sn])]*np.tan(ac[sn]/2))*np.sin(theta0[sn]))
		'''
		Need to stroke the path... how best to do this?
		hmmm...
		
		Generate new objects with as x,y coordinates, append interpolated points?
		'''
		
		# Move to next blanket segment
		sn=sn+1
		si=si+pv # define new starting point
		pv=1 # define new positional variable
	'''
	Need to modify for straight sections and areas of constant curvature
	The aim is to group areas which have common gradients or curvatures
	Several segments with identical gradients/curvatures are then given an averaged length
	The number of averaged sections, and the number of modules which have been averaged,
		is stored in 'navgd', which can be accessed by calling the script with
		outtypes='averaged'
	N.B. no further checks are made on length or on curvature
		this is intentional
			as if segments i:j all have equal curvature and have all been checked,
			or segments i:j all have L=<L_max, presumably with most L=L_max and one L<L_max
		then taking segments with averaged lengths,
			theta_av=theta_i(=theta_j)
			L_av = sum(L_n, for n=i:j)/(j-i) < L_max
	'''
	index=list(range(0,len(sind)))
	#Start from the first segment
	i=0
	#end when the final segment is reached
	while i<len(eind)-1:
		#define starting index
		sd=sind[i]
		#check to determine if sections have been averaged
		check=0
		#work backwards from the end segment to tell if a section can be averaged
		for j,ed in zip(index[i:][::-1],eind[i:][::-1]):
			#first averaging criteria, defined from the gradient
			#either the cradient is constant (to an approximation) or there is no change
			#in x OR no change in y, but not both (i.e. xor)
			gradcrit=np.std(dydx[sd:ed])/np.mean(dydx[sd:ed])<limit or (all(nody[sd:ed]) != all(nodx[sd:ed]))
			#second averaging criteria, defined from the curvature
			curvcrit=np.std(fullcurve[sd:ed])/np.mean(fullcurve[sd:ed])<limit
			#if either is met, then begin averaging
			if gradcrit or curvcrit:
				#calculating some useful info on an average section
				laveraged.append(np.sqrt((x[ed]-x[sd])**2+(y[ed]-y[sd])**2))
				avstartpoint.append(sd)
				avendpoint.append(ed)
				navgd.append(j-i)
				#iterate over the segments to be averaged
				#as j-i is the number of units from i to j,
				#the number of nodes from i to j is (j-i)+1
				for p in [int(q) for q in np.linspace(i,j,j-i+1)]:
					#define start and end points
					sind[p]=round((ed-sd)*((p)/(j-i+1))+sd)
					eind[p]=round((ed-sd)*((p+1)/(j-i+1))+sd)
					#finding new points, only roughly found
					xii[p]=x[sind[p]]
					xie[p]=x[eind[p]]
					yii[p]=y[sind[p]]
					yie[p]=y[eind[p]]
					#defining new curvature, midpoint and length
					theta0[p]=np.arctan((yie[p]-yii[p])/(xie[p]-xii[p]))
					dl[p]=np.sqrt((xie[p]-xii[p])**2+(yie[p]-yii[p])**2)
					mp[p]=int(round(np.mean([sind[p],eind[p]])))
					#defining outboard faces from thickness, length, curvature and/or previous points
					if p==0:
						xoi[p]=xii[p]+dr[int(mp[p])]*np.cos(theta0[p]+np.pi/2+ac[p]/2)/np.cos(ac[p]/2)
						yoi[p]=yii[p]+dr[int(mp[p])]*np.sin(theta0[p]+np.pi/2+ac[p]/2)/np.cos(ac[p]/2)
					else:
						xoi[p]=xoe[p-1]
						yoi[p]=yoe[p-1]
					xoe[p]=xoi[p]+(dl[p]+2*dr[int(mp[p])]*np.tan(ac[p]/2))*np.cos(theta0[p])
					yoe[p]=yoi[p]+(dl[p]+2*dr[int(mp[p])]*np.tan(ac[p]/2))*np.sin(theta0[p])
				i=j
				check=1
				break
				
'''
this was the working code for the curved section, held separate to the gradient section
The maths is the same, but the criteria is different
So the criteria were paired and the maths unified to simplify code.
Check the script works without this section and just the criteria sorted

'''
#			if np.std(fullcurve[sd:ed])/np.mean(fullcurve[sd:ed])<limit:
#				laveraged.append(np.sqrt((x[ed]-x[sd])**2+(y[ed]-y[sd])**2))
#				avstartpoint.append(sd)
#				avendpoint.append(ed)
#				navgd.append(j-i)
#				for p in [int(q) for q in np.linspace(i,j,j-i+1)]:
#					sind[p]=round((ed-sd)*((p)/(j-i+1))+sd)
#					eind[p]=round((ed-sd)*((p+1)/(j-i+1))+sd)
#					xii[p]=x[sind[p]]
#					xie[p]=x[eind[p]]
#					yii[p]=y[sind[p]]
#					yie[p]=y[eind[p]]
#					theta0[p]=np.arctan((yie[p]-yii[p])/(xie[p]-xii[p]))
#					dl[p]=np.sqrt((xie[p]-xii[p])**2+(yie[p]-yii[p])**2)
#					mp[p]=int(round(np.mean([sind[p],eind[p]])))
#					if p==0:
#						xoi[p]=xii[p]+dr[int(mp[p])]*np.cos(theta0[p]+np.pi/2+ac[p]/2)/np.cos(ac[p]/2)
#						yoi[p]=yii[p]+dr[int(mp[p])]*np.sin(theta0[p]+np.pi/2+ac[p]/2)/np.cos(ac[p]/2)
#					else:
#						xoi[p]=xoe[p-1]
#						yoi[p]=yoe[p-1]
#					xoe[p]=xoi[p]+(dl[p]+2*dr[int(mp[p])]*np.tan(ac[p]/2))*np.cos(theta0[p])
#					yoe[p]=yoi[p]+(dl[p]+2*dr[int(mp[p])]*np.tan(ac[p]/2))*np.sin(theta0[p])
#				i=j
#				check=1
#				break
		if check==0:
			i=i+1
	#this section is largely for nice presentation; the points are joined together
	for sn in index:
		xpath.append([])
		xpath[sn]= np.linspace(xii[sn],xoi[sn],30).tolist()
		xpath[sn].extend(np.linspace(xoi[sn],xoe[sn],90).tolist())
		xpath[sn].extend(np.linspace(xoe[sn],xie[sn],30).tolist())
		xpath[sn].extend(np.linspace(xie[sn],xii[sn],90).tolist())
		
		ypath.append([])
		ypath[sn]= np.linspace(yii[sn],yoi[sn],30).tolist()
		ypath[sn].extend(np.linspace(yoi[sn],yoe[sn],90).tolist())
		ypath[sn].extend(np.linspace(yoe[sn],yie[sn],30).tolist())
		ypath[sn].extend(np.linspace(yie[sn],yii[sn],90).tolist())
		
	'''
	this is used for bugchecking:
	if you want to look at individual numbers used in the function, add an entry
	and call the function with that entry
	'''
	if 'outtype' in kwargs:
		outtype=kwargs.get('outtype')
		if type(outtype)==str:
			if outtype.lower()=='normal':
				return xpath,ypath
			elif outtype.lower()=='points':
				return xoi,xoe,yoi,yoe
			elif outtype.lower()=='ind':
				return sind,eind
			elif outtype.lower()=='angle':
				return ac,theta0
			elif outtype.lower()=='lengths':
				return dl
			elif outtype.lower()=='angles':
				return fullcurve
			elif outtype.lower()=='averaged':
				return navgd
			elif outtype.lower()=='ind':
				return sind,eind
			elif outtype.lower()=='gradchecks':
				return nody,nodx
			elif outtype.lower()=='ij':
				return i,j
			elif outtype.lower()=='midpoints':
				return mp
		else:
			return xpath,ypath
	#without any input, it defaults to giving the blanket module segment perimeters
	else:
		return xpath,ypath