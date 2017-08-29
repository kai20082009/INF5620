def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    ):
	
	"""
	Problem 1)
	The scaled equation is given as
	x'' = -beta/(1-beta)(1 - beta/|L|)x	
	y'' = -beta/(1-beta)(1 - beta/|L|)(y-1) -beta
	
	|L| = ( x**2 + ( y-1)**2)**0.5
	with inital conditions
	x(0) = (1 +epsilon)sin(Theta) 	, x'(0)=0
	y(0) = 1- (1 +epsilon)cos(Theta), y'(0)=0
	The length of the wire, denoted as L, is depenedent on timestep n.
	|L|= L --> L(t) --> L(x(t(n),y(t(n))) = L(n)
	The  discretized expression of second derivatives is
	u'' = [u(n+1) -2*u(n) + u(n-1) ]/dt**2
	
	This leads 
	[x(n+1) -2*x(n) + x(n-1) ]/dt**2 = -beta/(1-beta)(1-beta/L(n)x(n) 
	[y(n+1) -2*y(n) + y(n-1) ]/dt**2 = -beta/(1-beta)(1-beta/L(n))(y(n)-1) - beta 
	Rearranging these expressions yields
	x(n+1) =  [dt**2*beta/(1-beta)]*(1-beta/L(n))*x(n) +2*x(n) - x(n-1)
	y(n+1) =  [dt**2*beta/(1-beta)]*(1-beta/L(n))*(y(n)-1) - beta*dt**2 +2*y(n) - y(n-1)
	We use of the first derivative of x and y to find the values of x(-1) and y(-1). The first 	   derivatives are given as 
	D_2t = [u(n+1) - u(n-1)]/2*dt = 0  
	Thus 
		x(-1) = x(1)  and y(-1) ) y(1)
	Hence the equtions for x(1) and y(1)  becomes
	x(1) = -0.5*(dt**2)*beta/(1-beta)(1-beta/L(0))x(0) +x(0)
	y(1) = -0.5*(dt**2)*beta/(1-beta)(1-beta/L(0))(y(0)-1) +y(0) - 0.5*beta*dt**2
	
	"""
	
	P= np.pi*2
	time_steps = int(P*time_steps_per_period)
	time = num_periods*P
	dt = float(time / time_steps)
	t = [ i*dt for i in range(time_steps)] 



	
	x = [ 0.0 for i in t]
	y = [ 0.0 for i in t]
	x[0] = (1. +epsilon)*math.sin(math.radians(Theta))
	y[0] = 1.- (1. +epsilon)*math.cos(math.radians(Theta))

	dt2=dt**2
	
	if beta ==1 :
	  	beta -= 0.001

	if Theta == 0:
		Theta-= 0.0001

	
 	aux = beta/(1-beta)
	
	L = lambda x,y : (x**2+(y-1.0)**2)**0.5	
	
	dt2=dt**2
	Linv = beta/L(x[0],y[0])

	x[1] = -0.5*dt2*aux*(1- Linv)*x[0]      +x[0]
	y[1] = -0.5*dt2*aux*(1- Linv)*(y[0]-1)  +y[0] - 0.5*beta*dt2

	for n in range(1,time_steps-1) : 
		Linv = beta/L(x[n],y[n])
		x[n+1] = -dt2*aux*(1.- Linv)*x[n] + 2.0*x[n] - x[n-1]
		y[n+1] = -dt2*aux*(1.- Linv)*(y[n]-1)  + 2.0*y[n] - y[n-1] - beta*dt2

	
	theta_func = lambda x,y : math.degrees(math.atan(x/(1-y)))  
	theta = map(theta_func,x,y)
	

	if plot : 
		visulasize_pendulum(x,y,save_id="beta%f_theta%f_epsilon%f"%(beta,Theta,epsilon))
		visulasize_pendulum_angle(t,theta,save_id="beta%.3f_theta%.3f_epsilon%.3f"%(beta,Theta,epsilon))
		if ( Theta < 10 ) :
			x_non,y_non,theta_non=simulate(1,Theta,epsilon, num_periods,  time_steps_per_period, plot=False)
			visulasize_compared_pendulums(theta,theta_non,t,save_id="beta%.3f_theta%.3f_epsilon%.3f"%(beta,Theta,epsilon))
		
			
	return x,y,theta

	
	
	
	def demo ( beta,Theta):
		"""
		Function demo takes 2 arguments 
		---------------------------------
		beta :
		Theta: The intial angle of the pendelum. 
		"""
	
		x,y,theta = simulate(beta,Theta,num_periods=3,time_steps_per_period=600)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	