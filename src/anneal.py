import sys, os, string, math, random
aas = 'ACDEFGHIKLMNPQRSTVWY'

def prAccept(E, E_new, T):
	if E_new < E:
		p = 1.0
	else:
		p = math.exp((float(E)-E_new)/T)
	return p

def getTemperature(T_0, time, max_time):
	return T_0/math.log(time+1,2)

# class Evaluator:
# 	def setInitialCandidate(self, initial_value)
#	def newCandidate(self)
# 	def acceptCandidate(self)
#	def rejectCandidate(self)
#	def onImprovement(self)
#	def onFinish(self)

# The central algorithm implementing simulated annealing.
def design(initial_value, evaluator, temperatureFxn, initial_temperature, max_time):
	energy_info = evaluator.setInitialCandidate(initial_value)
	E = energy_info.energy
	E_max = 0.0
	time = 1
	last_acceptance = time
	print "# Set initial candidate, t=%d, E=%1.5f" % (time, E)
	print time, "%1.4E" % temperatureFxn(initial_temperature, time, max_time), "%1.4f" % 0.0, E
	while time < max_time and E > E_max:
		time += 1
		(candidate, energy_info) = evaluator.newCandidate()
		E_new = energy_info.energy
		T = temperatureFxn(initial_temperature, time, max_time)
		p_move = prAccept(E, E_new, T)
		#print p_move
		# Accept move?
		rand_p = random.random()
		if rand_p < p_move:
			evaluator.acceptCandidate()
			# Status report; save progress
			if E != E_new:
				print "#", energy_info
				print time, "%1.4E" % T, "%1.4f" % p_move, E_new
				sys.stdout.flush()
				if E_new < E and time - last_acceptance > 1e4: # Last bit prevents huge amounts of writing early on when gains are easy
					evaluator.onImprovement()
					last_acceptance = time
			E = E_new
		else:
			# Reject move
			#print "# -- %s rand.p=%1.4f p.move=%1.4f" % (energy_info, rand_p, p_move)
			evaluator.rejectCandidate()

	final_value = evaluator.onFinish()
	return final_value

