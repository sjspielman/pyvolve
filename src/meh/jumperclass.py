










#########################################################################

class JumperEvolver(Evolver):
	''' Simulation strategy: use exponential waiting times and a jump chain to evolve along a branch; does not require re-calculation of P(t) '''
	def __init__(self, *args):
		super(JumperEvolver, self).__init__(*args)
		
		######## Initialization should include creating a jump matrix and saving the Q diagonal for each partition ########
		for i in range(self.NUMPARTS):
		
			self.PARTS[i][1].jumpMat = np.empty([61,61])
			
			# Extract -1*diagonal from instantaneous rate matrix and save
			self.PARTS[i][1].Qdiag = -1. * (np.diag( self.PARTS[i][1].Q ))
			
			for n in range(61):
				self.PARTS[i][1].jumpMat[n] = np.divide( self.PARTS[i][1].Q[n] , self.PARTS[i][1].Qdiag[n] )
				self.PARTS[i][1].jumpMat[n][n] = 0. 
				assert (abs (np.sum( self.PARTS[i][1].jumpMat[n] ) - 1. ) < self.ZERO), "Row in jump matrix does not sum to 1."
		
	
	def getSubRate(self, i, seq): 
		''' Retrieve substitution rate for a given codon, seq, in partition i  '''
		return self.PARTS[i][1].Qdiag[seq]



	def selectSite(self, normrates):
		''' Given a vector of normalized substitution rates, select one to evolve. This is basically the same as generateCodon, but we repeat it for readability. '''
		r = rn.uniform(0,1)
		i=0
		sum=normrates[i]
		while sum < r:
			i+=1
			sum+=normrates[i]
		return i		


	def evolve_branch(self, node, baseSeq):
		
		bl = self.checkBranch(node, baseSeq)
		if bl < self.ZERO:
			print bl, "branch length of 0 detected"
			node.seq = baseSeq
		
		else:		
		
			newSeq = np.empty(self.SEQLEN, dtype=int)
			start = 0 

			for i in range(self.NUMPARTS):

				# Starting sequence for this partition, again beginning with baseSeq
				tempSeq = baseSeq[ start: self.PARTS[i][0] ]
						
				# Retrieve site-wise and total substitution rates from tempSeq
				siteRates = np.empty( len(tempSeq) )
				for site_index in range( len(siteRates) ):
					siteRates[site_index] = self.getSubRate( i, tempSeq[site_index] )
				totalSubRate = np.sum(siteRates)
				
				total_time = rn.expovariate(totalSubRate)

				# We have an event
				while total_time < bl:
					
					# Randomly select a position to evolve and get its codon
					norm_siteRates = np.divide(siteRates, totalSubRate )
					site_index = self.selectSite(norm_siteRates)
					site_value = tempSeq[site_index]
					
					# Evolve the site					
					tempSeq[site_index] = self.generateCodon( self.PARTS[i][1].jumpMat[site_value] )						
	
					# Replace that index in siteRates to represent the position's new codon state			
					siteRates[site_index] = self.getSubRate(i, tempSeq[site_index])
					
					# Recalculate totalSubRate and from that, get a new waiting time
					totalSubRate = np.sum(siteRates)
					total_time += rn.expovariate(totalSubRate)		
					
					 
				newSeq[start : self.PARTS[i][0]] = tempSeq
				start += self.PARTS[i][0]
		
			# Attach final sequence to node
			node.seq = newSeq
	
	
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
				
			
			
			
			
			
		