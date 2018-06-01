import sys, os

file_input = open("cancella.txt", "r")
file_ouput = open("output.txt", "w")

curr_run = curr_lumi = curr_event = curr_pt = 0

for each_line in file_input:
	try:
		(bhu, row, run, lumi, event, mass, eta, pt, bhu2) = each_line.split("*")
		if(curr_run != run or curr_event != event or curr_lumi != lumi):
			curr_run = run
			curr_event = event
			curr_lumi = lumi
			curr_pt = pt
# 			if(mass > 120 and mass < 400):
			file_ouput.write(mass + '\t' + pt + '\t' + run + ':' + lumi + ':' + event + '\n')
		else:
			if(pt < curr_pt):
# 				print 'cancella'
				file_ouput.write('cancella secondo = ' + mass + '\t' + pt + '\t' + run + ':' + lumi + ':' + event + '\n')
			else:
				file_ouput.write('cancella prima = ' + mass + '\t' + pt + '\t' + run + ':' + lumi + ':' + event + '\n')
			
			
			
			
	except:
		pass