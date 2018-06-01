file = open('micro_ntuple.temp_Our2012_SORT.txt')
noScale = open('micro_ntuple.temp_Our2012_NoScale_SORT.txt')
for each_line, each_No in zip(file, noScale):
	run,lumi,event,mass = each_line.split("\t",3)
	run_no,lumi_no,event_no,mass_no = each_No.split("\t",3)
	if(float(mass) - float(mass_no) > 0.00001 or float(mass_no) - float(mass) > 0.00001):
		print '%s | %s' % (mass_no,mass)
