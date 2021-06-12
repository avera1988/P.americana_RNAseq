#!/usr/bin/env python3
##################################################################################################
#	Scritp to Retrieve csv tables and treemap.R from REVIGO 
#	Author Arturo Vera
#	May 2021
###################################################################################################
import requests
import time
import sys

if len(sys.argv) != 4:
	print("Usage: python scrape_revigo.mod.py <Go_file> <Output file> <Semantical_Cluster_value (0.9,0.7,0.5 or 0.4)>")
	sys.exit (1)

Go=open(sys.argv[1],'r').read()
fBP=open(sys.argv[2]+".BP.csv",'w')
rBP=open(sys.argv[2]+".BP.R",'w')



# Submit job to Revigo
print('submit job to Revigo ...')
payload = {'cutoff':sys.argv[3], 'valueType':'pvalue', 'speciesTaxon':'0', 'measure':'SIMREL', 'goList':Go}
r = requests.post("http://revigo.irb.hr/StartJob.aspx", data=payload)

jobid = r.json()['jobid']

# Check job status
print('checking job ...')
running = 1
while (running!=0):
    r = requests.post("http://revigo.irb.hr/QueryJobStatus.aspx",data={'jobid':jobid})
    running = r.json()['running']
    time.sleep(1)

# Fetch results
rBPtable = requests.post("http://revigo.irb.hr/ExportJob.aspx", data={'jobid':jobid, 'namespace':'1', 'type':'csvtable'})
rRtree = requests.post("http://revigo.irb.hr/ExportJob.aspx", data={'jobid':jobid, 'namespace':'1', 'type':'rtree'})

print('Writing results ...')
# Write CSV table
fBP.write(rBPtable.text)
fBP.close()

#Write R treemap

rBP.write(rRtree.text)
rBP.close()

