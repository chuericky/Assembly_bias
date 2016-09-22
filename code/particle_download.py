### This program downloads the particle files after the queues are done.
### Input file: job id file (argv[1])
### Updated: Aug 29, 2015

import os, re, sys

user = "xxx"   ### You have to enter your login_id and password here.
pw = "xxx"
web = "https://www.cosmosim.org/uws/query"
#link = "/Users/rickyccy/Documents/Research_assemblybias/mock/particles/all/"
link = "/Volumes/Seagate Expansion Drive/Bolshoi_sim/"
link2 = "/Users/rickyccy/Documents/Research_assemblybias/mock/particles/"
#link = "/Users/rickyccy/Documents/Research_assemblybias/mock/Neal_request/particles/"
#link2 = "/Users/rickyccy/Documents/Research_assemblybias/mock/Neal_request/particles/par/"

os.chdir(link2)

infile = open(sys.argv[1], "r")

id = []

for line in infile:
    sline = line.split()
    id.append(sline[0])


infile.close()

os.chdir(link2)


for i in range(len(id)):
    os.system("uws --host %s --user %s --password %s job results %s csv" % (web, user, pw, id[i]))
    ### Delete the job to clear the space for the account.
    #os.system("uws --host %s --user %s --password %s job delete %s" % (web, user, pw, id[i]))
    print str(i + 1) + " files out of %d have been downloaded" % len(id)

del id
